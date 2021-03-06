/*
 * Copyright (C) 2015 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <iostream>
#include <fstream>
#include <getopt.h>
#include <chrono>
#include <ctime>

#include "src/vg.hpp"

using namespace vg;
using namespace std;


// time as string
string timestring()
{
    std::chrono::time_point<std::chrono::system_clock> t;
    t = std::chrono::system_clock::now();
    std::time_t ct = std::chrono::system_clock::to_time_t(t);
    return std::ctime(&ct);
}

// get total length, total nodes in graph
pair<int, int> graph_stats(VG& graph)
{
    int total_length = 0;
    int total_nodes = 0;
    graph.for_each_node([&](Node* node) {
            total_length += node->sequence().length();
            ++total_nodes;
        });
    return make_pair(total_length, total_nodes);
}

// aggregate some bubble stats (count and length)
struct BubbleStats {

    struct Tally {
        int count; // number of bubbles counted
        int total_length; //length stats
        int min_length;
        int max_length;
        map<int, int> length_hist;        
        int total_nc; // node-count stats
        int min_nc;
        int max_nc;
        map<int, int> nc_hist;

        Tally() : count(0), total_length(0), min_length(999999), max_length(0),
                  total_nc(0), min_nc(999999), max_nc(0) {}

        // add a bubble's length and node count to our histogram (which doesn't bin)
        void update_hist(int length, int nc) {
            if (length_hist.count(length)) {
                ++length_hist[length];
            } else {
                length_hist[length] = 1;
            }
            if (nc_hist.count(nc)) {
                ++nc_hist[nc];
            } else {
                nc_hist[nc] = 1;
            }
        }
        
        // not sure why, but need to add std:: and vg:: here to compile despite using the namespaces...
        void add_bubble(VG& graph, std::pair<vg::NodeSide, vg::NodeSide> ends, const vector<vg::id_t>& bubble) {
            int nc = 0;
            int length = 0;
            for (auto nid : bubble) {
                ++nc;
                length += graph.get_node(nid)->sequence().length();
            }
            if (length > 0) {
                ++count;
                total_length += length;
                min_length = min(min_length, length);
                max_length = max(max_length, length);
                total_nc += nc;
                min_nc = min(min_nc, nc);
                max_nc = max(max_nc, nc);
                update_hist(length, nc);
            }
        }
    };
   
    vector<Tally> tally_map; // stats sfor given depth of nestedness
    std::map<vg::id_t, int> visit_count;  // keep track of depth (make sure we visit bubbles in order of size)
    pair<int, int> gs; // stats for whole graph (total length, total nodes)
    // dont include in tally (but keep for depth computations), used to break out cyclic/dag
    set<pair<vg::NodeSide, vg::NodeSide>> ignore_ends; 

    void add_bubble(VG& graph, std::pair<vg::NodeSide, vg::NodeSide> ends, const vector<vg::id_t>& bubble) {

        // figure out how many times we've seen the nodes in the bubble
        int depth = 1;
        Node* node = NULL;
        if (bubble.size() > 0) {
            node = graph.get_node(bubble[0]);
            if (visit_count.count(node->id())) {
                depth = visit_count[node->id()];
            }
        }
        // for empty bubble, we infer depth from endpoints
        else {
            node = graph.get_node(ends.first.node);
            if (visit_count.count(node->id())) {
                depth = visit_count[node->id()];
                assert (depth == visit_count[ends.second.node]);
            }
        }

        // sanity checks to make sure all nodes in bubble have same depth
        for (auto i : bubble) {
            if (visit_count.count(node->id())) {
                assert(visit_count[i] == depth);
            }
            else {
                assert(visit_count.count(i) == false);
            }
        }
     
        // increase by one
        for (auto i : bubble) {
            visit_count[i] = depth + 1;
        }

        add_bubble(graph, ends, bubble, depth);
    }

    void add_bubble(VG& graph, std::pair<vg::NodeSide, vg::NodeSide> ends, const vector<vg::id_t>& bubble, int depth) {
        if (!ignore_ends.count(ends) && !ignore_ends.count(make_pair(ends.second, ends.first))) {   
            // add bubble stats to tracker for the given depth
            while (tally_map.size() < depth + 1) {
                tally_map.push_back(Tally());
            }
            tally_map[depth].add_bubble(graph, ends, bubble);
            // add bubble stats for overall count
            tally_map[0].add_bubble(graph, ends, bubble);
        }
    }
   
    struct bub_less {
        bool operator()(const std::vector<vg::id_t>& b1, const std::vector<vg::id_t>& b2) const {
            return b1.size() < b2.size();
        }
    };
    
    void compute_stats(VG& graph, const std::map<std::pair<vg::NodeSide, vg::NodeSide>, std::vector<vg::id_t>>& bubbles) {
        // we want to visit bubbles in decreasing size (measured in number of nodes)
        // this will guarantee that if bubble A is nested inside B, then B is visited
        // first, and our lookup for computing depth will make sense
        std::multimap<std::vector<vg::id_t>, std::pair<vg::NodeSide, vg::NodeSide>, bub_less> sorted_bubbles;
        for (auto b : bubbles) {
            // strip off endpoints , we don't want to count them because they aren't *in* the bubble
            vector<vg::id_t> bub;
            for (auto u : b.second) {
                if (u != b.first.first.node && u != b.first.second.node) {
                    bub.push_back(u);
                }
            }
            assert(bub.size() == b.second.size() - 2);
       
            sorted_bubbles.insert(make_pair(bub, b.first));
        }
        assert(sorted_bubbles.size() == bubbles.size());
        tally_map.clear();
        visit_count.clear();
        for (auto i = sorted_bubbles.rbegin(); i != sorted_bubbles.rend(); ++i) {
            add_bubble(graph, i->second, i->first);
        }
        // this really should be moved outside to avoid doing over and over
        gs = graph_stats(graph);
    }

    // hack to let superbubbles run through same interface even though we dont get sides out of it. 
    void compute_stats(VG& graph, const std::map<std::pair<vg::id_t, vg::id_t>, std::vector<vg::id_t>>& bubbles) {
        std::map<std::pair<vg::NodeSide, vg::NodeSide>, std::vector<vg::id_t>> ns_bubbles;
        for (auto b : bubbles) {
            ns_bubbles[make_pair(NodeSide(b.first.first), NodeSide(b.first.second))] = b.second;
        }
        compute_stats(graph, ns_bubbles);
    }
};

ostream& operator<<(ostream& os, BubbleStats& bs) {
    os << "Depth\t" << "Count\t"
       << "Total-Length\t" << "Min-Length\t" << "Max-Length\t" << "Avg-Length\t"
       << "Frac-Length (/" << bs.gs.first << ")\t"
       << "Total-Nodes\t" << "Min-Nodes\t" << "Max-Nodes\t" << "Avg-Nodes\t"
       << "Frac-Nodes (/" << bs.gs.second << ")\t"
       << endl;
    for (int i = 0; i < bs.tally_map.size(); ++i) {
        const BubbleStats::Tally& t = bs.tally_map[i];
        os << (i == 0 ? "all" : to_string(i)) << "\t"

           << t.count << "\t"
           << t.total_length << "\t"
           << t.min_length << "\t"
           << t.max_length << "\t"
           << (t.count ? (double)t.total_length / t.count : 0) << "\t"
           << (bs.gs.first ? (double)t.total_length / bs.gs.first : 0) << "\t"
           << t.total_nc << "\t" 
           << t.min_nc << "\t"
           << t.max_nc << "\t"
           << (t.count ? (double)t.total_nc / t.count : 0) << "\t"
           << (bs.gs.second ? (double)t.total_nc / bs.gs.second : 0)

           << endl;
    }
    return os;
}

void histograms(ostream* hist, BubbleStats& bs, int bin_size) {
    if (hist == NULL){
        return;
    }
    ostream& os = *hist;

    function<map<int, int>(map<int, int>)> make_bins = [&](map<int, int> hist) -> map<int, int>{
        map<int, int> binned_hist;
        for (auto i : hist) {
            int bin = i.first / bin_size;
            if (binned_hist.count(bin)) {
                binned_hist[bin] += i.second;
            } else {
                binned_hist[bin] = i.second;
            }
        }
        return binned_hist;
    };
    
    for (int i = 0; i < bs.tally_map.size(); ++i) {
        const BubbleStats::Tally& t = bs.tally_map[i];

        os << "Length Histogram for depth " << (i == 0 ? "all" : to_string(i)) << "\n";
        map<int, int> binned_length = make_bins(t.length_hist);
        for (auto i : binned_length) {
            os << i.first << "\t" << i.second << "\n";
        }
        os << endl;
        os << "Node Histogram for depth " << (i == 0 ? "all" : to_string(i)) << "\n";
        map<int, int> binned_nc = make_bins(t.nc_hist);
        for (auto i : binned_nc) {
            os << i.first << "\t" << i.second << "\n";
        }
        os << endl;
    }
}

void cycle_stats(VG& graph, BubbleTree* bubble_tree, int size_cap, ostream* hist, int bin_size)
{
    cerr << "Computing cycle stats" << endl;

    // start by splitting bubbles into two groups by cyclicity
    std::map<std::pair<vg::id_t, vg::id_t>, std::vector<vg::id_t> > acyclic_bubbles;
    std::map<std::pair<vg::id_t, vg::id_t>, std::vector<vg::id_t> > cyclic_bubbles;

    BubbleStats acyclic_bs;
    BubbleStats cyclic_bs;

    // length of nodes that aren't endpoints in bubble
    function<size_t(Bubble&)> bub_length = [&](Bubble& bubble) {
        size_t len = 0;
        for (auto bi : bubble.contents) {
            if (bi != bubble.start.node && bi != bubble.end.node) {
                Node* node = graph.get_node(bi);
                len += node->sequence().length();
            }
        }
        return len;
    };

    bubble_tree->for_each_preorder([&](BubbleTree::Node* node) {
            // cut root to be consistent with superbubbles()
            if (node != bubble_tree->root) {
                Bubble& bubble = node->v;
                // sort nodes to be consistent with superbubbles
                sort(bubble.contents.begin(), bubble.contents.end());
                // all bubbles get added to each stats (so we can compute depth)
                // but we use the ignore_ends member to make sure we conly count ones we want
                acyclic_bubbles[make_pair(bubble.start.node, bubble.end.node)] = bubble.contents;
                cyclic_bubbles[make_pair(bubble.start.node, bubble.end.node)] = bubble.contents;
                if (bubble.dag) {
                    cyclic_bs.ignore_ends.insert(make_pair(bubble.start.node, bubble.end.node));
                    
                    // optionally ignore acyclic bubbles bigger than size_cap
                    if (size_cap >= 0 && bub_length(bubble) > size_cap) {
                        acyclic_bs.ignore_ends.insert(make_pair(bubble.start.node, bubble.end.node));
                    }
                } else {
                    acyclic_bs.ignore_ends.insert(make_pair(bubble.start.node, bubble.end.node));
                }
            }
        });

    // now just rerun our stats separately on each group
    acyclic_bs.compute_stats(graph, acyclic_bubbles);
    if (size_cap >= 0) {
        cout << "Acyclic ultra bubbles with length <= " << size_cap << " stats" << endl << acyclic_bs << endl;
    } else {
        cout << "Acyclic ultra bubble stats" << endl << acyclic_bs << endl;
        cyclic_bs.compute_stats(graph, cyclic_bubbles);
        cout << "Snarl stats" << endl << cyclic_bs << endl;

        if (hist) {
          *hist << "Acyclic ultra bubbles histograms" << endl;
        }
        histograms(hist, acyclic_bs, bin_size);

        if (hist) {
          *hist << "Snarls histograms" << endl;
        }
        histograms(hist, cyclic_bs, bin_size);
    }
}

BubbleStats chain_stats(VG& graph, BubbleTree* bubble_tree, ostream* hist, int bin_size) {
    cerr << "Computing chain stats" << endl;

    auto get_depth = [&](BubbleTree::Node* node) {
        int depth = 0;
        for (; node != bubble_tree->root; ++depth, node = node->parent);
        return depth;
    };
  
    BubbleStats chains_bs;
    chains_bs.gs = graph_stats(graph);

    bubble_tree->for_each_preorder([&](BubbleTree::Node* node) {

            Bubble& bubble = node->v;
            for (int i = 0; i < bubble.chain_offsets.size(); ++i)
            {
                int last = i == bubble.chain_offsets.size() - 1 ? node->children.size() - 1 : bubble.chain_offsets[i+1] - 1;
                vector<vg::id_t> chain_nodes;
                for (int j = bubble.chain_offsets[i] + 1; j <= last; ++j) {
                    // only look at edges flanked on both sides by bubbles in the chain
                    vg::id_t chain_node = 0;
                    if (node->children[j-1]->v.end.node == node->children[j]->v.start.node ||
                        node->children[j-1]->v.end.node == node->children[j]->v.end.node) {
                        chain_node = node->children[j-1]->v.end.node;
                    } else if (node->children[j-1]->v.start.node == node->children[j]->v.start.node ||
                               node->children[j-1]->v.start.node == node->children[j]->v.end.node) {
                        chain_node = node->children[j-1]->v.start.node;
                    }
                    assert(chain_node != 0);
                    chain_nodes.push_back(chain_node);
                }
                if (chain_nodes.size() > 0) {
                    chains_bs.add_bubble(graph, std::pair<vg::NodeSide, vg::NodeSide>(NodeSide(0), NodeSide(0)),
                                         chain_nodes, get_depth(node) + 1);
                }
            }      
        });

    cout << "Chains stats" << endl << chains_bs << endl;
    if (hist) {
        *hist << "Chains Histograms" << endl;
    }
    histograms(hist, chains_bs, bin_size);
    return chains_bs;
}

BubbleStats missing_stats(VG& graph, BubbleTree* bubble_tree, ostream* hist, int bin_size, ostream* missing) {
    cerr << "Computing missing stats" << endl;
    std::set<vg::id_t> found_nodes;
      
    BubbleStats missing_bs;
    missing_bs.gs = graph_stats(graph);

    for (auto tn : bubble_tree->root->children) {
        Bubble& bubble = tn->v;
        // insert every node inside a bubble
        for (auto bn : bubble.contents) {
            if (bn != bubble.start.node && bn != bubble.end.node) {
                found_nodes.insert(bn);
            }
        }
    }

    // todo: move this code to bubbles.hpp?  or otherwise merge with logic in chains_stats here.
    // insert every chain node
    BubbleTree::Node* node = bubble_tree->root;
    Bubble& bubble = node->v;
    for (int i = 0; i < bubble.chain_offsets.size(); ++i)
    {
        int last = i == bubble.chain_offsets.size() - 1 ? node->children.size() - 1 : bubble.chain_offsets[i+1] - 1;
        vector<vg::id_t> chain_nodes;
        for (int j = bubble.chain_offsets[i] + 1; j <= last; ++j) {
            // only look at edges flanked on both sides by bubbles in the chain
            vg::id_t chain_node = 0;
            if (node->children[j-1]->v.end.node == node->children[j]->v.start.node ||
                node->children[j-1]->v.end.node == node->children[j]->v.end.node) {
                chain_node = node->children[j-1]->v.end.node;
            } else if (node->children[j-1]->v.start.node == node->children[j]->v.start.node ||
                       node->children[j-1]->v.start.node == node->children[j]->v.end.node) {
                chain_node = node->children[j-1]->v.start.node;
            }
            assert(chain_node != 0);      
            found_nodes.insert(chain_node);
        }
    }

    vector<vg::id_t> missing_nodes;
    graph.for_each_node([&](Node* node) {
            if (found_nodes.count(node->id()) == false) {
                missing_nodes.push_back(node->id());
            }
        });
  
    for (auto node : missing_nodes) {
        missing_bs.add_bubble(graph, std::pair<vg::NodeSide, vg::NodeSide>(NodeSide(0), NodeSide(0)),
                              vector<vg::id_t>(1, node), 1);
    }

    cout << "Missing stats" << endl << missing_bs << endl;
    if (hist) {
        *hist << "Missing Histograms" << endl;
    }
    histograms(hist, missing_bs, bin_size);

    if (missing != NULL) {
        for (auto n : missing_nodes) {
            *missing << n << endl;
        }
    }
    return missing_bs;
}

void overall_stats(VG& graph, BubbleStats& bs, BubbleStats& cs)
{
    int len = bs.tally_map[1].total_length + cs.tally_map[1].total_length;
    int nodes = bs.tally_map[1].total_nc + cs.tally_map[1].total_nc;
    cout << "Summary Stats" << endl;
    cout << "Top Level Chains + Bubbles (Length)" << "\t"
         << "Top Level Chains + Bubbles (Length Frac/" << bs.gs.first << ")\t"
         << "Top Level Chains + Bubbles (Nodes)" << "\t"
         << "Top Level Chains + Bubbles (Nodes Frac/" << bs.gs.second << endl;
    cout << len  << "\t"
         << ((double)len / bs.gs.first) << "\t"
         << nodes  << "\t"
         << ((double)nodes / bs.gs.second) << endl;
    cout << "Missing Nodes" << "\t" << "Missing Length" << endl;
    cout << (bs.gs.second - nodes) << "\t" << (bs.gs.first - len) << endl;
}

void ultra_stats(VG& graph, ostream* hist, int bin_size, ostream* missing, int size_cap)
{
    cerr << "Computing ultrabubbles " << timestring() << endl;
    auto bubbles = vg::ultrabubbles(graph);
    cerr << "Running ultra stats on " << bubbles.size() << " bubbles " << timestring() << endl;
    BubbleStats bs;
    bs.compute_stats(graph, bubbles);
    cout << "Ultra bubbles stats" << endl << bs << endl;
    if (hist) {
        *hist << "Bubble Histograms" << endl;
    }
    histograms(hist, bs, bin_size);
    
    cerr << "Computing ultrabubbles" << endl;
    BubbleTree* bubble_tree = ultrabubble_tree(graph);

    cycle_stats(graph, bubble_tree, -1, hist, bin_size);
    cycle_stats(graph, bubble_tree, size_cap, hist, bin_size); 
    BubbleStats cs = chain_stats(graph, bubble_tree, hist, bin_size);
    missing_stats(graph, bubble_tree, hist, bin_size, missing);  
    overall_stats(graph, bs, cs);
    
    delete bubble_tree;
}

void super_stats(VG& graph, ostream* hist, int bin_size)
{
    cerr << "Computing superbubbles" << endl;
    auto bubbles = vg::superbubbles(graph);
    cerr << "Running super stats on " << bubbles.size() << " bubbles" << endl;
    BubbleStats bs;
    bs.compute_stats(graph, bubbles);
    cout << "Super bubbles stats" << endl << bs << endl;
  
}

void help_main(char** argv)
{
    cerr << "usage: " << argv[0] << " [options] VGFILE" << endl
         << "Generate some reports on bubble decomposition of given graph" << endl
         << "    -h, --help          print this help message" << endl
         << "    -b, --bin N         bin size for histograms [1]" << endl
         << "    -i, --hist FILE     file name to write all histograms to in tsv format" << endl
         << "    -m, --missing FILE  file name to write all node ids that are not in covered by chain or bubble" << endl
         << "    -l, --length N      print separate statistics of bubbles <= N bases [default=100]" << endl;
}

int main(int argc, char** argv)
{
    
    if(argc == 1) {
        // Print the help
        help_main(argv);
        return 1;
    }

    int bin_size = 1;
    string hist_path;
    string missing_path;
    int size_cap = 100;
    
    optind = 1; // Start at first real argument
    bool optionsRemaining = true;
    while(optionsRemaining) {
        static struct option longOptions[] = {
            {"help", no_argument, 0, 'h'},
            {"bin", required_argument, 0, 'b'},
            {"hist", required_argument, 0, 'i'},
            {"missing", required_argument, 0, 'm'},
            {"length", required_argument, 0, 'l'},
            {0, 0, 0, 0}
        };

        int optionIndex = 0;

        switch(getopt_long(argc, argv, "w:o:hb:i:m:l:", longOptions, &optionIndex)) {
            // Option value is in global optarg
        case -1:
            optionsRemaining = false;
            break;
        case 'h': // When the user asks for help
            help_main(argv);
            exit(1);
            break;
        case 'b':
            bin_size = atoi(optarg);
            break;
        case 'i':
            hist_path = optarg;
            break;
        case 'm':
            missing_path = optarg;
            break;
        case 'l':
            size_cap = atoi(optarg);
            break;            
        default:
            cerr << "Illegal option" << endl;
            exit(1);
        }
    }

    if(argc - optind < 1) {
        // We don't have one positional arguments
        // Print the help
        help_main(argv);
        return 1;
    }

    string vg_file = argv[optind++];

    ostream* hist = NULL;
    ofstream hist_file;
    if (!hist_path.empty()) {
        hist_file.open(hist_path.c_str());
        hist = &hist_file;
    }

    ostream* missing = NULL;
    ofstream missing_file;
    if (!missing_path.empty()) {
        missing_file.open(missing_path.c_str());
        missing = &missing_file;
    }
    
    // Open the vg file
    cerr << "Loading vg " << timestring() << endl;
    ifstream vgStream(vg_file);
    if(!vgStream.good())
    {
        cerr << "Could not read " << vg_file << endl;
        exit(1);
    }
    VG graph(vgStream);
    cerr << "Finished loading v " << timestring() << endl;

    ultra_stats(graph, hist, bin_size, missing, size_cap);
    super_stats(graph, hist, bin_size);
  
    return 0;
}

