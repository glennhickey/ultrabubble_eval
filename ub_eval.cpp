/*
 * Copyright (C) 2015 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <iostream>
#include <fstream>
#include <getopt.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "src/vg.hpp"

using namespace vg;
using namespace std;


void help_main(char** argv)
{
  cerr << "usage: " << argv[0] << " [options] VGFILE OUTPUT_DIR" << endl
       << "Generate some reports on bubble decomposition of given graph" << endl
       << "    -h, --help          print this help message" << endl;
}

// aggregate some bubble stats (count and length)
struct BubbleStats {
   
   struct Tally {
      int count; // number of bubbles counted
      int total_length; //length stats
      int min_length;
      int max_length;
      int total_nc; // node-count stats
      int min_nc;
      int max_nc;


      Tally() : count(0), total_length(0), min_length(999999), max_length(0),
        total_nc(0), min_nc(999999), max_nc(0) {}

      // not sure why, but need to add std:: and vg:: here to compile despite using the namespaces...
      void add_bubble(VG& graph, std::pair<vg::id_t, vg::id_t> ends, const vector<vg::id_t>& bubble) {
        int nc = 0;
        int length = 0;
        for (auto nid : bubble) {
          ++nc;
          length += graph.get_node(nid)->sequence().length();
        }
        ++count;
        total_length += length;
        min_length = min(min_length, length);
        max_length = max(max_length, length);
        total_nc += nc;
        min_nc = min(min_nc, nc);
        max_nc = max(max_nc, nc);
      }
   };
   
   vector<Tally> tally_map; // stats sfor given depth of nestedness
   std::map<vg::id_t, int> visit_count;  // keep track of depth (make sure we visit bubbles in order of size)

   void add_bubble(VG& graph, std::pair<vg::id_t, vg::id_t> ends, const vector<vg::id_t>& bubble) {

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
       node = graph.get_node(ends.first);
       if (visit_count.count(node->id())) {
         depth = visit_count[node->id()];
         assert (depth == visit_count[ends.second]);
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

   void add_bubble(VG& graph, std::pair<vg::id_t, vg::id_t> ends, const vector<vg::id_t>& bubble, int depth) {
     // add bubble stats to tracker for the given depth
     while (tally_map.size() < depth + 1) {
       tally_map.push_back(Tally());
     }
     tally_map[depth].add_bubble(graph, ends, bubble);
     // add bubble stats for overall count
     tally_map[0].add_bubble(graph, ends, bubble);            
   }
   
   struct bub_less {
      bool operator()(const std::vector<vg::id_t>& b1, const std::vector<vg::id_t>& b2) const {
          return b1.size() < b2.size();
        }
   };
   
   void compute_stats(VG& graph, const std::map<std::pair<vg::id_t, vg::id_t>, std::vector<vg::id_t>>& bubbles) {
     // we want to visit bubbles in decreasing size (measured in number of nodes)
     // this will guarantee that if bubble A is nested inside B, then B is visited
     // first, and our lookup for computing depth will make sense
     std::multimap<std::vector<vg::id_t>, std::pair<vg::id_t, vg::id_t>, bub_less> sorted_bubbles;
     for (auto b : bubbles) {
       // strip off endpoints , we don't want to count them because they aren't *in* the bubble
       vector<vg::id_t> bub;
       for (auto u : b.second) {
         if (u != b.first.first && u != b.first.second) {
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
   }
};

ostream& operator<<(ostream& os, const BubbleStats::Tally& t) {
  os << t.count << "\t"
     << t.total_length << "\t"
     << t.min_length << "\t"
     << t.max_length << "\t"
     << (t.count ? (double)t.total_length / t.count : 0) << "\t"
     << t.total_nc << "\t" 
     << t.min_nc << "\t"
     << t.max_nc << "\t"
     << (t.count ? (double)t.total_nc / t.count : 0);
  return os;
}

ostream& operator<<(ostream& os, BubbleStats& t) {
  os << "Depth\t" << "Count\t"
     << "Total-Length\t" << "Min-Length\t" << "Max-Length\t" << "Avg-Length\t"
     << "Total-Nodes\t" << "Min-Nodes\t" << "Max-Nodes\t" << "Avg-Nodes"
     << endl;
  for (int i = 0; i < t.tally_map.size(); ++i) {
    os << i << "\t" << t.tally_map[i] << endl;
  }
  return os;
}

void cycle_stats(VG& graph, BubbleTree* bubble_tree)
{
  cerr << "Computing cycle stats" << endl;

  // start by splitting bubbles into two groups by cyclicity
  std::map<std::pair<vg::id_t, vg::id_t>, std::vector<vg::id_t> > acyclic_bubbles;
  std::map<std::pair<vg::id_t, vg::id_t>, std::vector<vg::id_t> > cyclic_bubbles;
  
  bubble_tree->for_each_preorder([&](BubbleTree::Node* node) {
      // cut root to be consistent with superbubbles()
      if (node != bubble_tree->root) {
        Bubble& bubble = node->v;
        // sort nodes to be consistent with superbubbles
        sort(bubble.contents.begin(), bubble.contents.end());
        if (bubble.acyclic) {
          acyclic_bubbles[make_pair(bubble.start.node, bubble.end.node)] = bubble.contents;
        } else {
          cyclic_bubbles[make_pair(bubble.start.node, bubble.end.node)] = bubble.contents;
        }
      }
    });

  // now just rerun our stats separately on each group
  BubbleStats acyclic_bs;
  acyclic_bs.compute_stats(graph, acyclic_bubbles);
  BubbleStats cyclic_bs;
  cyclic_bs.compute_stats(graph, cyclic_bubbles);

  cout << "Acyclic ultra bubbles stats" << endl << acyclic_bs << endl;  
  cout << "Cyclic ultra bubbles stats" << endl << cyclic_bs << endl;
}

void chain_stats(VG& graph, BubbleTree* bubble_tree) {
  cerr << "Computing chain stats" << endl;

  auto get_depth = [&](BubbleTree::Node* node) {
    int depth = 0;
    for (; node != bubble_tree->root; ++depth, node = node->parent);
    return depth;
  };
  
  BubbleStats chains_bs;

  bubble_tree->for_each_preorder([&](BubbleTree::Node* node) {

      Bubble& bubble = node->v;
      for (int i = 0; i < bubble.chain_offsets.size(); ++i)
      {
        int last = i == bubble.chain_offsets.size() - 1 ? node->children.size() - 1 : bubble.chain_offsets[i+1] - 1;
        vector<vg::id_t> chain_nodes;
        for (int j = bubble.chain_offsets[i]; j <= last; ++j) {
          chain_nodes.push_back(node->children[i]->v.start.node);
          if (j == last) {
            chain_nodes.push_back(node->children[i]->v.end.node);
          }
        }

        chains_bs.add_bubble(graph, std::pair<vg::id_t, vg::id_t>(0, 0), chain_nodes, get_depth(node));
      }      
    });

  cout << "Chains stats" << endl << chains_bs << endl;
}


void ultra_stats(VG& graph, const string& out_dir)
{
  cerr << "Computing ultrabubbles" << endl;
  auto bubbles = vg::ultrabubbles(graph);
  cerr << "Running ultra stats on " << bubbles.size() << " bubbles" << endl;
  BubbleStats bs;
  bs.compute_stats(graph, bubbles);
  cout << "Ultra bubbles stats" << endl << bs << endl;
  
  cerr << "Computing ultrabubbles" << endl;
  BubbleTree* bubble_tree = ultrabubble_tree(graph);

  cycle_stats(graph, bubble_tree);
  chain_stats(graph, bubble_tree);
  delete bubble_tree;
}

void super_stats(VG& graph, const string& out_dir)
{
  cerr << "Computing superbubbles" << endl;
  auto bubbles = vg::superbubbles(graph);
  cerr << "Running super stats on " << bubbles.size() << " bubbles" << endl;
  BubbleStats bs;
  bs.compute_stats(graph, bubbles);
  cout << "Super bubbles stats" << endl << bs << endl;
  
}

int main(int argc, char** argv)
{
    
  if(argc == 1) {
    // Print the help
    help_main(argv);
    return 1;
  }
    
  optind = 1; // Start at first real argument
  bool optionsRemaining = true;
  while(optionsRemaining) {
    static struct option longOptions[] = {
      {"help", no_argument, 0, 'h'},
      {0, 0, 0, 0}
    };

    int optionIndex = 0;

    switch(getopt_long(argc, argv, "w:o:h", longOptions, &optionIndex)) {
      // Option value is in global optarg
    case -1:
      optionsRemaining = false;
      break;
    case 'h': // When the user asks for help
      help_main(argv);
      exit(1);
      break;
    default:
      cerr << "Illegal option" << endl;
      exit(1);
    }
  }

  if(argc - optind < 2) {
    // We don't have two positional arguments
    // Print the help
    help_main(argv);
    return 1;
  }

  string vg_file = argv[optind++];
  string out_dir = argv[optind++]; 
    
  // Open the vg file
  cerr << "Loading vg" << endl;
  ifstream vgStream(vg_file);
  if(!vgStream.good())
  {
    cerr << "Could not read " << vg_file << endl;
    exit(1);
  }
  VG graph(vgStream);
  cerr << "Sorting vg" << endl;
  graph.sort();

  // Make the output directory
  mkdir(out_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

  ultra_stats(graph, out_dir);
  super_stats(graph, out_dir);
  
  return 0;
}

