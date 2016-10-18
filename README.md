# ultrabubble_eval
Extract detailed ultra bubble stats from a vg graph.  Statistics are broken down into ultrabubbles and snarls (cyclclic ultrabubbles), nesting level, and size in bases and nodes.  Length distributions are optinally written as histograms.  

## Generating 1000 Genomes Results

Note: loading vg graphs into memory is RAM intensive.  Expect to use about 40G for chr1 and 700G for whole genome.  In practice, the graphs are chunked into smaller pieces when using the ultrabubble algorithm for defining sites for variant calling. 

### Download data

     wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
     wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz
	  wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz.tbi
     gunzip hs37d5.fa.gz

### Chromosome 1 Results

     vg construct -f -R 1 -r hs37d5.fa -v ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz -t 15 -m 100 > 1.vg
	  ./ub_eval 1.vg -i hist_1.tsv -m missing_1.tsv -l lengths_1.tsv > stats_1.tsv

### Whole Genome Results

This takes approximately ?? days and ??G of memory
    for i in $(seq 1 22; echo X; echo Y);
    do
         vg construct -f -R $i -r hs37d5.fa -v ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz -t 15 -m 100 >$i.vg 2>$i.err
    done
    vg ids -j $(for i in $(seq 22; echo X; echo Y); do echo $i.vg; done)
	 cat $(for i in $(seq 22; echo X; echo Y); do echo $i.vg; done) > wg.vg
	 ./ub_eval wg.vg -i hist_wg.tsv -m missing_wg.tsv -l lengths_wg.tsv > stats_wg.tsv

### Histograms

Desired histogram (for bubble type and depth) manually cut from hist_[1|wg].tsv and saved to separate file (including just the two data columns) then `histogram.py` used to plot.  Ex:

     ./histogram.py chr1_chop_hist_top_level_bub.txt --x_max 100 --bins 100 --save ub_sizes.pdf --log_counts --x_label "Top-level bubb
le size (bp)" --y_label "Log frequency" --no_n --title "" --width 3 --height 2.2

     ./scripts/histogram.py chr1_chop_hist_top_level_bub.txt --x_min 101 --bins 100 --save ub_sizes_tail.pdf  --x_label "Top-level bubble size
 (bp)" --y_label "Frequency" --no_n --title "" --sparse_ticks --width 3 --height 2.2

