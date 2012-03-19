#!/usr/bin/env ruby
#
#  use_indels
#
#  Created by Dan MacLean (TSL) on 2012-01-17.
#  Copyright (c)  . All rights reserved.
###################################################

### An example script to specify a region, get small deletion positions, make density threads for different kernels, 
### cluster for different values of k and then draw the threads, bands and signal. Generates plots for each new set of parameters

require 'bio-gngm'



  g = Bio::Util::Gngm.new(:file => "aln.sorted.bam", 
               :format => :bam, 
               :fasta => "reference.fasta", 
               :samtools => {:r => "Chr1:1-3000000",
                             :q => 20,
                             :Q => 50
               }
    )
    g.snp_positions(:deletions_only => true)
    g.collect_threads 
    g.calculate_clusters(:k => 9, :adjust => 0.5, :control_chd => 0.5, :expected_chd => 1.0)
    filename = "indels_all_threads.png" 
    g.draw_threads(filename)
    filename = "indels_clustered_bands.png"
    g.draw_bands(filename)
    filename = "indels_signal.png"
    g.draw_signal(filename)
    filename = "indels_peaks.png"
    g.draw_peaks(filename)
    g.close

