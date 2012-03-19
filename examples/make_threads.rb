#!/usr/bin/env ruby
#
#  make_threads
#
#  Created by Dan MacLean (TSL) on 2012-01-17.
#  Copyright (c)  . All rights reserved.
###################################################

### An example script to specify a region, get SNP positions, make density threads for different kernels, 
### cluster for different values of k and then draw the threads, bands and signal. Generates plots for each new set of parameters 


$LOAD_PATH.unshift(File.join(File.dirname(__FILE__), '..', 'lib'))
$LOAD_PATH.unshift(File.dirname(__FILE__))
require 'bio-gngm'


  g = Bio::Util::Gngm.new(:file => "/Users/macleand/Desktop/ngm/fli.sorted.bam", 
               :format => :bam, 
               :fasta => "/Users/macleand/Desktop/ngm/TAIR9_chr_all.fas", 
               :samtools => {:r => "Chr1:1000000-2000000",
                             :q => 20,
                             :Q => 50
               },
    )
    g.snp_positions( :min_depth => 10, :mapping_quality => 40.0, :min_non_ref_count => 5)
    g.collect_threads(:start => 0.3, :stop => 0.8, :slide => 0.01, :size => 0.2 )
    [0.25, 0.5, 1.0].each do |kernel_adjust|
      [4, 9, 11].each do | k |
        filename = "#{name}_#{k}_#{kernel_adjust}_all_threads.png"  
        g.calculate_clusters(:k => k, :adjust => kernel_adjust, :control_chd => 1.0, :expected_chd => 0.5)
        g.draw_threads(filename)
        filename = "#{name}_#{k}_#{kernel_adjust}_clustered_bands.png"
        g.draw_bands(filename)
        filename = "#{name}_#{k}_#{kernel_adjust}_signal.png"
        g.draw_signal(filename)
      end
    end
    g.close
end

sam.close
