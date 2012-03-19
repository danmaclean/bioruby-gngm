#!/usr/bin/env ruby
#
#  make_threads_isize
#
#  Created by Dan MacLean (TSL) on 2012-01-17.
#  Copyright (c)  . All rights reserved.
###################################################

### An example script to luse th insert size between paired end reads to infer large deletions. The script, make density threads based on the proportion of reads with insert size above 0.5 for different kernels, 
### cluster for different values of k and then draw the threads, bands and signal. Generates plots for each new set of parameters. Contains a rescue clause in case something goes wrong, these sorts of calculations
### can be very data-sparse and density curves can be hard to make and cluster when lots of window are empty 

require 'bio-gngm'


  g = Bio::Util::Gngm.new(:file => "aln.sort.bam", 
               :format => :bam, 
               :fasta => "reference.fasta", 
               :samtools => {:r => "1:1000000-10000000",
                             :q => 20,
                             :Q => 50
               }
    )
    g.get_insert_size_frequency(:ref_window_size => 152, :ref_window_slide => 152, :isize => 150)
    g.collect_threads
    [0.25, 0.5, 1.0].each do |kernel_adjust|
      [4, 9, 11].each do | k |
        begin
            g.calculate_clusters(:k => k, :adjust => kernel_adjust, :control_chd => 0.5, :expected_chd => 0.9)
            filename = "isize_#{k}_#{kernel_adjust}_all_threads.png" 
            g.draw_threads(filename, :draw_legend => "isize_#{k}_#{kernel_adjust}_legend.png")
            filename = "isize_#{k}_#{kernel_adjust}_bands.png"
            g.draw_bands(filename)
            filename = "isize_#{k}_#{kernel_adjust}_signal.png"
            g.draw_signal(filename)
          rescue  Exception => e 
              puts "failed on #{k} #{kernel_adjust}"
              puts e.message, e.backtrace
          end
        end
      end
