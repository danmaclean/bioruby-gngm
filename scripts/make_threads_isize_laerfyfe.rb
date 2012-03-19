#!/usr/bin/env ruby
#
#  make_bands
#
#  Created by Dan MacLean (TSL) on 2012-01-17.
#  Copyright (c)  . All rights reserved.
###################################################

### An example script to loop over each reference in the BAM file, get SNP positions, make density threads for different kernels, 
### cluster for different values of k and then draw the threads, bands and signal. Generates plots for each new set of parameters 


$LOAD_PATH.unshift(File.join(File.dirname(__FILE__), '..', 'lib'))
$LOAD_PATH.unshift(File.dirname(__FILE__))
require 'bio-gngm'
require 'bio'
length = 0
chr_name = "" 
file = Bio::FastaFormat.open("/Users/macleand/Desktop/laerfyve_vs_stitched_ler/ler_contigs_stitched.fa")
file.each do |entry|
  length = entry.length
  chr_name = entry.entry_id
end

puts chr_name
  g = Bio::Util::Gngm.new(:file => "/Users/macleand/Desktop/laerfyve_vs_stitched_ler/aln.sort.bam", 
               :format => :bam, 
               :fasta => "/Users/macleand/Desktop/laerfyve_vs_stitched_ler/ler_contigs_stitched.fa",
               :samtools => {:q => 20, :Q => 50, :r => "#{chr_name}:1-#{length}"
               }
    )
      
      g.get_insert_size_frequency(:ref_window_size => 750, :ref_window_slide => 750, :isize => 184)
      g.collect_threads
      [0.25, 0.5, 1.0].each do |kernel_adjust|
        [4, 9, 11].each do | k |
          begin
              g.calculate_clusters(:k => k, :adjust => kernel_adjust, :control_chd => 0.5, :expected_chd => 0.2)
              filename = "isize_#{k}_#{kernel_adjust}_all_threads.png" 
              g.draw_threads(filename, :draw_legend => "isize_#{k}_#{kernel_adjust}_legend.png")
              filename = "isize_#{k}_#{kernel_adjust}_bands.png"
              g.draw_bands(filename)
              filename = "isize_#{k}_#{kernel_adjust}_signal.png"
              g.draw_signal(filename)
              filename = "isize_#{k}_#{kernel_adjust}_hits.png"
              g.draw_hit_count(filename)
            rescue  Exception => e
                puts "failed on #{k} #{kernel_adjust}"
                puts e.message, e.backtrace
            end
          end
        end