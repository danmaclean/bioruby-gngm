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

ctrl_thread_values = []
expected_thread_values = []
interval_width = 10000000
#puts chr_name
(1..length).step(interval_width) do |start|
  stop = start + interval_width
  region = "#{chr_name}:#{start}-#{stop}"
  
  puts "analyzing - #{region}"

  g = Bio::Util::Gngm.new(:file => "/Users/macleand/Desktop/laerfyve_vs_stitched_ler/aln.sort.bam", 
               :format => :bam, 
               :fasta => "/Users/macleand/Desktop/laerfyve_vs_stitched_ler/ler_contigs_stitched.fa",
               :samtools => {:q => 20, :Q => 50, :r => region}
    )
      
      g.get_unmapped_mate_frequency(:ref_window_size => 76, :ref_window_slide => 76)
      g.collect_threads(:start => 0.0, :stop => 0.5, :slide => 0.1, :size => 0.1)
      #puts g.threads
      [0.5].each do |kernel_adjust|
        [4].each do | k |
          begin  
              g.calculate_clusters(:k => k, :adjust => kernel_adjust, :control_chd => 0.0, :expected_chd => 0.4, :pseudo => false)
              #filename = "unmapped_2_#{region}_#{k}_#{kernel_adjust}_all_threads.png" 
              #g.draw_threads(filename, :draw_legend => "unmapped_#{region}_#{k}_#{kernel_adjust}_legend.png")
              ctrl_threads = g.get_band(0.0)
              expected_threads = g.get_band(0.4)
              ctrl_thread_values += g.threads.select {|x| ctrl_threads.include?(x.first) }.last
              expected_thread_values += g.threads.select {|x| expected_threads.include?(x.first) }.last
              filename = "unmapped_2_#{region}_#{k}_#{kernel_adjust}_bands.png"
              g.draw_bands(filename)
              #filename = "unmapped_#{region}_#{k}_#{kernel_adjust}_signal.png"
              #g.draw_signal(filename)
              #filename = "unmapped_#{region}_#{k}_#{kernel_adjust}_hits.png"
              #g.draw_hit_count(filename)

          rescue  Exception => e
                puts "failed on #{k} #{kernel_adjust}"
                puts e.message, e.backtrace
          end
        end
      end
    g.close      

    
end
File.open("ctrl_thread.txt", 'w') {|f| f.write(ctrl_thread_values.join("\n")) }
File.open("epxected_thread.txt", 'w') {|f| f.write(expected_thread_values.join("\n")) }