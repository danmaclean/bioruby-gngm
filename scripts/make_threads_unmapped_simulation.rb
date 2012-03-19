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
file = Bio::FastaFormat.open("/Users/macleand/Desktop/deletion_simulation/NC_000962.fna")
file.each do |entry|
  length = entry.length
  chr_name = entry.entry_id
end


  region = "gi|57116681|ref|NC_000962.2|:1-#{length}"
  
  puts "analyzing - #{region}"

  g = Bio::Util::Gngm.new(:file => "/Users/macleand/Desktop/deletion_simulation/aln.sort.bam", 
               :format => :bam, 
               :fasta => "/Users/macleand/Desktop/deletion_simulation/NC_000962.fna",
               :samtools => {:q => 20, :Q => 50, :r => region
               }
    )
      
      g.get_unmapped_mate_frequency(:ref_window_size => 76, :ref_window_slide => 76)
      g.collect_threads(:start => 0.0, :stop => 0.5, :slide => 0.1, :size => 0.1)
      puts g.threads

      begin
        g.calculate_clusters(:pseudo => true)
        filename = "sim_2_#{region}_all_threads.png" 
        g.draw_threads(filename, :draw_legend => "sim_#{region}_legend.png")
        ##no bands or signal to draw without clustering... 
        filename = "sim_#{region}_hits.png"
        g.draw_hit_count(filename)
      rescue  Exception => e
        puts e.message, e.backtrace
      end



        
