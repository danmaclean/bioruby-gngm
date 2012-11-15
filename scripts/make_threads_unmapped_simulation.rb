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
require 'pp'
length = 0
chr_name = "" 
#file = Bio::FastaFormat.open("/Users/macleand/Desktop/deletion_simulation/NC_000962.fna")
#file.each do |entry|
#  length = entry.length
#  chr_name = entry.entry_id
#end


  region = "1:2000000-3000000"
  
  puts "analyzing - #{region}"

  g = Bio::Util::Gngm.new(:file => "/Users/macleand/Desktop/insertion_finding/athal/aln.sort.bam", 
               :format => :bam, 
               :fasta => "/Users/macleand/Desktop/insertion_finding/athal/chr1.fa",
               :samtools => {:q => 50, :Q => 13, :r => region
               }
    )
      
      g.get_unmapped_mate_frequency(:ref_window_size => 152, :ref_window_slide => 10)
      
      
      
      g.collect_threads(:start => 0.0, :stop => 1.0, :slide => 0.1, :size => 0.1)
      pp g.threads
      g.threads.delete_if {|x| x.last.length <= 3 }

      begin
        #g.calculate_clusters(:pseudo => true)
        g.calculate_clusters(:k => 4, :adjust => 0.5, :control_chd => 0.0, :expected_chd => 0.3, :pseudo => false)
        filename = "deletion_real_data#{region}_all_threads.png" 
        g.draw_threads(filename, :draw_legend => "deletion_real_data#{region}_legend.png")
        ##no bands or signal to draw without clustering... 
        filename = "deletion_real_data#{region}_bands.png"
        g.draw_bands(filename)
        filename = "deletion_real_data#{region}_signal.png"
        g.draw_signal(filename)
        
        
        filename = "deletion_real_data_#{region}_hits.png"
        g.draw_hit_count(filename)
      rescue  Exception => e
        puts e.message, e.backtrace
      end



        
