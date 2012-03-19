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
require 'bio-samtools'
require 'rinruby'
require 'bio'

length = 0
chr_name = "" 
file = Bio::FastaFormat.open("/Users/macleand/Desktop/laerfyve_vs_stitched_ler/ler_contigs_stitched.fa")
file.each do |entry|
  length = entry.length
  chr_name = entry.entry_id
end

[250000,7000000,100000].each do |start|
  b = Bio::DB::Sam.new(:bam => "/Users/macleand/Desktop/laerfyve_vs_stitched_ler/aln.sort.bam",
              :fasta => "/Users/macleand/Desktop/laerfyve_vs_stitched_ler/ler_contigs_stitched.fa",
              :q => 20, :Q => 50)
              b.open
              arr = []              
              b.fetch(chr_name, start, start + 100000).each do |aln|
                puts aln.mate_unmapped
              end
       #puts arr
       #r = RinRuby.new(echo = false, interactive = false)
       #r.eval "options(warn=-1)"
       #r.a = arr
       #r.eval "med = median(a)"
       #a = r.med
       #puts a
       #r.eval "png(\"#{start}.png\")"
       #r.eval "hist(a)"
       #r.eval "dev.off()"
       #r.quit       
end


=begin
sample_area_length = 100000
[250000,7000000,100000].each do |start|
  [184, 336, 500, 1000].each do |window|
    [50, 184, 336, 250].each do |isize|
      begin
        stop = start + sample_area_length
         g = Bio::Util::Gngm.new(:file => "/Users/macleand/Desktop/laerfyve_vs_stitched_ler/aln.sort.bam", 
                      :format => :bam, 
                      :fasta => "/Users/macleand/Desktop/laerfyve_vs_stitched_ler/ler_contigs_stitched.fa",
                      :samtools => {:q => 20, :Q => 50, :r => "#{chr_name}:#{start}-#{stop}"
                      }
         )
        g.get_insert_size_frequency(:ref_window_size => window, :ref_window_slide => window, :isize => isize)
        $stderr.puts "collecting"
        g.collect_threads
        #pp g.threads
        #puts g.snp_positions.length
        #pp g.snp_positions
        filename = "#{start}-#{stop}_#{window}_#{isize}.png"
        g.draw_hit_count(filename)
        rescue  Exception => e
            puts "failed on #{start} #{window} #{isize}"
            puts e.message, e.backtrace
        ensure
            g.close
        end
        
    end
  end  
end
=end
