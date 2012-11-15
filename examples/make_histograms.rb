#!/usr/bin/env ruby
#
#  make_histograms
#
#  Created by Dan MacLean (TSL) on 2012-01-17.
#  Copyright (c)  . All rights reserved.
###################################################

### An example script to get SNP positions and make histograms
### of the frequncy of discordant SNPs. Generates plots for each. 

$LOAD_PATH.unshift(File.join(File.dirname(__FILE__), '..', 'lib'))
$LOAD_PATH.unshift(File.dirname(__FILE__))
require 'bio-gngm'
require 'bio-samtools'

#def make_snp_array(id, file)
#  a = []
#  File.open(file, "r").each do |line|
#    arr = line.split(/\t/)
#    next if arr[0] !~ /#{id}/
#    a << [a]
#  end
#end

sequences = Bio::DB::FastaLengthDB.new(:file => ARGV[0])
$stderr.puts "Loaded sequences..."


sequences.each do |id, length|
  g = Bio::Util::Gngm.new(:file => ARGV[1], 
               #:format => :pileup, 
               :format => :bam,
               :fasta => ARGV[0],
               :chromosome => id,
               :start => 1,
               :stop => length,
               :samtools => {
                             :q => 20,
                             :Q => 20
                             },
              :variant_call => {
                 :indels => false, 
                 :deletions_only => false, 
                 :insertions_only => false, 
                 :min_depth => 6, 
                 :max_depth => 250, 
                 :mapping_quality => 20.0, 
                 :min_non_ref_count => 2, 
                 :ignore_reference_n => true,
                 :min_snp_quality => 20,
                 :min_consensus_quality => 20
                 }               
  )
  $stderr.puts "getting #{id}.."
  
  g.snp_positions
  puts "Found #{g.snp_positions.length} SNPs"
  [100000, 250000, 500000].each do |bin_width|
    $stderr.puts "working on #{bin_width} windows"
    file_name = "test_#{id}_#{bin_width}.png"
    g.frequency_histogram("#{file_name}",bin_width)
  end
g.close
end

