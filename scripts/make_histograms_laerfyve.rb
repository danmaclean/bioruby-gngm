#!/usr/bin/env ruby
#
#  make_histograms
#
#  Created by Dan MacLean (TSL) on 2012-01-17.
#  Copyright (c)  . All rights reserved.
###################################################

### An example script to loop over each reference in the BAM file, get SNP positions and make histograms
### of the frequncy of discordant SNPs. Generates plots for each. 

$LOAD_PATH.unshift(File.join(File.dirname(__FILE__), '..', 'lib'))
$LOAD_PATH.unshift(File.dirname(__FILE__))
require 'bio-gngm'
require 'bio-samtools'
require 'bio'
=begin
module Bio
  class DB
    class Sam
      def each_reference
        index_stats.each_pair do |k, v|
          yield k, v[:length].to_i
        end
      end
    end
  end
end
=end

#ff = Bio::FlatFile.new(Bio::FastaFormat, "/Users/macleand/Desktop/laerfyfe_vs_ler/Ler-1.SHORE.scaffolds.2010-09-30.500bp.fa")
file = Bio::FastaFormat.open("/Users/macleand/Desktop/laerfyfe_vs_ler/Ler-1.SHORE.scaffolds.2010-09-30.500bp.fa")
file.each do |entry|
    next if entry.length < 10000
    $stderr.puts "doing #{entry.entry_id} - #{entry.length}"                      
    g = Bio::Util::Gngm.new(:file => "/Users/macleand/Desktop/laerfyfe_vs_ler/aln.sort.bam", 
                 :format => :bam,
                 :samtools => {:q => 20, :Q => 50, :r => "#{entry.entry_id}:1-#{entry.length}"}, 
                 :fasta => "/Users/macleand/Desktop/laerfyfe_vs_ler/Ler-1.SHORE.scaffolds.2010-09-30.500bp.fa"

      )
       g.snp_positions(:min_depth => 10, :mapping_quality => 40.0, :min_non_ref_count => 5)
          puts "found #{g.snp_positions.length} SNPs .."
          begin
          [1000, 2500, 5000, 10000, 25000, 50000].each do |bin_width|
            next if bin_width > entry.length
            file_name = "ler_contigs_#{entry.entry_id}_#{bin_width}.png"
            g.frequency_histogram("#{file_name}",bin_width)
          end
          rescue Exception => e
            puts "failed #{e}"
          ensure 
            g.close
          end
end
#sam = Bio::DB::Sam.new(:bam => "/Users/macleand/Desktop/laerfyfe_vs_ler/aln.sort.bam", 
#                      :fasta => "/Users/macleand/Desktop/laerfyfe_vs_ler/Ler-1.SHORE.scaffolds.2010-09-30.500bp.fa")
#sam.open                      
#sam.fetch("Scaffold_2", 10000, 10500).each do |a|
#  puts a.qname
#end

#sam.close
#sam.each_reference do |name, length|
#  $stderr.puts "skipping..."
#  next if length < 10000
#  $stderr.puts "doing #{name}"                      
#  g = Bio::Util::Gngm.new(:file => "/Users/macleand/Desktop/laerfyfe_vs_ler/aln.sort.bam", 
#               :format => :bam,
#               :samtools => {:q => 20, :Q => 50}, 
#               :fasta => "/Users/macleand/Desktop/laerfyfe/Ler-1.SHORE.scaffolds2010-09-30.bp.fa"
  
#    )
#    g.snp_positions(:min_depth => 10, :mapping_quality => 40.0, :min_non_ref_count => 5)
    #puts g.snp_positions.length
    
#    [10000, 25000, 50000, 100000, 250000, 500000].each do |bin_width|
#      file_name = "ler_contigs_#{name}_#{bin_width}.png"
#      g.frequency_histogram("#{file_name}",bin_width)
#    end
#    g.close
#end
#sam.close