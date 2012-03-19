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

length = 0
chr_name = "" 
file = Bio::FastaFormat.open("/Users/macleand/Desktop/laerfyve_vs_stitched_ler/ler_contigs_stitched.fa"
file.each do |entry|
  length = entry.length
  chr_name = entry.entry_id
end



  g = Bio::Util::Gngm.new(:file => "/Users/macleand/Desktop/laerfyve_vs_stitched_ler/aln.sort.bam", 
               :format => :bam,
               :samtools => {:q => 20, :Q => 50, :r => "#{chr_name}:1-#{length}"}, 
               :fasta => "/Users/macleand/Desktop/laerfyve_vs_stitched_ler/ler_contigs_stitched.fa"
  
    )
    g.snp_positions(:min_depth => 10, :mapping_quality => 40.0, :min_non_ref_count => 5)
    puts g.snp_positions.length
    
    [10000, 25000, 50000, 100000, 250000, 500000].each do |bin_width|
      file_name = "stitched_contigs_snps_q_filt_#{bin_width}.png"
      g.frequency_histogram("#{file_name}",bin_width)
    end
    g.close


