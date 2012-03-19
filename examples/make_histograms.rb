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




#open the BAM file and specify the region of interest
g = Bio::Util::Gngm.new(:file => "aln.bam", 
               :format => :bam, 
               :fasta => "reference.fasta", 
               :samtools => {:r => "Chr1:1-6000000",
                             :q => 20,
                             :Q => 50
               }
    )
#retrieve the SNPs from the BAM file
g.snp_positions

#plot a frequency histogram for different bin sizes
[100000, 250000, 500000].each do |bin_width|
    file_name = "#{bin_width}.png"
    g.frequency_histogram("#{file_name}",bin_width)
end

#close the BAM file
g.close

