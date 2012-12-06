#!/usr/bin/env ruby

# Author : Naveed Ishaque (inspired by Dan Maclean) edited again by Dan to include ChD value setting and deletion of SNPs file
# naveed.ishaque@tsl.ac.uk; naveed.ishaque@hotmail.co.uk
# Date: 20th June 2012 and 1st November 2012

# This scripts produces a HTML with embedded images showing the SNP density and chastity plots for a given BAM file
# It will automatically iterate over all contigs, form begining to end
# NOTE - run using ruby executable: /home/programs/gngm/ruby/bin/ruby
# NOTE - this will only run on the cluster via a bsub command


require 'rubygems'
require 'bio'
require 'bio-gngm'
require 'base64'
require 'getoptions'

usage = "\n#{$PROGRAM_NAME} reads in a fasta and bam files and produces a html file indicating QTL locations as peaks\n\n\n\t #{$PROGRAM_NAME}\n\n\n -f [reference fasta file]\n -b [bam file]\n -e expected ChD (allele freq) default 1\n -c control ChD (allele freq) default 0.5\n -s List of known SNPS [tab delimited file]\n\n"

# PARSE INPUTS

opt = GetOptions.new(%w(h help f=@s b=@s e=@s c=@s s=@s))

puts "#{usage}" if opt[:h]
exit if opt[:h]
puts "#{usage}" if opt[:help]
exit if opt[:help]

puts "ERROR - no fasta file provided (-f)\n#{usage}" unless opt[:f]
exit unless opt[:f]
puts "ERROR - fasta file '#{opt[:f][0]}' does not exist\n#{$usage}" unless FileTest.exist?("#{opt[:f][0]}")
exit unless FileTest.exist?("#{opt[:f][0]}")
warn "\nUsing FASTA file #{opt[:f][0]}"

puts "ERROR - no bam file provided (-b)\n#{usage}" unless opt[:b]
exit unless opt[:b]
puts "ERROR - BAM file '#{opt[:b][0]}' does not exist\n#{usage}" unless FileTest.exist?("#{opt[:b][0]}")
exit unless FileTest.exist?("#{opt[:b][0]}")
warn "Using BAM file #{opt[:b][0]}"


expected_chd = 1.0
expected_chd =opt[:e][0].to_f if expected_chd and opt[:e]

control_chd = 0.5
control_chd = opt[:c][0].to_f if control_chd and opt[:c]

known_snps = Hash.new { |h,k| h[k] = Array.new }
if opt[:s].first
  File.open(opt[:s].first).each  do |line|
    chr,pos = line.split("\t")
    known_snps[chr] << pos.to_i
  end
end


$stderr.puts "using expected ChD: #{expected_chd} and control ChD: #{control_chd}"

hist_bins = [100000, 250000, 500000]
ks = [5, 7, 9, 11]
kadjusts = [0.5, 0.25, 0.1, 0.05, 0.01]
warn "Histogram bin sizes: #{hist_bins}\nThread clusters (K): #{ks}\nKernal adjusts: #{kadjusts}\n"

# LOAD FASTA and find contigs
sequences = Bio::DB::FastaLengthDB.new(:file => "#{opt[:f][0]}")


# For Each contig in the fasta file analyse...

sequences.each do |id,length|
  warn "\nProcessing #{id}:1 - #{length}..."

  warn "Skipping #{id} as too short ..." if length < (4 * hist_bins.max)
  next if length < (4 * hist_bins.max)

  g = Bio::Util::Gngm.new(:file => "#{opt[:b][0]}", 
                 :format => :bam, 
                 :fasta => "#{opt[:f][0]}",
                 :start => 1,
                 :stop => 10000,
                 :chromosome => id, 
                 :samtools => {
                               :q => 20,
                               :Q => 20
                              },
                 :ignore_file => "#{opt[:s][0]}",
                 :write_pileup => "pileup.txt",
                 :write_vcf => "snps.vcf"
    )

  # predict SNPs

  warn "  Prediciting SNPs for #{id}:1-#{length}..."
  g.snp_positions

          #delete SNPs from known snp_list
          #a = g.snp_positions.dup
          #known_snps[seq.entry_id].each {|snp_pos| a.delete_if{|x| x.first == snp_pos} }
          #$stderr.puts "deleted #{g.snp_positions.length - a.length} snps appearing in #{opt[:s]}"
          #g.snp_positions = a



  # produce SNP density histograms

  warn "  Iterating over different histogram bin sizes..."
  hist_bins.each do |bin_width|
    warn "  Makings PNG for bin size #{bin_width}..."
    file_name = "#{id}_SNP_histogram_bin#{bin_width}.png"
    g.frequency_histogram("#{file_name}",bin_width, :title => "#{id}: SNP density histogram (bin width - #{bin_width})", :width => 1066, :height => 300)
  end

  # Write to embedded HTML  
  
  htmlout = File.open("#{id}.html", 'w') 
  htmlout.puts "<html>\n"
  htmlout.puts "  <head>\n"
  htmlout.puts "    <title>GNGM #{id} - QTL mapping</title>\n"
  htmlout.puts "    <style type=\"text/css\">\n"
  htmlout.puts "      table,\n"
  htmlout.puts "      td,\n"
  htmlout.puts "      tbody,\n"
  htmlout.puts "      thead,\n"
  htmlout.puts "      thead th,\n"
  htmlout.puts "      tr.even,\n"
  htmlout.puts "      tr.odd {\n"
  htmlout.puts "        border: 0;\n"
  htmlout.puts "      }\n"
  htmlout.puts "    </style>\n"
  htmlout.puts "  </head>\n"
  htmlout.puts "  <body>\n\t\t"
  htmlout.puts "    <table>\n"
  htmlout.puts "      <tr>\n"
  hist_bins.each do |bin_width|
    htmlout.puts "        <td>\n"
    htmlout.puts "<img src=\"data:image/gif;base64,"  
    htmlout.puts [open("#{id}_SNP_histogram_bin#{bin_width}.png").read].pack("m")
    File.delete("#{id}_SNP_histogram_bin#{bin_width}.png")
    htmlout.puts "\" width=\"533\" height=\"150\"/>\n"
    htmlout.puts "        </td>\n"
  end
  htmlout.puts "      </tr>\n"
  htmlout.puts "    </table>\n"

  # Perform chastity calculations 

  warn "  Collecting threads..."
  g.collect_threads
  warn "  Iterating over k and kernel adjusts..."
  ks.each do | k |
    begin
      warn "  Makings PNG for k = #{k} ..."
      warn "    Calculating threads ..."
      g.calculate_clusters(:k => k, :adjust => 0.5, :control_chd => control_chd, :expected_chd => expected_chd)
      warn "    Drawing threads ..."
      filename = "#{id}_k#{k}_threads.png"  
      g.draw_threads(filename, :title => "#{id}: Chastity bands - all phases (k=#{k})", :width => 700, :height => 300)
      warn "    Clustering bands ..."
      filename = "#{id}_k#{k}_clustered_bands.png"
      g.draw_bands(filename, :title => "#{id}: Homozygous and heterozygous chastity belts (k=#{k})", :width => 800, :height => 300)
      kadjusts.each do |kernel_adjust|
        begin
          warn "    Calculating threads (with kernal adjust #{kernel_adjust}) ..."
          g.calculate_clusters(:k => k, :adjust => kernel_adjust, :control_chd => control_chd, :expected_chd => expected_chd)
          warn "    Calculating signal ..."
          filename = "#{id}_k#{k}_kadjust#{kernel_adjust}_signal.png"
          g.draw_signal(filename, :title => "#{id}: Homo/Het signal ratio (k=#{k}, kernal=#{kernel_adjust})", :width => 800, :height => 300)
          warn "    Estimating peaks ..."
          filename = "#{id}_k#{k}_kadjust#{kernel_adjust}_peaks.png"
          g.draw_peaks(filename, :title => "#{id}: Signal peaks (k=#{k}, kernal=#{kernel_adjust})", :width => 800, :height => 300)
        rescue => e
          $stderr.puts "skipping #{k} #{kernel_adjust} => #{e}"
        end
      end
    rescue => e
      $stderr.puts "Skipping #{k} => #{e}"
    end
  end

  g.close

  # Write to embedded HTML  

  htmlout.puts "    <table>\n"

  # all bands
  htmlout.puts "      <tr>\n"
  ks.each do | k | 
    htmlout.puts "        <td>\n"
    htmlout.puts "<img src=\"data:image/gif;base64,"  
    htmlout.puts [open("#{id}_k#{k}_threads.png").read].pack("m")
    File.delete("#{id}_k#{k}_threads.png")
    htmlout.puts "\" width=\"400\" height=\"150\"/>"
    htmlout.puts "        </td>\n"
  end
  htmlout.puts "      </tr>\n"

  # homo/het bands
  htmlout.puts "      <tr>\n"
  ks.each do | k | 
    htmlout.puts "        <td>\n"
    htmlout.puts "<img src=\"data:image/gif;base64,"  
    htmlout.puts [open("#{id}_k#{k}_clustered_bands.png").read].pack("m")
    File.delete("#{id}_k#{k}_clustered_bands.png")
    htmlout.puts "\" width=\"400\" height=\"150\"/>"
    htmlout.puts "        </td>\n"
  end
  htmlout.puts "      </tr>\n"

  # k/adjusts
  kadjusts.each do |kernel_adjust| 
    htmlout.puts "      <tr>\n"
    ks.each do | k |
      htmlout.puts "        <td>\n"
      htmlout.puts "<img src=\"data:image/gif;base64,"  
      htmlout.puts [open("#{id}_k#{k}_kadjust#{kernel_adjust}_signal.png").read].pack("m")
      File.delete("#{id}_k#{k}_kadjust#{kernel_adjust}_signal.png")
      htmlout.puts "\" width=\"400\" height=\"150\"/>"
      htmlout.puts "        </td>\n"
    end
    htmlout.puts "      </tr>\n"
    htmlout.puts "      <tr>\n"
    ks.each do | k |
      htmlout.puts "        <td>\n"
      htmlout.puts "<img src=\"data:image/gif;base64,"  
      htmlout.puts [open("#{id}_k#{k}_kadjust#{kernel_adjust}_peaks.png").read].pack("m")
      File.delete("#{id}_k#{k}_kadjust#{kernel_adjust}_peaks.png")
      htmlout.puts "\" width=\"400\" height=\"150\"/>"
      htmlout.puts "        </td>\n"
    end
    htmlout.puts "      </tr>\n"
  end

  htmlout.puts "    </table>\n"
  htmlout.puts "\n  </body>\n</html>\n"

  htmlout.close

end



