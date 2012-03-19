#!/usr/bin/env ruby
#
#  untitled
#
#  Created by Dan MacLean (TSL) on 2012-03-01.
#  Copyright (c)  . All rights reserved.
###################################################

require 'bio'

file = Bio::FastaFormat.open(ARGV[0])
file.each do |entry|
  section = entry.seq[ARGV[1].to_i..ARGV[2].to_i]
  puts ">#{entry.entry_id}:#{ARGV[1]}..#{ARGV[2]}"
  puts section
end