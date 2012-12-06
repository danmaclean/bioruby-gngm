require 'helper'

class TestBioGngm < Test::Unit::TestCase
  
  context "Classes in bio/db/sam" do
    should "contain Pileup" do
      assert_instance_of Bio::DB::Pileup, Bio::DB::Pileup.new("seq1\t272\tT\t24\t,.$.....,,.,.,...,,,.,..^+.\t<<<+;<<<<<<<<<<<=<;<;7<&")
    end
  
    should "contain Vcf " do
      assert_instance_of Bio::DB::Vcf, Bio::DB::Vcf.new("20	14370	rs6054257	G	A	29	0	NS=3;DP=14;AF=0.5;DB;H2	GT:GQ:DP:HQ	0|0:48:1:51,51	1|0:48:8:51,51	1/1:43:5:-1,-1")
    end
  end
  
  context "A Pileup instance" do
    setup do
      @p1 = Bio::DB::Pileup.new("seq1\t272\tT\t10\t,.$.....,,.^+.\t<<<+;<<<<<<<<<<<=<;<;7<&")
      @p2 = Bio::DB::Pileup.new("seq1\t272\tT\t24\taaaaaaaaaaaa............\t<<<+;<<<<<<<<<<<=<;<;7<&")
      @p3 = Bio::DB::Pileup.new("seq1\t272\tT\t24\taaaaaaaaaaaaaaaaaaaaaaaa\t<<<+;<<<<<<<<<<<=<;<;7<&")
      @p4 = Bio::DB::Pileup.new("seq1\t272\tT\t9\taaaaaaaGGC\t<<<+;<<<<<<<<<<<=<;<;7<&")
    end
    
    should "calculate discordant chastity" do
      assert_equal @p1.discordant_chastity, 0.0
      assert_equal @p2.discordant_chastity, 0.5
      assert_equal @p3.discordant_chastity, 1.0
      assert_equal @p4.discordant_chastity, 0.7777777777777778
    end
    
    should "know it is a snp" do
      assert @p2.is_snp?(:min_depth => 5, :min_non_ref => 5)
      assert !@p1.is_snp?(:min_depth => 5, :min_non_ref => 5)
    end
  end
  
  context "A Gngm object" do
    setup do
      @g1 = Bio::Util::Gngm.new(:file => "test/data/testu.bam", 
                     :format => :bam, 
                     :fasta => "test/data/test_chr.fasta", 
                     :samtools => {:r => "chr_1:100-300"}
                     )
      @g1.snp_positions
      
      @g2 = Bio::Util::Gngm.new(:file => "/Users/macleand/Desktop/ngm/fli.sorted.bam", 
                     :format => :bam, 
                     :fasta => "/Users/macleand/Desktop/ngm/TAIR9_chr_all.fas", 
                     :samtools => {:r => "1:100-10000"},
                     :histo_bin_width => 100

      )
    end
    
    
    should "open BAM files" do
      assert_instance_of Bio::DB::Sam, @g1.file
    end
    
    should "return array of snps annotated with discordant chastity" do
      assert_instance_of Array, @g1.snp_positions
      assert_equal [[185, 0.3333333333333333], [203, 0.3333333333333333], [279, 0.5]], @g1.snp_positions
    end
    
    should "draw a picture of a SNP histogram" do
      @g1.frequency_histogram("test_histo.png")
      assert File.exists?("test_histo.png")
      File.delete("test_histo.png")
    end
    
    should "get chastity threads" do
      @g1.snp_positions
      assert_equal [[0.3, [185, 203]],[0.4, []],[0.5, [279]], [0.6000000000000001, []], [0.7, []],[0.8, []],[0.9000000000000001, []],[1.0, []]], @g1.threads(:start => 0.3, :stop => 1, :slide => 0.1, :size => 0.1)
    end
    
    should "draw chastity threads" do
      @g2.snp_positions
      @g2.calculate_clusters(:k =>4, :set_seed => 1)
      @g2.draw_threads("test_threads.png")
      assert File.exists?("test_threads.png")
      #File.delete("test_threads.png")
    end


    should "return lists of densities for each window from R" do
      @g1.snp_positions
      #puts @g1.densities[0]
      assert_instance_of Array, @g1.densities
      assert_instance_of Array, @g1.densities.first
      assert_equal 0.24, @g1.densities[0].first
    end
    
    should "find index of window in the @densities array" do
      @g2.snp_positions
      assert_equal 65, @g2.find_index(1.0)
    end
    
    should "find bands ie cluster and return list of threads in that cluster together" do
      @g2.snp_positions
      @g2.calculate_clusters(:k => 4, :set_seed => 1)  ##due to randomness of clustering, its currently hard to make sure this is exactly tested.. 
                                                        #setting the seed ensures consistent clustering ...
      assert_equal [0.21000000000000002, 0.22, 0.8200000000000001, 0.8300000000000001, 0.8400000000000001, 0.8500000000000001, 0.8600000000000001, 0.97, 0.98, 0.99, 1.0], @g2.get_band(0.22)  
    end
    
    should "draw found bands" do
      @g2.snp_positions
      @g2.calculate_clusters(:k => 4, :set_seed => 1)
      @g2.draw_bands("test_bands.png")
      assert File.exists?("test_bands.png")
      #File.delete("test_bands.png")
    end
    
    should "draw_signal_plots" do
      @g2.snp_positions
      @g2.calculate_clusters(:k => 4, :set_seed => 1)
      @g2.draw_signal("signal.png")
      File.delete("signal.png")
    end
    
    should "calculate signal" do
      @g2.snp_positions
      @g2.calculate_clusters(:k => 4, :set_seed => 1)
      assert_equal 1.157971267690456, @g2.signal.last
    end
  end
  
end
