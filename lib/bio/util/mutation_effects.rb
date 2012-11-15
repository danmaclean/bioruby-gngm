require 'bio'


module Bio
  class Util
    def self.read_gff3(fn)
      genearray=Array.new
      mRNAhash=Hash.new
      exonhash=Hash.new
      tehash=Hash.new
      lastid = ''
      lastrecord = nil
      gff3 = Bio::GFF::GFF3.new(File.read(fn))
      gff3.records.each do | record |
       feature_type = record.feature_type
       if(feature_type == 'gene')
         genearray << record.id
       elsif(feature_type == 'mRNA')
         parent = record.get_attribute('Parent')
        if mRNAhash[parent] == nil
          mRNAhash[parent] = [record]
        else 
          mRNAhash[parent] << record
        end
       elsif(feature_type == 'transposable_element')
        #--- not yet implemented
       elsif(feature_type == 'exon')
        parents = record.get_attributes('Parent')
        parents.each do |parent|
         if exonhash[parent] == nil
          exonhash[parent] = Array.new
         end
         exonhash[parent] << record
         end
       end
      end
    end
  end
end