
Neanderthal_DB_ADAPTER = 'postgresql'
Neanderthal_DATABASE = "neanderthal"
Neanderthal_DB_HOST = 'localhost'
Neanderthal_DB_USERNAME = 'tools'
Neanderthal_DB_PASSWORD = ''

module Toolbox 
  module Neanderthal
    
    include ActiveRecord
    
    class DBConnection < ActiveRecord::Base
      self.abstract_class = true
    
      def self.connect(version="")

        establish_connection(
                              :adapter => Neanderthal_DB_ADAPTER,
                              :host => Neanderthal_DB_HOST,
                              :database => "#{Neanderthal_DATABASE}#{version}",
                              :username => Neanderthal_DB_USERNAME,
                              :password => Neanderthal_DB_PASSWORD
                              #:port => port
                            )
      end
    
    end
    
    class Contig < DBConnection
      
      has_many :blasts
      has_many :organisms, :through => :blasts
      
      def to_fasta
        return Bio::FastaFormat.new(Bio::Sequence::NA.new(self.sequence).to_fasta(self.name))
      end
      
    end
    
    class Blast < DBConnection
    
      belongs_to :contig, :foreign_key => "contig_id"
      belongs_to :organism, :foreign_key => "organism_id"
      
      def ensembl_gene
        elements = self.target_def.split("|")
        return elements.shift
      end
      
      def ensembl_seq_region
        elements = self.target_def.split("|")
        return elements[1]
      end
      
      def ensembl_start
        elements = self.target_def.split("|")
        return elements[2].to_i + self.query_start
      end
      
      def ensembl_stop
        elements = self.target_def.split("|")
        return elements[2].to_i + self.query_end
      end
      
      def ensembl_strand
        elements = self.target_def.split("|")
        return elements[4].to_i
      end
      
      def length
        return self.target_seq.length
      end
      
    end
    
    class Organism < DBConnection
      
      has_many :blasts
      
    end
    
  end
  
end
