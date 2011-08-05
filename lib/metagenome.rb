
Metagenome_DB_ADAPTER = 'postgresql'
Metagenome_DATABASE = "metagenomics"
Metagenome_DB_HOST = 'localhost'
Metagenome_DB_USERNAME = 'tools'
Metagenome_DB_PASSWORD = 'analysis'

module Toolbox 
  module MetagenomeDB
    
    include ActiveRecord
    
    class DBConnection < ActiveRecord::Base
      self.abstract_class = true
    
      def self.connect(version="")

        establish_connection(
                              :adapter => Metagenome_DB_ADAPTER,
                              :host => Metagenome_DB_HOST,
                              :database => "#{Metagenome_DATABASE}#{version}",
                              :username => Metagenome_DB_USERNAME,
                              :password => Metagenome_DB_PASSWORD
                              #:port => port
                            )
      end
    
    end
    
    class Project < DBConnection
    	has_many :dnafrags
    end
    
    class Dnafrag < DBConnection
    	belongs_to :project, :foreign_key => "project_id"
    	has_many :blast_hits
    	has_many :infernal_hits
    	
    	def to_fasta
    		return Bio::FastaFormat.new(Bio::Sequence::NA.new(self.sequence).to_fasta(self.accession))
    	end
    	
    end
    
    class BlastHit < DBConnection
    	belongs_to :dnafrag, :foreign_key => "dnafrag_id"
		end
		
		class InfernalHit < DBConnection
		  belongs_to :dnafrag, :foreign_key => "dnafrag_id"
		end    
   	
  end
  
end