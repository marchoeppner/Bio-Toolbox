Leca_DB_ADAPTER = 'postgresql'
Leca_DATABASE = "leca"
Leca_DB_HOST = 'localhost'
Leca_DB_USERNAME = 'tools'
Leca_DB_PASSWORD = 'analysis'

module Toolbox 
  module LecaDB
    
    include ActiveRecord
    
    class DBConnection < ActiveRecord::Base
      self.abstract_class = true
    	self.pluralize_table_names = false
    	
      def self.connect(version="")

        establish_connection(
                              :adapter => Leca_DB_ADAPTER,
                              :host => Leca_DB_HOST,
                              :database => "#{Leca_DATABASE}#{version}",
                              :username => Leca_DB_USERNAME,
                              :password => Leca_DB_PASSWORD
                              #:port => port
                            )
      end
    
    end
       
    class Supergroup < DBConnection
      has_many :organisms
    end
    
    class Organism < DBConnection
      has_many :genes
      belongs_to :supergroup, :foreign_key => "supergroup_id"
      
      def db_args
      	args = {}
      	self.locator.split(";").each do |element|
      		e = element.split("=")
      		args[e[0].strip] = e[1].strip
      		if e[0] == "dbname"
      			args["version"] = e[1].slice(/5[0-9]/).to_i
      		end
      	end
      	return args
     	end
     	
      def connect_to_genome(version=56)
        if self.db_args["host"] == "mysql.ebi.ac.uk"
        	Ensembl::Core::DBConnection.connect("#{self.snake_case}",self.db_args["version"],{:ensembl_genomes => true, :database => self.db_args["dbname"]})
        elsif self.locator.include?("eukaryotes")
          EukaryoteDB::DBConnection.connect(1)
        else
        	Ensembl::Core::DBConnection.connect("#{self.snake_case}",self.db_args["version"])
        end
      end
      
      def snake_case
      	return self.name.gsub(/\s/, '_').downcase
      end
      
    end
    
    class Gene < DBConnection
      belongs_to :organism, :foreign_key => "organism_id"
      has_many :transcripts
      has_many :xref_member_pfams
      has_many :pfams, :through => :xref_member_pfams
    end
    
    class Transcript < DBConnection
      belongs_to :gene, :foreign_key => "gene_id"
      has_one :translation
    end
    
    class Translation < DBConnection
      belongs_to :transcript, :foreign_key => "transcript_id"
      has_many :xref_translation_pfams
      has_many :pfams, :through => :xref_translation_pfams
    end
    
    class Pfam < DBConnection
      has_many :xref_gene_pfams
      has_many :xref_translation_pfams
    end
    
    class XrefTranslationPfam < DBConnection
      belongs_to :pfam, :foreign_key => "pfam_id"
      belongs_to :translation, :foreign_key => "translation_id"
    end
    
    class XrefGenePfam < DBConnection
      belongs_to :pfam, :foreign_key => "pfam_id"
      belongs_to :gene, :foreign_key => "gene_id"
    end
    
    
  end
  
end
