GenomeComparaDB_ADAPTER = 'postgresql'
GenomeCompara_DATABASE = "genome_compara"
GenomeComparaDB_HOST = 'localhost'
GenomeComparaDB_USERNAME = 'tools'
GenomeComparaDB_PASSWORD = 'analysis'

module Toolbox 
  module GenomeComparaDB
    
    require 'composite_primary_keys'
    
    include ActiveRecord
    
    class DBConnection < ActiveRecord::Base
      self.abstract_class = true
    	self.pluralize_table_names = false
    	
      def self.connect(version="")

        establish_connection(
                              :adapter => GenomeComparaDB_ADAPTER,
                              :host => GenomeComparaDB_HOST,
                              :database => "#{GenomeCompara_DATABASE}_#{version}",
                              :username => GenomeComparaDB_USERNAME,
                              :password => GenomeComparaDB_PASSWORD
                              #:port => port
                            )
      end
    
    end
       
    class Organism < DBConnection
       has_many :members
    		
    	def connect_to_genome(version=61)
         if self.database.include?("EnsEMBL")
            Ensembl::Compara::DBConnection.connect(version,{:ensembl_genomes => true, :local => true})
            organism = Ensembl::Compara::GenomeDb.find_by_name(self.name)
            organism.connect_to_genome_new(version)
         else        
            GenomeDB::DBConnection.connect(self.name)
         end  
      
      end
    end
    
    class Member < DBConnection
        belongs_to :organism, :foreign_key => "organism_id"
        has_many :homology_members
        has_one :homology_group
    end
    
    class PeptideMember < DBConnection
        belongs_to :member, :foreign_key => "member_id"
        has_many :homology_members
    end
    
    class HomologyGroup < DBConnection
        belongs_to :member, :foreign_key => "member_id"
        has_many :homology_members
    end
    
    class HomologyMember < DBConnection
        set_primary_keys :member_id,:homology_group_id
        belongs_to :member, :foreign_key => "member_id"
        belongs_to :homology_group, :foreign_key => "homology_group_id"
        belongs_to :peptide_member, :foreign_key => "peptide_member_id"
    end
    
  end
  
end