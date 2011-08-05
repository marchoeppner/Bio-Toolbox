
Pfam_DB_ADAPTER = 'mysql'
Pfam_DATABASE = "pfam"
Pfam_DB_HOST = 'localhost'
Pfam_DB_USERNAME = 'tools'
Pfam_DB_PASSWORD = 'analysis'

module Toolbox 
  
  module PfamDB
    
    require 'composite_primary_keys'
    
    class DBConnection < ActiveRecord::Base
      self.abstract_class = true
    
      def self.connect(version="24")

        establish_connection(
                              :adapter => Pfam_DB_ADAPTER,
                              :host => Pfam_DB_HOST,
                              :database => "#{Pfam_DATABASE}#{version}",
                              :username => Pfam_DB_USERNAME,
                              :password => Pfam_DB_PASSWORD
                              #:port => port
                            )
      end
    
    end
    
    
    class NcbiTaxonomy < DBConnection
    	set_table_name "ncbi_taxonomy"
    	set_primary_key "ncbi_taxid"
    	has_many :pfamseqs, :foreign_key => "ncbi_taxid"
    	
    	def adl_taxon
    	  
    	  # All exceptions....  	
    	  if self.taxonomy.include?("Bacteria")
    	    return "Bacteria"
    	  elsif self.taxonomy.include?("Archaea")
    	    return "Archaea"
    	  elsif self.taxonomy.downcase.include?("virus")
          return "Virus"
        elsif self.taxonomy.include?("unclassified") or self.taxonomy.length == 0
          return "unclassified"
        elsif self.taxonomy.include?("other sequences")
          return "other"
    	  end
    	  
    		levels = self.taxonomy.split(";").collect{|l| l.strip.gsub(/\./, '')}
    		lookup = PfamDB::TaxLookup.find_by_name(levels[1])
    		
    		return nil if lookup.nil?
    		return lookup.tax_supergroup.name   		
    	end
    	
    	
    	def self.get_taxa_at_level(level)   		
    		answer = []  		
    		self.find(:all).each do |taxon|    			
    			this_level = taxon.tax_string.split(";")
    			next if this_level.nitems-1 < level
    			answer << this_level[level-1] unless answer.include?(this_level[level-1]) or this_level[level-1].include?("candidate")   			
    		end   		
    		return answer.compact.uniq.sort    		
    	end
    	
    	def self.get_eukaryote_taxa_at_level(level)  		
    		answer = []
    		self.find(:all).each do |taxon|    		
    			levels = taxon.tax_string.split(";").collect{|a| a.gsub(/\./, '').strip}
    			next unless levels.include?("Eukaryota")
    			answer << levels[level] unless answer.include?(levels[level])    	
    		end
				return answer    					
    	end
    	
    end
    
    class Pfamseq < DBConnection
    	set_table_name "pfamseq"
    	set_primary_key "auto_pfamseq"
    	belongs_to :ncbi_taxonomy, :foreign_key => "ncbi_taxid"
    	has_many :pfam_reg_full_significants, :foreign_key => "auto_pfamseq"
    	
    	def ncbi_id
    		return self.taxonomy.ncbi_id
    	end
    	
    	def tax_string
    		self.taxonomy.tax_string
    	end
    	
    	def tax_entries
    		return self.taxonomy.tax_string.split(";").collect{|e| e.strip}
    	end
    end
    
    class PfamA < DBConnection
    	set_primary_key "auto_pfamA"
    	set_table_name 'pfamA'
      has_many :pfamA_reg_full_significants, :foreign_key => 'auto_pfamA'
      has_many :pfamA_to_gos, :foreign_key => "auto_pfamA"
      has_many :pfamA_supergroups, :foreign_key => "auto_pfamA"
      has_one :interpro, :foreign_key => "auto_pfamA"
      has_one :clan_membership, :foreign_key => "auto_pfamA"
      has_one :clan, :through => :clan_membership
      
      def go_terms
        return self.pfamA_to_gos.collect{|g| g.go_term }.join(",")
      end
    	
    end
    
    class Interpro < DBConnection
      set_primary_key 'auto_pfamA'
      belongs_to :pfamA, :foreign_key => "auto_pfamA"
    end
    
    class PfamARegFullSignificant < DBConnection
			set_table_name 'pfamA_reg_full_significant'
			set_primary_keys :auto_pfamA, :auto_pfamseq
      belongs_to :pfam, :foreign_key => 'auto_pfamA'
      belongs_to :pfamseq, :foreign_key => "auto_pfamseq"
      has_one :ncbi_taxonomy, :through => :pfamseq
      has_many :pfamA_tax_supergroups
      
      def adl_taxon
      	tax_string = self.pfamseq.ncbi_taxonomy.adl_taxon
      end
      
    end
    
    class TaxLookup < DBConnection
    	belongs_to :tax_supergroup, :foreign_key => "supergroup_id"
    end
    
    class TaxSupergroup < DBConnection
    	has_many :tax_lookups
    	has_many :pfamA_supergroups
    end 
    
    class PfamASupergroup < DBConnection
      set_table_name "pfamA_supergroup"
      belongs_to :pfamA, :foreign_key => "auto_pfamA"
      belongs_to :tax_supergroup, :foreign_key => "tax_supergroup_id"
    end
    
    class PfamAToGo < DBConnection
      set_table_name "pfamA_to_go"
      set_primary_key "id"
      belongs_to :pfamA, :foreign_key => "auto_pfamA"
    end
    
    class ClanMembership< DBConnection
      set_primary_keys :auto_pfamA,:auto_clan
      belongs_to :pfamA, :foreign_key => "auto_pfamA"
      belongs_to :clan, :foreign_key => "auto_clan"
    end
    
    class Clan < DBConnection
      set_primary_key 'auto_clan'
      has_many :clan_memberships, :foreign_key => "auto_clan"
      has_many :pfamAs, :through => :clan_memberships
    end
    
  end
  
end
