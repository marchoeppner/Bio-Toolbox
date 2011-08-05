
Rfam_DB_ADAPTER = 'mysql'
Rfam_DATABASE = "rfam"
Rfam_DB_HOST = 'localhost'
Rfam_DB_USERNAME = 'tools'
Rfam_DB_PASSWORD = 'analysis'

module Toolbox 
  
  module RfamDB
    
    class DBConnection < ActiveRecord::Base
      self.abstract_class = true
    
      def self.connect(version="_91")

        establish_connection(
                              :adapter => Rfam_DB_ADAPTER,
                              :host => Rfam_DB_HOST,
                              :database => "#{Rfam_DATABASE}#{version}",
                              :username => Rfam_DB_USERNAME,
                              :password => Rfam_DB_PASSWORD
                              #:port => port
                            )
      end
    
    end
    
    
    class GenomeEntry < DBConnection
    	
    	has_many :chromosome_builds, :foreign_key => "auto_genome"
    	has_many :rfam_reg_fulls, :foreign_key => "auto_genome"
    	has_many :rfam_reg_seeds, :foreign_key => "auto_genome"
    	belongs_to :taxonomy, :foreign_key => "ncbi_id"
    	
    end
    
    class ChromosomeBuild < DBConnection
    	
    	belongs_to :rfamseq, :foreign_key => "auto_rfamseq"
    	belongs_to :genome_entry, :foreign_key => "auto_genome"
    	set_primary_key "auto_genome"
    	
    end
    
    class Taxonomy < DBConnection
    	set_table_name "taxonomy"
    	set_primary_key "ncbi_id"
    	has_many :rfamseqs, :foreign_key => "ncbi_id"
    	has_many :genome_entries, :foreign_key => "ncbi_id"
    	
    	def adl_taxon  	
    		tax_string = self.tax_string
				return "unclassified" unless self.tax_string.include?(";")
    		levels = self.tax_string.split(";").collect{|l| l.strip.gsub(/\./, '')}
    		lookup = RfamDB::TaxLookup.find_by_name(levels[1].gsub(/\./, ''))
    		if lookup.nil?
    			return nil
    		else
    			return lookup.tax_supergroup.name   		
    		end
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
    
    class Rfamseq < DBConnection
    	set_table_name "rfamseq"
    	set_primary_key "auto_rfamseq"
    	belongs_to :taxonomy, :foreign_key => "ncbi_id"
    	has_many :rfam_reg_seeds, :foreign_key => "auto_rfamseq"
    	has_many :rfam_reg_fulls, :foreign_key => "auto_rfamseq"
    	has_one :chromosome_build, :foreign_key => "auto_rfamseq"
    	
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
    
    class Rfam < DBConnection
    	set_primary_key "auto_rfam"
    	set_table_name 'rfam'
    	belongs_to :wikitext, :foreign_key => "auto_wiki"
      has_many :rfam_reg_fulls, :foreign_key => 'auto_rfam'
      has_many :rfam_reg_seeds, :foreign_key => 'auto_rfam'
      has_many :rfam_counts
      has_one :clan_membership, :foreign_key => "auto_rfam"
      has_one :clan, :through => :clan_membership
      has_many :rfam_database_links, :foreign_key => "auto_rfam"
      has_many :rfam_taxonomies,:foreign_key => "auto_rfam"
      has_many :rfam_xref_seqs, :foreign_key => "auto_rfam"
      
      
      def clan_name
        if self.has_clan?
          return self.clan.clan_acc
        else
          return ""
        end
      end
      
      def biotype
      	if self.typ.include?("snoRNA")
      		return "snoRNA"
      	else
      		return self.typ
      	end
      end
      
      def count_species
      	return self.rfam_counts.select{|c| c.typ == "species"}
      end
      
      def count_total
      	return self.rfam_counts.select{|c| c.typ == "total"}
      end
      
      def present_in_leca?
        answer = false
        present = 0
        self.rfam_counts.each do |count|
          present += 1 if count.species_count > 1
        end
        present > 1 ? answer = true : answer = false
        return answer
      end
      
      def has_clan?
        self.clan_membership ? true : false
      end
      
      def sum
      	counter = 0
	      self.rfam_counts.each {|c| counter += c.count}
      	return counter
      end
      
      def self.find_all_snornas    
      	answer = []
      	self.find(:all).select{|r| r.typ.include?("snoRNA")}.each{|s| answer << s}
      	self.find(:all).select{|r| r.typ.include?("scaRNA")}.each{|s| answer << s}
      	return answer
      end
    	
    end
    
    class RfamRegFull < DBConnection
			set_table_name 'rfam_reg_full'
			set_primary_keys :auto_rfam, :auto_rfamseq
      belongs_to :rfam, :foreign_key => 'auto_rfam', :primary_key => 'auto_rfam'
      belongs_to :rfamseq, :foreign_key => "auto_rfamseq"
      has_one :taxonomy, :through => :rfamseq
      belongs_to :genome_entry, :foreign_key => "auto_genome"
      
      def species
        return self.rfamseq.taxonomy.species
      end
      
      def tax_string
        return self.rfamseq.taxonomy.tax_string
      end
      
      def adl_taxon
      	answer = nil
      	tax_string = self.rfamseq.taxonomy.tax_string.clone.split(";")
      	if tax_string.include?("Bacteria")
      	  answer = "Bacteria"
      	elsif tax_string.include?("Archaea")
      	  answer = "Archaea"
      	elsif tax_string.include?("Eukaryota")
      	  lookup = TaxLookup.find_by_name(tax_string[1].strip.gsub(/\./, ''))
      	  raise "Encountered unknown eukaryotic taxon #{tax_string}" if lookup.nil?
      	  answer = lookup.tax_supergroup.name
      	 elsif tax_string.include?("Viruses")
      	 	tax_supergroup = TaxSupergroup.find_by_name("Virus")
      	 	answer = "Virus"
      	end
      	tax_string = nil
      	return answer
      end
      
    end
    
    class RfamTaxonomy < DBConnection
    	belongs_to :rfam, :foreign_key => "auto_rfam"
    end
    
    class RfamRegSeed < DBConnection
			set_table_name 'rfam_reg_seed'
			set_primary_keys :auto_rfam, :auto_rfamseq
      belongs_to :rfam, :foreign_key => "auto_rfam", :primary_key => "auto_rfam"
      belongs_to :rfamseq, :foreign_key => "auto_rfamseq"
      belongs_to :genome_entry, :foreign_key => "auto_genome"
      
      def adl_taxon
      	tax_string = self.rfamseq.taxonomy.tax_string.split(";")
      	return nil unless tax_string.include?("Eukaryota")
      	lookup = TaxLookup.find_by_name(tax_string[1].strip)
      	return lookup.tax_supergroup.name
      end
    end
    
    class Wikitext < DBConnection
    	set_table_name 'wikitext'
    	set_primary_key "auto_wiki"
    end	
    
    class GenomeEntry < DBConnection
    	set_table_name 'genome_entry'
    	set_primary_key "auto_genome"
    end
    
    class RfamDatabaseLink < DBConnection
        belongs_to :rfam, :foreign_key => "auto_rfam"
    end
    
    # A view
    class RfamXrefSeq < DBConnection
    	belongs_to :rfam, :foreign_key => "auto_rfam"
    
    end
        
    class TaxLookup < DBConnection
    	belongs_to :tax_supergroup, :foreign_key => "supergroup_id"
    end
    
    class TaxSupergroup < DBConnection
    	has_many :tax_lookups, :foreign_key => "supergroup_id"
    	has_many :rfam_counts, :foreign_key => "tax_supergroup_id"
    end 	
    
    class RfamCount < DBConnection
    	belongs_to :rfam, :foreign_key => "rfam_id"
    	belongs_to :tax_supergroup, :foreign_key => "tax_supergroup_id"
    	
    	def lineage_specific?
    		unique = true
    		return false if self.species_count == 0
    		self.rfam.rfam_counts.each do |c_count|
					unique = false if c_count.tax_supergroup.name != self.tax_supergroup.name and c_count.species_count > 0
				end
				return unique
    	end
    	
    end
    
    class Clan < DBConnection
      set_primary_key "auto_clan"
    	has_many :clan_memberships, :foreign_key => "auto_clan"
    	has_many :rfams, :through => :clan_memberships
    	has_many :clan_literature_references, :foreign_key => "auto_clan"
    	has_many :literature_references, :through => "clan_literature_references"
    end
    
    class ClanMembership < DBConnection
      belongs_to :rfam, :foreign_key => "auto_rfam"
      belongs_to :clan, :foreign_key => "auto_clan"
    end
    
    class ClanLiteratureReference < DBConnection
    	belongs_to :clan, :foreign_key => "auto_clan"
    	belongs_to :literature_reference, :foreign_key => "auto_lit"
    end
    
    class LiteratureReference < DBConnection
			set_primary_key 'auto_lit'
    	has_many :clan_literature_references
    	
    end
    
  end
  
end
