
Ncrna_DB_ADAPTER = 'postgresql'
Ncrna_DATABASE = "ncrna"
Ncrna_DB_HOST = 'localhost'
Ncrna_DB_USERNAME = 'tools'
Ncrna_DB_PASSWORD = 'analysis'

module Toolbox  
  module NcrnaDB
    
    include ActiveRecord
    
    class DBConnection < ActiveRecord::Base
      self.abstract_class = true
    	self.pluralize_table_names = false
    	
      def self.connect(version="56")

        establish_connection(
                              :adapter => Ncrna_DB_ADAPTER,
                              :host => Ncrna_DB_HOST,
                              :database => "#{Ncrna_DATABASE}#{version}",
                              :username => Ncrna_DB_USERNAME,
                              :password => Ncrna_DB_PASSWORD
                              #:port => port
                            )
      end
      
    end
     
    class Organism < DBConnection     
       has_many :genes
       has_many :snornas
       
       def connect_to_genome
         Ensembl::Core::DBConnection.connect(self.database,version)     
       end    
        
    end
     
    class Gene < DBConnection    
      belongs_to :organism, :foreign_key => "organism_id"
      has_many :expressions
      has_many :ncrnas, :through => :xref_ncrna_genes

      def display_details
	return "#{self.stable_id}|#{self.start}>#{self.stop}|#{self.strand}|#{self.description[0..30]}"
      end

      def expression_data
        return nil if expressions.empty?
	answer = {}
	self.expressions.each do |e|
	  answer[e.expression_tissue.name] = e.value
	end
	return answer
      end
      
    end
     
    class Ncrna < DBConnection 
      belongs_to :organism, :foreign_key => "organism_id"
      has_one :xref_ncrna_gene
      has_one :gene, :through => :xref_ncrna_gene    
      
      def is_intronic?
	 if self.xref_ncrna_gene
	   return true
	 else
	   return false
	 end
      end
      
    end
    
    class HomologyGroup < DBConnection
      belongs_to :gene, :foreign_key => "gene_id"
      has_many :homology_members
    end
    
    class HomologyMember < DBConnection
      belongs_to :gene, :foreign_key => "gene_id"
      belongs_to :homology_group, :foreign_key => "homology_group_id"
    end
    
    class GoTerm < DBConnection
      belongs_to :gene, :foreign_key => "gene_id"
    end
    
    class ExpressionTissue < DBConnection
      has_many :expressions
    end
    
    class Expression < DBConnection
      belongs_to :expression_tissue, :foreign_key => "expression_tissue_id"
      belongs_to :gene, :foreign_key => "gene_id"

      def tissue
	return self.expression_tissue.name
      end
    end
    
    class Grouping < DBConnection
      has_many :xref_ncrna_groupings
      has_many :ncrnas, :through => :xref_ncrna_groupings
      belongs_to :genomic_align_block, :foreign_key => "genomic_align_block_id"
      has_many :xref_grouping_dollos
      has_many :dollos, :through => :xref_grouping_dollos
      
      def dollo_state_at_node(node)      
        node = NcrnaDB::Dollo.find_by_node(node) unless node.kind_of?(NcrnaDB::Dollo)    
        state = self.xref_grouping_dollos.select{|x| x.dollo_id == node.id}.shift.state.strip
	#puts "Checking numeric state at node #{node.node} (is #{state})"
        if state == "."
          if node.node == 1
	    return 0
	  else
	     self.dollo_state_at_node(node.get_parent)
	  end
	else
	  return state.to_i
        end             
      end
      
      def external_acc
	return self.ncrnas.collect{|n| n.external_acc if n.external_acc }.shift
      end
      
      def is_intronic?
	self.ncrnas.each do |ncrna|
	  return true if ncrna.is_intronic?
	end
      end
      
      def deepest_node
	nodes = NcrnaDB::Dollo.find(:all, :order => "id DESC")
	answer = 0
	nodes.each do |node|
	  answer = node.node if self.dollo_state_at_node(node) == 1
	end
	return answer
      end
     
    end
    
    class XrefNcrnaGrouping < DBConnection
      belongs_to :ncrna, :foreign_key => "ncrna_id"
      belongs_to :grouping, :foreign_key => "grouping_id"
    end
    
    class XrefNcrnaGene < DBConnection
      belongs_to :ncrna, :foreign_key => "ncrna_id"
      belongs_to :gene, :foreign_key => "gene_id"
    end
    
    class Dollo < DBConnection
      has_many :xref_grouping_dollos
      has_many :groupings, :through => :xref_grouping_dollos
      
      def get_parent
        return nil if self.child_of.nil?
        return NcrnaDB::Dollo.find_by_node(self.child_of)
      end
      
      def get_children
        return NcrnaDB::Dollo.find_all_by_child_of(self.id)
      end
    
    end
    
    class XrefGroupingDollo < DBConnection
      belongs_to :dollo, :foreign_key => "dollo_id"
      belongs_to :grouping, :foreign_key => "grouping_id"
    end
    
    class XrefGroupingGab < DBConnection
      belongs_to :grouping, :foreign_key => "grouping_id"
      belongs_to :genomic_align_block, :foreign_key => "genomic_align_block_id"
    end
    
    class GenomicAlignBlock < DBConnection
      has_many :groupings
    end
    
  end   
end
  

