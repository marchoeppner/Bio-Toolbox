
Compara_DB_ADAPTER = 'postgresql'
Compara_DATABASE = "compara"
Compara_DB_HOST = 'localhost'
Compara_DB_USERNAME = 'tools'
Compara_DB_PASSWORD = 'analysis'

module Toolbox 
  module ComparaDB
    
    require 'composite_primary_keys'
    
    include ActiveRecord
    
    class DBConnection < ActiveRecord::Base
      self.abstract_class = true
    
      def self.connect(version="")

        establish_connection(
                              :adapter => Compara_DB_ADAPTER,
                              :host => Compara_DB_HOST,
                              :database => "#{Compara_DATABASE}#{version}",
                              :username => Compara_DB_USERNAME,
                              :password => Compara_DB_PASSWORD
                              #:port => port
                            )
      end
    
    end
  
    class Organism < DBConnection      
        has_many :ncrnas
        has_many :genes
        
        def connect_to_genome(version)
          Ensembl::Core::DBConnection.connect("#{self.name.gsub(/\s/, '_').downcase}",version.to_i)
        end
        
    end
    
    class Dataset < DBConnection
      has_many :groupings, :foreign_key => "dataset_id"
      has_many :genomic_align_blocks
      has_many :xref_dataset_organisms
      has_many :organisms, :through => :xref_dataset_organisms
      has_one :dollo_tree
      has_many :dollo_nodes, :through => :dollo_tree
    end
    
    class XrefDatasetOrganism < DBConnection
      set_primary_keys :dataset_id,:organism_id
      belongs_to :dataset, :foreign_key => "dataset_id"
      belongs_to :organism, :foreign_key => "organism_id"
    end
     
    class Ncrna < DBConnection
      belongs_to :organism, :foreign_key => "organism_id"
      has_one :xref_ncrna_gene
      has_one :gene, :through => :xref_ncrna_gene
      has_many :xref_ncrna_groupings
      has_many :groupings, :through => :xref_ncrna_groupings
      
      def to_fasta        
        return Bio::FastaFormat.new(Bio::Sequence::NA.new(self.sequence).to_fasta(self.stable_id))     
      end
      
      def deepest_node_name(dataset)
      	
      	grouping = self.groupings.select{|g| g.dataset_id == dataset.id}[0]
      	return "Species-specific" if grouping.nil?
      	
      	return grouping.deepest_node_name
      
      end
      
    end
    
    
    class Gene < DBConnection
      belongs_to :organism, :foreign_key => "organism_id"
      has_many :xref_ncrna_genes
      has_many :ncrnas, :through => :xref_ncrna_genes
    end
    
    class Grouping < DBConnection
      belongs_to :dataset, :foreign_key => "dataset_id"
      belongs_to :genomic_align_block, :foreign_key => "genomic_align_block_id"
      has_many :xref_ncrna_groupings
      has_many :xref_grouping_nodes
      has_many :ncrnas, :through => :xref_ncrna_groupings
      
      def rfam_acc
        self.ncrnas.each do |ncrna|
          return ncrna.rfam_acc unless ncrna.rfam_acc.nil?
        end
        return nil
      end
      
      def state_at_node?(node) 				
        node = ComparaDB::DolloNode.find_by_node(node) unless node.kind_of?(ComparaDB::DolloNode)            
        state = self.xref_grouping_nodes.select{|x| x.dollo_node_id == node.id }.shift.state 				
        if state == "." and node.node.include?("root")
          return "0"
        elsif state == "."
          self.state_at_node?(node.get_parent)
        else
          return state
        end 			
      end

      def deepest_node
        nodes = self.xref_grouping_nodes.collect{|x| x.dollo_node}.sort_by{|n| n.id}
        nodes.each do |node|
          return node if self.state_at_node?(node) == "1"
        end
        return nil
      end

      def deepest_node_name
        deepest_node = self.deepest_node
        if deepest_node.nil?
          return "Species-specific"
        else
          return deepest_node.description
        end
      end     
      
      def align
        
        ff = Bio::MultiFastaFormat.new(self.ncrnas.collect{|n| n.to_fasta})
        factory = Bio::Clustalw.new
        aln = ff.do_align(factory)
        return aln
        
      end
      
    end
    
    class XrefNcrnaGrouping < DBConnection
      set_primary_keys :ncrna_id,:grouping_id
      belongs_to :ncrna, :foreign_key => "ncrna_id"
      belongs_to :grouping, :foreign_key => "grouping_id"
      
    end  
    
    class XrefNcrnaGene < DBConnection
      set_primary_keys :ncrna_id,:gene_id
      belongs_to :ncrna, :foreign_key => "ncrna_id"
      belongs_to :gene, :foreign_key => "gene_id"
    end
  
    
    class DolloTree < DBConnection
      has_many :dollo_nodes
      belongs_to :dataset, :foreign_key => "dataset_id"
      
      def output_phylip
        tree_string = self.tree_string
        tree_organisms = self.tree.leaves.collect{|t| t.name.strip }
        tree_organisms.each do |t_o|
          phylip_name = t_o.gsub(/\s/, '_').phylip_name
          t_o = t_o.gsub(/\s/ , '_')
          tree_string.gsub!(/#{t_o}/, "#{phylip_name}")
        end
        return tree_string
      end
      
      def tree
         return Bio::Newick.new(self.tree_string).tree
      end
      
    end
    
    class DolloNode < DBConnection
      
      has_many :xref_grouping_nodes     
      belongs_to :dollo_tree, :foreign_key => "dollo_tree_id"
      has_one :dollo_node, :foreign_key => "child_of"
         
      def get_parent
        return nil if self.child_of.nil?
        return ComparaDB::DolloNode.find(self.child_of)
      end
      
      def get_children
        return ComparaDB::DolloNode.find_all_by_child_of(self.id)
      end
      
    end
    
    class XrefGroupingNode < DBConnection
      set_primary_keys :grouping_id,:dollo_node_id
      belongs_to :grouping, :foreign_key => "grouping_id"
      belongs_to :dollo_node, :foreign_key => "dollo_node_id"      
    end
    
    class GenomicAlignBlock < DBConnection
      set_primary_key "id"
      has_many :groupings
      belongs_to :dataset, :foreign_key => "dataset_id"
    end
   
  end
end
