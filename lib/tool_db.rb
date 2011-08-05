
Tool_DB_ADAPTER = 'postgresql'
Tool_DATABASE = "tooldb"
Tool_DB_HOST = 'localhost'
Tool_DB_USERNAME = 'tools'
Tool_DB_PASSWORD = 'analysis'

module Toolbox 
  module ToolDB
    
    require 'composite_primary_keys'
    
    include ActiveRecord
    
    class DBConnection < ActiveRecord::Base
      self.abstract_class = true
    
      def self.connect(version="")

        establish_connection(
                              :adapter => Tool_DB_ADAPTER,
                              :host => Tool_DB_HOST,
                              :database => "#{Tool_DATABASE}#{version}",
                              :username => Tool_DB_USERNAME,
                              :password => Tool_DB_PASSWORD
                              #:port => port
                            )
      end
    
    end
  
    class Organism < DBConnection
      has_many :mirbase_targets
      has_many :ensembl_ncrnas
      has_many :ensembl_hostgenes
      has_many :xref_dataset_organisms
      has_many :datasets, :through => :xref_dataset_organisms
      
      def fetch_ncrnas(type)
        return ToolDB::EnsemblNcrna.find(:all, :conditions => ["organism_id = ? and biotype = ?", self.id,"#{type}"])
      end
      
      def connect_to_genome(version)
        Ensembl::Core::DBConnection.connect("#{self.name.gsub(/\s/, '_').downcase}",version.to_i)
      end
      
      def snake_case
        return self.name.downcase.gsub(/\s/, '_')
      end
      
      # = DESCRIPTION
      # creates a 10 character version of the organism name
      def phylip_name
        
        return self.name.gsub(/\s/, '_')[0..9]
        
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
    
    class GoTerm < DBConnection
      
      has_many :go_to_gos
      has_one :go_description
      has_one :go_slim
      
      set_primary_key 'id'
      
    end  
     
    class EnsemblNcrna < DBConnection
      has_one :ensembl_hostgene, :through => :xref_ncrna_hostgene
      has_one :xref_ncrna_hostgene
      has_one :ncrna_group
      has_one :ncrna_sequence
      has_one :xref_ncrna_mirbase
      has_one :mirna_accession, :through => :xref_ncrna_mirbase
      belongs_to :organism
      has_one :grouping, :through => :ncrna_group

      # = DESCRIPTION
      # Returns only npcRNAs located in the genome alignment.
      def self.fetch_conserved(biotype)
        
        answer = Array.new
        ncrnas = ToolDB::EnsemblNcrna.find_all_by_biotype(biotype)
        ncrnas.each do |ncrna|
          answer.push(ncrna) if ncrna.has_group? 
        end
        
        return answer
        
      end
      
      def deepest_node_name
        self.has_group? ? self.grouping.deepest_node_name : self.organism.name
      end
      
      def deepest_node_nr
        self.has_group? ? self.grouping.deepest_node_nr : 0
      end
        
      # = DESCRIPTION
      # Returns the corresponding sequencing from ncrna_sequence
      def seq
        return "" if self.ncrna_sequence.nil?
        return self.ncrna_sequence.seq
      end
      
      def find_external_id
        
        return self.external_id if self.biotype == "snoRNA"
        return self.mirbase_family if self.biotype == "miRNA"
        
      end 
      
      def naseq
        return Bio::Sequence::NA.new(self.seq)
      end
      
      # = DESCRIPTION
      # Checks whether this ncRNA is sense or
      # anti sense to it's host gene
      def is_sense?        
        return false unless self.is_intronic?       
        "#{self.strand}" == "#{self.fetch_hostgene.strand}" ? answer = true : answer = false        
        return answer        
      end
      
      # = DESCRIPTION
      # Returns the Ensembl::Core::Gene object
      # for this ncRNA.  
      def fetch_ensembl_gene        
        self.organism.ensembl_db_connection
        return Ensembl::Core::Gene.find_by_stable_id("#{self.stable_id}")      
      end
      
      # = DESCRIPTION
      # Returns the Ensembl::Compara::Member object
      # for a protein-coding host gene. Used to 
      # retrieve homologs.
      def fetch_ensembl_member(version=54)       
        Ensembl::Compara::DBConnection.connect(version)
        return Ensembl::Compara::Member.find_by_stable_id(self.stable_id)        
      end
      
      def to_fasta       
        return Bio::FastaFormat.new(Bio::Sequence::NA.new(self.seq).to_fasta(self.stable_id))       
      end
            
      def present_at?(node)
        
        return false unless self.has_group?
        return self.group.present_at?(node)
        
      end
        
      def has_group?
        answer = false    
        self.ncrna_group ? answer = true : answer = false
        return answer        
      end
        
      def group
        
        return self.ncrna_group.grouping
        
      end
      
      def fetch_group_members
        answer = Array.new
        
        group = self.ncrna_group.grouping
        group.ncrna_groups.each do |ngroup|
          answer.push(ngroup.ensembl_ncrna)
        end
        
        return answer
        
      end
      
      def fetch_node_presence
        
        answer = []
        self.fetch_nodes.each do |node|
          if "#{node.val}" == "2"
            answer.push("#{node.tree_node.nr}")
          end
        end
        return answer  
        
      end
        
      def fetch_nodes
        
        return self.ncrna_group.grouping.group_nodes
        
      end
      
      def overlaps_exon?
        
        return false unless self.has_hostgene?
        
        gene = self.fetch_hostgene.fetch_ensembl_gene
        exons = []
        gene.transcripts.each do |t|
          t.exons.each do |exon|
            return true if self.start >= exon.start and self.stop <= exon.stop
          end
          
        end
        
        return false
        
      end
      
      def is_intronic?
        self.has_hostgene?
      end
      
      def has_hostgene?
        self.ensembl_hostgene.nil? ? answer = false : answer = true
        return answer
      end
      
      def targets_hostgene?
        
        return false unless self.has_hostgene?
        
        hostgene = self.fetch_hostgene
        
        return false if self.best_mirbase_match.nil?
        
        Ensembl::Core::DBConnection.connect(self.organism.name.downcase.gsub(/\s/, '_'))
        
        self.best_mirbase_match.mirbase_products.each do |product|
          
          product.mirbase_targets.collect{|t| t.target }.each do |target|
            transcript = Ensembl::Core::Transcript.find_by_stable_id(target)
            next if transcript.nil?
            return true if transcript.gene.stable_id == hostgene.stable_id
            
          end
          
        end
        
        return false
        
      end
      
      def best_mirbase_match        
        raise "Must be a microRNA!" unless self.biotype == "miRNA"  
        return self.mirna_accession.name  
      end
        
    end
    
    class NcrnaSequence < DBConnection
      belongs_to :ensembl_ncrna, :foreign_key => "ensembl_ncrna_id"
      
    end
    
    class EnsemblHostgene < DBConnection
      has_many :xref_ncrna_hostgenes
      has_many :ensembl_ncrnas, :through => :xref_ncrna_hostgenes
      has_many :hostgene_goterms
      belongs_to :organism, :foreign_key => "organism_id"
      
      set_primary_key 'id'
      
      def hosts_ncrna?(type)
        answer = false
        self.ensembl_ncrnas.select{|n| n.biotype == "#{type}"}.empty? ? answer = false : answer = true
        return answer
      end
      
      def get_expression_data
        ExpressionDB::DBConnection.connect
        expressions = ExpressionDB::Expression.find_all_by_ensembl_gene_id(self.stable_id)
        answer = { :median => nil, :specificity => nil, :entropy => nil, :mean => nil }
        unless expressions.empty?
          answer[:median] = expressions.select{|e| e.tissue.name == "median" }[0].value
          answer[:specificity] = expressions.select{|e| e.tissue.name == "specificity" }[0].value
          answer[:entropy] = expressions.select{|e| e.tissue.name == "entropy" }[0].value
          answer[:mean] = expressions.select{|e| e.tissue.name == "mean" }[0].value
        end
        return answer
      end
      
      def present_at?(type,node)
        
        raise "Node must be either ToolDB::TreeNode or ToolDB::DolloNode" unless node.kind_of?(ToolDB::TreeNode) or node.kind_of?(ToolDB::DolloNode)
        
        self.fetch_ncrnas(type).each do |ncrna|
          return true if ncrna.present_at?(node)
        end
        
        return false
      
      end
      
      def fetch_all_ncrnas
        
        answer = []
        self.xref_ncrna_hostgenes.each do |xref|
          answer.push(xref.ensembl_ncrna)
        end
        
        return answer
        
      end
      
      def fetch_ncrnas(type)
        
        answer = []
        
        self.xref_ncrna_hostgenes.each do |xref|
          answer.push(xref.ensembl_ncrna) if xref.ensembl_ncrna.biotype == type
        end
        
        return answer 
        
      end
      
      def belongs_to_group?(type)
        
        ncrnas = self.fetch_ncrnas(type)
        
        return false if ncrnas.empty?
        
        ncrnas.each do |ncrna|
          return true if ncrna.has_group?
        end
        
        return false
        
      end
      
      def has_ncrna?(type)
        
        if self.fetch_ncrnas(type).empty?
          return false
        else
          return true
        end
      
      end
      
      def fetch_ensembl_gene(ensembl_version=54)
        
        self.organism.ensembl_db_connection(ensembl_version)
        return gene = Ensembl::Core::Gene.find_by_stable_id(self.stable_id)
      
      end  
      
      def fetch_ensembl_member(ensembl_version=54)
        
        Ensembl::Compara::DBConnection.connect(ensembl_version)
        return Ensembl::Compara::Member.find_by_stable_id(self.stable_id)
        
      end
      
      def fetch_pan_member(ensembl_version=52)
        
        return Ensembl::PanCompara::Member.find_by_stable_id(self.stable_id)
        
      end
        
      def fetch_goterms
        
        answer = Array.new
        self.hostgene_goterms.each do |go|
          answer.push(go)
        end
        
        return answer
        
      end
      
      def probe
        xref = Toolbox::ToolDB::XrefProbeEnsembl.find_by_ensembl_stable_id(self.stable_id)
        return xref.probe unless xref.nil?
        return nil
      end
      
      def expression_stats
        
        filter = [ "mean", "median", "sum_of_all", "specificity", "entropy"]
        return self.expressions.select{|e| filter.include?("#{e.expression_tissue.name}")}.sort_by{|e| e.expression_tissue.name}
        
      end
        
      def expression_levels
        
        filter = [ "mean", "median", "sum_of_all", "specificity", "entropy"]
        return self.expressions.select{|e| filter.include?("#{e.expression_tissue.name}") == false }
        
      end
       
    end
    
    class NcrnaGroup < DBConnection
      belongs_to :ensembl_ncrna, :foreign_key => "ensembl_ncrna_id"
      belongs_to :grouping, :foreign_key => 'grouping_id'
      
      def nodes
        
        answer = Array.new
        
        self.group_nodes.each do |node|
          answer.push(node)
        end
        
        return answer
        
      end
      
    end
    
    class GroupNode < DBConnection
      belongs_to :grouping, :foreign_key => "grouping_id"
      belongs_to :tree_node, :foreign_key => "tree_node_id"
      
    end  
    
    class XrefNcrnaHostgene < DBConnection
      belongs_to :ensembl_ncrna, :foreign_key => "ensembl_ncrna_id"
      belongs_to :ensembl_hostgene, :foreign_key => "ensembl_hostgene_id"
    end
    
    class Grouping < DBConnection
      has_many :group_nodes
      has_many :ncrna_groups
      has_many :ensembl_ncrnas, :through => :ncrna_groups, :order => "organism_id"
      belongs_to :genomic_align_block, :foreign_key => "genomic_align_block_id"
      has_many :tree_nodes, :through => :group_nodes
      has_many :gerps
      has_many :xref_grouping_dollo_nodes
      has_many :dollo_nodes, :through => :xref_grouping_dollo_nodes
      belongs_to :dataset, :foreign_key => "dataset_id"
      
      
      def deepest_node_name
        node = self.deepest_node
        if node.nil?
          return self.ensembl_ncrnas[0].organism.name
        else
          return node.description
        end
      end
      
      def deepest_node_nr
        node = self.deepest_node
        if node.nil?
          return 0
        else
          return self.deepest_node.node.to_i
        end
      end
      
      def deepest_node
        nodes = ToolDB::DolloNode.find_all_by_dollo_tree_id(self.dataset.dollo_tree.id).sort_by{|n| n.id}
        nodes.each do |node|
          return node if self.dollo_state_at_node(node) == "1"
        end
        return nil
      end
        
      # = DESCRIPTION
      # Returns the ancestral state from the dollo reconstruction
      def dollo_state_at_node(node)
        
        node = ToolDB::DolloNode.find_by_node(node) unless node.kind_of?(ToolDB::DolloNode)
        
        state = self.xref_grouping_dollo_nodes.select{|x| x.dollo_node_id == node.id}.shift.state.strip
        
        if state == "."
          return "0" if node.node == "1"
          parent = node.get_parent
          self.dollo_state_at_node(parent)
        else
          return state
        end
                
      end
      
      # = DESCRIPTION
      # Collects all ncRNAs belonging to a given
      # group and returns them as an Array of 
      # ToolDB::EnsemblNcrna objects.
      def fetch_ncrnas
        
        answer = Array.new
        self.ncrna_groups.each do |ncrnagroup|
          answer.push(ncrnagroup.ensembl_ncrna)
        end
        
        return answer
        
      end
      
      # = DESCRIPTION
      # FInds the host genes of all member ncRNAs
      def fetch_hostgenes(organism_id=nil)
        
        answer = []
        
        if organism_id.nil?
          self.ensembl_ncrnas.each do |ncrna|
            answer.push(ncrna.fetch_hostgene) unless ncrna.fetch_hostgene.nil?
          end
        else
          self.ensembl_ncrnas.each do |ncrna|
            next if ncrna.organism_id != organism_id
            answer.push(ncrna.fetch_hostgene) unless ncrna.fetch_hostgene.nil?
          end
        end
        
        return answer
        
      end  
      
      # = Description
      # Checks if the group can be classified as intronic,
      # based on any of the member ncRNAs
      def is_intronic?
        
        self.ensembl_ncrnas.each do |ncrna|
          return true if ncrna.is_intronic?
        end
        
        return false
        
      end
      
      # = DESCRIPTION
      # Checks whether this ncRNA is transcribed sense or anti-sense
      # to the host gene
      def is_sense?
        
        return "Not intronic!" unless self.is_intronic?
        
        self.ensembl_ncrnas.each do |ncrna|
          return true if ncrna.is_sense?
        end
        
        return false
        
      end
        
      # DESCRIPTION
      # Returns the gerp scores for this section of the alignment
      def gerp_scores(window_size)
        
        return [] if self.gerps.empty?
        return self.gerps.select{|g| g.window_size == window_size}.shift.scores.split(" ")
      
      end  
      # = DESCRIPTION 
      # finds the corresponding RFam accession for a group
      # of ncRNAs (or else miRBase)
      
      # = DESCRIPTION
      # Calculates the mean of the gerp scores for 
      # a given window size
      def gerp_mean(window_size)
        
        answer = 0.0
        
        gerp_scores = self.gerp_scores(window_size)
        
        return "N/A" if gerp_scores.empty?
        
        gerp_scores(window_size).each{|s| answer += s.to_f }
                
        return answer/gerp_scores.nitems
        
      end
      
      # = DESCRIPTION
      # Calculates the median of gerp scores
      # for a given window size
      def gerp_median(window_size)
        
        answer = 0.0
        
        gerp_scores = self.gerp_scores(window_size).sort.delete_if{|s| s == 0.0}.collect{|s| s.to_f}
        
        n = (gerp_scores.length - 1) / 2 
        n2 = (gerp_scores.length) / 2 
        
        if gerp_scores.length % 2 == 0
          return (gerp_scores[n] + gerp_scores[n2]) / 2
        else
          return gerp_scores[n]
        end
         
      end  
      
      def fetch_ncrna_plain_name
        
        if self.biotype == "miRNA"
          return self.fetch_ncrna_group
        else
          return ToolDB::Rfam.find_by_acc("#{self.fetch_ncrna_group}").name
        end
        
      end
      
      def rfam_clan
      
      	rfam_id = self.fetch_ncrna_group
      	rfam = ToolDB::Rfam.find_by_accession(rfam_id).rfam_clan
      	
      end
      
      def rfam_clan_name
      	if self.rfam_clan
      		return self.rfam_clan.name
      	else
      		return ""
      	end
      end
      
      def external_id 
        return self.fetch_ncrna_group
      end
      
      def fetch_ncrna_group
        
        answer = nil
        self.ensembl_ncrnas.each do |ncrna|
          if self.biotype == "miRNA"
            answer = ncrna.mirbase_family unless ncrna.mirbase_family.nil?
          else
            answer = "#{ncrna.external_id.strip}" unless ncrna.external_id.nil?
          end
        end

        warn "Group #{self.nr} has no ncrna group" if answer.nil?

        return answer
        
      end
      
      # = DESCRIPTION
      # Creates and returns a hash of nodes and the
      # corresponding values for a given group.
      def node_values
        
        answer = Hash.new
        
        self.group_nodes.each do |gr_node|
          answer["#{gr_node.tree_node.nr}"] = "#{gr_node.val}"
        end
        
        return answer 
        
      end
      
      # = DESCRIPTION
      # Checks whether a group (ToolDB::Grouping object)
      # is present_at? at a given node (ToolDB::TreeNode object)
      # returns true or false
      def present_at?(node)
        
        if node.kind_of?(ToolDB::TreeNode)
          
          self.group_nodes.each do |gr_node|
            if gr_node.tree_node == node and "#{gr_node.val}" == "1"
              return true
            end
          end
        
          return false
          
        elsif node.kind_of?(ToolDB::DolloNode)
          
          state = self.dollo_state_at_node(node)
          return false if state == "0"
          return true if state == "1"
          
        else
          
          raise "Node must be either ToolDB::TreeNode or ToolDB::DolloNode"
          
        end
        
      end
      
      # = DESCRIPTTION
      # Finds all GO terms associated with any host genes
      # belonging to all ncRNAs that are children to this
      # group
      def fetch_unique_goterms
        
        answer = Array.new
        genes = self.fetch_hostgenes
        genes.each do |gene|
          gene.fetch_goterms.each do |go|
            answer.push(go.go_term) unless answer.include?(go.go_term)
          end
        end
        
        return answer
        
      end
      
      # = DESCRIPTION
      # Returns all organisms that are part of this group through their ncRNAs
      def fetch_organisms       
        return self.ensembl_ncrnas.collect{|ncrna| ncrna.organism }      
      end
      
      # = DESCRIPTION
      # Returns the genomic alignment from Ensembl for this ncRNA group
      def fetch_genomic_alignment
        gal_block = self.fetch_ensembl_genomic_align_block(57)
        aln = gal_block.get_clustalw(self.start,self.stop)
        return aln
      end
      
      # = DESCRIPTION
      # Returns the correspoing Ensembl::Compara::GenomicAlignBlock object
      def fetch_ensembl_genomic_align_block(version=54)
        
        Ensembl::Compara::DBConnection.connect(version)
        return Ensembl::Compara::GenomicAlignBlock.find(self.genomic_align_block.id)
        
      end
      
      # = DESCRIPTION
      # Returns a Bio::Alignment object of aligned ncRNA sequences
      # Can be "fasta" or "clustalw"
      def do_align(type="fasta")
      
        raise "not enough sequences" if self.ensembl_ncrnas.nitems == 1
        
        a = Bio::Alignment.new(self.ensembl_ncrnas.collect{|n| n.to_fasta })
        factory = Bio::ClustalW.new

        return a.do_align(factory) if type == "fasta"
        return factory.query(a) if type == "clustalw"
        
      end
            
    end
    
    class Gerp < DBConnection
      belongs_to :grouping, :foreign_key => "grouping_id"
    end
    
    class HostgeneGoterm < DBConnection
      belongs_to :ensembl_hostgene, :foreign_key => "ensembl_hostgene_id"
      
      def is_slim?
        return ToolDB::GoTerm.find_by_term(self.go_term).is_slim?
      end
      
      def name
        ToolDB::GoTerm.find_by_term(self.go_term).name
      end
      
      def description
        ToolDB::GoTerm.find_by_term(self.go_term).go_description.description
      end
      
      def namespace
        ToolDB::GoTerm.find_by_term("#{self.go_term}").namespace
      end
        
    end
  
    class TreeNode < DBConnection
      
      has_many :group_nodes
      has_many :groupings, :through => :group_nodes, :order => "nr", :conditions => ["\"group_node\".val = 2"]
      
      # = DESCRIPTION
      # Returns all groups that are present at this particular node
      def fetch_present_groups(biotype)
        
        answer = Array.new
        self.groupings.each do |group|
          answer.push(group) if group.biotype == "#{biotype}"
        end
        return answer
        
      end
      
      def fetch_hostgenes(organism_id,type)
        
        answer = []
        
        self.fetch_present_groups(type).each do |group|
          hostgene = group.fetch_hostgenes(organism_id).shift
          answer.push(hostgene) unless answer.include?(hostgene) or hostgene.nil?
        end
        
        return answer
        
      end
      
      def fetch_all_hostgenes
        
        answer = []
        self.groupings.each do |group|
          group.ensembl_ncrnas.each do |ncrna|
            answer.push(ncrna.ensembl_hostgene) unless ncrna.ensembl_hostgene.nil?
          end
          
        end
        
        return answer
      end
      
    end
    
    class DolloTree < DBConnection
      has_many :dollo_nodes
      belongs_to :dataset, :foreign_key => "dataset_id"
    end
    
    class DolloNode < DBConnection
      
      has_many :xref_grouping_dollo_nodes
      has_many :xref_family_nodes
      belongs_to :dollo_tree, :foreign_key => "dollo_tree_id"
      has_one :dollo_node, :foreign_key => "child_of"
      
      def get_parent
        return nil if self.child_of.nil?
        return ToolDB::DolloNode.find(self.child_of)
      end
      
      def get_children
        return ToolDB::DolloNode.find_all_by_child_of(self.id)
      end
      
    end
    
    class XrefGroupingDolloNode < DBConnection
      belongs_to :grouping, :foreign_key => "grouping_id"
      belongs_to :dollo_node, :foreign_key => "dollo_node_id"
      
    end
    
    class GenomicAlignBlock < DBConnection
      set_primary_key "id"
      has_many :groupings
      belongs_to :dataset, :foreign_key => "dataset_id"
    end
    
    class Family < DBConnection
      
      has_many :xref_family_nodes
      
      # = DESCRIPTION
      # Returns the ancestral state from the dollo reconstruction
      def dollo_state_at_node(node)

        node = ToolDB::DolloNode.find_by_node(node) unless node.kind_of?(ToolDB::DolloNode)

        state = self.xref_family_nodes.select{|x| x.dollo_node_id == node.id}.shift.state.strip

        if state == "."
          return "0" if node.node == "1"
          parent = node.get_parent
          self.dollo_state_at_node(parent)
        else
          return state
        end

      end
      
      def present_at?(node)
        
        answer = nil
        self.dollo_state_at_node(node) == "1" ? answer = true : answer = false
        return answer
        
      end
    
    end
    
    class XrefFamilyNode < DBConnection
      
      belongs_to :family, :foreign_key => "family_id"
      belongs_to :dollo_node, :foreign_key => "dollo_node_id"
      
    end
    
    
  
    
  end
end
