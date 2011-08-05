Homolog_DB_ADAPTER = 'postgresql'
Homolog_DATABASE = "homolog"
Homolog_DB_HOST = 'localhost'
Homolog_DB_USERNAME = 'tools'
Homolog_DB_PASSWORD = 'analysis'

module Toolbox 
  module HomologDB
    
    require 'composite_primary_keys'
    
    include ActiveRecord
    
    class DBConnection < ActiveRecord::Base
      self.abstract_class = true
    	self.pluralize_table_names = false
    	
      def self.connect(version="")

        establish_connection(
                              :adapter => Homolog_DB_ADAPTER,
                              :host => Homolog_DB_HOST,
                              :database => "#{Homolog_DATABASE}#{version}",
                              :username => Homolog_DB_USERNAME,
                              :password => Homolog_DB_PASSWORD
                              #:port => port
                            )
      end
    
    end
       
    class Organism < DBConnection
      has_many :members
      has_many :infernal_hits
      has_many :homolog_groups, :through => :members
      has_many :ncrna_members
      has_many :ncrnas
      has_many :rrnas, :foreign_key => "organism_id"
      has_many :xref_ncrna_rrnas, :through => :rrnas
      has_many :infernals
      has_many :xref_infernal_rrnas, :through => :rrnas
      belongs_to :supergroup, :foreign_key => "supergroup_id"
      
      def db_args
      	args = {}
      	self.database.split(";").each do |element|
      		e = element.split("=")
      		args[e[0].strip] = e[1].strip
      		if e[0] == "dbname"
      			args["version"] = e[1].slice(/[56][0-9]/).to_i
      		end
      	end
      	return args
     	end
     	
      def connect_to_genome(version=58)
        if self.database.include?("eukaryotes")
        	version = self.database.slice(/[0-9]*$/).to_i
          EukaryoteDB::DBConnection.connect(version)
        elsif self.db_args["host"] == "mysql.ebi.ac.uk"
        	Ensembl::Core::DBConnection.connect("#{self.snake_case}",version,{:ensembl_genomes => true })#, :database => self.db_args["dbname"]})  
        else
        	Ensembl::Core::DBConnection.connect("#{self.snake_case}",version)
        end
      end
      
      def connect_to_genome_new(version=60)
         if self.database.include?("eukaryote")
            version = self.database.slice(/[0-9]*$/).to_i
            EukaryoteDB::DBConnection.connect(version) 
         else
            Ensembl::Compara::DBConnection.connect(version,{:ensembl_genomes => true, :local => true})
            organism = Ensembl::Compara::GenomeDb.find_by_name(self.name.to_snake_case)
            organism.connect_to_genome_new(version)
         end  
      
      end
      
      def snake_case
      	return self.name.gsub(/\s/, '_').downcase
      end
      
      def infernal_hits
        answer = []
        self.members.each do |m|
          m.infernal_hits.each do |hit|
            answer << hit
          end
        end
        return answer
      end
      
      def lsu
      	return self.rrnas.select{|r| r.typ == "lsu"}.shift
      end
      
      def ssu
      	return self.rrnas.select{|r| r.typ == "ssu"}.shift
			end
			      
    end
    
    class Rrna < DBConnection
    
    	belongs_to :organism, :foreign_key => "organism_id"
    	has_many :xref_ncrna_rrnas
    	has_many :xref_infernal_rrnas
    	
    	def aligned_seq
    		return AlignSeq.new(self.sequence,self.cigar_line).align
    	end
    	
    end
    
    class RrnaLocus < DBConnection
    	has_many :rrna_locus_xref_ncrnas, :foreign_key => "rrna_locus_id"
    	has_many :rrna_locus_xref_infernals , :foreign_key => "rrna_locus_id"
    	has_many :xref_ncrna_rrnas, :through => :rrna_locus_xref_ncrnas
    	has_many :xref_infernal_rrnas, :through => :rrna_locus_xref_infernals
    	has_many :xref_rrna_locus_nodes
    	has_many :dollo_nodes, :through => :xref_rrna_locus_nodes
    	
    	def deepest_node				
 				nodes = self.xref_rrna_locus_nodes.collect{|x| x.dollo_node}.sort_by{|n| n.id}
 				nodes.each do |node|
 					return node if self.state_at_node?(node) == "1"
 				end				
	 			return nil			
 			end
 			
 			def deepest_node_name
 				deepest_node = self.deepest_node
 				if deepest_node.nil?
 					return ""
 				else
 					return deepest_node.description
 				end
 			end
 			
 			def state_at_node?(node) 				
 				node = HomologDB::DolloNode.find_by_name(node) unless node.kind_of?(HomologDB::DolloNode) 				
 				state = self.xref_rrna_locus_nodes.select{|x| x.dollo_node_id == node.id }.shift.state				
 				if state == "." and node.name.include?("root")
 					return "0"
 				elsif state == "0"
 					return "0"
 				elsif state == "."
 					self.state_at_node?(node.get_parent)
 				else
 					return state
 				end 			
 			end
 			
    end
    
    class XrefRrnaLocusNode < DBConnection
    	belongs_to :rrna_locus, :foreign_key => "rrna_locus_id"
    	belongs_to :dollo_node, :foreign_key => "dollo_node_id"
    end
    
    class RrnaLocusXrefNcrna < DBConnection
    	belongs_to :rrna_locus, :foreign_key => "rrna_locus_id"
    	belongs_to :xref_ncrna_rrna, :foreign_key => "xref_ncrna_rrna_id"
    	has_one :ncrna, :through => :xref_ncrna_rrna
    end
    
    class RrnaLocusXrefInfernal < DBConnection
    	belongs_to :rrna_locus, :foreign_key => "rrna_locus_id"
    	belongs_to :xref_infernal_rrna, :foreign_key => "xref_infernal_rrna_id"
    	has_one :infernal, :through => :xref_infernal_rrna
    end
    
    class Ncrna < DBConnection
      belongs_to :organism, :foreign_key => "organism_id"  
      has_many :xref_ncrna_rrnas
      
      def is_intronic?
        
        organism.connect_to_genome(60)
        
        if organism.database.include?("eukaryote")
          gene = EukaryoteDB::Gene.find_by_stable_id(self.stable_id)
          return gene.is_intronic?
        else
          gene = Ensembl::Core::Gene.find_by_stable_id(self.stable_id)
          return false if gene.nil?
          warn "#{self.stable_id} not in DB, setting false" if gene.nil?
          return gene.is_intronic?
        end      
        
      end
      
      def get_hostgene
      
      	organism.connect_to_genome(60)
      	if organism.database.include?("eukaryote")
      		gene = EukaryoteDB::Gene.find_by_stable_id(self.stable_id)
      		return gene.slice.genes(true).select{|g| g!= gene}[0]
      	else
      		gene = Ensembl::Core::Gene.find_by_stable_id(self.stable_id)
      		return nil if gene.nil?
      		return gene.slice.genes(true).select{|g| g != gene }[0]
      	end
      
      end
      
      def gene
      	organism.connect_to_genome(60)
        
        if organism.database.include?("eukaryote")
          gene = EukaryoteDB::Gene.find_by_stable_id(self.stable_id)
        else
          gene = Ensembl::Core::Gene.find_by_stable_id(self.stable_id)
        end  
        return gene
      end    
      
      
      def to_fasta
      	return Bio::FastaFormat.new(Bio::Sequence::NA.new(self.sequence).to_fasta(self.stable_id))
      end
      
    end
    
    class XrefNcrnaRrna < DBConnection
	    set_primary_key 'id'
    	belongs_to :ncrna, :foreign_key => "ncrna_id"
    	belongs_to :rrna, :foreign_key => "rrna_id"
    	has_one :rrna_locus_xref_ncrna, :foreign_key => "xref_ncrna_rrna_id"
    	has_one :rrna_locus, :through => :rrna_locus
    	
    	def center_position
    		return self.rrna_from+((self.rrna_to-self.rrna_from)/2)
    	end
    	
    	def aligned_center_position
    		
    		pos = self.center_position
    		cigar = self.rrna.cigar_line
    		mapped_pos = Bio::Mapper.new(pos,cigar).map
    		return mapped_pos
    	
    	end
    	
    end
    
    class Infernal < DBConnection
    	belongs_to :organism, :foreign_key => "organism_id"
    	has_many :xref_infernal_rrnas
    	
    	def is_intronic?
    	  
    	  organism.connect_to_genome(60)
    	  
    	  if organism.database.include?("eukaryote")
    	    organism = EukaryoteDB::Organism.find_by_name(self.organism.name)
    	    dna = EukaryoteDB::Dna.find_by_organism_id_and_accession(organism.id,self.seq_region)
          slice = EukaryoteDB::Slice.new(dna.id,self.start,self.stop,self.strand)
          slice.genes(true).empty? ? answer = false : answer = true
          return answer
        else
          seq_region = Ensembl::Core::SeqRegion.find_by_name(self.seq_region)
          warn "Seq region not found (#{self.seq_region})" if seq_region.nil?   
          return false if seq_region.nil?     
          slice = Ensembl::Core::Slice.new(seq_region,self.start,self.stop,self.strand)
          slice.genes(true).empty? ? answer = false : answer = true
          return answer
        end
    	  
    	end
    	
    	def genomic_context
    	
    		organism.connect_to_genome(60)
    	  
    	  if organism.database.include?("eukaryote")
    	    organism = EukaryoteDB::Organism.find_by_name(self.organism.name)
    	    dna = EukaryoteDB::Dna.find_by_organism_id_and_accession(organism.id,self.seq_region)
          slice = EukaryoteDB::Slice.new(dna.id,self.start,self.stop,self.strand)         
         	genes = slice.genes(true).select{|g| g.biotype == "protein_coding"}
         	if genes.empty? == false
         		hostgene = genes[0]
         		return "INTRONIC:#{hostgene.stable_id}|#{self.rfam_acc}[#{self.strand}]|#{hostgene.stable_id}[#{hostgene.strand}]"
         	else
         		upstream = EukaryoteDB::Gene.find(:all, :conditions => ["dna_id = ? and start > ? and biotype = ?",dna.id,self.stop,"protein_coding"],:order => "start ASC")[0]
         		downstream = EukaryoteDB::Gene.find(:all, :conditions => ["dna_id = ? and stop < ? and biotype = ?",dna.id,self.start,"protein_coding"],:order => "start DESC")[0]
         		if upstream and downstream
  						updist = upstream.start-self.stop
  						downdist = self.start-downstream.stop
  					elsif upstream
  						updist = upstream.start-self.stop
  						downdist = "NA"
  					elsif downstream
  						updist = "NA"
  						downdist = self.start-downstream.stop
  					else
  						updist,downdist = "NA","NA"
  					end
  					upstream ? upstream_string = "#{upstream.stable_id}[#{upstream.strand}]" : upstream_string = "NA[NA]"
  					downstream ? downstream_string = "#{downstream.stable_id}[#{downstream.strand}]" : downstream_string = "NA[NA]"
  			
  					genomic_context= "#{downstream_string}|#{downdist}|#{self.rfam_acc}[#{self.strand}]|#{updist}|#{upstream_string}"

          	return "INTERGENIC:#{genomic_context}"
  					
         	end
         	
        else
          seq_region = Ensembl::Core::SeqRegion.find_by_name(self.seq_region)
          warn "Seq region not found (#{self.seq_region})" if seq_region.nil?   
          return false if seq_region.nil?     
          slice = Ensembl::Core::Slice.new(seq_region,self.start,self.stop,self.strand)
          genes = slice.genes(true).select{|g| g.biotype == "protein_coding"}
          genes.empty? ? answer = false : answer = true
          if answer        
          	hostgene = genes[0]
          	return "INTRONIC:#{hostgene.stable_id}|#{self.rfam_acc}[#{self.strand}]|#{hostgene.stable_id}[#{hostgene.strand}]"
          else
          	downstream = Ensembl::Core::Gene.find(:all, :conditions => ["seq_region_id = ? and seq_region_end < ? and biotype = ?","#{slice.seq_region.id}","#{slice.start}","protein_coding"], :order => "seq_region_start DESC")[0]
  					upstream = Ensembl::Core::Gene.find(:all, :conditions => ["seq_region_id = ? and seq_region_start > ? and biotype = ?","#{slice.seq_region.id}","#{slice.stop}","protein_coding"], :order => "seq_region_end ASC")[0]  
  		
  					if upstream and downstream
  						updist = upstream.start-self.stop
  						downdist = self.start-downstream.stop
  					elsif upstream
  						updist = upstream.start-self.stop
  						downdist = "NA"
  					elsif downstream
  						updist = "NA"
  						downdist = self.start-downstream.stop
  					else
  						updist,downdist = "NA","NA"
  					end
  			
  					upstream ? upstream_string = "#{upstream.stable_id}[#{upstream.strand}]" : upstream_string = "NA[NA]"
  					downstream ? downstream_string = "#{downstream.stable_id}[#{downstream.strand}]" : downstream_string = "NA[NA]"
  			
  					genomic_context= "#{downstream_string}|#{downdist}|#{self.rfam_acc}[#{self.strand}]|#{updist}|#{upstream_string}"

          	return "INTERGENIC:#{genomic_context}"
          end
        end
    	
    	end
    	
    	
    end
    
    class XrefInfernalRrna < DBConnection
    	belongs_to :infernal, :foreign_key => "infernal_id"
    	belongs_to :rrna, :foreign_key => "rrna_id"
    	has_one :rrna_locus_xref_infernal #, :foreign_key => "xref_infernal_rrna_id"
    	has_one :rrna_locus, :through => :rrna_locus_xref_infernal
    	
    	def center_position
    		return self.rrna_from+((self.rrna_to-self.rrna_from)/2)
     	end
     
     	def aligned_center_position
     		pos = self.center_position
     		cigar = self.rrna.cigar_line
     		mapped_pos = Bio::Mapper.new(pos,cigar).map
     		return mapped_pos
     	end
     
    end
     
    class XrefPeptideNcrna < DBConnection
      set_primary_keys :peptide_member_id, :ncrna_member_id
      belongs_to :peptide_member, :foreign_key => "peptide_member_id"
      belongs_to :ncrna_member, :foreign_key => "ncrna_member_id"
    end
    
    class XrefPeptideInfernal < DBConnection
      set_primary_keys :peptide_member_id, :infernal_hit_id
      belongs_to :peptide_member, :foreign_key => "peptide_member_id"
      belongs_to :infernal_hit, :foreign_key => "infernal_hit_id"
    end
    
    class Member < DBConnection
      
      belongs_to :organism, :foreign_key => "organism_id"
      has_one :homology_group, :foreign_key => "member_id"
      has_many :homology_members
      has_many :ncrna_members
      has_many :peptide_members
      has_many :transcripts
      has_many :infernal_hits
      has_many :blast_hits
      has_many :xref_member_pfams
      has_many :pfams, :through => :xref_member_pfams
      has_many :member_go_terms
      
      
      def display_details
      	return "#{self.stable_id}|#{self.start}>#{self.stop}|#{self.strand} (#{self.description[0..20]}...)"
      end
      
      def group_has_tag?
        answer = false
        self.homology_members.each do |hmember|
          return true if hmember.homology_group.snorna == true
        end
        return answer
      end
      
      def has_snornas? 	
      	self.ncrna_members.empty? ? answer = false : answer = true
      	return answer    	
      end
      
      def has_infernal_hits?(infernal_cutoff=30.0)      
      	self.infernal_hits.select{|h| h.score >= infernal_cutoff }.empty? ? answer = false : answer = true
      	return answer      	
      end
      	
      def infernal_hits_by_biotype(biotype)     	
				return self.infernal_hits.select{|h| h.biotype.include?("#{biotype}")}
      end	
      
	  	def get_genomic_region(range)
		
				self.organism.connect_to_genome
		
				gene = Ensembl::Core::Gene.find_by_stable_id(self.stable_id)
				slice = Ensembl::Core::Slice.new(gene.seq_region,gene.start-range,gene.stop+range,gene.strand)
		
				return slice.seq
	  
	  	end
	  
	  	def get_alignment_transcript
	  		return Ensembl::Core::Translation.find_by_stable_id("#{self.peptide_sequence.stable_id}").transcript
	  	end
	  	
      def get_longest_transcript
        
        answer = nil
        self.transcripts.each do |transcript|
          if answer.nil?
            answer = transcript if transcript.has_translation?
          else
            answer = transcript if transcript.cds_length > answer.cds_length and transcript.has_translation?
          end
        end
        return answer
        
      end
      
      def get_infernal_hits_by_score(score)
        raise "Must provide a bit score cut_off" if score.nil?
        return self.infernal_hits.select{|h| h.score >= score }.sort_by{|h| h.score }
      end
      
      def to_fasta
        self.organism.connect_to_genome
        if self.organism.database.include?("eukaryotes")
        	return EukaryoteDB::Gene.find_by_stable_id(self.stable_id).to_fasta
        else
        	return Ensembl::Core::Gene.find_by_stable_id(self.stable_id).to_fasta	
        end#return Bio::FastaFormat.open(self.organism.nt_fasta).find{|e| e.definition == self.stable_id }
      end
      
      def longest_peptide_member
        answer = nil
        self.peptide_members.each do |pm|
          answer = pm if answer.nil? or answer.sequence.length < pm.sequence.length
        end
        return answer
      end
      
      def peptide_to_fasta
        seq = Bio::FastaFormat.new(Bio::Sequence::AA.new(self.peptide_seq).to_fasta("#{self.organism.supergroup_id}|#{self.stable_id}"))
      end
      
      def nucleotide_to_fasta
        return Bio::FastaFormat.open(self.organism.nt_fasta).select{|e| e.definition == self.stable_id }[0]
        #return Bio::FastaFormat.new(Bio::Sequence::NA.new(self.nucleotide_seq).to_fasta(self.stable_id))
      end
      
      def peptide_seq
        return "#{self.longest_peptide_member.sequence}"
      end
      
      
      def get_snornas
        return self.ncrna_members.select{|s| s.biotype == "snoRNA"}
      end
      
      def get_pan_member(version=55)
        Ensembl::PanCompara::DBConnection.connect(version)
        return Ensembl::PanCompara::Member.find_by_stable_id(self.stable_id)
      end
  
      def to_svg
     		
     	  length = self.stop-self.start
     	  graph_width = 800
     	  x = 0
     	  y = 0
	      factor=graph_width.to_f/self.nucleotide_seq.length.to_f
     	  this_y = y
     	  this_x = x
     	  track_step = 80
     	  graph = SVGWriter::Graph.new(graph_width)
     		
     	  graph.add_rectangle({:x => 0, :y => 0, :fill => "#ffffff"},graph_width,400)
	      graph.add_label("#{self.stable_id} (#{self.organism.name})",{:x => x, :y => this_y+30},20)
		
	      this_y += track_step
	
	      graph.add_rectangle({ :x => x, :y => this_y, :fill => "#dcdcdc"},graph_width)
				
	      self.transcripts.each do |transcript|
	        #puts "#{transcript.stable_id} (#{transcript.start}>#{transcript.stop})"
	        this_y += track_step
	        graph.add_path({:x => x, :y => this_y+10 },graph_width)
	        transcript.exons.each do |exon|
	          #puts "\t#{exon.stable_id} (#{exon.start}>#{exon.stop})"
	          this_x = factor*(exon.start-self.start)
	          graph.add_rectangle({:x => this_x, :y => this_y, :fill => "#ffffff"},factor*(exon.stop-exon.start))
	        end				
	      end
				
	      this_y += track_step
				
	      self.infernal_hits.each do |hit|
	        next if hit.score < 20.0
	        this_x = factor*hit.member_start
	        this_width = factor*(hit.member_stop-hit.member_start)
	        graph.add_rectangle({:x => this_x, :y => this_y, :fill => "#CC0000"},this_width)
	        graph.add_label("#{hit.rfam.accession} (#{hit.score})", {:x => this_x+this_width+5, :y => this_y},10)
	      end

	      graph.height = this_y+50
				
	      return graph
				
      end
     	
    end

    class PeptideMember < DBConnection
      set_primary_key "id"
      belongs_to :member, :foreign_key => "member_id"
      has_many :xref_peptide_ncrnas
      has_many :ncrna_members, :through => :xref_peptide_ncrnas
      has_many :xref_peptide_infernals
      has_many :infernal_hits, :through => :xref_peptide_infernals
      has_many :peptide_intron_positions
      has_many :homology_members
    
      def to_fasta
        return Bio::FastaFormat.new(Bio::Sequence::AA.new(self.sequence).to_fasta("#{self.member.organism.supergroup_id}|#{self.stable_id}"))
      end
      
      def has_homology_members?
      	self.homology_members.empty? ? false : true
      end
      
      def snorna_positions
        answer = []
        self.xref_peptide_ncrnas.each{|x| answer << x.position}
        return answer
      end
      
      def intron_positions    
      	return self.peptide_intron_positions.collect{|i| i.position}
      end
      
      def _intron_positions
        answer = []
        self.member.organism.connect_to_genome
        if self.member.organism.database.include?("eukaryotes")
          transcript = EukaryoteDB::Translation.find_by_stable_id(self.stable_id).transcript
          warn "#{self.stable_id} has no transcript in the eukaryote DB" if transcript.nil?
          transcript.introns.each do |i|
            pos = transcript.peptide_position(i.previous_exon.stop-1)
            answer << pos unless pos.nil?
          end
        else
          transcript = Ensembl::Core::Translation.find_by_stable_id(self.stable_id).transcript
          warn "#{self.stable_id} has no transcript in Ensembl?" if transcript.nil?
          transcript.introns.each do |intron|
            next if intron.stop < transcript.coding_region_genomic_start or intron.start > transcript.coding_region_genomic_end
            answer << transcript.peptide_position(intron.previous_exon.stop)
          end
        end        
        return answer
      end
      
      def infernal_positions
        answer = []
        self.xref_peptide_infernals.each{|x| answer << x.position }
        return answer
      end
      
    end
    
    class PeptideIntronPosition < DBConnection
    	belongs_to :peptide_member, :foreign_key => "peptide_id"
    end
    	
    class NcrnaMember < DBConnection
      
      belongs_to :organism, :foreign_key => "organism_id"
      belongs_to :member, :foreign_key => "member_id"
      has_many :xref_peptide_ncrnas
      has_many :peptide_members, :through => :xref_peptide_ncrnas
      has_one :homology_group_position_snorna_member, :foreign_key => "ncrna_member_id"
      
      def to_fasta     
        return Bio::FastaFormat.new(Bio::Sequence::NA.new(self.seq).to_fasta(self.stable_id))        
      end
      
      def supergroup
        self.organism.supergroup
      end
      
      def rfam      
      	return HomologDB::Rfam.find_by_accession(self.external_id)     	
      end
      
      def member_start
      	
      	if self.member.strand
      		return self.start - self.member.start
      	else
      		return self.member.stop-self.stop
      	end
      
      end
      
      def rfam_accession
      	if self.rfam
      		return self.rfam.accession
      	else
      		return ""
      	end
      end
      
      def rfam_clan 
      	if self.rfam
      		return self.rfam.rfam_clan
      	else
      		return false
      	end
      end
      
      def rfam_clan_name      
      	if self.rfam_clan
      		return self.rfam_clan.name
      	else
      		return ""
      	end      	
      end
      
      def seq
        self.sequence
      end
      
      def present_at?(node)       
        ToolDB::DBConnection.connect("_54")      
        ncrna = ToolDB::EnsemblNcrna.find_by_stable_id(self.stable_id)
        node = ToolDB::DolloNode.find_by_node(node)        
        return ncrna.present_at?(node)       
      end
      
      def to_fasta       
        return Bio::FastaFormat.new(Bio::Sequence::NA.new(self.sequence).to_fasta("#{self.stable_id}"))       
      end
      
      def hostgene_start
        if self.member.strand == false
          return self.member.stop - self.stop
        else
          return self.start-self.member.start
        end
      end
      
      def hostgene_stop
        if self.member.strand == false
          return self.member.stop - self.start
        else  
          return self.stop-self.member.start
        end
      end
      
      def display_details
     		if self.rfam 
      		return "#{self.stable_id}| #{self.rfam.accession} (#{self.rfam_clan_name}) | #{self.rfam.typ} | #{self.strand} | #{self.hostgene_start} -> #{self.hostgene_stop}"
      	else
      		return "#{self.stable_id}| N/A | N/A | N/A | #{self.strand} | #{self.hostgene_start} -> #{self.hostgene_stop}"
      	end
      end
      
    end
    
    class HomologyGroup < DBConnection
      belongs_to :member, :foreign_key => "member_id"
      has_many :homology_members
      has_many :homology_group_positions 
      has_many :homology_group_intron_positions
      has_many :xref_homology_group_nodes, :foreign_key => "homology_group_id"
      
      def find_member_by_stable_id(stable_id)      
      	self.homology_members.select{|h| h.stable_id == stable_id }[0]      
      end
   
      def state_at_node?(node)
         node = HomologDB::DolloNode.find_by_name(node) unless node.kind_of?(HomologDB::DolloNode)

         state = self.xref_homology_group_nodes.select{|x| x.dollo_node_id == node.id }.shift.state

         if state == "." and node.name.include?("root") or state == "0"
            return "0"
         elsif state == "."
            self.state_at_node?(node.get_parent)
         else
            return state
         end
  
      end
      
      def deepest_node
      
         xrefs = self.xref_homology_group_nodes.sort_by{|x| x.dollo_node_id }
         xrefs.each do |x|
            return x.dollo_node if self.state_at_node?(x.dollo_node) == "1"
         end       
         return nil        
      end

      def deepest_node_name
         
         node = self.deepest_node
         node.nil? ? answer = "species-specific" : answer = node.description
         return answer
         
      end
      
      def position_deepest_node
      	nodes = []
      	self.homology_group_positions.each do |position|
      		nodes << position.deepest_node
      	end
         
         nodes.compact.empty? ? answer = nil : answer = nodes.compact.sort_by{|n| n.id}[0]
      	return answer
         
      end
      
      def position_deepest_node_name
      	node = self.position_deepest_node
      	node.nil? ? answer = "" : answer = node.description
      	return answer
      end
      
      def active_members
      	
      	return self.homology_members.select{|h| h.member.organism.active }
      
      end
      
      def has_snornas?
        answer = false
        self.homology_members.select{|h| h.member.organism.active }.each do |hmember|
          answer = true if hmember.has_snornas?
        end
        return answer
      end
      
      def has_infernals?
      	answer = false
      	self.homology_members.select{|h| h.member.organism.active}.each do |hmember|
      		answer = true if hmember.has_infernals?
      	end
      	return answer
      end
      
      def get_all_peptides
      	answer = self.homology_members.collect{|hmember| hmember.peptide_to_fasta }
      	return answer
      end
      
      
      def get_clustalw(type="peptide")    
      	aln = []      	
      	self.homology_members.select{|hm| hm.member.organism.active}.each do |hmember|     		
      		if type == "peptide" or type == "protein"
      			aln << hmember.aligned_peptide_to_fasta
      		else     		
      			aln << hmember.aligned_nucleotide_to_fasta     			
      		end      		
      	end

      	return Bio::Alignment::MultiFastaFormat.new(aln.join)
      
      end
      
      def do_align      
        bucket = []        
        self.homology_members.each do |hm|
          bucket << hm.member.to_fasta
        end       
        factory = Bio::ClustalW.new
        aln = Bio::Alignment::OriginalAlignment.new(bucket)        
        a = aln.do_align(factory)        
        return a        
      end
      
      def fetch_overlapping_alignments(ncrna)       
        start = ncrna.hostgene_start
        stop = ncrna.hostgene_stop        
        answer = []        
        self.homology_members.each do |member|          
          m_aligns = member.fetch_overlapping_alignments(start,stop)
          m_aligns.each {|a| answer << a }         
        end        
        return answer       
      end
            
    end
    
    class XrefHomologyGroupNode < DBConnection
      set_primary_keys :homology_group_id, :dollo_node_id
      belongs_to :homology_group, :foreign_key => "homology_group_id"
      belongs_to :dollo_node, :foreign_key => "dollo_node_id"
      
    end
    
    class HomologyMember < DBConnection
      
      set_primary_keys :homology_group_id, :member_id
      
      belongs_to :homology_group, :foreign_key => "homology_group_id"
      belongs_to :member, :foreign_key => "member_id"
      belongs_to :peptide_member, :foreign_key => "peptide_member_id"
      
      def has_snornas?    
      	self.get_snornas.empty? ? answer = false : answer = true
      	return answer      
      end
      
      def has_infernals?
      	peptide = self.peptide_member
      	raise "No peptide member for #{self.member_id}/#{self.homology_group_id} (#{self.member.stable_id}/#{self.member.organism.name})" if peptide.nil?
      	self.peptide_member.infernal_hits.empty? ? answer = false : answer = true
      	return answer
      end
      
      def peptide_member
        return HomologDB::PeptideMember.find(self.peptide_member_id)
      end
      
      def infernal_hits
      	return self.peptide_member.infernal_hits
      end
      
      def has_infernal_hits?
      	self.infernal_hits.empty? ? false : true
      end
      
      def display_details
      	return "#{self.id}|Member: #{self.stable_id}|#{self.member.start}->#{self.member.stop} (#{self.member.strand})|Group: #{self.homology_group.id}|Species: #{self.member.organism.name}"
      end
      
      def supergroup
        return self.member.organism.supergroup
      end
			
			def blast_hits(threshold=20.0)
				return self.member.blast_hits.select{|b| b.score >= threshold}.sort_by{|b| b.score}.reverse
			end
			
			def organism_name
				return self.member.organism.name
			end
			
			def infernal_hits_simple
				return self.pmember.infernal_hits
			end
			
			def stable_id
				return self.member.stable_id
			end
			
      def aligned_seq(type="peptide")
        
        if type == "peptide" 
          seq = self.peptide_member.sequence
          cigar = "#{self.peptide_cigar_line}"
        else  
          seq = self.member.nucleotide_seq
          cigar = "#{self.nt_cigar_line}"
        end
        
        warn "No cigar line (#{type}) found for #{self.member.stable_id}" if cigar.length == 0
        aligned_seq = ""
        
        until cigar.length == 0         
          x = cigar.slice!(/^\d+/).to_i
          char = cigar.slice!(/^[A-Z]/)         
          if char == "D"
            x.times do 
              aligned_seq += "-"
            end
          else
            x.times do 
              aligned_seq += seq.slice!(0..0)
            end
          end
          
        end
        
        return aligned_seq += seq
        
      end
      
      def aligned_nucleotide_to_fasta     
      	return Bio::FastaFormat.new(Bio::Sequence::NA.new(self.aligned_seq(type="nt")).to_fasta(self.member.stable_id))       
      end
      
      def aligned_peptide_to_fasta     
      	return Bio::FastaFormat.new(Bio::Sequence::AA.new(self.aligned_seq(type="peptide")).to_fasta("#{self.member.organism.supergroup_id}|#{self.peptide_member.stable_id}"))     
      end
      
      def map_peptide_position(pos,verbose=false)
      
      	return nil unless self.peptide_cigar_line      	
      	cig = "#{self.peptide_cigar_line.clone}"     
      	gapped_counter = 0
      	counter = 0     	
      	until cig.length == 0
      		
      		i = cig.slice!(/^\d+/).to_i
      		char = cig.slice!(/^\w/)
      		
      		i.times do 
      		
      			if char == "D"
      				gapped_counter += 1
      			else
      				gapped_counter += 1
      				counter += 1
      				return gapped_counter if counter == pos
      			end
    	
      		end
      	
      	end
      	
      	puts "Didn't meet the end... #{gapped_counter},#{counter} where position was #{pos}" if verbose
      	return gapped_counter
      
      end
      
      def peptide_to_fasta        
        return self.member.peptide_to_fasta        
      end
      
      def nucleotide_to_fasta      
      	return self.member.nucleotide_to_fasta     	
      end
      
      def get_rnaz(hit)        
        seq = self.member.sequence[hit.target_start-50..hit.target_end+50]
        rnaz = Bio::Vienna::RNAz.new()        
      end
  
      def get_snornas
        return self.member.ncrna_members.select{|n| n.biotype=="snoRNA"}
      end
  
      def get_aligned_position(pos,type="nt")      	
      	if type == "nt"
      		cig = self.nt_cigar_line.clone    		
      	else
      		cig = self.peptide_cigar_line.clone
      	end      	
      	curr_pos = 0
      	aligned_pos = 0      
      	while curr_pos < pos
      		x = cig.slice!(/^[0-9]*/).to_i
      		char = cig.slice!(/^[A-Z]/).strip	
      		puts "#{x} / #{char} (#{cig[0..5]})" if x == 0
      		raise if x == 0
      		if char == "M"
      			x.times do
      				curr_pos += 1
      				aligned_pos += 1
      				return aligned_pos if curr_pos == pos
      			end
      		else
      			aligned_pos += x
      		end
      	end      	
      	return aligned_pos      	
      end      
    end
    
    class InfernalHit < DBConnection
    
    	belongs_to :member, :foreign_key => "member_id"
    	belongs_to :organism, :foreign_key => "organism_id"
    	has_many :xref_peptide_infernals
    	has_many :peptide_members, :through => :xref_peptide_infernals
    	has_one :homology_group_position_infernal_member
    	
    	def genomic_sequence(offset=0)
        warn "#{self.id} infernal hit has no member?" if self.member.nil?
        return "" if self.member.nil?
    	  self.member.organism.connect_to_genome
    	  seq = self.member.to_fasta.naseq   	  
    	  return "" if seq.nil?
      	if self.genomic_strand == self.member.strand
      	  return seq[self.member_start-offset..self.member_stop-1+offset]
      	else
      	  seq = seq[self.member_start-1-offset..self.member_stop-1+offset]
      	  return "" if seq.nil?
      	  return Bio::Sequence::NA.new(seq).complement
      	end
    	end
    	    	  
    	def start 
    	  #puts "#{self.genomic_strand}/#{self.member.strand}"
    	  if self.genomic_strand == false
    	    if self.strand == 1
    	      return self.member.stop-self.member_stop
    	    else
    	      return self.member.start+self.member_start-1
    	    end
    	  end
    	  self.genomic_strand == self.member.strand ? self.member.start+self.member_start-1 : self.member.start+((self.member.stop-self.member.start)-self.member_stop-1)
    	end
    	
    	def stop 
    	  if self.genomic_strand == false
    	    if self.strand == 1
    	      return self.member.stop-self.member_start
    	    else
    	      return self.member.start+self.member_stop-1
    	    end    	    
    	  end
    	  self.genomic_strand == self.member.strand ? self.member.start+self.member_stop-1 : self.member.start+((self.member.stop-self.member.start)-self.member_start-1) 
    	end
    	
    	def genomic_strand
    	  if self.member.strand == false 
    	    self.strand == 1 ? answer = false : answer = true
    	  else
    	    self.strand == 1 ? answer = true : answer = false
    	  end
    	  return answer
    	end
    	
		  def display_details
			  return "#{self.rfam_id} (#{self.rfam_clan_name}) | #{self.rfam_type} | #{self.score} | #{self.start}->#{self.stop} | #{self.peptide_position} |#{self.strand} (#{self.member.strand})"
		  end
		
		  def rfam
		    RfamDB::DBConnection.connect("10")
		    return RfamDB::Rfam.find_by_rfam_acc(self.rfam_id)
      end
    
		  def rfam_clan
		  	return nil unless self.rfam
		  	if self.rfam.clan
			  	return self.rfam.clan
			  else
			  	return nil
			  end
		  end
		
		  def rfam_type
		  	return nil unless self.rfam
		    self.rfam.typ
      end

		  def rfam_clan_name
			  if self.rfam_clan 
				  return self.rfam_clan.name
			  else
				  return ""
			  end	
		  end
		
		  def to_fasta
			  return Bio::FastaFormat.new(Bio::Sequence::NA.new(self.genomic_sequence).to_fasta("#{self.member.stable_id}|INFERAL|#{self.rfam_id}|#{self.score}"))
		  end
			
	  	
    end
    
    # = DESCRIPTION
 		# Groups snoRNAs and Infernal loci belonging to the Homology Group
 		class HomologyGroupPosition < DBConnection
 			belongs_to :homology_group, :foreign_key => "homology_group_id"
 			has_many :homology_group_position_snorna_members
      has_many :ncrna_members, :through => :homology_group_position_snorna_members
 			has_many :homology_group_position_infernal_members
      has_many :infernal_hits, :through => :homology_group_position_infernal_members
 			has_many :xref_position_nodes
 			
 			def self.fetch_all_by_rfam_acc(rfam_acc) 				
 				answer = []
 				HomologDB::HomologyGroupPosition.find(:all).each do |pos|
 					answer << pos if pos.homology_group_position_snorna_members.select{|sm| sm.ncrna_member.rfam_acc == rfam_acc}.nitems > 0 or pos.homology_group_position_infernal_members.select{|im| im.infernal_hit.rfam_id == rfam_acc}.nitems > 0
 				end
				return answer 			
 			end
 			
 			def state_at_node?(node)
 				
 				node = HomologDB::DolloNode.find_by_name(node) unless node.kind_of?(HomologDB::DolloNode)
 				
 				state = self.xref_position_nodes.select{|x| x.dollo_node_id == node.id }.shift.state
 				
 				if state == "." and node.name.include?("root")
 					return "0"
 				elsif state == "0"
 					return "0"
 				elsif state == "."
 					self.state_at_node?(node.get_parent)
 				else
 					return state
 				end
 			
 			end
 			
 			def deepest_node
 				
 				nodes = self.xref_position_nodes.collect{|x| x.dollo_node}.sort_by{|n| n.id}
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
 			
 		end
 		
 		# = DESCRIPTION
 		# snoRNAs that belong to a particular aligned locus
 		class HomologyGroupPositionSnornaMember < DBConnection
 			belongs_to :homology_group_position, :foreign_key => "homology_group_position_id"
 			belongs_to :peptide_member, :foreign_key => "peptide_member_id"
 			belongs_to :ncrna_member, :foreign_key => "ncrna_member_id"
 		end
    
    class HomologyGroupPositionInfernalMember < DBConnection
 			belongs_to :homology_group_position, :foreign_key => "homology_group_position_id"
 			belongs_to :peptide_member, :foreign_key => "peptide_member_id"
 			belongs_to :infernal_hit, :foreign_key => "infernal_hit_id"
 		end
 		
 		class HomologyGroupIntronPosition < DBConnection
 			belongs_to :homology_group, :foreign_key => "homology_group_id"
 			has_many :homology_group_intron_position_members
 			has_many :xref_intron_position_nodes
 			
 			def state_at_node?(node) 				
 				node = HomologDB::DolloNode.find_by_name(node) unless node.kind_of?(HomologDB::DolloNode)            
 				state = self.xref_intron_position_nodes.select{|x| x.dollo_node_id == node.id }.shift.state 				
 				if state == "." and node.name.include?("root")
 					return "0"
 				elsif state == "."
 					self.state_at_node?(node.get_parent)
 				else
 					return state
 				end 			
 			end
         
         def deepest_node
            nodes = self.xref_intron_position_nodes.collect{|x| x.dollo_node}.sort_by{|n| n.id}
            nodes.each do |node|
               return node if self.state_at_node?(node) == "1"
            end
            return nil
         end

         def deepest_node_name
            deepest_node = self.deepest_node
            if deepest_node.nil?
               return ""
            else
               return deepest_node.description
            end
         end     
              			
 		end
 		
 		class HomologyGroupIntronPositionMember < DBConnection
 			belongs_to :homology_group_intron_position, :foreign_key => "homology_group_intron_position_id"
 			belongs_to :peptide_member, :foreign_key => "peptide_member_id"
 		end
 		
 		class XrefIntronPositionNode < DBConnection
 			belongs_to :homology_group_intron_position, :foreign_key => "homology_group_intron_position_id"
 			belongs_to :dollo_node, :foreign_key => "dollo_node_id"
 		end
 		
    # = DESCRIPTION
    # Phylogenetic tree
    class DolloTree < DBConnection   	
    	has_many :dollo_nodes
    	has_many :dollo_positions
    	belongs_to :supergroup, :foreign_key => "supergroup_id"
      
      def tree
         return Bio::Newick.new(self.tree_string).tree
      end
      
      def output_phylip
				tree_string = self.tree_string
				tree_organisms = self.tree.leaves.collect{|t| t.name.strip }
				tree_organisms.each do |t_o|
					phylip_name = t_o.phylip_name
					t_o = t_o.gsub(/\s/ , '_')
					tree_string.gsub!(/#{t_o}/, "#{phylip_name}")
				end
				return tree_string
      end
      
    end
    
    class DolloNode < DBConnection
    	belongs_to :dollo_tree, :foreign_key => "dollo_tree_id"
    	has_many :xref_position_nodes
    	has_many :xref_intron_position_nodes
    	has_many :xref_rrna_locus_nodes
      has_many :xref_homology_group_nodes
    	
    	def get_parent

    		if self.child_of.nil?
    			return self
    		else	
    			return DolloNode.find(self.child_of)
    		end
    	end
    	
    	def get_children
    		return DolloNode.find_all_by_child_of(self.id)
    	end
    	
    end
    
    
    class XrefPositionNode < DBConnection
    	belongs_to :homology_group_position, :foreign_key => "homology_group_position_id"
    	belongs_to :dollo_node, :foreign_key => "dollo_node_id"
    end	
    
    class Supergroup < DBConnection
     has_many :organisms   
     has_one :dollo_tree
    end
     
    class XrefOrganismSupergroup < DBConnection
    
      belongs_to :organism, :foreign_key => "organism_id"
      belongs_to :supergroup, :foreign_key => "supergroup_id"
    
    end
    
    class MemberGoTerm < DBConnection
    	belongs_to :member, :foreign_key => "member_id"
    end
  end
  
end