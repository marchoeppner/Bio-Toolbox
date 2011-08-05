Eukaryote_DB_ADAPTER = 'postgresql' 
Eukaryote_DATABASE  = "eukaryotes" 
Eukaryote_DB_HOST = 'localhost' 
Eukaryote_DB_USERNAME =  'tools' 
Eukaryote_DB_PASSWORD = 'analysis'

require 'zlib'

module Toolbox 
  
  module EukaryoteDB
    
    include ActiveRecord
    
    class DBConnection < ActiveRecord::Base
      self.abstract_class = true
    	self.pluralize_table_names = false
    	
      def self.connect(version="1")

        establish_connection(
                              :adapter => Eukaryote_DB_ADAPTER,
                              :host => Eukaryote_DB_HOST,
                              :database => "#{Eukaryote_DATABASE}#{version}",
                              :username => Eukaryote_DB_USERNAME,
                              :password => Eukaryote_DB_PASSWORD
                              #:port => port
                            )
      end
      
    end
      
    class Supergroup < DBConnection
     	has_many :organisms
     end
      
    class Organism < DBConnection	
     	belongs_to :supergroup, :foreign_key => "supergroup_id"
     	has_many :genes
     	has_many :dnas, :order => "accession ASC"
     	has_many :inparanoid_clusters
     	
     	def snake_case
     		return name.gsub(/\s/, '_').downcase
     	end
     	
     	def accessions
     		accs = self.dnas.select{|d| d.accession.length > 1}.collect{|d| d.accession}.sort
     		if accs.nitems > 1
     		  accs = accs.compact
	     		return "#{accs[0]}-#{accs[-1]}"
  			else
  				return accs[0]
  			end
     	end	
     	
     	def ncrnas 
     	  answer = []
     	  self.dnas.each do |dna|
     	    dna.infernal_hits.each{|i| answer << i}
     	  end
     	  return answer
     	end
     	
     end
     
    class Dna < DBConnection
     	belongs_to :organism, :foreign_key => "organism_id"
     	has_many :genes
     	has_many :infernal_hits
     	
     	def slice
     	  return EukaryoteDB::Slice.new(self.id,0,self.seq.length,1)
     	end
     	
     	# = DESCRIPTION
     	# Retrieves unzipped sequence
     	def seq
     		ff = Bio::FastaFormat.open("/home/marc/databases/eukaryotes/genome_files/#{self.organism.snake_case}/#{self.organism.snake_case}.dna.fasta")
     		return ff.find{|e| e.definition == self.accession }.naseq
     	end
     	
     	def to_fasta
     	  return Bio::FastaFormat.new(Bio::Sequence::NA.new(self.seq).to_fasta(self.accession))
     	end
     	
     end
      
    class Gene < DBConnection
     	has_one :nt_seq
     	has_many :transcripts
     	belongs_to :organism, :foreign_key => "organism_id"
     	belongs_to :dna, :foreign_key => "dna_id"
     	has_many :homologs
     	has_many :ensembl_homologs
      has_many :xref_gene_pfams, :order => "score"
      has_many :pfams, :through => :xref_gene_pfams
      
	    def display_details
	      return "#{self.stable_id}|#{self.start}>#{self.stop}|#{self.strand}|#{self.description[0..30]}..."
	    end
      
      def slice
        return EukaryoteDB::Slice.new(self.dna_id,self.start,self.stop,self.strand)
      end
        
      def genomic_context
      
      	puts "Self: #{self.inspect}"
      	
      	downstream = EukaryoteDB::Gene.find(:all, :conditions => ["dna_id = ? AND start < self.start"], :order => "start DESC")[0]
      	upstream = EukaryoteDB::Gene.find(:all, :conditions => ["dna_id = ? AND start > self.start"], :order => "start ASC")[0]
      	
      	if upstream and downstream
      		down_dist = self.start - downstream.stop
      		updist = upstream.start -self.stop
      	elsif upstream
      		down_dist = "NA"
      		updist = upstream.start -self.stop
      	elsif downstream
      		down_dist = self.start - downstream.stop
      		updist = "NA"
      	end
      	
      	downstream ? downstream_string = "#{downstream.stable_id}[#{downstream.strand}]" : downstream_string = "NA[NA]"
      	upstream ? upstream_string = "#{upstream.stable_id}[#{upstream.strand}]" : upstream_string = "NA[NA]"
      	
      	return "#{downstream_string}|#{downdist}>|#{self.stable_id}[#{self.strand}]|<#{updist}|#{upstream_string}"
      	
      end
      
      def fetch_hostgene
        if self.is_intronic?
          return self.slice.genes(true)[0]
        else
          return nil
        end
      end
        
      def is_intronic?
        return false if self.biotype == "protein_coding"
        slice = self.slice
        slice.genes(true).select{|g| g != self }.empty? ? answer = false : answer = true
        return answer
      end
        
	    def infernal_hits_by_threshold(threshold=15.0)
	      return self.infernal_hits.select{|h| h.score >= threshold}	  
	    end
	    
	    def infernal_ncrnas(infernal_version)
	      return EukaryoteDB::InfernalHit.find(:all, :conditions => ["dna_id = ? and start > ? and stop < ? and infernal_version = ?", self.dna_id,self.start,self.stop,infernal_version])
	    end
	    
	    def infernal_ncrnas_for_homologdb
	      ncrnas = self.infernal_ncrnas
	      answer = []
	      ncrnas.each do |ncrna|
	        this_ncrna = {}
	        this_ncrna[:gene_start] = ncrna.start-self.start
	        this_ncrna[:gene_stop] = ncrna.stop-self.start
	        this_ncrna[:strand] = ncrna.strand
	        this_ncrna[:score] = ncrna.score
	        this_ncrna[:rfam_id] = ncrna.rfam_acc
	        this_ncrna[:query_seq] = ncrna.query_seq
	        this_ncrna[:target_seq] = ncrna.target_seq
	        this_ncrna[:midline] = ncrna.midline
	        answer << this_ncrna
	      end
	      return answer
	    end
	        
	    def infernal_snornas
	      answer = []
	      hits = self.infernal_ncrnas
	      return answer if hits.empty?
     	  RfamDB::DBConnection.connect("_91")
     	  rfams = RfamDB::Rfam.find_all_snornas.collect{|s| s.rfam_acc}
     	  hits.each do |hit|
     	    answer << hit if rfams.include?(hit.rfam_acc)
     	  end
     	  return answer
     	end
     	
     	def seq
     	  return self.sequence
     	end
     	
     	def sequence
     		seq =self.dna.seq
     		return nil if seq.nil? 
     		seq = Bio::Sequence::NA.new(seq)
     		if self.strand == 1
     			return seq.subseq(self.start,self.stop)
     		else
     			return seq.subseq(self.start,self.stop).complement
     		end
     	end
     	
     	def snornas
     		return EukaryoteDB::Gene.find_all_by_biotype_and_organism_id_and_dna_id("snoRNA",self.organism_id,self.dna_id).select{|s| s.start >= self.start and s.stop <= self.stop }
     	end
     	
     	def to_fasta
     		return Bio::FastaFormat.new(Bio::Sequence::NA.new(self.sequence).to_fasta(self.stable_id))
     	end
     	
     	def get_longest_transcript

     		answer = nil
     		self.transcripts.each do |transcript|
     			if answer.nil?
     				 answer = transcript
     			elsif answer.product.length < transcript.product.length
     				answer = transcript
     			end
     		end
     		return answer
     	end

      def peptide_to_fasta
     		transcript = self.get_longest_transcript
     		if transcript.nil?
     			return ""
     		else
     			return Bio::FastaFormat.new(Bio::Sequence::AA.new(transcript.translation_seq).to_fasta(transcript.translation.stable_id))
     		end
     	end
     	
    end
      
    class Transcript < DBConnection
     	belongs_to :gene, :foreign_key => "gene_id"
     	has_many :xref_transcript_exons
     	has_many :exons, :through => :xref_transcript_exons, :order => "start"
     	has_one :translation
     	
     	def introns

			  answer = []
			  exons = self.exons.clone
			  this_exon = exons.shift

			  exons.each_with_index do |exon,index|
				  answer << Intron.new(this_exon,exon,self)
				  this_exon = exon
			  end

			  return answer

		  end
		
		  def protein_seq 
		  	if self.translation.nil?
		  		warn "#{self.stable_id} has no translation product"
		  		return ""
		  	end 
			  return self.translation_seq
		  end
		
		  def seq
		    gene_seq = self.gene.seq
		    return gene_seq[self.start-self.gene.start..self.stop-self.gene.start]
		  end
		    
		  def translation_seq
     		return self.translation.sequence
     	end
     	
     	def product
     		
     		seq = Bio::Sequence::NA.new(self.gene.dna.seq)
     		cds = ""
     		self.exons.each do |exon|
					cds += seq.subseq(exon.start,exon.stop)   			
     		end
     		
     		if self.strand == 1
     			return Bio::Sequence::NA.new(cds).translate
     		else
     			return Bio::Sequence::NA.new(cds).complement.translate
     		end
     	
     	end
     	
     	def product_stable_id
     		return self.translation.stable_id
     	end

		  def transcript_to_peptide(position)

  			pos = 0

  			if self.strand == 1

  				self.exons.each do |exon|
  					if exon.stop < position
  						pos += exon.length
  					elsif exon.start < position and exon.stop >= position
  						pos += position-exon.start
  					end
  				end

  				return pos
  			else
  				self.exons.sort_by{|e| e.start}.reverse.each do |exon|
  					if exon.start >= position
  						pos += exon.length
            end
  				end
  			end
			
  			return pos
  		end

		  def peptide_position(position,verbose=false)
  			answer = nil
  			self.introns.each do |intron|
  				if position >= intron.start and position <= intron.stop
  					puts "position (#{position}) is located in this intron (#{intron.start}>#{intron.stop})" if verbose
  					answer = intron.previous_exon 
  				end
  			end
  			if answer.nil?
  				introns = self.introns.collect{|i| "#{i.start}>#{i.stop}" }.join(",") if verbose
  				puts "position (#{position}) not in any of the introns of this gene (#{introns})" if verbose
  				return nil
  			else
  				return transcript_to_peptide(answer.stop)
  			end
  		end

     	
    end
     
    class XrefTranscriptExon < DBConnection
     	belongs_to :transcript, :foreign_key => "transcript_id"
     	belongs_to :exon, :foreign_key => "exon_id"
     end
     
    class Exon < DBConnection
     	 has_many :xref_transcript_exons
		   has_many :transcripts, :through => :xref_transcript_exons
     	
     	 def details
     		 return "#{self.stable_id}|#{self.start}|#{self.stop}"
     	 end

		   def seq
			   nt_seq = Bio::Sequence::NA.new(self.transcripts.first.gene.dna.seq).subseq(self.start,self.stop)
			   if self.strand == -1
				   return nt_seq.complement
			   else
				   return nt_seq
			   end
		   end

		   def length
			   return self.stop-self.start
		   end	  

     end
     
	  class Intron
		
  		attr_reader :start, :stop, :exon_1, :exon_2, :transcript
  		def initialize(exon_1,exon_2,transcript)
  			@exon_1 = exon_1
  			@exon_2 = exon_2
        @transcript = transcript
  			@start  = exon_1.stop+1
  			@stop = exon_2.start-1
  		end

  		def previous_exon
  			return self.exon_1
  		end
		
  		def next_exon
  			return self.exon_2
  		end

  	end
	
    class Translation < DBConnection
    	belongs_to :transcript, :foreign_key => "transcript_id"
    	has_many :xref_translation_pfams
    	has_many :pfams, :through => :xref_translation_pfams
  	
  	  def seq
  	  	ff = Bio::FastaFormat.open("/home/marc/databases/eukaryotes/genome_files/#{self.transcript.gene.organism.snake_case}/#{self.transcript.gene.organism.snake_case}.aa.fasta")
  	  	s = ff.find{|e| e.definition == self.stable_id}.naseq.upcase
  	  	return s
  	  end
  	  
  	  def sequence
  	  	self.seq
  	  end
  	  
    	def seq_start
    		return 0
    	end
  	
    	def seq_end
    		return self.end_exon.length
    	end
  	
    	def start_exon
    	   if self.transcript.strand == 1
    		  return self.transcript.exons[0]
    		 elsif self.transcript.strand == -1
    		  return self.transcript.exons[-1]
    		 else
    		  raise "something's wrong with the transcript strand..."
    		 end
    	end
  	
    	def end_exon
    	  if self.transcript.strand == 1
    		  return self.transcript.exons[-1]
    		else
    		  return self.transcript.exons[0]
    		end
    	end
    end
    	
    class InparanoidCluster < DBConnection
      has_many :inparanoid_members
      belongs_to :organism, :foreign_key => "organism_id"
    
      def self.fetch_all_by_gene_stable_id(stable_id)
        return EukaryoteDB::InparanoidMember.find_all_by_gene_stable_id(stable_id).collect{|im| im.inparanoid_cluster}
      end
    
      def self.fetch_all_homologues_by_gene_stable_id(stable_id,cutoff=0.5)
        answer = []
        self.fetch_all_by_gene_stable_id(stable_id).each do |cluster|
          cluster.inparanoid_members.select{|im| im.gene_stable_id.include?("ENSG") == false }.each{|im| answer << im.gene}
        end
        return answer.compact
      end
    
    end
  
    class InparanoidMember < DBConnection
      belongs_to :inparanoid_cluster, :foreign_key => "inparanoid_cluster_id"
      def gene
        return EukaryoteDB::Gene.find_by_stable_id(self.gene_stable_id)
      end
      
      def protein 
      	transcript = EukaryoteDB::Transcript.find_by_stable_id(self.protein_stable_id)    	
        return transcript.translation
      end
    end
		
		class HumanHomolog < DBConnection
			belongs_to :organism, :foreign_key => "homolog_organism_id"
			belongs_to :gene, :foreign_key => "homolog_gene_stable_id"
			belongs_to :transcript, :foreign_key => "homolog_transcript_stable_id"
			
			def gene
				return EukaryoteDB::Gene.find_by_stable_id(self.homolog_gene_stable_id)
			end
			
			def protein
				transcript = EukaryoteDB::Transcript.find_by_stable_id(self.homolog_transcript_stable_id)
				return transcript.translation
			end
			
		end
		
  	class XrefGenePfam < DBConnection
  	  belongs_to :gene, :foreign_key => "gene_id"
  	  belongs_to :pfam, :foreign_key => "pfam_id"
  	end
	
  	class XrefTranslationPfam < DBConnection
  	  belongs_to :translation, :foreign_key => "translation_id"
  	  belongs_to :pfam, :foreign_key => "pfam_id"
  	end
	
  	class Pfam < DBConnection
  	  has_many :xref_translation_pfams
  	  has_many :xref_gene_pfams
  	end
  
    class InfernalHit < DBConnection
    	belongs_to :dna, :foreign_key => "dna_id"
  	
    	def display_details
    		return "#{self.rfam_acc}|#{self.score}|#{self.start}>#{self.stop}|#{self.strand}"
    	end
  	
  	  def slice
  	    return EukaryoteDB::Slice.new(self.dna_id,self.start,self.stop,1)
  	  end
  	  
  	  def is_intronic?
  	    self.slice.genes(true).empty? ? answer = false : answer = true
  	    return answer
  	  end
  	  
    	def genomic_seq
    	  if self.strand == 1
    	    return self.dna.seq.clone[self.start..self.stop]
    	  else
    	    return Bio::Sequence::NA.new(self.dna.seq[self.start..self.stop]).complement
    	  end
    	end
    	
    end
    
    class Slice
      
      attr_reader :dna_id, :start, :stop, :strand
      def initialize(dna_id,start,stop,strand)
        @dna_id, @start, @stop, @strand = dna_id,start,stop,strand
      end
      
      def dna
        return EukaryoteDB::Dna.find(self.dna_id)
      end
      
      def genes(overlap=false)
        if overlap
          return EukaryoteDB::Gene.find(:all, :conditions => ["dna_id = ? AND start < ? AND stop > ?", self.dna_id, self.stop,self.start])
        else
  	      return EukaryoteDB::Gene.find(:all, :conditions => ["dna_id = ? AND start > ? and stop < ?", self.dna_id, self.start,self.stop])
        end
      end
      
      def seq
        if self.strand == 1
          return self.dna.seq[self.start-1..self.stop-1]
        else
          return Bio::Sequence::NA.new(self.dna.seq[self.start-1..self.stop-1]).complement
        end
      end
      
    end
  
  end   
end
