Genome_DB_ADAPTER = 'postgresql'
Genome_DB_HOST = 'localhost'

Genome_DB_USERNAME = 'tools'
Genome_DB_PASSWORD = 'analysis'

module Toolbox 
  module GenomeDB
    
    require 'composite_primary_keys'
    
    include ActiveRecord
    
    class DBConnection < ActiveRecord::Base
      self.abstract_class = true
    	self.pluralize_table_names = false
    	
      def self.connect(species)

        establish_connection(
                              :adapter => Genome_DB_ADAPTER,
                              :host => Genome_DB_HOST,
                              :database => species,
                              :username => Genome_DB_USERNAME,
                              :password => Genome_DB_PASSWORD
                              #:port => port
                            )
      end
      
    end
      
    class Dna < DBConnection
      has_many :genes
      has_many :transcripts, :through => :genes
      def sequence
        return Bio::FastaFormat.open("/home/marc/databases/genomes/fasta/#{self.connection.current_database}.fasta").find{|e| e.definition == self.accession }.naseq
      end
      
    end

    class Gene < DBConnection
      set_primary_key 'id'
      belongs_to :dna, :foreign_key => "dna_id"
      has_many :transcripts
      
      def slice
        return Slice.new(self.dna,self.start,self.stop,self.strand)
      end
      
      def seq
        return self.slice.seq
      end
      
      def longest_transcript        
         transcripts = self.transcripts
         if transcripts.length == 1
            return transcripts.shift
         else
            answer = nil
            self.transcripts.each do |t|
               if answer.nil?
                  answer = t 
               else
                  answer = t if answer.cds.length < t.cds.length
               end
            end
            return answer
         end   
            
      end
      
    end

    class Transcript < DBConnection
      belongs_to :gene, :foreign_key => "gene_id"
      has_one :dna, :through => :gene
      has_many :xref_transcript_exons
      has_many :exons, :through => :xref_transcript_exons, :order => "start ASC"
      has_one :translation
      
      def slice
        return Slice.new(self.gene.dna,self.start,self.stop,self.strand)
      end

      def seq
        return self.slice.seq
      end
      
      def translation_product
             
        if self.strand == 1
          return Bio::Sequence::NA.new(self.cds).translate(self.translation.codon_start,self.translation.transl_table)
        else
          return Bio::Sequence::NA.new(self.cds).complement.translate(self.translation.codon_start,self.translation.transl_table)
        end
       
      end
      
      def cds
        seq = self.gene.dna.sequence
        cds = ""
        coding_start,coding_end = self.translation.coding_start,self.translation.coding_end
        self.exons.each do |exon|
            if coding_start > exon.start and coding_start < exon.stop and coding_end > exon.start and coding_end < exon.stop
               cds += seq.subseq(coding_start,coding_end)
            elsif coding_start > exon.start and coding_start < exon.stop
               cds += seq.subseq(coding_start,exon.stop)
            elsif coding_end > exon.start and coding_end < exon.stop
               cds += seq.subseq(exon.start,coding_end)
            else
               cds += seq.subseq(exon.start,exon.stop)
            end
        end
        return cds
      end
    
    end

    class XrefTranscriptExon < DBConnection
      set_primary_keys :exon_id,:transcript_id
      belongs_to :transcript, :foreign_key => "transcript_id"
      belongs_to :exon, :foreign_key => "exon_id"
    end

    class Exon < DBConnection
      has_many :xref_transcript_exons
      has_many :transcripts, :through => :xref_transcript_exons
    end

    class Translation < DBConnection
      belongs_to :transcript, :foreign_key => "transcript_id"
      
      def seq
        
        self.transcript.translation_product
      
      end
      
    end 
    
    class Slice
      
      attr_reader :dna, :start, :stop, :strand
      def initialize(dna,start,stop,strand)
        @dna,@start,@stop,@strand = dna,start,stop,strand
      end
      
      def seq
        answer = nil
        self.strand == 1 ? answer = self.dna.sequence.subseq(start,stop) : answer = self.dna.sequence.subseq(start,stop).complement
        return answer
      end
      
    end

  end # GenomeDB

end # Toolbox
