  
  
  def Math.log2( x )
    Math.log( x ) / Math.log( 2 ) 
  end
  
# = DESCRIPTION
# This is an internal class to reconstruct aligned
# sequences from cigar lines.
# Can take a range argument (start,stop) to return a sub-region.
# The indent argument ignores indent gaps for 2X 
# coverage genomes (several genomic_align objects for one 
# genomic align block are normally padded via "X")
class AlignSeq
  
  attr_reader :seq, :cigar_line, :start, :stop, :noindent
  
  def initialize(seq,cigar_line,start=0,stop=nil,noindent=false)
    @seq = seq
    @cigar_line = cigar_line
    @start = start
    @stop = stop
    @noindent = noindent
  end
  
  def align
    
    cigar_line = "#{self.cigar_line}"   # The encoded gap structure
    seq = "#{self.seq}"                 # the ungapped sequence
    gapped_seq = ""                     # the reconstructed, gapped sequence
    
    until cigar_line.length == 0
      
      x = cigar_line.slice!(/^[0-9]*/)
      x = 1 if x.nil? or x.to_i == 0 or x.length == 0
      char = cigar_line.slice!(/^[A-Z]/)
      x = x.to_i
      
      if char == "M"
        gapped_seq = gapped_seq + seq.slice!(0..x-1)
      else
        gapped_seq = gapped_seq + ("-"*x) unless self.noindent==true and char == "X"
      end  
      
    end
    
    start = self.start
    stop = self.stop unless self.stop.nil?
    stop = gapped_seq.length if self.stop.nil?
    
    return gapped_seq[start..stop]
    
  end
     
end

class String

  def phylip_name
    string = self
    answer = ""
    if string.include?("_")
      return string.capitalize[0..9]
    else
      string.split(" ").each do |s|
        answer += "_#{s[0..4]}"
      end
      return answer[0..9].gsub(/^_/, '')
    end
    
  end	

	def strip_stockholm
		return self.gsub(/^\#=G[A-Z]\s[A-Z]*(\s+)*/, '')
	end
	
  def mask(start,stop,char=".")
    return char * (stop-start+1) + self[stop+1..-1] if start == 0 # if the mask starts at the beginning, do not slice off the first nt
    return self[0..start-1] + char * (stop-start+1) + self[stop+1..-1]
  end
  
  
  def transform_date
  	
  	elements = self.split("/")
  	
  	return "#{elements[2]}-#{elements[0]}-#{elements[1]}"
  	
  end
  
  def sql_compliant
    return self.gsub(/\'/, '').gsub(/\"/, '').gsub(/\/*/, '').gsub(/\´/, '').gsub(/;/, ",").gsub(/\n/, '')
  end
  
	def to_snake_case
		return self.gsub(/\s+/, '_').downcase
	end
		
  def binary_to_ascii
    
    return self.unpack('d')
  end
  
  def to_width(width)
    answer = []
    this_seq = self
    until this_seq.length <= width
      answer << this_seq.slice!(0..width-1)
    end
    answer << this_seq
    return answer.join("\n")
  end
  
  def softmasked_seq_to_html(window_size)
    seq = self.to_s
    exon = "red"
    intron = "grey"

    answer = []

    answer << "<table width='#{window_size}px' cellspacing='1px' height='10px'>\n<tr>\n"
    answer << "<td width='#{window_size/2-1}px'></td><td width='2px' bgcolor='blue'></td><td width='#{window_size/2-1}px'></td>\n</tr>\n<tr>\n</table>"
    answer << "<table width='#{window_size}px' cellspacing='0px' height='10px'>\n<tr>\n"

    until seq.length == 0
      if seq[0..0].match(/[A-Z]/)
        this_seq = seq.slice!(/^[A-Z]*/)
        answer << "<td width='#{this_seq.length}px' bgcolor='#{exon}'></td>\n"
      else
        this_seq = seq.slice!(/^[a-z]*/)
        answer << "<td width='#{this_seq.length}px' bgcolor='#{intron}'></td>\n"
      end
    end
    answer << "\n</tr>\n</table>"

    return answer.join

  end
  
  
end
 
def mask_seq(seq,coordinates, char=".")
  last = 0
  coordinates.each do |start,stop|
    next if stop == 0
    seq = seq.mask(last,start)
    last = stop
  end
  seq = seq[0..last-1] + char*(seq.length-last) unless seq.length-last <= 0
  return seq
end

def parse_baseml(output)
  
  answer = {}
  
  result = output.split("\n")
  result.each do |line|
    
    answer["homogeneity"] = line.gsub(/Homogeneity\sstatistic\:\s/, '') if line.include?("Homogeneity")
    answer["constant_sites"] = line.gsub(/\#\sconstant\ssites:\s/, '') if line.include?("constant sites")
    
  end
  
  return answer
    
end
    


class String
  
  def mask(last,start,char=".")
    if last == 0 and start > 0
      return char*start + self[start..-1]
    elsif last > 0  
      return self[0..last-1] + char*(start-last) + self[start..-1] 
    else
      return self
    end
  end
  
end
  
class Array
  
  def coordinates_overlap?(coordinates)
    
    self.sort.each do |start,stop|
      
      coordinates.sort.each do |c_start,c_stop|
        
        if start < c_stop and stop > c_start
          return true
        end
        
      end
      
      return false
      
    end
    
  end
  
end

module Bio

  module Phylip

    class Dollop

      attr_accessor :infile, :outfile, :intreefile, :ancseq, :filename
      def initialize(infile,intreefile,ancseq=true)
        @infile = infile
        @filename = $$
	      @outfile = "#{@filename}.fdollop"
	      @ancseq = ancseq
	      @intreefile = intreefile
      end

      def run
	      system("fdollop -infile #{self.infile} -intreefile #{self.intreefile} -outfile /tmp/#{self.outfile} -ancseq -progress N")
        # Parse the output
	      valid = false
	      node1 = nil
	      node2 = nil
	      dollo = {}
				FileUtils.cd("/tmp/") { |dir|
	      	IO.readlines("#{self.outfile}").compact.each do |line|
	      		valid = true if line.include?("same as in the node below")
	      		next unless valid
	        	if line.include?("yes") or line.include?("\sno\s")  # a new nodecat 
	        		node1 = line.slice!(/^(\s*)[A-Za-z0-9_]+/).strip
	        		node2 = line.slice!(/^(\s*)[A-Za-z0-9_]+/).strip
	        		steps = line.slice!(/^(\s*)[A-Za-z0-9_]+/).strip
	        		states = line.strip.gsub(/\s/, '')
	        	  dollo["#{node1}_#{node2}"] = states.split(//).compact
	        	elsif node1 # node1 must be specified already, or skip
	        		line.strip.gsub(/\s/, '').split(//).compact.each {|e| dollo["#{node1}_#{node2}"] << e}
	        	end
	      	end
					cleanup
				}
	      return dollo
	
      end
      
      def cleanup
      	system("rm /tmp/#{self.filename}*")
      end
    end

  end

  class Feature
    def fetch_translation
	raise "not a mRNA/CDS" unless self.feature == "CDS" or self.feature == "mRNA"
	hash = self.assoc
	entry = Bio::Fetch.query('refseq', "#{hash['protein_id'].gsub(/\.(\d+)$/, '')}","raw","fasta").split("\n")
	entry.shift
	return entry.join
    end
	
  end
	
  class Mapper
	  
    attr_accessor :position, :cigar_line, :verbose

    def initialize(position,cigar_line,verbose=false)
      @position = position+1
      @verbose = verbose
      @cigar_line = cigar_line
    end
    
    def map
      
      cig = self.cigar_line.clone
      puts cig if verbose
      ungapped_pos = 0
      gapped_pos = 0
      
      until cig.length == 0
        
        count = cig.slice!(/^(\d+)/).to_i
        char = cig.slice!(/^[A-Z]/).strip
        
        puts "#{count} #{char}" if verbose
        if char == "D"
          gapped_pos += count
          puts "\t#{ungapped_pos}/#{gapped_pos} | #{self.position}" if verbose
        elsif char == "M"
          count.times do 
            ungapped_pos += 1
            gapped_pos += 1
            puts "\t#{ungapped_pos}/#{gapped_pos} | #{self.position}" if verbose
            return gapped_pos if ungapped_pos == self.position
          end
        elsif count == nil
          raise "Count illegal #{count} #{char}"
        else
          raise "Unknown operation! #{count},#{char}"
        end
         
      end
      #warn "Reached end of alignment before match"
      return gapped_pos
    end
    
    def ungapped
    
    	cig = self.cigar_line.clone
    	gapped_pos = 0
    	ungapped_pos = 0
    	
    	until cig.length == 0
    	
				count = cig.slice!(/^(\d+)/).to_i
        char = cig.slice!(/^[A-Z]/).strip
        if char == "D"
        	gapped_pos += count
    		else
    			count.times do   				
    				ungapped_pos += 1
    				gapped_pos += 1
    				return ungapped_pos if gapped_pos >= self.position  			
    			end
    		end
    	
    	end
    	
    end
    	
	end
  
  class FastaFormat
    
    def length
      return self.naseq.length
    end
    
    def cigar_line
      
      answer = ""
      
      seq = "#{self.data.upcase.strip.gsub(/\n/, '')}"
      
      until seq.strip.length == 0
        if seq.match(/^[-]/)
          this_seq = seq.slice!(/^[-]+/)
          answer += "#{this_seq.length}D"
        else
          this_seq = seq.slice!(/^[A-Z\*]+/)
          answer += "#{this_seq.length}M"
        end
      end
      
      return answer
      
    end
      
    
  end
  
  class Tree
    
    def get_taxa
      
      answer = []
      
      self.each_node do |node|
        answer.push(node.name) unless node.name.nil? or node.name.length == 0
      end
      
      return answer
      
    end
    
  end
  
  module Alignment
    
    class OriginalAlignment
      
      def length
        return self.consensus_string.length
      end
      
      def run_baseml(tree)
        
        aln = Bio::Alignment::OriginalAlignment.new
        self.each_pair { |k,s| aln.store(k.gsub(/\s/,'')[0..6],s)}
        
        baseml = Bio::PAML::Baseml.new
        report = baseml.query(aln,tree)

        result = parse_baseml(baseml.output)
        return result
        
      end
      
      def select_taxa(taxa)
        
        answer = Bio::Alignment::OriginalAlignment.new
        self.each_pair do |taxon,seq|
          taxa.each do |t|
            answer.add_seq(seq,taxon) if taxon.include?(t)
          end
        end
        
        return answer
        
      end
        
      def print_for_paml(tree)
        
        answer = []
        answer.push("\t#{self.number_of_sequences}\t#{self.alignment_length}")
        
        self.each_pair do |k,seq|
          answer.push("#{k.gsub(/\s/, '_')[0..9]}")
          while seq.length > 0
            answer.push("#{seq.slice!(0..29)}")
          end
        end
        
        return answer.join("\n")
      end
        
      def ffd_codons
        
        ffd = [ "TC", "CT", "CC", "CG", "AC", "GT", "GC", "GG" ]
        
        keys = self.sequence_names
        
        count = 0
        sites = []    # groups 3 sites 
        codons = []    # ... into one codon
        ffd_codons = []   # holds all ffd codons
        
        self.each_site do |site|
          
          count += 1
          sites.push(site)
          
          if count == 3
            codons = Array.new
            keys.each do |k|
              one = sites[0].shift
              two = sites[1].shift
              three = sites[2].shift
              codons.push("#{one}#{two}#{three}")
            end
            add = true
            codons.each do |codon|
              add = false unless ffd.include?(codon[0..1].upcase) or codon.include?("-")
            end
            ffd_codons.push(codons) if add == true
            count = 0 
            sites.clear
          end
          
        end
        
        answer = []
        keys.each do |k|
          seq = ""
          ffd_codons.each do |site|
            seq += site.shift.slice(2..2)
          end
          answer.push(Bio::FastaFormat.new(Bio::Sequence::NA.new(seq).to_fasta(k)))
        end
        
        return Bio::Alignment::OriginalAlignment.new(answer)
        
      end  
        
      def strip_columns
        
        aln = {}
        keys = []
        
        self.each_pair do |k,v|
          aln[k] = ""
          keys.push(k)
        end
        
        self.each_site do |site|
          unless site.include?(".") or site.join.count("-") == site.nitems # removed masked sites and gap-only sites
            keys.each do |k|
              aln[k] = aln.fetch(k) + site.shift
            end
          end
        end
        
        seqs = []
        aln.each do |k,seq|
          seqs.push(Bio::FastaFormat.new(Bio::Sequence::NA.new(seq).to_fasta(k)))
        end
        
        return Bio::Alignment::OriginalAlignment.new(seqs)
        
      end
    
      def reverse_alignment
        
        new_alignment = Bio::Alignment::OriginalAlignment.new
        
        self.each_pair do |k,s|
          new_alignment.add_seq(Bio::Sequence::NA.new(s).complement.to_s,k)
        end
        
        return new_alignment
        
      end
      
      def mask_and_strip(coordinates,organism)
        
        new_alignment = Bio::Alignment::OriginalAlignment.new
        
        self.each_pair do |k,seq|
          if k.include?("#{organism}")
            seq = mask_seq(seq.to_s,coordinates.sort)
            seq.gsub!(/-/, '.')
          end
          new_alignment.add_seq(seq,k)
        end
        return new_alignment.strip_columns
      end
      
      def percent_similarity
        
        return 0 if self.number_of_sequences == 1

        site_values = []

        self.each_site do |site|

          site_pairwise = []

          valid_sites = site.delete_if {|s| s == "-" }
          
          until valid_sites.empty?    # all sites of the column that are nucleotides/amino acids (gaps are ignored)
            ref = valid_sites.shift   # the site against which to compare (removed from the search space)
            valid_sites.each do |s|
              ref == s ? site_pairwise.push(1.0) : site_pairwise.push(0.0)
            end
          end
          next if site_pairwise.empty?
          matches = site_pairwise.select { |v| v == 1.0 }         # the number of matches
          site_values.push(matches.nitems/site_pairwise.nitems) # versus the number of valid sites
        end
        total = 0.0
        site_values.each { |v| total += v }
        return (total/self.alignment_length)*100
        
      end
      
      def percent_similarity_by_column
        
        answer = []
        
        return answer if self.number_of_sequences == 1
        
        self.each_site do |site|
          
          valid_sites = site.delete_if{|s| s == "-" or s == "n"}
          next if valid_sites.nitems < 2
                    
          site_pairwise = []
          
          until valid_sites.empty?
            ref = valid_sites.shift
            valid_sites.each do |s|
              "#{ref}" == "#{s}" ? site_pairwise.push(1.0) : site_pairwise.push(0.0)
            end
          end
          
          matches = site_pairwise.select{|v| v == 1.0 }
          answer.push(matches.nitems.to_f/site_pairwise.nitems.to_f)
          
        end 
        
        return answer
        
      end
        
      def print_nexus(title,tree)
        
        taxa = {}
        
        self.each_pair do |key,seq|
          taxa["#{key}"] = "#{seq}"
        end
        
        nex = Toolbox::Nexus::Writer.new(title,taxa,tree)
        
        tree_string = Toolbox::Converter::Trees.new(tree).plain_tree
        puts nex.header
        puts nex.taxon_list
        puts nex.characters
        
        puts "BEGIN TREES;"
        puts "Tree"
        puts "#{title}_tree=#{tree_string}"
        puts "END;"
        
      end
      
      def gap_score
        
        return 0 if self.number_of_sequences == 1
        
        gaps = []
        
        self.each_site do |site|
          gaps.push((site.select{|s| "#{s}" == "-"}.nitems.to_f/site.nitems.to_f))
        end
        
        score = 0.0
        gaps.each {|g| score += g}
        return (score/self.alignment_length)*100
        
      end
      
      def check_taxa(tree)
        
        new_alignment = Bio::Alignment::OriginalAlignment.new
        
        aln_taxa = {}
        self.each_pair { |k,s| aln_taxa["#{k.gsub(/\s/, '_')}"] = s }
        
        l = self.alignment_length
        
        self.each_pair { |k,s| raise "Sequence not of the same size as alignment (#{k})" if s.length != l }
        
        tree_taxa = tree.get_taxa.collect! { |t| t.gsub(/\s/, '_')}
        
        tree_taxa.each do |tree_taxon|
          new_alignment.add_seq("-"*self.alignment_length,tree_taxon) unless aln_taxa.include?(tree_taxon)
        end
        
        aln_taxa.each_pair { |k,s| new_alignment.add_seq(s,k)}
        
        return new_alignment
        
      end
      
      def map_position(pos,ref,ref_seq=nil,homolog_seqs=[])
      
        copy = self
        
        answer = {}
      
        raise "Key (#{ref}) NOT part of the alignment" unless copy.has_key?(ref)
        copy.each_pair do |name,seq|
          if name == ref
            ref_seq = seq
          end 
        end
        
        aln_counter = 0
        internal_counter = 0
                
        ref_seq = ref_seq.split(/[A-Za-z-]/)
        
        until internal_counter == pos
          char = ref_seq.shift
          if char == "-"
            aln_counter += 1
          else
            internal_counter += 1
            aln_counter += 1
          end
        end
        
        copy.delete(ref)
        
        copy.each_pair do |name,seq|
          
          seq = seq.split(//)
          homolog_internal_counter = 0
          homolog_aln_counter = 0
          
          until homolog_aln_counter == aln_counter
            char = seq.shift
            if char == "-"
              homolog_aln_counter += 1
            else
              homolog_aln_counter += 1
              homolog_internal_counter += 1
            end
          end
          
          answer[name] = homolog_internal_counter
          
        end
          
        return answer
          
      end
      
    end
    
  end
  
  module Rfam
  
  	class Hit
  		attr_accessor :query_def,:query_from,:query_to,:bitscore,:strand,:evalue,:model_start,:model_end,:rfam_acc,:rfam_id
  		def initialize(query_def,query_from,query_to,bitscore,strand,evalue,model_start,model_end,rfam_acc,rfam_id)
  			@query_def = query_def
  			@query_from = query_from.to_i
  			@query_to = query_to.to_i
  			@bitscore = bitscore
  			@strand = strand
  			@evalue = evalue
  			@model_start = model_start
  			@model_end = model_end
  			@rfam_acc = rfam_acc
  			@rfam_id = rfam_id
  		end
  	
  	end
    
    class Report
      attr_accessor :file, :hits
      def initialize(file,hits)
        @file,@hits = file,hits
      end
    end
  	
  	class Parser
  	
  		attr_accessor :file
  		def initialize(file)
  			@file = file
  		end
  		
  		def parse

  			hits = []
  			IO.readlines(self.file).each do |line|
  				next if line.match(/^\#.*$/)
  				query_def,bla1,bla2,query_from,query_to,bitscore,strand,bla3,rfam_string = line.split("\t")
  				strand == "+" ? strand = 1 : strand = -1
  				evalue,gc_content,ncrna_id,model_end,model_start,rfam_acc,rfam_id,score = rfam_string.split(";").collect{|e| e.gsub(/^[a-z_-]*\=/, '')}
  				
  				hits << Bio::Rfam::Hit.new(query_def,query_from,query_to,bitscore,strand,evalue,model_start,model_end,rfam_acc,rfam_id)
  				
  			end
  			return Bio::Rfam::Report.new(self.file,hits)
  		end # end parser
  		
  	end # end Rfam
  	
  end
  
  module HMMER3
    
    class HMMreport
      
      attr_accessor :target, :target_accession, :tlen, :query_name, :pfam_accession, :qlen, :evalue, :full_score, :full_bias, :no, :of, :cvalue, :ivalue, :domain_score, :domain_bias, :hmm_from, :hmm_to, :align_from, :align_to, :env_from, :env_to, :acc, :target_description
      
      def initialize(parameters)
        raise "Error, parameters array invalid" if parameters.nitems == 1
        @target = parameters.shift
        @target_accession = parameters.shift
        @tlen = parameters.shift
        @query_name = parameters.shift
        @pfam_accession = parameters.shift
        @qlen = parameters.shift
        @evalue = parameters.shift
        @full_score = parameters.shift
        @full_bias = parameters.shift
        @no = parameters.shift
        @of = parameters.shift
        @cvalue = parameters.shift
        @ivalue = parameters.shift
        @domain_score = parameters.shift
        @domain_bias = parameters.shift
        @hmm_from = parameters.shift
        @hmm_to = parameters.shift
        @hmm_to = parameters.shift
        @align_from = parameters.shift
        @align_to = parameters.shift
        @env_from = parameters.shift
        @env_to = parameters.shift
        @acc = parameters.shift
        @target_description = parameters.shift  
      end
    end
    
    class HMMsearch
      
      attr_accessor :hmmfile, :infile, :outfile, :tbl_outfile
      def initialize(infile,hmmfile="/home/marc/databases/pfam/Pfam-A.hmm")
        @hmmfile = hmmfile
        @infile = infile
        @outfile = "#{infile}.hmm.out"
        @tbl_outfile = "#{infile}.hmm.domtblout"
      end
      
      def run        
        system("hmmsearch --domtblout #{self.tbl_outfile} --cpu 3 -o #{self.outfile} #{self.hmmfile} #{self.infile}")
        parse(self.tbl_outfile)
      end
      
      def parse(hmm_file)
        reports = []
        lines = IO.readlines(hmm_file).delete_if{|l| l.match(/^\#.*/)}
        lines.each do |line|
          elements = line.split(/(\s+)/).compact.collect{|l| l.strip }.delete_if{|e| e.length == 0 }
          raise "Too few elements #{elements.inspect}" if elements.nitems < 2
          reports << Bio::HMMER3::HMMreport.new(elements)
        end
        return reports 
      end
      
    end
    
  end

  class Exonerate

    attr_accessor :model, :query, :target, :params

    def initialize(model,query_seqs,target_seq,params="",name=nil)
  	  
  	  @model = model
  	  @query = query_seqs
  	  @target = target_seq
  	  @params = params
  	  
  	  @output = []
  	  
    end

    def run
      
      IO.popen("exonerate -m #{self.model} #{self.params} #{self.query} #{self.target}") do |f|
        while line = f.gets do
          @output.push(line)
        end
      end
      
      return @output

    end
    
    def parse(output)
      
    end
    
  end
  
  # DESCRIPTION
  # Wrapper for the Vienna RNA package
  module Vienna
    
    class RNAfold
      
      def initialize(sequence)
        @sequence = sequence
      end
      
      def run
        
      end
      
    end
    
    class RNAalifold
      
      attr_accessor :infile, :alignment, :options, :remove
      attr_reader :structure, :consensus
      
      def initialize(alignment,options="",remove=true)
        
        @infile = "#{rand(100000)}_alifold.aln"
        
        f = File.new(@infile,"a")
        f.puts "#{alignment}"
        f.close
        
        @alignment = alignment
        @options = options
        @output = []
        @structure = nil
        @consensus = nil
        @remove = remove
        
      end
      
      def run
        
        IO.popen("RNAalifold #{options} < #{self.infile} 2>&1") do |f|
          while line = f.gets do
            @output.push(line)
          end
        end
        
        @consensus = @output[1]
        @structure = @output[2]
        
        system("rm #{self.infile}") if remove == true
        
      end
    
    end
    
    class RNAz
      
      attr_accessor :infile, :options, :remove, :forward, :reverse
      
      def initialize(alignment,options="",remove=true)
        
        @infile = "#{rand(100000)}_alifold.aln"
        
        f = File.new(@infile,"a")
        f.puts "#{alignment.raw}"
        f.close
        
        @alignment = alignment
        @options = options
        @remove = remove
        
        @output = []
        @forward = nil
        @reverse = nil
        
      end
      
      def run
        
        IO.popen("RNAz -b #{self.infile}") do |f|
          while line = f.gets do
            @output.push(line)
          end
        end
        
        @this_strand = nil
        
        @output.each do |line|
          
          if line.match("Reading direction")
            @this_strand = line.gsub(/Reading\sdirection\:\s/, '').strip
          elsif line.match("Prediction")
            @status = line.gsub(/Prediction\:\s/, '')
            @this_strand == "forward" ? self.forward = @status : self.reverse = @status
          end
          
        end
        
        system("rm #{self.infile}") if remove == true
        
      end
    end
    
  end
  

	module Stockholm
	
		class Entry
			
			attr_accessor :de, :au, :cc, :mb, :seqs, :ac, :lookup
			def initialize(de)
				@ac = nil
				@de = de
				@au = nil
				@cc = ""
				@mb = []
				@seqs = {}
				@lookup = {}
			end
		
		end
	
		class Parser
			
			attr_accessor :infile
			def initialize(infile)		
				@infile = infile
			end
			
			def parse

				lines = IO.readlines(infile).join
				entries = lines.split(/^\/\/$/)
				
				reports = []
				
				entries.each do |entry|

					lines = entry.split("\n")
					next if lines.empty?

					definition = lines.find{ |l| l.include?("#=GF DE") }.strip_stockholm
					accession = lines.find{ |l| l.include?("#=GF AC") }.strip_stockholm
					author = lines.find{|l| l.include?("#=GF AU")}.strip_stockholm
					
					report = Bio::Stockholm::Entry.new(definition)
										
					translations = lines.select{|l| l.include?("#=GS")}.collect{|l| l.strip_stockholm}
					translations.each do |t|
						name,accession = t.split("AC").collect{|e| e.strip}
						report.lookup[name] = accession
					end

					report.ac = accession
					report.au = author
					
					lines.each do |line|
						next if line.match(/^\#.*$/) or line.match(/^$/)
						seq_acc,seq = line.split(" ").collect{|e| e.strip}.compact
						raise "#{line}" if seq.length == 0
						report.seqs.has_key?(seq_acc) ? report.seqs[seq_acc] += seq : report.seqs[seq_acc] = seq					
					end
					
					reports << report
				
				end

				return reports
						
			end
			
		end
		
	end
	
	
  # = DESCRIPTION
  # Wrapper for the RNA prediction tool Infernal
  module Infernal
  
  	class Report
  	
  		attr_accessor :model, :strand, :query, :query_start, :query_stop, :target, :target_start, :target_stop, :host_start, :host_stop, :score, :query_seq, :midline, :target_seq
  		def initialize(model)
  			@model = model
  			@strand = nil
  			@query = nil
  			@query_start = nil
  			@query_stop = nil
  			@target = nil
  			@target_start = nil
  			@target_stop = nil
  			@host_start = nil
  			@host_stop = nil
  			@score = nil
  			@query_seq = ""
  			@midline = ""
  			@target_seq = ""
  		end
  	end
    
    class Cmsearch
      
      attr_accessor :sequence, :model, :modelpath, :options, :output, :verbose, :is_empty, :cutoff
      attr_reader :models
      
      def initialize(sequence,model,verbose=false,modelpath="/home/marc/databases/rfam/models")
        @sequence = sequence
        @model = model.gsub(/\.cm/, '')
        @modelpath = modelpath
        @options = { "T" => "0" , 
          "E" => nil,
          "format" => nil,
          "toponly" => nil,
          "local" => false,
          "noalign" => false          
          }
        
        @output = []
        @verbose = verbose
        @is_empty = true
        
        RfamDB::DBConnection.connect("_91")
        rfam = RfamDB::Rfam.find_by_rfam_acc(@model)
        rfam.nil? ? @cutoff = 15.0 : @cutoff = rfam.gathering_cutoff/2.0
        puts "Scanning #{sequence.definition} with #{@model} using cutoff: #{@cutoff}"
        
      end
      
      def run
        
        f = File.new("#{self.sequence.definition}.fasta", "w+")
        f.puts self.sequence.to_s
        f.close
        
 				IO.popen("cmsearch #{modelpath}/#{model}.cm #{self.sequence.definition}.fasta") do |f|
        	while line = f.gets do
          	self.output << line
          	self.is_empty = false if line.include?("results")
          end
        end
        system("rm #{self.sequence.definition}.fasta")
				parse_new 
      
      end
      
      
      def parse_new
      
      	reports = []
      	strands = output.join("\n").split("results")
				this_strand = strands.shift
				
				return reports if self.is_empty
				
				this_strand.include?("Plus") ? strand = 1 : strand = -1
      	
      	strands.each do |this_strand|      		
      		this_strand = this_strand.split("\n").collect{|e| e.strip}.select{|e| e.length > 1}
      		hits = this_strand.join("\n").split("Query")
      		hits.each do |hit|
						next unless hit.include?("Target")      	
      			report = Bio::Infernal::Report.new(self.model)     			
      			report.strand = strand      			
      			lines = hit.split("\n")
      			coords = lines.shift.scan(/[0-9]*\s-\s[0-9]*/)
      			report.query_start = coords[0].slice(/^[0-9]*/).to_i
      			report.query_stop = coords[0].slice(/[0-9]*$/).to_i
      			report.target_start = coords[1].slice(/^[0-9]*/).to_i
      			report.target_stop = coords[1].slice(/[0-9]*$/).to_i      			
      			score = lines.shift.slice(/[0-9]*\.[0-9*]/).to_f     			
      			report.score = score
      			next if score.to_f < self.cutoff.to_f      			
      			until lines.nitems < 4
      				query_structure = lines.shift.strip
      				query_seq = lines.shift.strip.gsub(/[0-9]*/, '')
      				midline = lines.shift.strip
      				target_seq = lines.shift.strip.gsub(/[0-9]*/, '')      				
      				report.query_seq += query_seq.strip
      				report.target_seq += target_seq.strip
      				report.midline += midline      				
      			end      			
      			reports << report      			
      		end     	
      		strand = -1
      	end
      	
      	return reports
      	
      end	
    
    end
    
  end
  
  module Pecan
  
  	class Align
  		
  		attr_accessor :sequences
  		def initialize(sequences)			
  			@sequences = sequences 			
  		end
  		
  		def do_align

        aln = nil
        FileUtils.cd("/tmp") {  		
  			  f = File.new("output.#{$$}.mfa", "w+")
  			  f.close
  				
  			  a = File.new("#{$$}.fasta", "a")
				  self.sequences.each do |seq|
			  		f = File.new("#{seq.definition}.fasta", "w+")
			  		f.puts seq
			  		a.puts seq
			  		f.close
			  	end
			  	a.close				 
              
          system("clustalw -infile=#{$$}.fasta -tree")
			  	tree_string = IO.readlines("#{$$}.ph").collect{|l| l.strip.gsub(/\:(\d+)\.(\d+)/, '')}.join
			  	seqs =  tree_string.gsub(/[\(\);]/, '').split(",").collect{|s| s+".fasta"} 
          
			  	system("java -cp /home/marc/scripts/pecan_v0.8/pecan_v0.8.jar bp.pecan.Pecan -E '#{tree_string}' -F #{seqs.join(" ")}")
            
			  	aln = Bio::FlatFile.open(Bio::FastaFormat, "output.mfa")

				}
        
				return aln
  		
  		end
  		
  		
  		
  	end
  
  end

end

module RfamScan
  
  class Report
  	
  	attr_accessor :sequence_name, :rfam_acc, :start, :stop, :score, :target_seq, :query_seq, :midline
 		def initialize(sequence_name,rfam_acc,blasthit)
 			@sequence = sequence_name
 			@rfam_acc = rfam_acc
 			@start = blasthit.query_start
 			@stop = blasthit.query_end
 			@score = blasthit.bit_score
 			@target_seq = blasthit.target_seq
 			@query_seq = blasthit.query_seq
 			@midline = blasthit.midline
 		end
  		
 	end
  
 	class Scan
  	
 		attr_accessor :sequence, :start, :stop, :processes, :count, :answers, :blast_factory, :verbose
  		
 		def initialize(sequence,verbose=false)
 			@sequence = sequence
 			@answers = []
			@count = 0
 			@verbose = verbose
 			options = "-W 7 -F F -G 12 -E 8 -r 4 -q -5"
  		@blast_factory = Bio::Blast.local("blastn","/home/marc/databases/rfam/seeds/seeds.fasta",options) 		
 		end
 		
 		def run
  		puts "Querying #{self.sequence.definition}" if self.verbose
  		results = self.blast_factory.query(self.sequence)			
  		puts "...parsing BLAST output (#{self.sequence.definition})..."	if self.verbose
  		results.each do |hit|
  			rfam_acc = hit.target_def.split("|")[0]
  			present = false # Check whether this hit as already been reported
  			self.answers.select{|a| a.rfam_acc == rfam_acc}.each do |h|	# throw out overlapping hits of the same RFam id
  				present = true if h.start <= hit.query_end and h.stop >= hit.query_start
  			end
 				self.answers << RfamScan::Report.new(self.sequence.definition,hit.target_def.split("|")[0],hit) unless present
 			end
  		puts "Running infernal (#{self.sequence.definition})" if self.verbose
  		cmsearch(self.answers)			
  			
 		end
 		
 		def cmsearch(candidates,padding=200)
 		  answer = []
 		  candidates.each do |hit|
      	# Infernal scans sub-sequences only, so lets remember where those come from...
      	hit.start < padding ? start = 1 : start = hit.start-padding
      	hit.stop > self.sequence.naseq.length-padding ? stop = self.sequence.naseq.length-1 : stop = hit.stop+padding
      	# Now slice out the putative hit +/- some margin (padding)
      	this_seq = Bio::FastaFormat.new(Bio::Sequence::NA.new(self.sequence.naseq[start..stop]).to_fasta(self.sequence.definition))
      	# And run it through infernal - a sequence and a model must be specified..
      	infernal = Bio::Infernal::Cmsearch.new(this_seq,hit.rfam_acc)
      	infernal_results = infernal.run
      	infernal_results.each{|i|
      	  i.target_start += start # re-align the coordinates to the original sequence
      	  i.target_stop += start
      	  answer << i
      	}
 		  end
 		  return answer	  
 		end
  		
 	end
  	
end
