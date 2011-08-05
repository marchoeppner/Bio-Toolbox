
# Mixins into Bioruby

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
  
  
#
# = bio/appl/hmmer/hmmer3report.rb - hmmscan/hmmsearch parser
#
# Copyright::   Copyright (C) 2010
#               Christian Zmasek <cmzmasek@yahoo.com>
# License::     The Ruby License
#
# $Id:$
#
# == Description
#
# Parser class for hmmsearch and hmmscan in the HMMER 3 package.
#
# == Examples
#
#    #for multiple reports in a single output file (example.hmmpfam)
#    Bio::HMMER.reports(File.read("example.hmmpfam")) do |report|
#      report.program['name']
#      report.parameter['HMM file']
#      report.query_info['Query sequence']
#      report.hits.each do |hit|
#        hit.accession
#        hit.description
#        hit.score
#        hit.evalue
#        hit.hsps.each do |hsp|
#          hsp.accession
#          hsp.domain
#          hsp.evalue
#          hsp.midline
#      end
#    end
#
# == References
#
# * HMMER
#   http://hmmer.janelia.org/
#

  class Hmmer3Report
    def initialize(hmmer_output_file)
      @hits = Array.new
      @line_number = 0
      @type = nil
      my_hmmer_output_file = File.new(hmmer_output_file.to_s)
      my_hmmer_output_file.each_line() { |line| parse_line(line) }
    end

    attr_reader :hits

    private

    def parse_line(line)
      @line_number += 1
      if line  =~ /^#.+this\s+domain/
        @type = :domtblout
      elsif line =~ /^#.+best\s+1\s+domain/
        @type = :tblout
      elsif line =~ /\S/ && line !~ /^#/
        if @type == :domtblout
          @hits << HmmerPerDomainHit.new(line, @line_number)
        elsif  @type == :tblout
          @hits << HmmerPerSequenceHit.new(line, @line_number)
        else
          raise ArgumentError, "attempt to parse hmmscan/hmmsearch output style other than \"domtblout\" or \"tblout\""
        end
      end
    end
  end # class Hmmer3Report

  class Hmmer3Hit
    def initialize
      # This is an abstract class. Prevents 'new' being called on this class
      # and force implementation of 'initialize' in inheriting classes.
      raise NotImplementedError
    end
    attr_reader :target_name
    attr_reader :target_accession
    attr_reader :target_description
    attr_reader :query_name
    attr_reader :query_accession
    attr_reader :full_sequence_e_value
    attr_reader :full_sequence_score
    attr_reader :full_sequence_bias

  end # class Hmmer3hit


  class HmmerPerSequenceHit < Hmmer3Hit

    # Sets hit data.
    def initialize(line, line_number)

      #               tn    tacc      qn    qacc fs_eval fs_scor fs_bias   bst_e bst_scor bst_bias   exp     reg     clu      ov     env      dom     rep     inc   desc
      #                1       2       3       4       5       6       7       8       9      10      11      12      13      14      15       16      17      18     19
      if  line =~ /^(\S*)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(.*)/
        @target_name = $1
        @target_accession = $2
        @query_name = $3
        @query_accession = $4
        @full_sequence_e_value = $5.to_f
        @full_sequence_score = $6.to_f
        @full_sequence_bias =  $7.to_f
        @best_1_domain_e_value = $8.to_f
        @best_1_domain_score = $9.to_f
        @best_1_domain_bias =  $10.to_f
        @domain_number_est_exp = $11.to_i
        @domain_number_est_reg = $12.to_i
        @domain_number_est_clu = $13.to_i
        @domain_number_est_ov  = $14.to_i
        @domain_number_est_env = $15.to_i
        @domain_number_est_dom = $16.to_i
        @domain_number_est_rep = $17.to_i
        @domain_number_est_inc = $18.to_i
        @target_description = $19
      else
        raise ArgumentError, "line "+ line_number.to_s + " is in a unrecognized format"
      end

    end # initialize

    attr_reader :best_1_domain_e_value
    attr_reader :best_1_domain_score
    attr_reader :best_1_domain_bias
    attr_reader :domain_number_est_exp
    attr_reader :domain_number_est_reg
    attr_reader :domain_number_est_clu
    attr_reader :domain_number_est_ov
    attr_reader :domain_number_est_env
    attr_reader :domain_number_est_dom
    attr_reader :domain_number_est_rep
    attr_reader :domain_number_est_inc

  end # class HmmerPerSequenceHit

  class HmmerPerDomainHit < Hmmer3Hit

    # Sets hit data.
    def initialize(line, line_number)

      #                tn     acc    tlen   query     acc    qlen  Evalue   score    bias      #       of     c-E    i-E     score   bias      hf      ht      af      at     ef      et     acc  desc
      #                 1       2       3       4       5       6       7       8       9      10      11      12     13      14      15       16      17      18      19     20      21      22     23
      if  line =~ /^(\S*)\s+(\S+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s*(.*)/
        @target_name = $1
        @target_accession = $2
        @target_length = $3.to_i
        @query_name = $4
        @query_accession = $5
        @query_length = $6.to_i
        @full_sequence_e_value = $7.to_f
        @full_sequence_score = $8.to_f
        @full_sequence_bias =  $9.to_f
        @domain_number = $10.to_i
        @domain_sum = $11.to_i
        @domain_c_e_value = $12.to_f
        @domain_i_e_value = $13.to_f
        @domain_score = $14.to_f
        @domain_bias = $15.to_f
        @hmm_coord_from = $16.to_i
        @hmm_coord_to = $17.to_i
        @ali_coord_from = $18.to_i
        @ali_coord_to = $19.to_i
        @env_coord_from = $20.to_i
        @env_coord_to = $21.to_i
        @acc = $22.to_f
        @target_description = $23
      else
        raise ArgumentError, "line "+ line_number.to_s + " is in a unrecognized format"
      end

    end # initialize

    attr_reader :target_length
    attr_reader :query_length
    attr_reader :domain_number
    attr_reader :domain_sum
    attr_reader :domain_c_e_value
    attr_reader :domain_i_e_value
    attr_reader :domain_score
    attr_reader :domain_bias
    attr_reader :hmm_coord_from
    attr_reader :hmm_coord_to
    attr_reader :ali_coord_from
    attr_reader :ali_coord_to
    attr_reader :env_coord_from
    attr_reader :env_coord_to
    attr_reader :acc

  end # class HmmerPerDomainHit

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
