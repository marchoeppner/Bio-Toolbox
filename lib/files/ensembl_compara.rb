module Ensembl
  
  # = DESCRIPTION
  # The Ensembl::Compara module covers the compara database from
  # ensembldb.ensembl.org and includes genome alignments and orthology information.
  # For a full description of the database (and therefore the classes that
  # are available), see http://www.ensembl.org/info/software/compara/schema/index.html
  # and http://www.ensembl.org/info/software/compara/schema/schema_description.html
  #
  # = USAGE
  # To connect to the compara database, simply add 'Ensembl::Compara::DBConnection.connect'
  # to the beginning of your script. You may also want to 'include Ensembl::Compara'.
  # Note that the connection has to be spelled out fully to avoid conflicts with other 
  # database connections (Core, Vara).  
  module Compara  	
    require 'composite_primary_keys'
    
    # = DESCRIPTION
    # Contains a description of the taxonomic relationships
    # between all the taxa used in this database. It is usuall
    # read together with Ensembl::Compara::NcbiTaxaName
    class NcbiTaxaNode < DBConnection
      
    end
    
    # = DESCRIPTION
    # Contains descriptions of the taxonomic nodes
    # defined in EnsembL::Compara::NcbiTaxaNode
    class NcbiTaxaName < DBConnection

    end
    
    # = DESCRIPTION
    # Dnafrag (ments) are stretches of DNA that have been used
    # for alignments in pairwise or multiple species genomic
    # alignments. Sub-sections are referred to by the genomic_align
    # table in which information on aligned sequences are stored.
    #
    # = USAGE
    # fragment = Dnafrag.find_by_dnafrag_id(1)
    # slice = fragment.get_slice
    class Dnafrag < DBConnection
      set_primary_key 'dnafrag_id'
      belongs_to :genome_db, :foreign_key => "genome_db_id"
      has_many :genomic_aligns
      
      # = DESCRIPTION
      # Identifies and returns the coordinate system from which
      # this particular DNA fragment was taken (e.g. chromosome X)
      def coord_system_name
        self.slice.seq_region.coord_system.name
      end
      
      # = DESCRIPTION
      # Returns the id of the current Ensembl::Compara::Dnafrag object
      def display_id
        return self.dnafrag_id
      end
      
      # = DESCRIPTION
      # Creates a Ensembl::Core::Slice object from the
      # aligned DNA fragment. Accepts internal start and stop
      # position if needed. Check Ensembl::Core information
      # for available methods. 
      # = USAGE
      # fragment = Ensembl::Compara::Dnafrag.find(1)
      # fragment_slice = fragment.slice
      def slice(start=0,stop=self.length,strand=1)        
        self.genome_db.connect_to_genome_new
        seq_region = Ensembl::Core::SeqRegion.find_by_name(self.name)
        return Ensembl::Core::Slice.new(seq_region, start, stop, strand)        
      end

      # = DESCRIPTION
      # Allows to find a DNA fragment based on a Ensembl::Core::Slice
      # object, rather then the other way around. Uses a slice and the 
      # corresponding genome db as paramaters.
      # = USAGE
      # gene_slice = Ensembl::Core::Slice.find_by_gene_stable_id(ENSGxxxx)
      # genome_db = Ensembl::Compara::GenomeDB.find_by_name('homo_sapiens')
      # fragment = Esnembl::Compara::Dnafag.fetch_by_slice_and_genome_db(gene_slice,genome_db)
      def self.fetch_by_slice_and_genome_db(slice,genome_db)        
        fragment = nil
        coord = slice.seq_region.coord_system.name
        location = slice.seq_region.name
        fragment = Ensembl::Compara::Dnafrag.find(:first, :conditions => ["name = ? and coord_system_name = ? and genome_db_id = ?", location,coord,genome_db.genome_db_id,])
        return fragment        
      end
      
    end
    
    # = DESCRIPTION
    # Contains the genomic regions corresponding to every 
    # synteny relationship found. There are two genomic
    # regions for every synteny relationship.
    class DnafragRegion < DBConnection
      
      def slice
        return self.dnafrag.slice(self.dnafrag_start,self.dnafrag_end,self.dnafrag_strand)
      end
      
    end
      
    # = DESCRIPTION
    # Dna fragments are associated with the respective
    # genome database. A connection is required to fetch
    # any type of sequence data (slices) from the core database.
    # (done automatically when using any of the slice methods in
    # this module).
    #
    # = USAGE
    # gdb = GenomeDb.find_by_name('homo_sapiens')
    # gdb.connect_to_genome
    # slice = Slice.fetch_by_gene_stable_id('...')
    class GenomeDb < DBConnection
      set_primary_key 'genome_db_id'
      
      has_many :dnafrags
      has_many :members
      has_many :species_sets, :foreign_key => "genome_db_id"

      def get_all_db_links        
        answer = []
        sets = Ensembl::Compara::SpeciesSet.find_all_by_genome_db_id(self.genome_db_id)
        sets.each do |set|
          answer.push(Ensembl::Compara::MethodLinkSpeciesSet.find_all_by_species_set_id(set.species_set_id))
        end        
        return answer.sort.uniq!        
      end
            
      def db_args
        args = {}
      	self.locator.split(";").each do |element|
      		e = element.split("=")
      		args[e[0].strip] = e[1].strip
      		if e[0] == "dbname"
      			args["version"] = e[1].slice(/5[0-9]/).to_i
      		end
      	end
      	return args
     	end
     	     
      def connect_to_genome

        # new rules...let's play
        
        # Figure out where this bug lives..,
        # 1. Ensembl
        release = self.connection.instance_variable_get(:@config)[:database].split('/').last.split("_")[-1]
        name = self.name.gsub(/\s/, '_').downcase
        dummy_db = Ensembl::DummyDBConnection.connect({:port => 5306})
        dummy_connection = dummy_db.connection
        db_name = dummy_connection.select_values("SHOW DATABASES LIKE '%#{name}_core_#{release}%'").select{|e| e.include?("#{release}")}[0]
        
        if db_name
          Ensembl::Core::DBConnection.connect(name,release.to_i)
        else # 2. Ok, not EnsEMBL - try EBI
          dummy_db.disconnect!
          dummy_db = Ensembl::DummyDBConnection.connect({:host => "mysql.ebi.ac.uk", :port => 4157})
          dummy_connection = dummy_db.connection
          if self.db_args.has_key?("dbname")
            name = self.db_args["dbname"].split("_")[0..1].join("_")
            Ensembl::Core::DBConnection.connect(name,release,{:ensembl_genomes => true})
          else
            #db_name = dummy_connection.select_values("SHOW DATABASES LIKE '%#{self.name}_core'").select{|e| e.include?("#{release}")}[0]
            Ensembl::Core::DBConnection.connect(name,release,{:ensembl_genomes => true})
          end
          
        end
        
      end
      
      # = DESCRIPTION
      # Returns all ENsembl::Compara::SpeciesSet containing a given species.
      # I.e. all aligments including Homo sapiens
      def find_all_sets
        return Ensembl::Compara::SpeciesSet.find_all_by_genome_db_id(self.genome_db_id)
      end
      
      # = DESCRIPTION
      # Returns all comparisons for a given genome based on the employed method.
      # The method id may vary between database releases! 
      # Use MethodLink.fetch_methods to get a full list.
      def find_all_pairwise(method_link_id)
        answer = Array.new
        self.find_all_sets.each do |set|
          answer.push(MethodLinkSpeciesSet.find(:all, :conditions => ["species_set_id = ? and method_link_id = ?", set.species_set_id,method_link_id]))
        end
        return answer.flatten
      end  
      
      def find_pairwise_set_id(method_link_id,species_2)       
        answer = nil
        method_species_sets = self.find_all_pairwise(1)
        method_species_sets.each do |set|
          set = SpeciesSet.find_all_by_species_set_id(set.species_set_id)
          set.each do |s|
            if s.genome_db_id == species_2.genome_db_id
              answer = MethodLinkSpeciesSet.find_by_species_set_id(s.species_set_id)
            end
          end
        end          
        return answer        
      end
        
      def linked_genomes_by_method_link_id(method_link_id)
        answer = Array.new
        self.find_all_pairwise(method_link_id).each do |pair|
          set_id = pair.species_set_id
          set = Ensembl::Compara::SpeciesSet.find_all_by_species_set_id(set_id)
          set.each do |s|
            answer.push(Ensembl::Compara::GenomeDb.find(s.genome_db_id)) unless s.genome_db_id == self.genome_db_id
          end
        end
        return answer
      end

    end
    
    # = DESCRIPTION
    # Genomic Align Blocks are used to group aligned
    # sequences of a given set of species and particular alignment method.
    # Have many genomic_aligns (= aligned DNA fragments)
    class GenomicAlignBlock < DBConnection
      set_primary_key "genomic_align_block_id"
      belongs_to :method_link_species_set, :foreign_key => "method_link_species_set_id"
      has_many :genomic_aligns
      has_many :genomic_align_groups
      
      # = DESCRIPTION
      # Returns the ancestral sequence
      def fetch_ancestral_sequence(start=0,stop=self.length)
        self.genomic_aligns.select{|c| c.find_organism == "Ancestral sequences"}.each do |contig|
          puts contig.find_organism
          #return contig.aligned_sequence[start..stop]
        end
      end
            
      # = DESCRIPTION
      # Returns all sequences for an alignment block without gaps as 
      # array of Bio::Fasta objects. Slow for long blocks!
      def fetch_unaligned_sequences        
        answer = Array.new        
        self.genomic_aligns.each do |piece|          
          sequence = piece.get_slice.seq
          fas = Bio::FastaFormat.new(Bio::Sequence::NA.new(sequence).to_fasta(piece.genomic_align_id))
          answer.push(fas)          
        end        
        return answer        
      end
      
      # = DESCRIPTION
      # Returns all sequences for an alignment block including gaps as 
      # array of Bio::Fasta objects. VERY slow for long blocks!
      # Sub-sections of an alignment can be retrieved by passing
      # the start and stop relative to the alignment in parentheses.
      # As a third parameter an Array of organisms can be passed 
      # to only retrieve those sequences that are of interest.
      def alignment_strings(start=0,stop=self.length,organisms=nil)        
        answer = Array.new        
        self.genomic_aligns.each do |contig|
          if organisms.nil? # if no organisms were specified to limit the results
            sequence = contig.aligned_sequence(start,stop)
            answer << Bio::FastaFormat.new(Bio::Sequence::NA.new(sequence).to_fasta(contig.find_organism.name)) unless sequence.nil?
          else
            if organisms.include?(contig.find_organism)
              sequence = contig.aligned_sequence(start,stop)
              answer << Bio::FastaFormat.new(Bio::Sequence::NA.new(sequence).to_fasta(contig.find_organism.name))
            end
          end          
        end
        return answer        
      end
      
      # = DESCRIPTION
      # Calls alignment_strings and transforms the result into a 
      # Bio::Alignment object. Please refer to the Bioruby docs
      # for further information on available methods.
      def get_clustalw(start=0,stop=self.length,organisms=nil)        
        aln = self.alignment_strings(start,stop,organisms)        
        return Bio::Alignment::MultiFastaFormat.new(aln.join("\n"))        
      end
      
      # = DESCRIPTION
      # Returns the expected conservation score corresponding
      # to this alignment  - subregion can be passed as arguments.
      # block.get_expected_scores(0,5000,10) returns to the scores
      # for position 0 to 5000 in a 10nt window. 
      # Available window sizes are 1,10,100 and 500
      def get_expected_scores(start=0,stop=self.length,window_size=1)        
        answer = []       
        warn "Must have a window size!" if window_size.nil?        
        self.conservation_scores.select{|s| s.position >= start and s.position <= stop and s.window_size == window_size }.each do |score|
          score.get_expected_scores.each {|e| answer.push(e)}
        end        
        return answer        
      end
        
      def get_diff_scores(start=0,stop=self.length,window_size=1)        
        answer = []        
        warn "Must have a window size!" if window_size.nil?        
        self.conservation_scores.select{|s| s.position >= start and s.position <= stop and s.window_size == window_size }.each do |score|
          score.get_diff_scores.each {|e| answer.push(e)}
        end        
        return answer        
      end
      
    end
    
    # = DESCRIPTION
    # Contains conservation scores calculated from the whole-genome 
    # multiple alignments stored in the genomic_align_block table.
    class ConservationScore < DBConnection
      
      # = DESCRIPTION
      # Returns the decoded binary scores
      # (diff between expected and observed)
      # as integers in an array
      def get_diff_scores        
        answer = []        
        scores = self.diff_score
        i = 0
        while i < scores.length
          s = scores[i, 4]
          answer.push(s.unpack('4f')[0])
          i += 4
        end        
        return answer        
      end
      
      # = DESCRIPTION
      # Returns the decoded binary scores in an array
      def get_expected_scores        
        answer = []        
        scores = self.expected_score
        i = 0
        while i < scores.length
          s = scores[i, 4]
          answer.push(s.unpack('4f')[0])
          i += 4
        end        
        return answer           
      end
      
      def _get_aligned_scores_from_cigar_line        
        align_start = 1
        align_end = self.genomic_align_block.length       
      end
      
      def fetch_all_by_genomic_align_block(start=1,stop=self.genomic_align_block.length,slice_length=self.genomic_align_block.length,display_size=700,display_type="AVERAGE",window_size=nil)
        
        scores = []       
        align_start = start
        align_end = stop        
        bucket_size = (slice_length/display_size)        
        window_sizes = [ 1, 10, 100, 500]        
        warn "Invalid window size!" unless window_sizes.include?(window_size)        
        bucket = { :diff_score => 0,
          :start_pos => 0,
          :end_pos => 0,
          :start_seq_region_pos => 0,
          :end_seq_region_pos => 0,
          :called => 0,
          :cnt => 0,
          :size => bucket_size,
          :current => 0 }        
      end
      
    end
    
		# = DESCRIPTION
		class ConstrainedElement < DBConnection
			
		end
		
    # = DESCRIPTION
    # Genomic Align objects refer to each sequence (dnafrag) 
    # in a genomic alignment, grouped by a 
    # genomic_align_block_id (referring to the species
    # and alignment method).
    class GenomicAlign < DBConnection
      set_primary_key "genomic_align_id"
      belongs_to :dnafrag, :foreign_key => "dnafrag_id"
      belongs_to :genomic_align_block, :foreign_key => "genomic_align_block_id"
      has_one :genome_db, :through => :dnafrag
      
      def _cigar_element        
      end 
      
      def count_cigar_elements        
      end
      
      # = DESCRIPTION
      # This is a conviencence method that collects information
      # on organism, genome_db_id,coordinate system,name of the fragment,
      # start and stop as well as orientation and writes it into one string.
      def display_id
        organism_name = self.dnafrag.genome_db.taxon_id
        genomedb = self.dnafrag.genome_db.genome_db_id
        coord = self.dnafrag.slice.seq_region.coord_system.name
        dna_name = self.dnafrag.name
        start = self.dnafrag_start
        stop = self.dnafrag_end
        
        return "#{organism_name}:#{genomedb}:#{coord}:#{dna_name}:#{start}:#{stop}:#{self.dnafrag_strand}"
        
      end  
      
      # = DESCRIPTION
      # Alignments in ENSEMBL are not stored with gaps, but
      # with reference to the DNA fragment (e.g. chromosome, contig) and compressed
      # information on gap positions (cigar line). This method returns the
      # gapped sequence for a given GenomicAlign object. Slow!!
      # Can take the relative positions within the alignment to retrieve a subsequence.
      # = USAGE
      # fragment = GenomicAlign.find(:first)
      # aln_seq = fragment.aligned_sequence(14,500)
      # Will retrieve the gapped sequence from position 14 to 500 of this particular member of the alignment.
      def aligned_sequence(start=0,stop = nil,noindent=false)       
        self._get_aligned_sequence_from_original_sequence_and_cigar_line
        #seq = AlignSeq.new(self.get_slice.seq,self.cigar_line,start,stop).align
        #return Bio::FastaFormat.new(Bio::Sequence::NA.new(seq).to_fasta("#{self.find_organism}"))
      end 
      
      def _get_aligned_sequence_from_original_sequence_and_cigar_line      
        sequence = self.get_slice.seq
        cigar_line = self.cigar_line
        cig_elements = cigar_line.scan(/\d*[A-Z]/)
        aln_seq = ""     
        cig_elements.each do |cig|
          count = cig.slice!(/^[0-9]*/)
          count.length == 0 ? count = 1 : count = count.to_i
          type = cig.strip
          if type == "M"
            aln_seq += sequence.slice!(0..count-1)
          else
            aln_seq += ("-"*count)
          end          
        end        
        return aln_seq         
      end
      
      # = DESCRIPTION
      # Returns the name of the organism a particular aligned 
      # sequence comes from.        
      def find_organism
        return self.dnafrag.genome_db
      end
      
      # = DESCRIPTION
      # Turns a DNA fragment into a Ensembl::Core::Slice object (see Core documentation for available methods).
      # Establishes a connection to a Ensembl::Core DB - cancels any previously established Core connections!
      def get_slice        
        start = self.dnafrag_start
        stop = self.dnafrag_end
        strand = self.dnafrag_strand        
        self.dnafrag.genome_db.connect_to_genome        
        return Ensembl::Core::Slice.new(Ensembl::Core::SeqRegion.find_by_name(self.dnafrag.name), start, stop, strand)                
      end    
      
      def get_mapper(slice,verbose=false)
      
        cigar_line = self.cigar_line.clone
        aln_length = self.genomic_align_block.length
  
        if self.dnafrag_strand == -1
          start = self.dnafrag_end-slice.stop
          stop = self.dnafrag_end-slice.start	
        else
          start = slice.start-self.dnafrag_start
          stop = slice.stop-self.dnafrag_start	
        end
  
        elements = cigar_line.scan(/\d*[A-Z]/)
  
        count_ungapped,count_gapped,gapped_start = 0,0,0
  
        elements.each do |element|
          element.match(/(\d+)[A-Z]/) ? x = element.slice!(/^[0-9]*/).to_i : x = 1
          char = element
  
          x.times do			
            if char == "M"
              count_ungapped += 1
              count_gapped += 1
              gapped_start = count_gapped if count_ungapped == start
              if count_ungapped == stop
                return [ aln_length-count_gapped,aln_length-gapped_start ] if self.dnafrag_strand == 1
                return [ aln_length-count_gapped , aln_length-gapped_start ] if self.dnafrag_strand == -1
              end
            else
              count_gapped += 1
            end		
          end
  
        end	
        raise "Finished cigar line (#{count_gapped},#{count_ungapped}) without reaching coordinates (#{start}>#{stop})"
      
      end
      
      # = DESCRIPTION
      # Not all sequences in a given alignment cover the full length.
      # To determine where an aligned DNA fragment is located, use this
      # method. Returns an array (start,stop) - relative to the alignment.
      def aligned_position 
        
        cigar_line = "#{self.cigar_line}"
        
        x = cigar_line.slice!(/^[0-9]*/)
        char = cigar_line.slice!(/^[A-Z]/)
        
        x.nil? ? x = 1 : x = x.to_i
        char == "X" ? start = x : start = 0
        
        char = cigar_line.slice!(/[A-Z]$/)
        y = cigar_line.slice!(/[0-9]*$/)
        
        if y.nil?
          y = 1
        else
          y = y.to_i
        end
        
        if char == "X"
          stop = self.genomic_align_block.length - y
        else
          stop = self.genomic_align_block.length
        end
        answer = Array.new
        answer.push(start)
        answer.push(stop)
        return answer
      
      end
      
      # = DESCRIPTION
      # Returns the genomic_align_group of a given type      
      # for this genomic_align object.
      def genomic_align_group_by_type(type)
        return self.genomic_align_groups.find_by_type(type)
      end
        
    end
    
    # = DESCRIPTION
    # This table is used to index tree alignments, e.g.
    # EPO alignments. These alignments include inferred ancestral sequences.
    # The tree required to index these sequences is stored in this
    # table. This table stores the structure if the tree. 
    # Each node links to an entry in the genomic_align_group
    # table, whih links to one or several entries in the 
    # genomic_align table. 
    class GenomicAlignTree < DBConnection
      
      def genomic_align_groups
        
        return Ensembl::Compara::GenomicAlignGroup.find_all_by_group_id(self.node_id)
        
      end
      
      def aligned_sequence
        
        gapped_sequence = ""
        
        groups = self.genomic_align_groups
        puts "#{groups.nitems}"
        return groups.shift.genomic_align.aligned_sequence if groups.nitems == 1
        
        groups.each do |group|
          puts "#{group.group_id}"
          aln = group.genomic_align.aligned_sequence(multiple=true)
          puts aln
          gapped_sequence = gapped_sequence + aln
        end
        return gapped_sequence
        
      end
        
    end
    
    # = DESCRIPTION
    # This table is used to group genomic_aligns. It is
    # used to support 2X genomes in the EPO alignments.
    # These 2X genomes are split into vey small regions.
    # We than can use one of these regions in one single EPO
    # alignment. This table allows us to link one single 
    # genomic_align_tree.node_id with several genomic_align
    # entries
    class GenomicAlignGroup < DBConnection
      
      belongs_to :genomic_align_block, :foreign_key => "genomic_align_block_id"
      has_one :genomic_align, :foreign_key => "genomic_align_id"
      
      def self.inheritance_column
        nil
      end
      
      def find_organism
        return self.genomic_align.find_organism
      end
      
    end
    
    # = DESCRIPTION
    # Contains all the syntenic relationships found in the
    # relative orientation of both syntenic regions.
    class SyntenyRegion < DBConnection

    end
    
    # = DESCRIPTION
    # Species sets group species using a common
    # identifier - used in the method_link_species_set
    # class (joining organisms and alingment method information). 
    class SpeciesSet < DBConnection
	set_primary_keys :genome_db_id, :species_set_id
	belongs_to :genome_db, :foreign_key => 'genome_db_id'
     
    end
    
    # = DESCRIPTION
    # Information on the methods used (e.g. BLASTZ-NET)
    # Note: Due to database naming coventions, methods of this class
    # are accessed through object[:method], not object.method
    class MethodLink < DBConnection     
      set_primary_key "method_link_id"
      has_many :method_link_species_sets
      def self.inheritance_column
        nil
      end
    end
    
    # = MethodLinkSpeciesSets group a particular alignment method
    # with a set of species (i.e. genomic alignment of all mammals).
    # Prior knowledge on the method is required (see MethodLink).
    # = USAGE
    # dataset = MethodLinkSpeciesSet.find_by_method_link_species_set_id(332)
    # For release 50, this will give access to the genomic alignment of 12 vertebrate 
    # species.
    class MethodLinkSpeciesSet < DBConnection
      set_primary_key 'method_link_species_set_id'
      has_many :genomic_align_blocks, :foreign_key => "method_link_species_set_id"
      belongs_to :method_link, :foreign_key => "method_link_id"
      has_many :species_set, :foreign_key => "species_set_id"
     
      def self.inheritance_column
        nil
      end
      
      # = DESCRIPTION
      # Returns an array of all GenomeDb objects (access to genomes) belonging
      # to this particular dataset. 
      # = USAGE
      # dataset = MethodLinkSpeciesSet.find_by_method_link_species_set_id(338)
      # organisms = dataset.fetch_organisms
      # organisms.each do |organism|
      #   puts organism.name
      # end
      def fetch_organisms
        return SpeciesSet.find_all_by_species_set_id(self.species_set_id).collect{|s| s.genome_db }
      end
          
    end
    
    # = DESCRIPTION
    # Links sequences to the Ensembl core DB or
    # to external databases.
    class Member < DBConnection
      
      belongs_to :genome_db, :foreign_key => "genome_db_id"
      def find_organism
        return self.genome_db.name
      end
        
      # = DESCRIPTION
      # Fetches a member based on its stable id and 
      # source
      # = USAGE
      # Member.fetch_by_source_stable_id('ENSEMBLGENE','ENSG00000004059')
      def self.fetch_by_source_stable_id(s,q)
        answer = Ensembl::Compara::Member.find(:first, :conditions => ["stable_id = ? and source_name = ?",q,s])
        return answer
      end
      
      # = DESCRIPTION
      # Fetches homologues of a given  object as Ensembl::Compara::Member objects
      # = USAGE
      # some_gene = Member.find(3)
      # orthologues = some_gene.fetch_homologues('ortholog')
      # orthologues.each do |o|
      #   puts o.member.stable_id
      # end
      def fetch_homologues(type=nil) # default: all type of homologues, or else orthologues or paralogues (ortho, para)      
        answer = self.homology_members.collect{|hm| hm.homology}  
        answer = answer.select{|h| h.description.include?("#{type.downcase}")} unless type.nil? # filter for a particular kind of homology  
        return answer.collect{|h| h.homology_members.select{|hm| hm.member_id != self.id and hm.member.genome_db_id != self.genome_db_id }}.flatten.collect{|hm| hm.member }        
      end
    
      # = DESCRIPTION
      # Returns a Ensembl::Core::Gene object, if member is of
      # type "ENSEMBLGENE" - else 'nil' is returned.
      def get_gene(ensembl_version=55)       
        if self.source_name == "ENSEMBLGENE"
          self.genome_db.connect_to_genome(ensembl_version)
          return Ensembl::Core::Gene.find_by_stable_id(self.stable_id)
        else
          return nil
        end        
      end
      
      # = DESCRIPTION
      # Returns a Ensembl::Core::Transcript object, if the member is of
      # type "ENSEMBLPEP", else 'nil' is returned.
      def get_transcript        
        if self.source_name == "ENSEMBLPEP"
          self.genome_db.connect_to_genome
          return Ensembl::Core::Transcript.find_by_stable_id(self.stable_id)
        else
          return nil
        end        
      end
      
      # = DESCRIPTION
      # Returns a Ensembl::Core::Translation object, if the member is of
      # type "ENSEMBLPEP", else nil is returned.
      def get_translation        
        if self.source_name == "ENSEMBLPEP"
          self.genome_db.connect_to_genome
          return Ensembl::Core::Translation.find_by_stable_id(self.stable_id)
        else
          return nil
        end
      end
      
      # = DESCRIPTION
      # Returns array of protein members for this 
      # gene member (if member is of type "ENSEMBLGENE")
      def get_all_peptide_members        
        raise "This Ensembl::Compara::Member is not a gene!" unless self.source_name == "ENSEMBLGENE"
        return Ensembl::Compara::Member.find_all_by_gene_member_id(self.member_id)        
      end
      
      # = DESCRIPTION
      # Returns the longest protein for a particular member gene
      # (IF member is of type "ENSEMBLGENE")
      def get_longest_peptide_member       
        answer = nil
        raise "Member is not a gene!" unless self.source_name == "ENSEMBLGENE"       
        self.get_all_peptide_members.each do |peptide|
        	next unless peptide.sequence
          answer = peptide if answer.nil? or answer.seq.length <= peptide.seq.length 
        end       
        return answer        
      end
      
      # = DESCRIPTION
      # Returns the sequence length of this member
      def seq_length        
        if self.source_name == "ENSEMBLGENE"
          return self.get_gene.slice.length
        elsif self.source_name == "ENSEMBLPEP"
          return self.sequence.length
        else
          return nil
        end        
      end
    
      # = DESCRIPTION
      # Returns a slice object for the current "member", IF the member is a gene.
      # For proteins, nil is returned (proteins don't have slices). 
      def slice        
        if self.source_name == "ENSEMBLGENE"
          self.connect_to_genome
          return Ensembl::Core::Slice.fetch_by_gene_stable_id(self.stable_id)
        else
          return nil
        end       
      end
      
      def seq
        return self.sequence.sequence
      end        
          
    end
    
    # = DESCRIPTION
    # Contains the protein sequences present in the member table, used
    # for the protein alignments in ENSEMBL compara.
    class Sequence < DBConnection     
      
    end
    
    class Analysis < DBConnection      
         
    end
    
    class AnalysisDescription < DBConnection     
     
    end
    
    # = DESCRIPTION
    # Stores the raw HSP local alignment results of peptide
    # to peptide alignments returned by a BLAST run.
    class PeptideAlignFeature < DBConnection      
      
    end
    
    # = DESCRIPTION
    # Contains all the genomic homologies found. There
    # are two homology_member entries for each homology entry.
    class Homology < DBConnection
      set_primary_key 'homology_id'
      has_many :homology_members, :foreign_key => "homology_id"
      # = DESCRIPTION
      # dN/dS ratio for a given pairwise homology
      def dn_ds_ration
        return (self.dn/self.ds)
      end
      
      # = DESCRIPTION
      # Returns Bio::Alignment object for this 
      # homology pair.
      def get_clustalw        
        answer = []
        self.homology_members.each do |hm|
          answer.push(hm.aligned_sequence)
        end        
        return Bio::Alignment::OriginalAlignment.new(answer)        
      end
        
      # = DESCRIPTION
      # return the pair of members for the homology
      def gene_list
        answer = []
        self.homology_members.each { |hm| answer.push(hm.member)}
        return answer
      end
      
      # Returns the 4fold degenerate sites for this alignment
      # FIX ME!
      def get_4d_simplealign
        raise "not implemented"
      end
      
      def get_all_peptide_align_features        
        members = self.gene_list        
        answer = []
        paf = Ensembl::Compara::PeptideAlignFeature.find(:first, :conditions => ["qmember_id = ? and hmember_id = ?", members[0].member_id, members[1].member_id])
        answer.push(paf) unless paf.nil?
        paf = Ensembl::Compara::PeptideAlignFeature.find(:first, :conditions => ["qmember_id = ? and hmember_id = ?", members[1].member_id, members[0].member_id])
        answer.push(paf) unless paf.nil?
        return answer        
      end
      
      def self.fetch_all_by_member(member)       
        raise "Must be a Ensembl::Compara::Member object!" unless member.kind_of?(Ensembl::Compara::Member)       
        return member.homology_members.collect{|hm| hm.homology }
        
      end
      
    end
    
    # = DESCRIPTION
    # Contains the sequences corresponding to every genomic relationship
    # found. There are two homology_member entries for each pairwise homology entry.
    # The original alignment is not stored, but can be retrieved using #fetch_aligned_sequence.
    class HomologyMember < DBConnection
      belongs_to :member, :foreign_key => "member_id"
      belongs_to :homology, :foreign_key => "homology_id"
      # = DESCRIPTION
      # Returns the aligned sequence of a member for this particular
      # homology relationship as Bio::FastaFormat object.
      # If no sequence is present for the member, nil is returned.
      def aligned_sequence        
        peptide_member = Ensembl::Compara::Member.find_by_member_id(self.peptide_member_id)
        seq = peptide_member.sequence.sequence
        return nil if seq.nil?
        aln = Ensembl::Compara::AlignSeq.new(seq,self.cigar_line).align
        return Bio::FastaFormat.new(Bio::Sequence::NA.new(aln).to_fasta("#{self.member.stable_id}|#{peptide_member.stable_id}"))       
      end
      
    end
    
    # = DESCRIPTION
    # Contains all the group homologies found. Thera are several 
    # family_member entries for each family entry.
    class Family < DBConnection
      
      # = DESCRIPTION
      # Returns the aligned family members as
      # Bio::Alignment object
      def read_clustalw        
        answer = []
        self.family_members.each do |fmember|
          answer.push(Bio::FastaFormat.new(Bio::Sequence::AA.new(fmember.aligned_sequence).to_fasta(fmember.member.stable_id)))
        end  
        return Bio::Alignment::OriginalAlignment.new(answer)       
      end
      
    end
    
    # = DESCRIPTION
    # Stores relationships between genes/proteins
    # and protein families, inlcuding alignments. 
    class FamilyMember < DBConnection
      
      def aligned_sequence       
        seq = self.member.sequence.sequence.upcase
        return nil if seq.nil? or seq == 0
        aln = Ensembl::Compara::AlignSeq.new(seq,self.cigar_line).align      
        return aln        
      end
      
    end
    
    # = Not used.
    class Domain < DBConnection      
 
    end
    
    # = Not used.
    class DomainMember < DBConnection
           
    end
    
    # = DESCRPTION
    # Contains the datastructure for each tree.
    class ProteinTreeNode < DBConnection     
      
    end
    
    # = DESCRIPTION
    # Contains the data about the leaf present in each tree.
    class ProteinTreeMember < DBConnection     
     
    end
    
    class ProteinTreeTag < DBConnection 

    end
    
  end
  
end
