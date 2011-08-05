
module Toolbox
  
  module NcRNA
   
    require 'fileutils'
    
    class Entry
      
      attr_reader :organism, :gene, :rel_start, :rel_stop, :contig, :biotype, :stable_id 
      
      def initialize(organism,gene,rel_start,rel_stop,contig,biotype,stable_id)
        @gene = gene
        @rel_start = rel_start
        @rel_stop = rel_stop
        @organism = organism
        @contig = contig
        @biotype = biotype
        @stable_id = stable_id
      end
      
      def connect_to_genome  
        Ensembl::Core::DBConnection.connect("#{@organism.gsub(/\s/, '_').downcase}")
      end
      
    end
    
    class Selection
      
      attr_reader :infile, :blocks, :dataset_id, :organisms, :type, :counter
      
      def initialize(infile)
        @infile = infile
        @blocks = IO.readlines(infile)
        @header = @blocks.shift.split("|")
        @blocks.sort!.uniq!
        @dataset_id = @header.shift
        @type = @header.shift
        @organisms = @header.shift.split(",")
        @group_counter = 0
      end
      
      def run_search
        
        @blocks.each do |block|
          
          Ensembl::Compara::DBConnection.connect
          gab = Ensembl::Compara::GenomicAlignBlock.find(block)
          
          @entries = Hash.new
          @entries["intron"] = []
          
          gab.genomic_aligns.each do |contig|
            Ensembl::Compara::DBConnection.connect
            organism = contig.find_organism.gsub(/\s/, '_')
            puts "Checking contig #{contig.genomic_align_id} (#{organism})"
            contig.get_slice.genes.each do |gene|
              
              puts "#{organism}:#{gene.stable_id}"
              
              if gene.biotype == "protein_coding" # introns and exons are only relevant for protein-coding genes...
                nc = gene.fetch_noncoding         # gets the non-coding parts of the gene  
                con_slice = contig.get_slice
                nc.each do |intron|
                  if intron.within?(con_slice)
                    pos = contig.aligned_positions(intron)
                    @entries["intron"] = @entries.fetch("intron").push(Toolbox::NcRNA::Entry.new("#{organism}",intron,pos[0],pos[1],contig,"intron","#{gene.stable_id}"))
                  end
                end
                                
              #else
              #  unless @entries.has_key?(gene.biotype)
              #    @entries[gene.biotype] = []
              #  end
              #  pos = contig.aligned_positions(gene)
              #  @entries[gene.biotype] = @entries.fetch(gene.biotype).push(Toolbox::NcRNA::Entry.new("#{organism}",gene.slice,pos[0],pos[1],contig,"#{gene.biotype}","#{gene.stable_id}"))
              
              end
              
            end
          end
          
          @entries.each do |type,entries|   # iterates over all types of features (introns, ncRNAs, etc)
            
            puts " ...creating groups of positionally conserved features (#{type})"
            
            @groups = Hash.new
            @ranges = Hash.new
            @counter = 0
            
            entries.each do |entry|         # all objects belonging to a given type of feature (introns, ncRNAs, etc)
            
              member = nil
            
              @ranges.each do |group,range|
            
                unless entry.rel_start > range[1] or entry.rel_stop < range[0]  # only if the features overlap
                  member = group
                end
              
              end
            
              if member.nil?                # features didn't overlap, creates a new group
                @counter += 1
                @ranges["#{@counter}"] = [ entry.rel_start, entry.rel_stop ]
                @groups["#{@counter}"] = [ entry ]
              else                          # features overlap, added to existing group
                @groups["#{member}"] = @groups.fetch("#{member}").push(entry) 
              end
          
            end
          
            puts "Summary for #{type}s"
            
            @groups.each do |group,members| # iterates over all groups within a given type of feature to get alignment
              
              @organisms = []
              @min_range = []
              puts " Group #{group}"
              members.each do |member|
                
                puts "  #{member.organism}|#{member.biotype}|#{member.rel_start}|#{member.rel_stop}"
                
                unless @organisms.include?(member.organism.gsub(/_/, ' '))
                  @organisms.push(member.organism.gsub(/_/, ' '))
                end
                
                if @min_range.empty?
                  @min_range = [ member.rel_start, member.rel_stop ]
                else
                  if member.rel_start > @min_range[0]
                    @min_range[0] = member.rel_start
                  elsif member.rel_stop < @min_range[1]
                    @min_range[1] = member.rel_stop
                  end
                end
                
              end
              
              puts "  minimal range covered by this group: #{@min_range[0]}/#{@min_range[1]}"
              
              f = File.new("#{gab.genomic_align_block_id}|#{group}_#{type}.fas", "a")
              
              gab.fetch_aligned_sequences(@min_range[0],@min_range[1],@organisms).each do |a|
                f.puts ">#{a.definition}"
                puts "#{a.definition}"
                f.puts "#{a.naseq}"
              end

              f.close
            end
            
          end
          
        end
        
      end
      
      def gene_alignment(start,stop,group,entries)
        
        answer = Array.new
        
        puts "A gene alignment..."
        
        @exon_groups = Hash.new
        @intron_groups = Hash.new
        @ranges = Hash.new
        
        entries.each do |entry|
          entry.connect_to_genome
          puts entry.gene.stable_id
          puts "Number of non-coding regions: #{entry.gene.introns.nitems}"
          puts "Number of transcripts: #{entry.gene.transcripts.nitems}"
          entry.gene.transcripts.each do |transcript|
            puts " #{transcript.stable_id},#{transcript.start},#{transcript.stop}"
            transcript.exons.each do |exon|
              puts "#{exon.stable_id}|#{exon.start}|#{exon.stop}"
              puts "Rel. position: #{exon.aligned_positions.join("<->")}"
            end
          end
          
        end
        
        return "Done.."
      end
      
    end  
      
    class Blocks
      
      attr_reader :outfile, :version
      
      def initialize(file,version)
        @outfile = "#{file}_blocks.txt"
        @version = version.to_i
      end
      
      def run_search
        
        Ensembl::Compara::DBConnection.connect(self.version)
        dataset = Ensembl::Compara::MethodLinkSpeciesSet.find(:first, :conditions => "name LIKE '%pecan%'")
        
        raise "\tDataset couldn't be found in the database (#{self.version}) :(" if dataset.nil?
        
        counter = 0
        @processed_blocks = Array.new
        
        organisms = dataset.fetch_organisms
        
        # The file header
        org = organisms.collect {|o| o.name}
        f = File.new("#{@outfile}", "a")
        f.puts "#{dataset.id}|#{org.join(",")}|#{self.version}|#{Date::today}|#{Time::now}"
        f.close
        
        #The blocks
        organisms.each do |organism|
          
          puts "#{organism.name}"
          Ensembl::Core::DBConnection.connect("#{organism.name.downcase.gsub(/\s/, '_')}",version)
          ncrnas = []
          
          Ensembl::Core::Gene.find_all_by_biotype("snoRNA").each {|s| ncrnas << s }
          Ensembl::Core::Gene.find_all_by_biotype("miRNA").each {|m| ncrnas << m }
          
          ncrnas.each do |ncrna|
            
            Ensembl::Compara::DBConnection.connect(self.version)
            
            Ensembl::Core::DBConnection.connect("#{organism.name.downcase.gsub(/\s/, '_')}",version)

            fragment = Ensembl::Compara::Dnafrag.fetch_by_slice_and_genome_db(ncrna.slice,organism)
            
            warn "No fragment found for #{organism.name},ncrna: #{ncrna.stable_id}" if fragment.nil?
            next if fragment.nil?
            gal = Ensembl::Compara::GenomicAlign.find(:first, :conditions => ["dnafrag_id = ? and dnafrag_start < ? and dnafrag_end > ? and method_link_species_set_id = ?", fragment.dnafrag_id, ncrna.start, ncrna.stop, dataset.id])

            next if gal.nil?
            
            gal_block = Ensembl::Compara::GenomicAlignBlock.find_by_genomic_align_block_id(gal.genomic_align_block_id)

            next if @processed_blocks.include?(gal_block.genomic_align_block_id)
            
            counter += 1
                
            puts "#{gal.genomic_align_block_id} (#{counter})"
            f = File.new("#{self.outfile}", "a")
            f.puts "#{gal.genomic_align_block_id}"
            f.close
                
            @processed_blocks.push(gal_block.genomic_align_block_id)

          end
          
        end
        
      end 
       
    end
    
    class ConservationEntry
      
      attr_reader :stable_id, :ncrna_id, :rel_pos, :length, :organism, :intronic, :hostgene, :hostgene_description, :goterms, :biotype
      
      def initialize(stable_id,ncrna_id,organism,length,rel_pos,intronic,hostgene,hostgene_description,goterms,biotype)
        @stable_id = stable_id
        @ncrna_id = ncrna_id
        @organism = organism
        @length = length
        @rel_pos = rel_pos
        @intronic = intronic
        @hostgene = hostgene
        @hostgene_description = hostgene_description
        @goterms = goterms
        @biotype = biotype
      end
      
    end
      
    class ConservationGroup
      
      attr_accessor :start, :stop, :members, :biotype
      
      def initialize(start,stop,member,biotype)
        @start = start
        @stop = stop
        @members = [ member ]
        @biotype = biotype
      end
      
    end
      
    class FindGroups
      
      attr_reader :infile, :dataset_id, :blocks, :external_db, :organisms, :base_name, :outfile, :version, :curr_block
      
      def initialize(base_name,curr_block=0,counter=0)
        @curr_block = curr_block
        @counter = counter
        @base_name = base_name
        @infile = "#{base_name}_blocks.txt"
        @outfile = "#{base_name}_conservation.xml"
        
        @blocks = IO.readlines("#{@infile}")
        
        @header = @blocks.shift
        
        @blocks.collect!{|b| b.strip.to_i }.sort!.uniq!
        
        elements = @header.split("|")
        @dataset_id = elements.shift
        @organisms = elements.shift.gsub(/_/, ' ').split(",")
        @version = elements.shift.to_i
        
        @missing = "-"  # define the characters used for the matrix
        @present = "1"
        @absent = "0"
        
        #FileUtils.mkdir_p "#{base_name}"
      end
      
      def output(string)
        
        f = File.new(self.outfile, "a")
        f.puts "#{string}"
        f.close
        
      end
      
      def print_alignment(group,aln,start,stop)
        puts "Printing alignment..."
        f = File.new("#{@base_name}/#{@type}_#{group}.fasta", "a")
        aln.each_pair do |k,seq|
          f.puts ">#{k}"
          f.puts "#{seq[start..stop]}"
        end
        f.close
      end
      
      def run_search(threshold)
        
        raise "Must have a threshold specified" if threshold.nil?
        
        round = 0
        @threshold = threshold 
        
        if self.curr_block == 0
          output("<Analysis id='#{self.dataset_id}' date='#{Date::today}'>")
          output("\s<SearchSpace nitems='#{@organisms.nitems}'>")
          @organisms.each do |o|
            output("\s\s<Organism name='#{o}' />")
          end
          output("\s</SearchSpace>")
        end
        
        group_counter = 0
        final_group_nr = 0
        
        @blocks.each do |block_id|
           
          next if block_id.to_i < @curr_block
          
          @curr_block = block_id.to_i
          round += 1
          puts "#{round}/#{@blocks.nitems}"
      
          Ensembl::Compara::DBConnection.connect(self.version)              
          
          gal_block = Ensembl::Compara::GenomicAlignBlock.find(block_id)
          
          #alignment = gal_block.get_clustalw
          
          contigs = gal_block.genomic_aligns

          @ncrnas = Array.new
          
          contigs.each do |contig|  # Find all the ncRNAs and determine the relative parameters...
            
            Ensembl::Compara::DBConnection.connect(self.version)              
            
            next if contig.find_organism == "Ancestral sequences"               # if, for whatever reason, an aligned contig is empty or does not belong to a species.....
            
            contig.get_slice.genes.select{|g| g.biotype == "snoRNA" or g.biotype == "miRNA"}.each do |gene|

              contig.dnafrag.genome_db.connect_to_genome
              
              if ! @ncrnas.include?("#{gene.stable_id}") # gene is of the right type and has not been processed before...
                
                s = File.new("#{base_name}_sequences.sql", "a")
                s.puts "INSERT INTO ncrna_sequence (ensembl_ncrna_id,seq,strand) SELECT id,'#{gene.slice.seq}','#{gene.strand}' FROM ensembl_ncrna WHERE stable_id = '#{gene.stable_id}';"
                s.close
                  
                organism = "#{contig.find_organism}"
                  
                stable_id = "#{gene.stable_id}"
                type = "#{gene.biotype}"
                ncrna_id = "#{gene.ncrna_id}"
                ncrna_id = "none" if ncrna_id.length == 0
                
                length = "#{gene.slice.seq.length}"
                rel_pos = contig.get_mapper(gene)                        # determines the position of the ncRNA gene in the alignment

                if gene.is_intronic?
                  intronic = "true"
                  hg = gene.fetch_hostgene
                  hostgene = "#{hg.stable_id}"
                  hg.description.nil? ? hostgene_description = "" : hostgene_description = "#{hg.description[0..100].gsub(/\'/, '').gsub(/\//, '').gsub(/\\/, '').gsub(/\./, '').gsub(/\,/, '-')}" 
                  goterms = hg.go_terms_with_evidence
                else
                  intronic = "false"
                  hostgene = "nil"
                  hostgene_description = "nil"
                  goterms = []
                end
                 
                @ncrnas.push(Toolbox::NcRNA::ConservationEntry.new(stable_id,ncrna_id,organism,length,rel_pos,intronic,hostgene,hostgene_description,goterms,gene.biotype))  
              
              end
            
            end
          
          end
          
          groups = Array.new
          
          until @ncrnas.empty?
            
            entry = @ncrnas.shift
            puts "processing #{entry.stable_id} | #{entry.rel_pos.join(",")}"
            start = entry.rel_pos[0].to_i
            stop = entry.rel_pos[1].to_i
            add = true
            biotype = entry.biotype
            groups.each do |group|
              break if add == false
              puts "\t checking group #{group.start}/#{group.stop}"
              if start >= group.start and stop <= group.stop
                group.stop = stop
                group.start = start
                group.members.push(entry)
                add = false
              elsif stop > group.start and stop < group.stop
                group.stop = stop
                group.members.push(entry)
                add = false
              elsif start > group.start and start < group.stop
                group.start = start
                group.members.push(entry)
                add = false
              elsif start <= group.start and stop >= group.stop
                add = false
                group.members.push(entry)
              end
              
            end
            puts "not part of above groups!" if add == true
            groups.push(Toolbox::NcRNA::ConservationGroup.new(start,stop,entry,biotype)) if add == true
            
          end
          
          groups.each do |group|
            
            group_counter += 1
                              
            organism_check = Hash.new
            self.organisms.collect { |o| organism_check[o] = @missing } # default assumption: no information for any of the genomes...
            
            output("\s<Group nr='#{group_counter}' genomic_align_block_id='#{gal_block.genomic_align_block_id}' start='#{group.start}' stop='#{group.stop}' biotype='#{group.biotype}'>")
            
            group.members.each do |member|
              organism_check.delete(member.organism)
              output("\s\s<NcRNA stable_id='#{member.stable_id}' ncrna_id='#{member.ncrna_id}' organism='#{member.organism}' intronic='#{member.intronic}' comment='#{@present}' >")  
              unless member.hostgene == "nil"
                output("\s\s\s<HostGene name='#{member.hostgene}' description='#{member.hostgene_description}'>")
                member.goterms.each do |go|
                  output("\s\s\s\s<GoId id='#{go.shift}' linkage='#{go.shift}' />")
                end
                output("\s\s\s</HostGene>")
              end
              output("\s\s</NcRNA>")
            end
            
            # check if the missing organisms have sequence data at that position,...
            
            organism_check.each do |org,val|
              contigs.each do |contig|
                if contig.find_organism == org
                  if contig.covers?(group.start)
                    organism_check.delete("#{org}")
                    output("\s\s<NcRNA stable_id='nil' ncrna_id='nil' rel_start='nil' organism='#{org}' intronic='false' comment='#{@absent}' />")
                  end
                end
              end
            end  
            
            # or else note that no sequence covers this part of the alignment (missing)
            organism_check.each do |org,val|
              output("\s\s<NcRNA stable_id='nil' ncrna_id='nil' rel_start='nil' organism='#{org}' intronic='false' comment='#{@missing}' />")
            end
            
            output("\s</Group>")
          end
          
        end  
        
        output("</Analysis>")
        
      end
     
    end  
    
    class Statistics
      
      attr_reader :type, :dataset_id, :outfile, :version
      
      def initialize(type,dataset_id,file,version=54)
        @type = type
        @dataset_id = dataset_id
        @outfile = "#{file}_all.xml"
        @version = version.to_i
        if @type == "snoRNA" or @type == "snRNA"
          @external_db = "RFAM"
        elsif @type == "miRNA"
          @external_db = "miRBase"
        end
      end
      
      def output(string)
        f = File.new("#{@outfile}", "a")
        f.puts "#{string}"
        f.close
      end
      
      def run_search
        
        Ensembl::Compara::DBConnection.connect(self.version)

        dataset = Ensembl::Compara::MethodLinkSpeciesSet.find_by_method_link_species_set_id(self.dataset_id)

        organisms = dataset.fetch_organisms
        
        output("<StatisticalOverview dataset_id='#{self.dataset_id}' type='#{self.type}'>")
        
        organisms.each do |organism|
          
          Ensembl::Compara::DBConnection.connect(self.version)
          
          output("\s<Organism name='#{organism.name}'>")
          
          organism_name = "#{organism.name.gsub(/\s/, '_').downcase}"
          Ensembl::Core::DBConnection.connect("#{organism_name}",self.version)
          
          ncrnas = Ensembl::Core::Gene.find_all_by_biotype("#{@type}")
         
          ncrnas.each do |ncrna|

            Ensembl::Core::DBConnection.connect("#{organism_name}",self.version)
            
            f = File.new("#{type}_#{self.version}_all_sequences.sql", "a")
            f.puts "INSERT INTO ncrna_sequence(ensembl_ncrna_id,seq,strand) SELECT id,'#{ncrna.slice.seq}','#{ncrna.strand}' FROM ensembl_ncrna WHERE stable_id = '#{ncrna.stable_id}';"
            f.close
            
            output("\s\s<NcRNA stable_id='#{ncrna.stable_id}' external_id='#{ncrna.ncrna_id(@external_db)}' intronic='#{ncrna.is_intronic?}'>")
            
            if ncrna.is_intronic?
              
              host_gene = ncrna.fetch_hostgene
              host_gene.description.nil? ? description = "" : description = host_gene.description.sql_compliant
             
              output("\s\s\s<HostGene stable_id='#{host_gene.stable_id}' description='#{description}' >")
              go_terms = host_gene.go_terms_with_evidence

              go_terms.each do |go|
                output("\s\s\s\s<GoId id='#{go.shift}' linkage='#{go.shift}' />")
              end
              output("\s\s\s</HostGene>")
            end

            output("\s\s</NcRNA>")

          end
          
        output("\s</Organism>")

        end

        output("</StatisticalOverview>")
      
      end
      
    end
    
    class ConstraintGroup
      
      attr_accessor :type, :block, :start, :stop, :members
      attr_writer :start, :stop, :members
      
      def initialize(block)
        @block = block
      end
      
    end
    
    # = DESCRIPTION
    # stores all annotated features within a Ensembl::Compara::GenomicAlign object
    class ContigMap
      
      attr_accessor :gal_id, :organism, :features
      
      def initalize(gal_id,organism)
        @gal_id = gal_id
        @organism = organism
        @features = {}
        get_features
      end
      
      def get_features
        
        contig = Ensembl::Compara::GenomicAlignBlock.find(self.gal_id)
        ncrs = contig.fetch_noncoding
        ncrs.each do |nc|
        end
        
        
      end
      
    end
    
    class Constraints
      
      attr_reader :blocks, :dataset_id, :type, :organisms
      
      def initialize(base_name,curr_block,counter)
        @curr_block = curr_block
        @counter = counter
        @infile = "#{base_name}_blocks.txt"
        @blocks = IO.readlines("#{@infile}")
        @header = @blocks.shift
        @blocks.sort!.uniq!
        elements = @header.split("|")
        @dataset_id = elements.shift
        @type = elements.shift
        @organisms = elements.shift.gsub(/_/, ' ').split(",")
      end
      
      def run_search
        
        Ensembl::Compara::DBConnection.connect
        #dataset = Ensembl::Compara::MethodLinkSpeciesSet.find(self.dataset_id)

        self.blocks.each do |block|
          
          gal_block = Ensembl::Compara::GenomicAlignBlock.find(block)
          alignment = gal_block.get_clustalw
          
          
          
        end  
        
      end  
      
     
    end
      
    class MIRBASE
      
      attr_reader :mirna, :mirbase_id, :organism
      attr_writer :mirbase_id
      
      def initialize(mirna,organism)
        @mirna = mirna
        @mirbase_id = nil
      end
      
      def map_family
        
        @options = "-m 7 -e 0.1"
        blastdb = File.dirname(__FILE__).gsub(/\/lib/, '') + "/data/hairpin.fa"
        factory = Bio::Blast.local('blastn',blastdb,@options)
        fas = Bio::FastaFormat.new(Bio::Sequence::NA.new(mirna.slice.seq).to_fasta(mirna.stable_id))
        report = factory.query(fas)
        count = 0
        family = "none"
        report.each_hit do |hit|
          if count == 0 # Return the best match as putative mirbase ID
            name = "#{hit.target_def.strip.split("\s")[0].downcase}"
            name.slice!(/-[a-z0-9]*$/) if name.match(/[a-z]*-[a-z0-9]*-[a-z0-9]*-[a-z0-9]/)
            @mirbase_id = name
            return name
            count += 1
          end
        end
        return "#{family}"
      
      end
      
      def targets(organism)
        
        return "none" if self.mirbase_id.nil?
        
        organism = organism.gsub(/\s/, '_').downcase
        
        target_file = File.dirname(__FILE__).gsub(/\/lib/, '') + "/data/targets_#{organism}.txt"
        
        IO.foreach(target_file) do |line|
          
          if line.include?(self.mirbase_id.gsub(/mir/, 'miR'))
            gene = line.scan(/ENS[A-Z]*[0-9]*/)
            puts "Checking #{gene} in #{organism}"
            Ensembl::Core::DBConnection.connect("#{organism}")
            ens_gene = Ensembl::Core::TranscriptStableId.find("#{gene}").transcript.gene
            puts "#{ens_gene.stable_id}"
            go_terms = ens_gene.go_terms
            puts "... #{go_terms.join(",")}"
          end
          
        end
        
      end
        
    end
      
  end

end