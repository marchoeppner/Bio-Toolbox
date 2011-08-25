module Ensembl
  
  module Core
    
    class Exon < DBConnection
      
      def self.find_by_stable_id(stable_id)
				return Ensembl::Core::ExonStableId.find_by_stable_id(stable_id).exon
      end
      
      def stable_id
	 			return self.exon_stable_id.stable_id
      end
            
    end

    class Gene < DBConnection
      
      def display_details
        return "#{self.stable_id} (#{self.biotype} |#{self.start}>#{self.stop}|#{self.strand})"
      end
      
      def self.fetch_by_embl_acc(acc)
        return Ensembl::Core::Xref.gene_by_embl_acc(acc)
      end
      
      def self.fetch_by_dbprimary_acc(acc)
        return Ensembl::Core::Xref.gene_by_dbprimary_acc(acc)
      end
      
      def slice_exception?
        aes = Ensembl::Core::AssemblyException.find_all_by_seq_region_id(self.seq_region_id)
        return false if aes.empty?
        this_slice = self.slice
        aes.each do |ae|
          return true if Slice.new(self.seq_region, ae.seq_region_start, ae.seq_region_end).overlaps?(this_slice)
        end
        return false
      end
      
      def structure_to_svg(x=10,y=150)

      	coordinates = []
      	introns = self.fetch_noncoding
      	
      	introns.each do |nc|
      		coordinates << nc.start
      		coordinates << nc.stop
      	end
      	
      	coordinates.sort!
      	
      	graph_width = 700
      	gene_length = self.stop-self.start.to_i
      	
      	graph = SVGWriter::Graph.new(graph_width+50)
      	
      	graph.add_path({:x => x, :y => y + 10},graph_width-20)

      	previous = self.start
      	
      	#iterate over non-coding regions and determine coordinates for coding regions
      	until coordinates.empty?      	
      		this_start = coordinates.shift.to_i
      		this_stop = coordinates.shift.to_i     		
      		this_width = (graph_width.to_f/gene_length.to_f)*(this_start-previous)
      		this_x = (graph_width.to_f/gene_length.to_f)*(previous-self.start)+x
      		this_y = y		
      		graph.add_rectangle({:x => this_x, :y => this_y},this_width) 		
      		if coordinates.empty? # coordinates for the final exon
      			previous = this_x+this_width+(graph_width.to_f/gene_length.to_f)*(this_stop-this_start)
      		else
      			previous = this_stop
      		end      		
      	end
      	
      	#draw the final exon
      	graph.add_rectangle({:x => previous, :y => y},graph_width-previous)      	
      	return graph
      	      
      end
      
      def snornas
        return self.slice.genes.select{|g| g.biotype == "snoRNA"}
      end
      
      def genomic_context
      	
      	downstream = Ensembl::Core::Gene.find(:all, :conditions => ["seq_region_id = ? AND seq_region_end < ? AND biotype = ?","#{self.seq_region_id}","#{self.seq_region_start}","protein_coding"],:order => "seq_region_start DESC")[0]
      	upstream = Ensembl::Core::Gene.find(:all, :conditions => ["seq_region_id = ? AND seq_region_start > ? AND biotype = ?","#{self.seq_region_id}","#{self.seq_region_end}","protein_coding"],:order => "seq_region_start ASC")[0]
      	
      	if downstream and upstream
      		downdist = self.start-downstream.stop
      		updist = upstream.start-self.stop 
      	elsif downstream
      		downdist = self.start-downstream.stop
      		updist = "NA"
      	elsif upstream
      		updist = upstream.start-self.stop 
      		downdist = "NA"
      	end
      	
      	downstream ? downstream_string = "#{downstream.stable_id}[#{downstream.strand}]" : downstream_string = "NA[NA]"
      	upstream ? upstream_string = "#{upstream.stable_id}[#{upstream.strand}]" : upstream_string = "NA[NA]"
      	
      	return "#{downstream_string}|#{downdist}>|#{self.stable_id}[#{self.strand}]|<#{updist}|#{upstream_string}"
      
      end
      
      def softmasked_seq        
        seq = self.seq.upcase        
        noncoding = self.fetch_noncoding        
        noncoding.each do |n|          
          seq.gsub!(/#{n.seq.upcase}/, "#{n.seq.downcase}")          
        end        
        return seq       
      end
      
      def genbank_id
        
        genbank_id = Ensembl::Core::ExternalDb.find_by_db_name('EMBL').id
        xref = self.all_xrefs.select{|x| x.external_db_id == genbank_id}[0]
        return nil if xref.nil?
        return xref.dbprimary_acc
        
      end
      
      def uniprot_id
        db_id = ExternalDb.find_by_db_name("Uniprot/SWISSPROT").external_db_id
        return self.all_xrefs.select{|r| r.external_db_id == db_id}.collect{|r| r.dbprimary_acc}.uniq.compact.shift
      end
      
      def ncrna_id
        if self.biotype == "snoRNA" or self.biotype == "snRNA" or self.biotype == "scaRNA"
          dbname = "RFAM"
        elsif self.biotype == "miRNA"
          dbname = "miRBase"
        end
        db = ExternalDb.find_by_db_name(dbname)
        if db.nil?
            return ""
        else
            db_id = db.external_db_id
            return self.all_xrefs.select{|r| r.external_db_id == db_id}.collect{|r| r.dbprimary_acc}.uniq.join
        end
      end
      
      def get_externaldb_accs(dbname)
        db_id = ExternalDb.find_by_db_name(dbname).external_db_id
        return self.all_xrefs.select{|r| r.external_db_id == db_id}.collect{|r| r.dbprimary_acc}.uniq
      end
      
      def to_fasta        
        return Bio::FastaFormat.new(Bio::Sequence::NA.new(self.slice.seq).to_fasta("#{self.stable_id}"))      
      end
      
      def get_pfam_features       
        answer = []
        check = []
        self.transcripts.each do |transcript|
          next if transcript.translation.nil?
          transcript.translation.get_pfam_features.each do |pfam|
            unless check.include?("#{pfam.hit_name}|#{pfam.hit_start}|#{pfam.hit_end}")
              answer << pfam
              check << "#{pfam.hit_name}|#{pfam.hit_start}|#{pfam.hit_end}"
            end
          end
        end        
        return answer.sort_by{|a| a.score }.reverse        
      end
      
      # = DESCRIPTION
      # This method is used to identify slices that are located
      # within the intron of protein coding genes.
      def is_intronic?
        genes = Ensembl::Core::Gene.find(:all, :conditions => ["biotype = ? AND seq_region_start < ? and seq_region_end > ? and seq_region_strand = ? and seq_region_id = ?","protein_coding",self.start,self.start,self.strand,self.seq_region_id])  
        genes.empty? ? answer = false : answer = true
        return answer
        #ncrna = self.slice
        #genes = ncrna.genes(true)
        #genes.select{|g| g.biotype == "protein_coding"}.each do |gene|
        #  unless gene.stable_id == self.stable_id
        #    gene_slice = Slice.fetch_by_gene_stable_id(gene.stable_id)
        #    return true if ncrna.within?(gene_slice)
        #  end
        #end
        #return false
      end
      
      def strictly_intronic?
        return false unless self.is_intronic?
        ncrna = self.slice
        host = self.fetch_hostgene
        return false if host.nil? or host.strand != self.strand
        host.transcripts.sort_by{|t| t.exons.nitems}.each do |t|
          t.introns.each do |i|
            return true if ncrna.overlaps?(i.slice)
          end
        end
        return false
      end
      
      def go_terms_with_evidence        
        answer = []
        go_db_id = Ensembl::Core::ExternalDb.find_by_db_name("GO").external_db_id
        self.transcripts.select{|t| t.translation.nil? == false }.each do |transcript|
          transcript.translation.object_xrefs.select{|o| o.xref.external_db_id == go_db_id}.each do |oxref|
            answer.push(["#{oxref.xref.dbprimary_acc}" , "#{oxref.go_xref.linkage_type}"])
          end
        end        
        return answer          
      end
      
      def fetch_hostgene
      
        genes = Ensembl::Core::Gene.find(:all, :conditions => ["biotype = ? AND seq_region_start < ? and seq_region_end > ? and seq_region_strand = ? and seq_region_id = ?","protein_coding",self.start,self.start,self.strand,self.seq_region_id])  
        genes.empty? ? answer = nil : answer = genes.shift
        
      end
      
      def get_longest_transcript        
        answer = nil
        self.transcripts.each do |transcript|
          next unless transcript.translation
          if answer.nil?
            answer = transcript
          elsif transcript.cds_seq.length > answer.cds_seq.length
            answer = transcript
          end 
        end
        return answer       
      end
      
      def fetch_noncoding        
        ranges = []        
        self.transcripts.each do |transcript|
          transcript.introns.each do |intron|
            coords = [ intron.seq_region_start,intron.seq_region_end].sort
            ranges << Range.new(coords.shift,coords.shift)            
          end          
        end        
        answer = []        
        until ranges.empty?
          this_range = ranges.shift
          add = true          
          ranges.each do |range|
            if this_range.eql?(range)   # the same intron is present in the remaining set, ignore the current one
              add = false
            elsif this_range.member?(range.first) and this_range.member?(range.end) # this intron contains another intron - ignore it
              add = false
            elsif this_range.member?(range.first) and this_range.member?(range.end) == false # overlapping upstream, take the overlap only 
              ranges.delete(range)
              this_range = Range.new(range.first,this_range.end)
            elsif this_range.member?(range.first) == false and this_range.member?(range.end) # overlapping downstream, take the overlap only
              this_range = Range.new(this_range.first,range.end)
              ranges.delete(range)
            end            
          end         
          answer << Ensembl::Core::Slice.new(self.seq_region,this_range.first,this_range.end,self.strand) if add == true         
        end       
       return answer        
      end      
    end
    
    class Transcript < DBConnection
      
      def to_fasta      	
      	return Bio::FastaFormat.new(Bio::Sequence::NA.new(self.slice.seq).to_fasta(self.stable_id))
      end
      
      def translation_exons
        answer = []
        self.exons.each do |exon|
          answer.push(exon) if exon.seq_region_end >= self.coding_region_genomic_start and exon.seq_region_start <= self.coding_region_genomic_end
        end
        return answer
      end
      
      def cdna_slices
        translation = self.translation
        answer = []
        self.translation_exons.each do |exon|
          if exon.exon_id == translation.start_exon_id
            slice = Ensembl::Core::Slice.new(exon.seq_region,self.coding_region_genomic_start,exon.seq_region_end,exon.seq_region_strand) if "#{self.strand}" == "1"
            slice = Ensembl::Core::Slice.new(exon.seq_region,exon.seq_region_start,self.coding_region_genomic_end,exon.seq_region_strand) if "#{self.strand}" == "-1"
            answer.push(slice)
          elsif exon.exon_id == translation.end_exon_id
            slice = Ensembl::Core::Slice.new(exon.seq_region,exon.seq_region_start,self.coding_region_genomic_end,exon.strand) if "#{self.strand}" == "1"
            slice = Ensembl::Core::Slice.new(exon.seq_region,self.coding_region_genomic_start,exon.seq_region_end, exon.seq_region_strand) if "#{self.strand}" == "-1"
            answer.push(slice)
          else
            answer.push(exon.slice)
          end
        end
        return answer.uniq
      end
     
      def all_xrefs
        return Ensembl::Core::ObjectXref.find_all_by_ensembl_id_and_ensembl_object_type(self.transcript_id,"Transcript").collect{|ox| ox.xref }
      end
      
      def get_externaldb_accs(dbname)
        db_id = ExternalDb.find_by_db_name(dbname).external_db_id
        return self.all_xrefs.select{|r| r.external_db_id == db_id}.collect{|r| r.dbprimary_acc}.uniq.join
      end
      
      def get_snorna_region(pos)        
        exon = self.exon_for_genomic_position(pos)       
      end
      
      def map_slice(sno_slice,verbose=false)
      	
      	genomic_start,genomic_end = self.coding_region_genomic_start,self.coding_region_genomic_end
      	      	
      	# First we check whether the slice is within the coding region - return 0 or cds+1 if not
      	
      	
      		if sno_slice.start < genomic_start
      			self.strand == 1 ? answer = 1 : answer = self.cds_seq.length+1
      			return answer
      		elsif sno_slice.start > genomic_end
      			self.strand == 1 ? answer = self.cds_seq.length+1 : answer = 1
      			return answer
      		end

      	pos = nil     	
      	self.introns.each do |intron|
      		if sno_slice.within?(intron)
						puts "snoRNA within intron: #{intron.previous_exon.stable_id}|#{intron.previous_exon.stop} <- #{intron.start}/#{intron.stop} -> #{intron.next_exon.start}|#{intron.next_exon.stable_id}" if verbose
					if intron.previous_exon.stop > genomic_end
						puts "exon (partly) non-coding, using genomic_coding_region... (#{self.coding_region_genomic_end} instead of #{intron.previous_exon.stop}), direction is #{self.strand}" if verbose
						pos = genomic_end
					elsif intron.previous_exon.stop < genomic_start
						puts "exon (partly) non-coding, using genomic_coding_region... (#{self.coding_region_genomic_start} instead of #{intron.previous_exon.stop}), direction is #{self.strand}" if verbose
						pos = genomic_start
					else
				 		pos = intron.previous_exon.stop
					end

					end
      	end
    	
      	if pos.nil?       		
      		#if snorna.slice.within?(self.slice)
      		if self.strand == -1
      			if sno_slice.start > genomic_start
      				pos = genomic_start
      			else
      				pos = genomic_end
      			end      				
      		else
      			if sno_slice.stop < genomic_start
      				pos = genomic_start
      			else
      				pos = genomic_end
      			end
      		end      		 
      	end
      	
      	
      	puts "\tUsing this position for the snoRNA: #{pos}" if verbose
      	raise "This snoRNA was not mapped to a neighbouring exon :( #{snorna.stable_id}" if pos.nil?
      	return peptide_position(pos,verbose)
      
      end
      
      def peptide_position(pos,verbose=false)
      
      	genomic_start,genomic_end = self.coding_region_genomic_start,self.coding_region_genomic_end
        
        this_exon = self.exon_for_genomic_position(pos)        
        raise "not within an exon (#{pos})" if self.exon_for_genomic_position(pos).nil?        
        this_pos = 0
        
        exons = self.exons        
        if self.strand == -1
          first_exon = self.exon_for_genomic_position(genomic_end)                    
          start = genomic_end
          stop = genomic_start          
          puts "\tGene is reverse, using start: #{start}, stop: #{stop}" if verbose          
          if first_exon.stop > start
            puts "\tFirst exon is #{first_exon.start}<>#{first_exon.stop}, start is #{start} (#{start-first_exon.start})"  if verbose
            this_pos += (start-first_exon.start)
          else
            this_pos += first_exon.length
          end          
          exons[1..-1].each do |exon|
          	next if exon == first_exon
            if exon.stop > first_exon.start
              puts "\tSkipping this exon (#{exon.start},#{exon.stop}) - outside of cDNA scope"  if verbose
            elsif exon == this_exon
              puts "\treached target exon #{exon.start}/#{exon.stop} - pos. was #{pos} (#{exon.stop-pos})"  if verbose
              this_pos += (exon.stop-pos)
              return this_pos
            elsif exon.start < this_exon.stop
            	puts "\tAlready found the correct exon, skipping" if verbose
            else
              puts "\tNot the target exon, adding length (#{exon.length} , now #{this_pos+exon.length})" if verbose
              this_pos += exon.length
            end            
          end          
        else     
          first_exon = self.exon_for_genomic_position(genomic_start)          
          start = genomic_start-1
          stop = genomic_end-1         
          puts "\tGene is forward, using start: #{start}, stop: #{stop}"  if verbose          
          if start > first_exon.start
            puts "\tFirst exon is #{first_exon.start}, start is #{start} (#{first_exon.stop-start})" if verbose
            this_pos += (first_exon.stop - start)
          else
            this_pos += first_exon.length
          end          
          if first_exon == this_exon
          	puts "\tThe first exon is the the target exon, returning position (#{this_pos})" if verbose
          	return this_pos
          end          
          exons[1..-1].each do |exon|
          	next if exon == first_exon
            puts "\tprocessing exon: #{exon.start}<->#{exon.stop}" if verbose
            if exon.stop < genomic_start
              puts "\tSkipping this exon (#{exon.start},#{exon.stop}) - outside of cDNA scope"  if verbose
            elsif exon == this_exon
              puts "\treached target exon #{exon.start}/#{exon.stop} - pos. was #{pos} (#{exon.length})"  if verbose
              this_pos += (pos-exon.start)
              return this_pos
            elsif exon.start > this_exon.stop
            	puts "\tAlready found the correct exon, skipping" if verbose
            else  
              puts "\tNot yet the target exon, adding length (#{exon.length})" if verbose
              this_pos += exon.length
            end            
          end              
        end
        return (this_pos)
      end
        
    end
    
    class Translation < DBConnection

      def all_xrefs
        return self.object_xrefs.collect{|ox| ox.xref }
      end

      def get_externaldb_accs(dbname)
        db_id = ExternalDb.find_by_db_name(dbname).external_db_id
        return self.all_xrefs.select{|r| r.external_db_id == db_id}.collect{|r| r.dbprimary_acc}.uniq.join
      end
      
      def get_pfam_features
        return self.protein_features.select{|pf| pf.hit_name.match(/^PF[0-9]*/) }.sort_by{|p| p.score}.reverse
      end

      def snorna_locations        
        answer = []        
        snornas = self.transcript.gene.snornas        
        self.transcript.introns.each do |intron|          
          snornas.each do |snorna|            
            if snorna.slice.within?(intron)
              #puts "#{intron.previous_exon.stable_id}|#{intron.previous_exon.stop} #{intron.start}/#{intron.stop} #{intron.next_exon.start}|#{intron.next_exon.stable_id}"
              if intron.previous_exon.stop > self.transcript.coding_region_genomic_end and self.transcript.strand == -1
                puts "exon partly non-coding, using genomic_coding_region... (#{self.transcript.coding_region_genomic_end} instead of #{intron.previous_exon.stop}), direction is #{self.transcript.strand}"
                answer << self.transcript.coding_region_genomic_end
              else
                answer << intron.previous_exon.stop
              end
            end            
          end          
        end        
        return answer        
      end
      
      def to_fasta
        return Bio::FastaFormat.new(Bio::Sequence::AA.new(self.transcript.protein_seq).to_fasta(self.stable_id))
      end
    end
    
  end
     
end
