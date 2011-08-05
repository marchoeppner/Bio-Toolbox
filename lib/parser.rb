module Toolbox
  
  #include Toolbox::Converter
  
  module Stockholm
  
  	class Report
  		
  		attr_accessor :accession, :rfam_id, :description, :type, :comment, :seeds
  		
  		def initialize(acc)
  			@accession = acc
  			@rfam_id = nil
  			@description = description
  			@type = type
  			@comment = ""
  			@seeds = {}
  		end	
  		
  		def seeds_to_fasta
  			answer = []
  			
  			self.seeds.each do |name,s|
  			
  				answer << Bio::FastaFormat.new(Bio::Sequence::NA.new(s.gsub(/\./, '').strip).to_fasta(name))
  			
  			end
  			
  			return answer
  		
  		end
  	end
  
  	class Parser
  	
  		attr_accessor :infile
  		def initialize(infile)
  				@infile = infile
  		end
  		
  		def run
  		
  			lines = IO.readlines(self.infile)
  			
  			answers = []
  			report = nil
  			
  			lines.each do |line|
  				
  				if line.match(/^#=GF\sAC/)

  					answers << report if report

  					accession = line.strip_stockholm.strip
  					report = Stockholm::Report.new(accession)
  					
  				elsif line.match(/^#=GF\sID/)
  					id = line.strip_stockholm
  					report.rfam_id = id
  					
  				elsif line.match(/^#=GF\sDE/)
  					description = line.strip_stockholm
  					report.description = description
  				
  				elsif line.match(/^#=GF\sTP/)
  					type = line.strip_stockholm
  					report.type = type
  					
  				elsif line.match(/^#=GF\sCC/)
  					comment = line.strip_stockholm.strip
  					report.comment += "#{comment} "
  					
  				elsif line.match(/^[A-Z].*/) 

  					elements = line.gsub(/\s+/, '|').split("|")

  					name = elements.shift.strip
  					seq = elements.shift.strip
  					
  					report.seeds.has_key?(name) ? report.seeds[name] += seq : report.seeds[name] = seq
  					
  				end
  				
  			
  			end
  			
  			return answers
  			
  		end
  		
  	end
  	
  	
  end
  
  class Parser
    
    attr_accessor :infile
    
    def initialize(infile)
      @infile = infile
    end
      
     
  end
    
  module NcRNA
    
    class Parser
      
      require 'rexml/document'
      include REXML
      
      attr_reader :infile, :version, :outfile
      
      def initialize(infile,version="56")
        @infile = infile
        @outfile = "#{infile}_groups_to_sql.txt"
        @version = version
      end
      
      # = DESCRIPTION
      # Converts the conservation XML file
      # into SQL inserts
      def xml_to_sql
       
        input = "#{self.infile}_conservation.xml"

        Toolbox::ToolDB::DBConnection.connect(self.version) 
				
        infile = File.new(input, "r")
        o = File.new("#{self.outfile}", "w+")
        doc = Document.new infile
        
        XPath.each( doc, "//Organism") do |organism|
          unless ToolDB::Organism.exists?(:name => organism.attributes["name"])
            o = ToolDB::Organism.new(:name => organism.attributes["name"])
            o.save
           end
        end

        XPath.each( doc, "//Group") do |group|

          type = group.attributes["biotype"]
          
          o.puts "INSERT INTO grouping (nr,biotype) VALUES (#{group.attributes["nr"]},'#{type}');"
          o.puts "INSERT INTO genomic_align_block (ensembl_galblock_id) VALUES (#{group.attributes["genomic_align_block_id"]});"
          o.puts "INSERT INTO xref_grouping_genomic_align_block (grouping_id,genomic_align_block_id,start,stop) VALUES((SELECT id FROM grouping WHERE nr = #{group.attributes["nr"]} and biotype = '#{type}'),(SELECT id from genomic_align_block WHERE ensembl_galblock_id = #{group.attributes["genomic_align_block_id"]}),#{group.attributes["start"]},#{group.attributes["stop"]});"

          group.elements.each("NcRNA") do |ncrna|

            next if ncrna.attributes["stable_id"] == "nil"

            organism_id = ToolDB::Organism.find_by_name("#{ncrna.attributes["organism"]}").id

            o.puts "INSERT INTO ensembl_ncrna (stable_id,organism_id,external_id,biotype) VALUES ('#{ncrna.attributes["stable_id"]}',#{organism_id},'#{ncrna.attributes["ncrna_id"]}','#{type}');"

            o.puts "INSERT INTO ncrna_group (grouping_id,ensembl_ncrna_id) VALUES((SELECT id FROM grouping WHERE nr = '#{group.attributes["nr"]}' AND biotype = '#{type}'),(SELECT id FROM ensembl_ncrna WHERE stable_id = '#{ncrna.attributes["stable_id"]}'));"

            next unless ncrna.attributes["intronic"] == "true"

            ncrna.elements.each("HostGene") do |gene|

              o.puts "INSERT INTO ensembl_hostgene (stable_id,organism_id,description) VALUES('#{gene.attributes["name"]}',#{organism_id},'#{gene.attributes["description"]}');"

              o.puts "INSERT INTO xref_ncrna_hostgene (ensembl_ncrna_id,ensembl_hostgene_id) VALUES((SELECT id FROM ensembl_ncrna WHERE stable_id = '#{ncrna.attributes["stable_id"]}'),(SELECT id FROM ensembl_hostgene WHERE stable_id = '#{gene.attributes["name"]}'));"

              gene.elements.each("GoId") do |go|
                o.puts "INSERT INTO hostgene_goterm (ensembl_hostgene_id,go_term,linkage) VALUES((SELECT id FROM ensembl_hostgene WHERE stable_id = '#{gene.attributes["name"]}'),'#{go.attributes["id"]}','#{go.attributes["linkage"]}');"
              end

            end

          end

        end
        
      end
      
      # = DESCRIPTION
      # Converts the remaining, non-conserved
      # ncRNAs to SQL inserts
      def all_xml_to_sql(type)
        
        raise "Must have a type" if type.nil?
        
        input = self.infile
        
        o = File.new("#{self.infile}_to_sql.txt", "a")
        
        ToolDB::DBConnection.connect(self.version) 

        infile = File.new(input, "r")
        doc = Document.new infile
        
        XPath.each( doc, "//Organism") do |organism|
          
          puts "#{organism.attributes["name"]}"
          
          organism_id = ToolDB::Organism.find_by_name("#{organism.attributes["name"]}").id
          raise "Organism not reckognized!" if organism_id.nil?
          
          organism.elements.each("NcRNA") do |ncrna|
            
            o.puts "INSERT INTO ensembl_ncrna(stable_id,organism_id,external_id,biotype) VALUES('#{ncrna.attributes["stable_id"]}',#{organism_id},'#{ncrna.attributes["external_id"]}','#{type}');"
            
            ncrna.elements.each("HostGene") do |hostgene|
              
              o.puts "INSERT INTO ensembl_hostgene(stable_id,organism_id,description) VALUES('#{hostgene.attributes["stable_id"]}',#{organism_id},'#{hostgene.attributes["description"]}');"
              
              o.puts "INSERT INTO xref_ncrna_hostgene (ensembl_ncrna_id,ensembl_hostgene_id) VALUES((SELECT id FROM ensembl_ncrna WHERE stable_id = '#{ncrna.attributes["stable_id"]}'),(SELECT id FROM ensembl_hostgene WHERE stable_id = '#{hostgene.attributes["stable_id"]}'));"
              
              hostgene.elements.each("GoId") do |go|
                
                o.puts "INSERT INTO hostgene_goterm (ensembl_hostgene_id,go_term,linkage) SELECT id,'#{go.attributes["id"]}','#{go.attributes["linkage"]}' FROM ensembl_hostgene WHERE stable_id = '#{hostgene.attributes["stable_id"]}';"
                
              end
              
            end
            
          end
          
        end
        
      end
        
      def map_mirnas
        
        ToolDB::DBConnection.connect(self.version)
        warn "Version is not set, using default DB" if self.version == ""
        
        blastdb = "/Users/marc/ncrna/blastdb/hairpin_mirbase_accs.fasta"
        options = "-m 7 -e 0.0001"
        factory = Bio::Blast.local("blastn",blastdb,options)
        
        mirnas = ToolDB::EnsemblNcrna.find_all_by_biotype("miRNA")
        
        mirnas.each do |mirna|
          fasta = mirna.to_fasta
          next if fasta.nil?
          results = factory.query(fasta)
          next if results.hits.empty?
          
          hit = results.hits.first
          
          puts "INSERT INTO xref_ncrna_mirbase(ensembl_ncrna_id,mirna_accession_id,evalue) VALUES((SELECT id FROM ensembl_ncrna WHERE stable_id = '#{mirna.stable_id}'),(SELECT id FROM mirna_accession WHERE acc = '#{hit.target_def}'),#{hit.evalue.to_f});"
          
        end
        
      end
      
      def parse_family_presence
        
        IO.readlines(self.infile).each do |line|
          
          if line.match(/^Char\./)
          
            elements = line.strip.split("\t")
            elements.shift
            @nodes = elements.collect{|n| n.strip}
            @node_counter = {}
            @nodes.each {|n| @node_counter["#{n}"] = 0}
          
          elsif line.match(/^character/)
            
            elements = line.split(/\t/)
            elements.shift
            
            @nodes.each do |node|
              val = elements.shift.strip
              @node_counter["#{node}"] += 1 if "#{val}" == "2"
            end
          
          end
        
        end
        
        @node_counter.sort_by {|n,c| n.to_i }.each do |node,count|
          
          puts "#{node}: #{count}"
        end
        
      end
      
      def parse_dollo_output_do_count(type,version="_54")
        
        valid = false

        @definition = nil
        @dollo = {}
        
        IO.readlines(self.infile).each do |line|
          
          if line.strip.match(/^[A-Za-z0-9_]+(\s+)[0-9]+(\s+)[yesno]/)  # a new node
            
            valid = true
            
            this_line = line.strip
            @definition = this_line.slice!(/^[A-Za-z0-9_]+(\s+)[0-9]+(\s+)[A-Za-z0-9_]+/).strip.split(" ")[1]
            @dollo[@definition] = []
            this_line.strip.gsub(/\s/, '').split(//).each {|c| @dollo[@definition] << c }
            
          elsif line.match(/^$/) and valid == true
            
          elsif line.strip.match(/^[A-Za-z0-9_]+(\s+)[A-Za-z_]+(\s+)[yes]/)
            
            valid = false
            
          elsif valid == true and line.include?("yes") == false
            
            line.strip.gsub(/\s/, '').split(//).each {|c| @dollo[@definition] << c }
            
          end
            
        end
        
        root = ToolDB::DolloNode.find_by_node(1)
        
        results = {}
        
        states = @dollo.fetch("#{root.node}").collect{|s| s.gsub(/\./, '0')}
        results["#{root.node}"] = states
        puts "#{root.node},#{states.select{|s| s == "1"}.nitems}"
        root.get_children.each do |child|
          _dollo_count(child,results,@dollo)
        end
        
      end
      
      def _dollo_count(node,results,dollo)
        
        states = dollo.fetch("#{node.node}")
        results["#{node.node}"] = []
        count = 0
        states.each do |state|
          "#{state}" == "." ? this_state = results.fetch("#{node.get_parent.node}")[count] : this_state = state
          results["#{node.node}"] << this_state
          count += 1
        end
        total = 0
        puts "#{node.node},#{results.fetch(node.node).select {|s| s == "1"}.nitems}"
        
        node.get_children.each do |c| 
          _dollo_count(c,results,dollo)
        end
        
      end
      
      def parse_dollo_output_for_sql(type)
        
        valid = false

        @definition = nil
        @dollo = {}
        
        IO.readlines(self.infile).each do |line|
          
          if line.strip.match(/^[A-Za-z0-9_]+(\s+)[0-9]+(\s+)[yes]/)  # a new node
            
            valid = true
            
            this_line = line.strip
            @definition = this_line.slice!(/^[A-Za-z0-9_]+(\s+)[0-9]+(\s+)[A-Za-z0-9_]+/).gsub(/(\s+)yes/, '')
            @dollo[@definition] = []
            this_line.strip.gsub(/\s/, '').split(//).each {|c| @dollo[@definition] << c.strip }
          elsif line.match(/^$/) and valid == true
            
          elsif line.strip.match(/^[A-Za-z0-9_]+(\s+)[A-Za-z_]+(\s+)[yes]/)
            
            valid = false
            
          elsif valid == true and line.include?("yes") == false
            
            line.strip.gsub(/\s/, '').split(//).each {|c| @dollo[@definition] << c.strip }
            
          end
            
        end
       
        #@dollo.each do |name,states|
        #  puts "#{name}\t#{states.join("<>")}"
        #end
        
        ToolDB::DBConnection.connect(self.version)
        groups = ToolDB::Grouping.find_all_by_biotype(type)
        
        groups.each do |group|
          
          @dollo.each do |node,states|
            this_node = node.strip.split(/\s+/)[1]
            puts "INSERT INTO xref_grouping_dollo_node (grouping_id,dollo_node_id,state) VALUES (#{group.id},(SELECT id FROM dollo_node WHERE node = '#{this_node.strip}'),'#{states.shift}');"
            #puts "\t#{node}: #{states.shift}"
          end
          
        end
      
      end
           
    end
    
  end
  
end