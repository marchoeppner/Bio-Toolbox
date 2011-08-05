class PairwiseCount
  
  attr_accessor :organisms, :count
  
  def initialize(organisms)
    @organisms = organisms
    @count = 0
  end
  
end

module Toolbox
  
  include ToolDB
  
  module DBParser
    
    class Summary
      
      attr_reader :type, :organisms, :base_name, :version
      
      def initialize(type,base_name,version="_54")
        @type = type
        ToolDB::DBConnection.connect(version)
        @organisms = Hash.new
        ToolDB::Organism.find(:all).each do |organism|
          @organisms[organism.name] = [ 0, 0 ,0 ,0]
        end
        @base_name = base_name
        @version = version
      end
      
      def summarize
        
        ncrnas = ToolDB::EnsemblNcrna.find_all_by_biotype(self.type)
        
        ncrnas.each do |ncrna|
          
          organism = ncrna.organism
          data = @organisms.fetch(organism.name)
          
          full = 1
          full_intronic = 0
          conserved = 0
          conserved_intronic = 0
          
          full_intronic = 1 if ncrna.has_hostgene?
          conserved = 1 if ncrna.has_group?
          conserved_intronic = full_intronic if ncrna.has_group?
          @organisms[organism.name] = [ data[0] += full , data[1] += full_intronic , data[2] += conserved , data[3] += conserved_intronic ]
          
        end
        
        o = File.new("ensembl_52_12_vertebrates_#{self.type}_summary.csv", "a")
        o.puts "Organism,Total,Intronic,Conserved,Conserved Intronic"
        @organisms.each do |organism,counts|
          o.puts "#{organism},#{counts[0]},#{counts[1]},#{counts[2]},#{counts[3]}"
        end
        o.close
      end
      
      def pairwise
        
        organisms = Organism.find(:all,:order=>"name")
        organism_list = organisms.collect{|o| o.name }.sort
        pairs = []
        
        until organism_list.empty?
          curr = organism_list.shift
          organism_list.each { |o| pairs.push(PairwiseCount.new([curr , o ])) }
        end
            
        groups = Grouping.find_all_by_biotype(type)
        
        groups.each do |group|
          
          species = group.fetch_organisms.collect{|o| o.name }
          
          pairs.each do |pair|
            pair.count += 1 if species.include?(pair.organisms[0]) and species.include?(pair.organisms[1])
          end
          
        end
      
      f = File.new("#{base_name}_#{type}_pairwise_counts.csv", "a")
       f.puts ",#{organisms.collect{|o| o.name}.join(",")}"
       organisms.each do |organism|
         f.print "#{organism.name},"
         organisms.each do |o|
           f.print "," if o.name == organism.name
           next if o.name == organism.name 
           pairs.each do |pair|
             f.print "#{pair.count}," if pair.organisms.include?(organism.name) and pair.organisms.include?(o.name)
           end
         end
         f.print "\n"
         
        end
      f.close
      end
      
      def summarize_detailed
        
        organisms = ToolDB::Organism.find(:all, :order => "name")
        
        organisms.each do |organism|
          
          external = Hash.new
          external_hostgenes = Hash.new
          
          ncrnas = organism.ensembl_ncrnas.select{|n| n.has_group? and n.biotype == self.type }
          
          ncrnas.each do |ncrna|
            
            ncrna.biotype == "miRNA" ? ext = ncrna.mirbase_family : ext = ncrna.external_id
            
            external.has_key?("#{ext}") ? external["#{ext}"] = external.fetch("#{ext}") + 1 : external["#{ext}"] = 1
            
            next unless ncrna.is_intronic?
            
            external_hostgenes.has_key?("#{ext}") ? external_hostgenes["#{ext}"] = external_hostgenes.fetch("#{ext}").push(ncrna.fetch_hostgene) : external_hostgenes["#{ext}"] = [ ncrna.fetch_hostgene ]
            
          end
          
          f = File.new("#{self.base_name}_#{self.type}_summary_detailed_#{organism.name}.csv", "a")
          f.puts "External ID,Total,Intronic"
          external.each do |acc,count|
            f.print "#{acc},#{count}"
            f.print ",#{external_hostgenes.fetch(acc).nitems}" if external_hostgenes.has_key?(acc)
            f.print "," 
            if self.type == "miRNA"
              mirbase = ToolDB::MirnaFamily.find_by_name("#{acc}")
              puts "#{acc}"
              f.print "#{mirbase.name}"
            else
              rfam = Toolbox::ToolDB::Rfam.find_by_acc(acc)
              f.print "#{rfam.name},#{rfam.family},#{rfam.targets}" unless rfam.nil?
            end
            f.print "\n"
          end
          f.close
        end
        
      end
      
      def create_mesquite_file(tree,version="_54")
        
        raise "No tree file specified" if tree.nil?
        
        ToolDB::DBConnection.connect(version)

        groups = ToolDB::Grouping.find_all_by_biotype("#{type}")
        organisms = ToolDB::Organism.find(:all, :order => "name")
        answer = {}
        organisms.each{ |o| answer[o.name] = "" }
        
        groups.each do |group|
          
          organisms.each do |organism|
            present = false
            group.ensembl_ncrnas.each do |ncrna|
              present = true if ncrna.organism.name == organism.name
            end
            present == false ? answer[organism.name] += "0" : answer[organism.name] += "1"
          end
          
        end
        
        nex = Toolbox::Nexus::Writer.new("#{type}",answer,tree)
        f = File.new("#{base_name}_#{type}_mesquite.nex", "a")
        f.puts nex.to_s
        f.close
        
      end
      
      def go_stat(organism_id,node,seqtype=nil)
        
        require 'ensembl'
        
        ToolDB::DBConnection.connect(self.version)
        
        raise "Wrong Node class" if node.kind_of?(Integer)
        
        if "#{node.node}" == "2"
          hostgenes = ToolDB::EnsemblHostgene.find_all_by_organism_id(organism_id).select{|h| h.present_at?(self.type,node)}
        else
          node2 = ToolDB::DolloNode.find_by_node(2)
          hostgenes = ToolDB::EnsemblHostgene.find_all_by_organism_id(organism_id).select{|h| h.present_at?(self.type,node) and h.present_at?(self.type,node2) == false}
        end
        ToolDB::Organism.find(organism_id).ensembl_db_connection(54) unless seqtype == "gene"
        
        answer = []
        
        hostgenes.each do |hgene|
          
          if seqtype == "gene"
            answer.push(hgene.stable_id)
          elsif seqtype == "transcript"
            gene = hgene.fetch_ensembl_gene(54)
            t = []
            gene.transcripts.each do |transcript|
              next if transcript.translation.nil?
              t << transcript.translation.stable_id
            end
            answer << t.uniq.shift
          elsif seqtype == "uniprot"
            gene = hgene.fetch_ensembl_gene(54)
            answer.push(gene.uniprot_id)
          end
      
        end
        
        answer.each{|a| puts a}
        
      end
      
      def html_gene_list(organism_id,node,node_2=nil)
        
        raise "Must have at least one node" if node.nil?
        
        ToolDB::DBConnection.connect(self.version)

        node = ToolDB::TreeNode.find_by_nr(node)
        node_2 = ToolDB::TreeNode.find_by_nr(node_2) unless node_2.nil?

        hostgenes = ToolDB::EnsemblHostgene.find_all_by_organism_id(organism_id).select{|h| h.present_at?(self.type,node)} if node_2.nil?
        hostgenes = ToolDB::EnsemblHostgene.find_all_by_organism_id(organism_id).select{|h| h.present_at?(self.type,node) and h.present_at?(self.type,node_2) == false } unless node_2.nil?
        
        ToolDB::Organism.find(organism_id).ensembl_db_connection(self.version)
        
        puts "Genes present at Node #{node.nr}/#{node_2.nr}"

        puts "<Table>"
          puts "<tr>"
          puts "<td>Gene Stable ID</td><td>Description</td><td>NcRNA(s)</td><td>PFAM domains</td>"
          puts "</tr>"

        hostgenes.each do |hgene|

          puts "<tr>"

          puts "<td><a href=\"http://www.ensembl.org/#{hgene.organism.name.gsub(/\s/,'_')}/Gene/Summary?g=#{hgene.stable_id}\">#{hgene.stable_id}<\/a></td>"
          
          gene = hgene.fetch_ensembl_gene
          
          if gene.description.nil?
            puts "<td></td>"
          else  
            uniprot = gene.description.slice(/Acc\:.*$/).gsub(/Acc\:/, '').gsub(/\]/, '') unless gene.description.nil?
            puts "<td>#{gene.description} (<a href=\"http://www.uniprot.org/uniprot/#{uniprot}\">#{uniprot}<\/a>)</td>"
          end
          puts "<td>"
          
          important_node = node if node_2.nil?
          important_node = node_2 unless node_2.nil?
          
          hgene.fetch_ncrnas(self.type).select{|n| n.present_at?(important_node)}.each do |ncrna|
            if self.type == "snoRNA"
              puts "<a href=\"http://www.ensembl.org/#{hgene.organism.name.gsub(/\s/,'_')}/Gene/Summary?g=#{ncrna.stable_id}\">#{ncrna.stable_id}<\/a> (<a href=\"http://rfam.janelia.org/cgi-bin/getdesc?acc=#{ncrna.external_id}\">#{ncrna.external_id}<\/a>),"
            else 
              print "<a href=\"http://www.ensembl.org/#{hgene.organism.name.gsub(/\s/,'_')}/Gene/Summary?g=#{ncrna.stable_id}\">#{ncrna.stable_id}<\/a>"
              print " (<a href=\"http://microrna.sanger.ac.uk/cgi-bin/sequences/mirna_entry.pl?acc=#{ncrna.mirna_accession.acc}\">#{ncrna.mirna_accession.acc}<\/a>)" unless ncrna.mirna_accession.nil?
              print ","
            end
          end
          print "\n"
          puts "</td>"

          gene.get_pfam_features.each do |pfam|
            puts "<td><a href=\"http://pfam.sanger.ac.uk//family/#{pfam}\">#{pfam}<\/a></td>"
          end
          puts "</tr>"

        end

        puts "</table>"
        
      end
      
      def create_dollo_file
        
        groups = ToolDB::Grouping.find_all_by_biotype(self.type)
        organisms = ToolDB::Organism.find(:all, :order => "name")
        
        species_states = {}
        organisms.each {|o| species_states["#{o.name}"] = "" }
        
        groups.each do |group|
          covered = group.ensembl_ncrnas.collect{|n| "#{n.organism.name}"}
          organisms.each do |organism|
            covered.include?("#{organism.name}") ? species_states["#{organism.name}"] += "1" : species_states["#{organism.name}"] += "0"
          end
        end
        puts "\t#{organisms.nitems}\t#{groups.nitems}"
        species_states.each do |org,states|
          puts "#{org.gsub(/\s/, '_')[0..9]} #{states}"
        end
        
      end
      
         
    end
    
    class Ancestral
      
      attr_reader :type, :base_name, :version
      
      def initialize(type,base_name,version="_54")
        raise "Format error" unless type == "miRNA" or type == "snoRNA"
        @type = type
        @base_name = base_name
        @version = version
        ToolDB::DBConnection.connect(version)
      end
      
      # = DESCRIPTION
      # Runs all currently implemented analyses
      def run_all
        
        puts "Running analysis for: #{self.base_name} | #{self.type}"
        puts "Summarzing nodes..."
        self.summarize_nodes
        puts "Finding intronic groups..."
        self.intronic_groups
        puts "Listing #{self.type} count at internal nodes..."
        self.ncrnas_per_node
        puts "Summarizing family distribution at internal nodes..."
        self.family_summary
        puts "Summarizing family presence/absence (local)"
        self.family_presence("local")
        puts "Summarizing host genes at internal nodes..."
        self.hostgene_at_nodes
        puts "...Done!"
        
      end
          
      def summarize_nodes
        
        nodes = ToolDB::TreeNode.find(:all, :order => "nr")
        
        nodes_counter = Hash.new
        nodes_counter_half = Hash.new
        nodes_stable_intronic = Hash.new
        
        nodes.each { |n| nodes_counter["#{n.nr}"] = 0 }
        nodes.each { |n| nodes_counter_half["#{n.nr}"] = 0 }
        nodes.each { |n| nodes_stable_intronic["#{n.nr}"] = 0 }
        
        groups = ToolDB::Grouping.find_all_by_biotype(self.type)
        
        groups.each do |group|
          
          group.group_nodes.each do |node|
            
            node_nr = "#{node.tree_node.nr}"
            nodes_counter["#{node_nr}"] = nodes_counter.fetch("#{node_nr}") + 1 if "#{node.val}" == "2"
            nodes_stable_intronic["#{node_nr}"] = nodes_stable_intronic.fetch("#{node_nr}") + 1 if "#{node.val.strip}" == "2" and group.is_intronic?
            nodes_counter_half["#{node_nr}"] = nodes_counter_half.fetch("#{node_nr}") +  1 if "#{node.val.strip}" == "1 2"
            
          end
          
        end
        
        f = File.new("#{@base_name}_node_count.csv", "a")
        f.puts "node,total(definitive),candidates,intronic (definitive)"
        nodes_counter_half.each do |node,count|
          f.puts "#{node},#{nodes_counter.fetch("#{node}")},#{count},#{nodes_stable_intronic.fetch("#{node}")}"
        end
        f.close
             
      end
    
      def summarize_nodes_intronic
        
        answer = Hash.new
        nodes = ToolDB::TreeNode.find(:all,:order => "nr")
        nodes.each { |n| answer["#{n.nr}"] = [] }
        groups = ToolDB::Grouping.find_all_by_biotype(self.type).select{|g| g.is_intronic?}
        
        groups.each do |group|
          
          nodes.each do |node|
            next unless group.present_at?(node)
            answer["#{node.nr}"] = answer.fetch("#{node.nr}").push(group)
          end
        
        end
        
        return nodes
        
      end
      
      # = DESCRIPTION
      # Finds all ncRNA groups where at least on
      # gene is located within an intron
      def intronic_groups
      
        answer = Hash.new
        
        groups = ToolDB::Grouping.find_all_by_biotype(self.type).select{|g| g.is_intronic? }
        nodes = ToolDB::TreeNode.find(:all, :order => "nr")
        nodes.each { |n| answer["#{n.nr}"] = [] }
        
        groups.each do |group|
          
          nodes.each do |node|
            
            next unless group.present_at?(node)
            answer["node.nr"] = answer.fetch("#{node.nr}").push(group)
          
          end
          
        end
        
        # RUN CVS compiler #
        _intronic_at_nodes(answer)
        
      end
      
      # = DESCRIPTION
      # Counts ncRNAs at internal node based on
      # data from Phylip Dollop
      def ncrnas_per_node_dollo

        groups = ToolDB::Grouping.find_all_by_biotype(self.type)

        states_total = {}
        states_intronic = {}

        nodes = ToolDB::DolloNode.find(:all, :order => "node")

        nodes.each {|n| states_total[n.node] = 0}
        nodes.each {|n| states_intronic[n.node] = 0}

        groups.each do |group|

          next if group.ensembl_ncrnas.nitems == 1

          nodes.each do |node|

            state = group.dollo_state_at_node(node.node).to_i
            states_total[node.node] += state
            states_intronic[node.node] += state if group.is_intronic?

          end

        end
        puts "Count of #{self.type}s at internal nodes (dollo)"
        puts "node,total,intronic"
        states_total.each do |node,count|

          puts "#{node},#{count},#{states_intronic.fetch(node)}"

        end
        
      end
      
      def ncrnas_per_node
        
        groups = ToolDB::Grouping.find_all_by_biotype(self.type)
        nodes = ToolDB::TreeNode.find(:all, :order => "nr")
        organisms = ToolDB::Organism.find(:all, :order => "name")
        
        nodes.each do |node|
          f = File.new("#{@base_name}_node_#{node.nr}.csv", "a")
          f.print ","
          organisms.each {|organism| f.print ",#{organism.name}"}
          f.print "\n"
          f.close
        end
          
        groups.each do |group|
          nodes.each do |node|
            if group.present_at?(node)
              ncrnas = group.fetch_ncrnas
              f = File.new("#{@base_name}_node_#{node.nr}.csv", "a")
              ncrna_group = "#{group.fetch_ncrna_group}"
              f.print "#{group.nr},#{ncrna_group}"
              organisms.each do |organism|
                state = "0"
                 ncrnas.each do |ncrna|
                   state = "X" if ncrna.organism.name == organism.name
                end
                f.print ",#{state}"
              end
              f.print ",intronic," if group.is_intronic?
              f.print ",intergenic," if group.is_intronic? == false
              if self.type == "miRNA"
                f.print "#{ncrna_group}"
              else
                rfam = ToolDB::Rfam.find_by_acc("#{ncrna_group}")
                f.print "#{rfam.name},#{rfam.family},#{rfam.targets},"
              end
              ncrnas.each do |ncrna|
                hostgene = ncrna.fetch_hostgene
                next if hostgene.nil?
                f.print "#{hostgene.stable_id}|#{hostgene.genbank_id}|#{hostgene.description.gsub(/\,/, ' ').strip}," 
              end
              f.print "\n"
              f.close
              
            end
          end  
        end  
          
      end
      
      # = DESCRIPTION
      # Checks each node for host genes and
      # lists their GO terms and expression profiles
      def nodes_expression
        
        ToolDB::DBConnection.connect

        nodes = ToolDB::TreeNode.find(:all, :order => "nr")

        nodes.each do |node|
          
          f = File.new("#{self.base_name}_node_#{node.nr}_expression.csv")
          node.fetch_present_groups("snoRNA").each do |group|

            puts "Group: #{group.nr}"
            group.fetch_hostgenes(1).each do |hostgene|
              puts ",#{hostgene.stable_id},#{hostgene.description}"
              puts
              puts ",GO terms"
              hostgene.fetch_goterms.each do |go_term|
                puts ",,#{go_term.name}"
              end
              puts
              puts ",Top tissues"
              hostgene.expression_profiles(5).each do |profile|
                puts ",,#{profile.sample_name}:#{profile.value}"
              end
            end

          end
          f.close
        end
        
      end
      
      # = DESCRIPTION
      # Returns stats for node-specific groups
      def dollo_nodes_stats
        
        nodes = ToolDB::DolloNode.find(:all, :order => "node")
        groups = ToolDB::Grouping.find_all_by_biotype(self.type)

        puts "Node,Group_nr,Seqs,perc_similarity,perc_gaps,gerp(mean)"

        nodes.each do |node|

          if node.get_parent.nil?
            g = groups.select{|g| g.present_at?(node,"dollo")}
          else
            g = groups.select{|g| g.present_at?(node,"dollo") and g.present_at?(node.get_parent,"dollo") == false }
          end

          g.each do |group|

            next if group.ensembl_ncrnas.nitems == 1
            aln = group.do_align
            puts "#{node.node},#{group.nr},#{group.ensembl_ncrnas.nitems},#{aln.percent_similarity},#{aln.gap_score},#{group.gerp_mean(1)}"

          end

        end
        
      end
      
      # = DESCRIPTION
      # internal
      def _intronic_at_nodes(organism_hash)
        
        organism_hash.each do |node,groups|

          f = File.new("#{@base_name}_#{node}_intronic.csv", "a")

          organisms = ToolDB::Organism.find(:all, :order => "name")
          f.print ","
          organisms.each { |o| f.print ",#{o.name}" }
          f.print "\n"

          groups.each do |group|

            f.print "#{group.nr}"
            external_id = group.fetch_ncrna_group
            f.print ",#{external_id}"
            organisms.each do |organism|
              state = ""
              group.ensembl_ncrnas.each do |ncrna|
                if ncrna.organism.name == organism.name
                  state = "X"
                end
              end
              f.print ",#{state}"
            end
            if "#{self.type}" == "snoRNA" or "#{self.type}" == "snRNA"
              rfam = ToolDB::Rfam.find_by_acc("#{external_id}")
              f.print ",#{rfam.name},#{rfam.family},#{rfam.targets},"
            else
              f.print "," + external_id + ","
            end
            group.fetch_hostgenes.each { |h| f.print "#{h.description.gsub(/\,/, ' ')}," }
            
            f.print "\n"

          end

        end
      end
      
      # = DESCRIPTION
      # Iterates over all npcRNA groups (within sno- and miRNAs) and lists counts per group and internal node
      # as well as the total count per genome for a given group.
      def family_summary
        
        @groups_total = Hash.new
        @groups_conserved = Hash.new
        
        organisms = ToolDB::Organism.find(:all, :order => "name")
        
        nodes = ToolDB::TreeNode.find(:all)
                
        groups = ToolDB::Grouping.find_all_by_biotype(self.type)
           
        groups.each do |group|
          
          ncrna_group = group.fetch_ncrna_group
          
          if @groups_conserved.has_key?("#{ncrna_group}")
            
            data = @groups_conserved.fetch("#{ncrna_group}")
            counter = 0
            nodes.each do |node|
              data[counter] += 1 if group.present_at?(node)
              counter += 1
            end
            
          else
            
            data = []
            nodes.each do |node|
              group.present_at?(node) ? data.push(1) : data.push(0)
            end
            @groups_conserved["#{ncrna_group}"] = data
            
          end
          
        end  
        
        f = File.new("#{@base_name}_#{self.type}_group_nodes.csv", "a")
        nodes.each { |n| f.print ",#{n.nr}"}
        f.print "," + organisms.collect{|o| o.name}.join(",")
        f.print "\n"
        @groups_conserved.each do |group,counts|
          if self.type == "miRNA"
            f.print "#{group},#{counts.join(",")},"
            organisms.each do |organism|
              f.print ToolDB::EnsemblNcrna.find_all_by_biotype_and_organism_id(self.type,organism.id).select{|n| n.mirbase_family == group }.nitems.to_s + ","
            end
            f.print "\n"
          else
            rfam = ToolDB::Rfam.find_by_acc("#{group}")
            f.print "#{group},#{counts.join(",")},#{rfam.name},#{rfam.family.strip},#{rfam.targets.strip},"
            organisms.each do |organism|
              f.print ToolDB::EnsemblNcrna.find_all_by_organism_id_and_biotype_and_external_id(organism.id,self.type,group).nitems.to_s + ","
            end
            f.print "\n"
          end
        end
        f.close  
      end
      
      # = DESCRIPTION
      # Identifies all host genes
      # and the number of ncRNAs at each internal node
      def hostgene_at_nodes
        
        genes = ToolDB::EnsemblHostgene.find(:all)
        nodes = ToolDB::TreeNode.find(:all)
        gene_nodes = Hash.new
        hostgene_ncrnas = Hash.new
        
        genes.each do |gene|
          
          data = []
          nodes.each { |n| data.push(0)}
          
          ncrnas = gene.fetch_ncrnas(self.type)
          
          unless ncrnas.nil?
            output = false
            ncrnas.each do |ncrna|
              if ncrna.has_group?     # Only consider ncRNAs that are part of the genome alignment
                counter = 0
                nodes.each do |node|
                  if ncrna.present_at?(node)
                    data[counter] += 1
                    output = true
                  end
                  counter += 1
                end
              end
            end
            gene_nodes["#{gene.stable_id}"] = data if output == true
          end
          
        end
        
        organisms = ToolDB::Organism.find(:all, :order => "name")
        
        organisms.each do |organism|
          f = File.new("#{@base_name}_hostgenes_through_time_#{organism.name}.html", "a")
          f.puts "<html>"
          f.puts "<body>"
          f.puts "<table>"
          f.puts "\t<tr>"
          f.print "<td></td><td></td>"
          nodes.each { |n| f.print "<td>Node #{n.nr}</td>"}
          f.print "\n"
          f.puts "\t</tr>"
          f.close
        end
          
        gene_nodes.each do |gene,node_list|
          hostgene = ToolDB::EnsemblHostgene.find_by_stable_id("#{gene.strip}")
          organism = hostgene.organism
          f = File.new("#{@base_name}_hostgenes_through_time_#{organism.name}.html", "a")
          f.puts "\t<tr bgcolor='lightgrey'>"
          f.print "\t\t<td>#{gene}</td><td></td>"
          node_list.each do |node| 
            f.print "<td>#{node}</td>"
          end
          f.print "<td>#{hostgene.description}...</td>"
          f.print "\n"
          f.puts "\t</tr>"
          hostgene.fetch_ncrnas(self.type).each do |ncrna|
            if ncrna.has_group?
              state = []
              nodes.each do |node|
                if ncrna.present_at?(node)
                  state.push("X")
                else
                  state.push("0")
                end
              end
              f.puts "\t<tr>"
              f.print "\t\t<td></td><td>#{ncrna.external_id}</td>"
              state.each do |s|
                f.print "<td>#{s}</td>"
              end
              f.print "\n"
              f.puts "\t</tr>"
            end
          end
          f.close
        end
        organisms.each do |organism|
          f = File.new("#{@base_name}_hostgenes_through_time_#{organism.name}.html", "a")
          f.puts "</table>"
          f.puts "</body>"
          f.puts "</html>"
          f.close
        end
      end
      
      # = DESCRIPTION
      # Creates a dollop compatible file of present families
      def dollo_family_presence(scope)
        
        organisms = ToolDB::Organism.find(:all, :order => "name")
        
        if self.type == "snoRNA"
          ncrna_families = ToolDB::EnsemblNcrna.find_all_by_biotype("snoRNA").collect{|n| "#{n.external_id}"}.uniq.select{|e| e.length > 1}.sort 
        elsif self.type == "miRNA"
          ncrna_families = ToolDB::EnsemblNcrna.find_all_by_biotype("miRNA").collect{|n| "#{n.mirbase_family.downcase}" }.uniq.select{|e| e.length > 1}.sort if self.type == "miRNA"
        end
        
        results = {}
        
        organisms.each do |organism|
          
          results["#{organism.phylip_name}"] = ""
          
          if self.type == "miRNA"
            families = ToolDB::EnsemblNcrna.find_all_by_biotype_and_organism_id("#{self.type}",organism.id).collect{|nc| "#{nc.mirbase_family.downcase}"}.uniq if scope == "global"
            families = ToolDB::EnsemblNcrna.find_all_by_biotype_and_organism_id("#{self.type}",organism.id).select{|n| n.has_group?}.collect{|nc| "#{nc.mirbase_family.downcase}"}.uniq if scope == "local"
          else
            families = ToolDB::EnsemblNcrna.find_all_by_biotype_and_organism_id("#{self.type}",organism.id).collect{|nc| "#{nc.external_id}"}.uniq if scope == "global"
            families = ToolDB::EnsemblNcrna.find_all_by_biotype_and_organism_id("#{self.type}",organism.id).select{|n| n.has_group?}.collect{|nc| "#{nc.external_id}"}.uniq if scope == "local"
          end
          
          ncrna_families.each do |fam|
            families.include?("#{fam}") ? state = 1 : state = 0
            results["#{organism.phylip_name}"] += "#{state}"
          end
          
        end
        
        puts "\t#{results.keys.nitems}\t#{results.to_a[0][1].length}"
        results.each do |o,states|
          puts "#{o}\s#{states}"
        end
        
      end
        
      # = DESCRIPTION
      # Creates a presence/absence matrix for each family (all ncRNAs!)
      # Scope can be "global" or "local" (conserved groups or all annotated
      # ncRNAs). Method can be 'dollo' or 'mesquite'.
      def family_presence(scope,tree)
        
        raise "Must have a scope (local/global)" if scope.nil?
        raise "Requires a tree file!" if tree.nil?
        
        organisms = ToolDB::Organism.find(:all, :order => "name")
        
        puts "#NEXUS"
        puts
        puts "BEGIN TAXA;"
        puts"\tTITLE #{self.type};"
        puts "\tDIMENSIONS NTAX=#{organisms.nitems};"
        puts "\tTAXLABELS"
        puts "\t#{organisms.collect{|o| o.name.gsub(/\s/, '_')}.join(" ")}"
        puts "\t;"
        puts "END;"
        
        ncrna_families = ToolDB::EnsemblNcrna.find_all_by_biotype("miRNA").collect{|n| "#{n.mirbase_family.downcase}" }.uniq if self.type == "miRNA"
        ncrna_families = ToolDB::EnsemblNcrna.find_all_by_biotype("snoRNA").collect{|n| "#{n.external_id}"}.uniq if self.type == "snoRNA"
        
        puts "BEGIN CHARACTERS;"
        puts "\tTITLE #{self.type};"
        puts "\tDIMENSIONS NCHAR=#{ncrna_families.nitems};"
        puts "\tFORMAT DATATYPE = STANDARD GAP = ? MISSING = - SYMBOLS = \" - 0 1 \";"
        puts "\tCHARSTATELABELS"
        counter = 1
        ncrna_families.each do |fam|
          print ", #{counter} #{fam.gsub(/\-/, '')}"
          counter += 1
        end
        print ";\n"
        puts "\tMATRIX"
        
        organisms.each do |organism|
          
          print "#{organism.name.gsub(/\s/, '_')}\t"
          
          if self.type == "miRNA"
            families = ToolDB::EnsemblNcrna.find_all_by_biotype_and_organism_id("#{self.type}",organism.id).collect{|nc| "#{nc.mirbase_family.downcase}"}.uniq if scope == "global"
            families = ToolDB::EnsemblNcrna.find_all_by_biotype_and_organism_id("#{self.type}",organism.id).select{|n| n.has_group?}.collect{|nc| "#{nc.mirbase_family.downcase}"}.uniq if scope == "local"
          else
            families = ToolDB::EnsemblNcrna.find_all_by_biotype_and_organism_id("#{self.type}",organism.id).collect{|nc| "#{nc.external_id}"}.uniq if scope == "global"
            families = ToolDB::EnsemblNcrna.find_all_by_biotype_and_organism_id("#{self.type}",organism.id).select{|n| n.has_group?}.collect{|nc| "#{nc.external_id}"}.uniq if scope == "local"
          end
          
          ncrna_families.each do |fam|
            families.include?("#{fam}") ? state = 1 : state = 0
            print "#{state}"
          end
          puts
        end
        puts ";"
        puts "END;"
        
        t = Converter::Trees.new(tree)
        puts t.to_nexus
        
      end
      
      # = DESRIPTION
      # Compares the copy number of snoRNA families
      # between nodes
      def family_copies(nodes)
        
        families = {}

        ToolDB::EnsemblNcrna.find_all_by_biotype(type).each { |n| families["#{n.external_id}"] = [ 0, 0, 0]} if type == "snoRNA"

        pos = 0

        nodes.each do |node|

          n = ToolDB::TreeNode.find_by_nr(node)
          puts n.inspect

          groups = ToolDB::Grouping.find_all_by_biotype(type).select{|g| g.present_at?(n)}

          groups.each do |group|

            families["#{group.fetch_ncrna_group}"][pos] += 1

          end

          pos += 1

        end

        puts "#{type} family,#{nodes.join(",")}"
        families.sort.each do |fam,counts|
          puts "#{fam},#{counts.join(",")}"
        end
        
      end
      
      def expression_at_nodes(node_nr,organism_id=1)
        
        tissues = ToolDB::ExpressionTissue.fetch_tissues
        tissue_hash = {}
        tissues.collect{|t| tissue_hash[t.name] = [] }

        basal_node = ToolDB::TreeNode.find_by_nr(2)
        node_nr == 2 ? node = basal_node : node = ToolDB::TreeNode.find_by_nr(node_nr)
        
        if node_nr == 2
          hostgenes = ToolDB::EnsemblNcrna.find_all_by_biotype_and_organism_id(self.type,organism_id).select{|n| n.present_at?(node) and n.is_intronic? }.collect{|n| n.fetch_hostgene}.uniq
        else
          hostgenes = ToolDB::EnsemblNcrna.find_all_by_biotype_and_organism_id(self.type,organism_id).select{|n| n.present_at?(node) and n.is_intronic? and n.present_at?(basal_node) == false}.collect{|n| n.fetch_hostgene}.uniq
        end
        
        hostgenes.each do |gene|

          expression = gene.expression_levels

          tissue_hash.each do |tissue,values|
            present = false
            expression.each do |exp|
              if exp.expression_tissue.name == tissue
                exp.value == 0.0 ? tissue_hash[tissue].push("") : tissue_hash[tissue].push(exp.value)
                present = true
              end
            end
            tissue_hash[tissue].push("") if present == false
          end


        end

        tissue_hash.each do |name,values|
          puts "#{name},#{values.join(",")}"
        end
        
      end
          
      def go_terms_at_nodes
        
        genes = ToolDB::EnsemblHostgene.find(:all)
        tree_nodes = ToolDB::TreeNode.find(:all)
        nodes = Array.new
        tree_nodes.each { |tn| nodes.push("#{tn.nr}")}
        
        organisms = ToolDB::Organism.find(:all)
        
        # Iterate over all organisms #
        organisms.each do |organism|
          
          node_genes = Hash.new
          nodes.each { |n| node_genes["#{n}"] = [] }
          
          genes = ToolDB::EnsemblHostgene.find_all_by_organism_id(organism.id)
          
          #Check all genes from that organism
          genes.each do |gene|
          
            if gene.has_ncrna?(self.type)
              
              gene_at_nodes = Hash.new
              
              nodes.each { |n| gene_at_nodes["#{n}"] = 0 }
              
              ncrnas = gene.fetch_ncrnas(self.type)
              
              ncrnas.each do |ncrna|
              
                if ncrna.has_group?
                  ncrna.group.group_nodes.each do |gr_node|
                    if "#{gr_node.val}" == "2"
                      gene_at_nodes["#{gr_node.tree_node.nr}"] = 1
                    end
                  end
                end
              end
              
              gene_at_nodes.each do |node,state|
                if state == 1
                  node_genes["#{node}"] = node_genes.fetch("#{node}").push("#{gene.stable_id}") unless node_genes.fetch("#{node}").include?("#{gene.stable_id}")        
                end
                
              end
              
            end
            
          end
          
          go_count = Hash.new
          
          node_genes.each do |node,genes|
            
            genes.each do |gene|
              
              ens_gene = ToolDB::EnsemblHostgene.find_by_stable_id("#{gene}")
              
              ens_gene.fetch_goterms.each do |term|
                
                if go_count.has_key?(term.go_term)
                  value = go_count.fetch(term.go_term)
                  value[nodes.index(node)] = value[nodes.index(node)] + 1
                  go_count[term.go_term] = value
                else
                  values = Array.new
                  nodes.each { |n| values.push(0) }
                  values[nodes.index("#{node}")] = values[nodes.index("#{node}")] + 1
                  go_count["#{term.go_term}"] = values
                end
                
              end
              
            end
            
          end
          
          hierarchy = self.go_hierarchy
          hierarchy.each do |entry|
            if go_count.has_key?("#{entry.strip}")
              puts "#{entry}\t#{go_count.fetch(entry.strip).join(",")}"
            end
          end
          
          f = File.new("#{@base_name}_goterms_#{organism.name}.csv", "a")
          f.puts ",#{nodes.join(",")}"
         
          go_count.each do |term,counts|

            go_term = ToolDB::GoTerm.find_by_term("#{term}")
            if go_term.nil?
              go_term = ""
            else
              go_term = go_term.name.gsub(/\,/, ' - ')
            end
            
            f.puts "#{term},#{counts.join(",")},#{go_term}"
          end
          f.close
        end
            
      end
      
      def node_diff(node1,node2)
        
        # Which groups are present
        
        node1_groups = Array.new
        node2_groups = Array.new
        
        node_1 = ToolDB::TreeNode.find_by_nr("#{node1}")
        node_1.fetch_present_groups("#{self.type}").each do |group|
          node1_groups.push(group.id)
        end
        
        node_2 = ToolDB::TreeNode.find_by_nr("#{node2}")
        node_2.fetch_present_groups("#{self.type}").each do |group|
          if node1_groups.include?(group.id)
            node1_groups.delete(group.id)
          else
            node2_groups.push(group.id)
          end
        end
        
        puts "node1 #{node1_groups.nitems}"
        puts "node2 #{node2_groups.nitems}"
        
        nodes_groups = [ node1_groups , node2_groups ]
        
        # Grouping by external IDs #
        ############################
        counter = 0
        @rfam_counter = Hash.new
        
        nodes_groups.each do |node_groups|
          node_groups.each do |group|
            group = ToolDB::Grouping.find(group)
            ncrna_group = group.fetch_ncrna_group
            if @rfam_counter.has_key?("#{ncrna_group}")
              values = @rfam_counter.fetch("#{ncrna_group}")
              values[counter] = values[counter] + 1
              @rfam_counter["#{ncrna_group}"] = values
            else
              values = [ 0 , 0 ]
              values[counter] = values[counter] + 1
              @rfam_counter["#{ncrna_group}"] = values
            end
          end
          counter += 1
        end
        
        f = File.new("#{@base_name}_comparison_accessions_node#{node1}_node#{node2}.csv","a")
        f.puts "Family,#{node1},#{node2}"
        @rfam_counter.each do |key,values|
          f.print "#{key},#{values.shift},#{values.shift}"
          rfam = ToolDB::Rfam.find_by_acc("#{key}")
          f.print ",#{rfam.name},#{rfam.family},#{rfam.targets}" unless rfam.nil?
          f.print "\n"
        end  
        f.close
        
        # Difference in host gene content between nodes #
        #################################################
        counter = 0
        @go_terms = Hash.new
        
        nodes_groups.each do |node_groups|
          
          node_groups.each do |group|
            
            group = ToolDB::Grouping.find(group)
            group.fetch_unique_goterms.each do |go|
              if @go_terms.has_key?(go)
                values = @go_terms.fetch(go)
                values[counter] = values[counter] += 1
              else
                values = [ 0 , 0 ]
                values[counter] = values[counter] += 1
              end
              @go_terms[go] = values
            end  
            
          end
          counter += 1
        end
        
        f = File.new("#{@base_name}_comparison_go_#{node1}_#{node2}.csv", "a")
        f.puts "Comparison of Nodes #{node1} and #{node2}"
        f.puts
        f.puts "GO term,#{node1},#{node2},name,namespace"
        @go_terms.each do |term,values|
          go = ToolDB::GoTerm.find_by_term(term)
          f.print "#{term},#{values.join(",")},"
          unless go.nil?
            f.print "#{go.name.gsub(/\,/, ' ')},#{go.namespace}"
          end
          f.print "\n"
        end
        
      end
      
      def go_hierarchy
        
        bp = ToolDB::GoTerm.find_by_name('biological_process')
        bp_hierarchy = bp.hierarchy
        return bp_hierarchy
        
      end
      
      # = DESCRIPTION
      # Checks whether a type of intronic npcRNA is
      # correlated to a host gene type
      def regulation_hypothesis(node)
        
        node = ToolDB::TreeNode.find_by_nr(node)
        raise "Node not present in dataset" if node.nil?
        list_of_terms = [ "regulation" ]
        group_counter = 0
        group_counter_intronic = 0
        translation_counter = 0
        hostgene_counter = 0
        
        groups = ToolDB::Grouping.find_all_by_biotype(self.type)
        
        hostgene_names = Array.new
        @matched_terms = Hash.new
        
        groups.each do |group|
          group_counter += 1
          if group.present_at?(node) and group.is_intronic?
            group_counter_intronic += 1
            hostgenes = group.fetch_hostgenes
            match_query = false
            hostgenes.each do |hostgene|
              hostgene.fetch_goterms.each do |go_term|
                unless match_query == true
                  go = ToolDB::GoTerm.find_by_term("#{go_term.go_term}")
                  unless go.nil?
                    match_query = true if go.check_terms(list_of_terms,@matched_terms) == true
                  end
                end
              end
            end
            translation_counter += 1 if match_query == true
          end
          
        end
        
        puts "#{@type} Groups in total #{group_counter} - of which #{group_counter_intronic} are intronic at node #{node.nr}"
        puts " of these #{group_counter_intronic} intronic groups, #{translation_counter} match the queries #{list_of_terms.join(",")}"
        
      end   
      
      # = DESCRIPTION
      # Analyzes the mode of evolution for a group
      # of ncRNAs across internal nodes
      def duplication
        
        nodes = ToolDB::TreeNode.find(:all, :order => "nr")
        organisms = ToolDB::Organism.find(:all, :order => "name")
        
        f = File.new("#{@base_name}_#{type}_duplications_summary.csv", "a")
        f.puts "Organism,trans duplication,cis duplication,intergenic duplications"
        f.close
        
        organisms.each do |organism|
          
          intronic_group_hash = Hash.new
          intergenic_group_hash = Hash.new
          
          hostgenes = organism.ensembl_hostgenes
          
          hostgenes.each do |hostgene|
            
            local_groups = Hash.new
            
            ncrnas = hostgene.fetch_ncrnas(self.type)
            
            ncrnas.each do |ncrna|
              
              if ncrna.has_group?
                
                if local_groups.has_key?(ncrna.external_id)
                  local_groups[ncrna.external_id] = local_groups.fetch(ncrna.external_id) + 1
                else
                  local_groups[ncrna.external_id] = 1
                end
              
              end
              
            end
            
            local_groups.each do |local_group,count|
              
              if intronic_group_hash.has_key?(local_group)
                intronic_group_hash[local_group] = intronic_group_hash.fetch(local_group).push(count)
              else
                intronic_group_hash[local_group] = [ count ]
              end  
              
            end
          
          end  
          
          
          organism.fetch_ncrnas(self.type).each do |ncrna|
            
            if ncrna.has_group? == true and ncrna.has_hostgene? == false
              
              if intergenic_group_hash.has_key?(ncrna.external_id)
                intergenic_group_hash[ncrna.external_id] = intergenic_group_hash.fetch(ncrna.external_id) + 1
              else
                intergenic_group_hash[ncrna.external_id] = 1
              end
              
            end
            
          end
          
          cis_duplications = 0
          trans_duplications = 0
          intergenic_duplications = 0
          
          f = File.new("#{@base_name}_#{type}_duplications_#{organism.name}.csv", "a")
          
          intronic_group_hash.each do |group,counts|
            
            if counts.nitems > 1
              trans_duplications = trans_duplications + (counts.nitems - 1)  # if a group of ncRNA is present at multiple locations, assume trans duplication and subract one (the original copy)
            end
            counts.each do |count|                                          # if a group of intronic ncRNAs has more than one copy within a gene, add them to the cis count, minus one (the original copy)  
              cis_duplications = cis_duplications + (count-1)
            end
            # for intergenic groups that also have intronic members, add the total count to the number of intergenically duplicated ncRNAs
            intergenic_duplications = intergenic_duplications + intergenic_group_hash.fetch(group) if intergenic_group_hash.has_key?(group)
            
            f.print "#{group},#{counts.join(";")}"
            f.print ", #{intergenic_group_hash.fetch(group)}" if intergenic_group_hash.has_key?(group)
            f.print "\n"
            intergenic_group_hash.delete(group)
          end
          
          # For the remaining groups of ncRNAs (without intronic members), add the count -1 to the number of intergenically duplicated ncRNAs
          intergenic_group_hash.each do |group,count|
            f.puts "#{group},no intronic members,#{count}"
            intergenic_duplications += (count-1)
          end
          
          f.close
          f = File.new("#{@base_name}_#{type}_duplications_summary.csv", "a")
          f.puts "#{organism.name},#{trans_duplications},#{cis_duplications},#{intergenic_duplications}"
          f.close
          
        end
       
      end  
      
      # = Cross-references miRNAs and hostgenes
      # with internal nodes
      def mirbase_targets_and_nodes
        
        nodes = Toolbox::ToolDB::TreeNode.find(:all, :order => "nr")
        mirnas = []
        Toolbox::ToolDB::EnsemblNcrna.find_all_by_biotype_and_organism_id("miRNA",1).each do |mirna|
          mirnas.push(mirna) if mirna.has_group?
        end
        
        nodes.each do |node|
          
          #f = File.new("#{base_name}_#{type}_mirbase_targets_#{node.nr}.csv", "a")
          
          mirnas.each do |mirna|
            next unless mirna.present_at?(node)
            
            puts "#{mirna.stable_id},"
            mirna.xref_mirbase_ensembls.each do |xref|
              puts "\t#{xref.mirbase_product.name},#{xref.score},"
              targets = xref.mirbase_product.mirbase_targets
              puts "\t\t#{targets.nitems} targets"
              Ensembl::Core::DBConnection.connect('homo_sapiens')
              go_list = {}
              targets.each do |target|
                transcript = Ensembl::Core::Transcript.find_by_stable_id(target.target)
                next if transcript.nil?
                transcript.gene.go_terms.each do |go|
                  if go_list.has_key?("#{go}")
                    go_list["#{go}"] = go_list.fetch("#{go}") + 1
                  else
                    go_list["#{go}"] = 1
                  end
                end
              end
              go_list.each do |term,count|
                puts "\t\t\t#{term}:#{count}"
              end
                  
            end
            
          end
          
        end
        
        
      end
      
      def hostgene_content
        
        ToolDB::DBConnection.connect
        self.type = "snoRNA"
        nodes = ToolDB::TreeNode.find(:all, :order => "nr")

        puts ",#{nodes.collect{|n| n.nr}.join(",")}"
        counter_total = { "1" => [] , "2" => [] , "3" => [] , "4" => [] , "5" => [] , "6" => [] , "7" => [] , "8" => [] }

        nodes.each do |node|

          node_counter = { "1" => 0 , "2" => 0 , "3" => 0 , "4" => 0 , "5" => 0 , "6" => 0 , "7" => 0 , "8" => 0 }

          node.fetch_hostgenes(1,"#{type}").each do |gene|
            count = 0
            gene.fetch_ncrnas(type).each do |ncrna|
              count += 1
            end
            node_counter["#{count}"] = node_counter.fetch("#{count}") + 1
          end

          node_counter.sort.each do |key,value|
            counter_total["#{key}"] = counter_total.fetch("#{key}").push(value)
          end

        end

        counter_total.sort.each do |key,values|
          puts "#{key} ncRNA(s),#{values.join(",")}"
        end
      end
      
      # = DESCRIPTION
      # Lists npcRNAs sorted by external ID, with node status
      # Requires an organism id as qualifier
      def list_by_rfam(organism_id)
        
        rfams = ToolDB::Rfam.find(:all, :order => "acc")
        nodes = ToolDB::TreeNode.find(:all, :order => "nr")

        puts ",#{nodes.collect{|n| n.nr}.join(",")}"

        rfams.each do |rfam|

          snornas = ToolDB::EnsemblNcrna.find_all_by_organism_id_and_external_id(organism_id,rfam.acc).select{|s| s.has_group? }

          next if snornas.empty?

          puts "#{rfam.acc}"

          snornas.each do |snorna|

            result = ""

            nodes.each do |node|
              snorna.present_at?(node) ? state = "X" : state = "0"
              result += ",#{state}"
            end  

            next unless result.include?("X")

            print ",#{snorna.stable_id},#{result}"

            if snorna.is_intronic?
              gene = snorna.fetch_hostgene
              print ",intronic,#{gene.stable_id},#{gene.description}"
            else
              print ",intergenic"
            end

            print "\n"

          end

        end
      
      end
      
      # = DESCRIPTION
      # Returns a list of Ensembl protein accession numbers
      # for host genes present at a given node
      def go_stat(organism_id,default_node,sub_node = nil,seqtype="gene",ensembl_version=54)
        
        def_node = ToolDB::TreeNode.find_by_nr(default_node)
        sub_node = ToolDB::TreeNode.find_by_nr(sub_node) unless sub_node.nil? or sub_node == default_node
        
        genes = ToolDB::EnsemblHostgene.find_all_by_organism_id(1)
        genes = genes.select{|g| g.present_at?(self.type,def_node)}
        genes = genes.select{|g| g.present_at?(self.type,sub_node) == false } unless sub_node.nil? or sub_node == default_node

        genes.each do |gene|
          
          if seqtype == "gene"
            
            #puts gene.stable_id
            puts gene.fetch_ensembl_gene.uniprot_id
            
          else
            
            e = gene.fetch_ensembl_gene

            answer = []
          
            e.transcripts.each do |t|
              next if t.translation.nil?
              answer.push(t.translation.stable_id)
            end
            answer.compact!
            next if answer.nil?
          
            puts "#{answer.shift}"
          end
          
        end
        
      end
      
    end
    
    class Process
      
      attr_reader :input
      
      def initialize(input)
        @input = Hash.new
      end
      
    end
       
  end
  
end