module Toolbox
  
  include Toolbox::ToolDB
  
  class SQLinject
    
    require 'rexml/document'
    include REXML
    
    attr_accessor :infile, :type, :version
    
    def initialize(infile,type,version="_54")
      
      @infile = infile
      @type = type
      @version = version
      
    end
    
    # = DESCRIPTION
    # Parses a project-compatible XML and stores the data in a
    # relational database
    def parse_xml
      
      Toolbox::ToolDB::DBConnection.connect
      
      infile = File.new(@infile, "r")
      doc = Document.new infile
      
      XPath.each( doc, "//Group") do |group|
        
        #IDENTIFY THE GROUP OF ORTHOLOGUES #
        g = Toolbox::ToolDB::Grouping.new(:nr => "#{group.attributes["nr"]}", :biotype => "#{self.type}", :genomic_align_block_id => "#{group.attributes["genomic_align_block_id"]}")
        g.save unless Toolbox::ToolDB::Grouping.exists?(:nr => "#{group.attributes["nr"]}", :biotype => "#{self.type}", :genomic_align_block_id => "#{group.attributes["genomic_align_block_id"]}")
        grouping_id = Toolbox::ToolDB::Grouping.find(:first, :conditions => ["nr = ? and biotype = ?", "#{group.attributes["nr"]}","#{self.type}"]).id
        
        group.elements.each("NcRNA") do |ncrna|
          
          unless ncrna.attributes["stable_id"] == "nil"
            organism_id = Toolbox::ToolDB::Organism.find_by_name("#{ncrna.attributes["organism"]}").id
            external_id = "#{ncrna.attributes["ncrna_id"]}"
            
            # STORE THE NCRNA #
            n = Toolbox::ToolDB::EnsemblNcrna.new(:stable_id => "#{ncrna.attributes["stable_id"]}", 
                            :external_id => "#{external_id}", 
                            :organism_id => "#{organism_id}", 
                            :biotype => "#{self.type}")
            n.save unless Toolbox::ToolDB::EnsemblNcrna.exists?(:stable_id => "#{ncrna.attributes["stable_id"]}")
            ncrna_id = Toolbox::ToolDB::EnsemblNcrna.find_by_stable_id("#{ncrna.attributes["stable_id"]}").id
            
            # STORE Xref NCRNAxGroup  = To which group does a NCRNA belong #
            ng = Toolbox::ToolDB::NcrnaGroup.new(:grouping_id => "#{grouping_id}", 
                            :ensembl_ncrna_id => "#{ncrna_id}")
            ng.save unless Toolbox::ToolDB::NcrnaGroup.exists?(:grouping_id => "#{grouping_id}", :ensembl_ncrna_id => "#{ncrna_id}")
            
            ncrna.elements.each("HostGene") do |hostgene|
              
              # STORE HOST GENE, if any #
              h = Toolbox::ToolDB::EnsemblHostgene.new(:stable_id => "#{hostgene.attributes["name"]}", 
                            :organism_id => "#{organism_id}", 
                            :description => "#{hostgene.attributes["description"]}")
              h.save unless Toolbox::ToolDB::EnsemblHostgene.exists?(:stable_id => "#{hostgene.attributes["name"]}")
              hostgene_id = Toolbox::ToolDB::EnsemblHostgene.find_by_stable_id("#{hostgene.attributes["name"]}").id
              
                # STORE HOST GENE AND NCRNA Xref = is the ncRNA intronic and what's the host gene #
              x = Toolbox::ToolDB::XrefNcrnaHostgene.new(:ensembl_hostgene_id => "#{hostgene_id}", :ensembl_ncrna_id => "#{ncrna_id}")
              x.save unless Toolbox::ToolDB::XrefNcrnaHostgene.exists?(:ensembl_hostgene_id => "#{hostgene_id}", 
                            :ensembl_ncrna_id => "#{ncrna_id}")
                            
                # Store the GO terms of the host gene #
              hostgene.elements.each("GoId") do |go_id|
                go = Toolbox::ToolDB::HostgeneGoterm.new(:ensembl_hostgene_id => "#{hostgene_id}",
                            :go_term => "#{go_id.attributes["id"]}")
                go.save unless Toolbox::ToolDB::HostgeneGoterm.exists?(:ensembl_hostgene_id => "#{hostgene_id}", :go_term => "#{go_id.attributes["id"]}")
              end
              
            end
              
          end
        end
      end
    
    end  
    
    def parse_mesquite
      
      ToolDB::DBConnection.connect(self.version)
            
      IO.foreach(@infile) do |line|
        
        if line.match(/^Char\./)
          
          elements = line.strip.split("\t")
          elements.shift
          @nodes = elements.collect{|n| n.strip}
         
        elsif line.match(/^character/)
          
          elements = line.split(/\t/)
          group = elements.shift.gsub(/^[a-z]*\s/, '')
          group_id = Toolbox::ToolDB::Grouping.find_by_nr_and_biotype("#{group}","#{self.type}").id
          
          @nodes.each do |node|
            node_id = ToolDB::TreeNode.find_by_nr("#{node}").id
            n = Toolbox::ToolDB::GroupNode.new(:tree_node_id => "#{node_id}", :grouping_id => "#{group_id}", :val => "#{elements.shift.strip}")
            n.save unless Toolbox::ToolDB::GroupNode.exists?(:tree_node_id => "#{node_id}", :grouping_id => "#{group_id}")
            
          end
          
        end
      
      end
      
    end
    
    # = DESCRIPTION
    # Parss the file of all ncRNAs (including non-syntenic)
    def parse_full
      
      Toolbox::ToolDB::DBConnection.connect
      
      infile = File.new(@infile, "r")
      doc = Document.new infile
      
      XPath.each( doc, "//Organism") do |organism|
        
        organism_id = ToolDB::Organism.find_by_name("#{organism.attributes["name"]}").id
        
        organism.elements.each("NcRNA") do |ncrna|
          
          external_id = "#{ncrna.attributes["external_id"]}"
          
          n = Toolbox::ToolDB::EnsemblNcrna.new(:stable_id => "#{ncrna.attributes["stable_id"]}", 
                          :external_id => "#{external_id}", 
                          :organism_id => "#{organism_id}", 
                          :biotype => "#{self.type}")
          n.save unless Toolbox::ToolDB::EnsemblNcrna.exists?(:stable_id => "#{ncrna.attributes["stable_id"]}")
          ncrna_id = Toolbox::ToolDB::EnsemblNcrna.find_by_stable_id("#{ncrna.attributes["stable_id"]}").id
          
          ncrna.elements.each("HostGene") do |hostgene|
            h = Toolbox::ToolDB::EnsemblHostgene.new(:stable_id => "#{hostgene.attributes["stable_id"]}", 
                          :organism_id => "#{organism_id}", 
                          :description => "#{hostgene.attributes["description"]}")
            h.save unless Toolbox::ToolDB::EnsemblHostgene.exists?(:stable_id => "#{hostgene.attributes["stable_id"]}")
            hostgene_id = Toolbox::ToolDB::EnsemblHostgene.find_by_stable_id("#{hostgene.attributes["stable_id"]}").id
          
              # STORE HOST GENE AND NCRNA Xref = is the ncRNA intronic and what's the host gene #
            x = Toolbox::ToolDB::XrefNcrnaHostgene.new(:ensembl_hostgene_id => "#{hostgene_id}", :ensembl_ncrna_id => "#{ncrna_id}")
            x.save unless Toolbox::ToolDB::XrefNcrnaHostgene.exists?(:ensembl_hostgene_id => "#{hostgene_id}", 
                          :ensembl_ncrna_id => "#{ncrna_id}")
                        
              # Store the GO terms of the host gene #
            hostgene.elements.each("GoId") do |go_id|
              go = Toolbox::ToolDB::HostgeneGoterm.new(:ensembl_hostgene_id => "#{hostgene_id}",
                          :go_term => "#{go_id.attributes["id"]}")
              go.save unless Toolbox::ToolDB::HostgeneGoterm.exists?(:ensembl_hostgene_id => "#{hostgene_id}", :go_term => "#{go_id.attributes["id"]}")
            end
          end
          
        end
          
      end
      
    end
    
    def run_gerp_scores(window_size=1,upstream=0,downstream=0,version="_54")
      
      ToolDB::DBConnection.connect(version)
      
      groups = ToolDB::Grouping.find_all_by_biotype(self.type)
      
      groups.each do |group|
        next if group.id < 2175
        scores = group.fetch_diff_scores(window_size,upstream,downstream)
        
        puts "INSERT INTO gerp (grouping_id,scores,window_size) VALUES(#{group.id},'#{scores.join(" ")}',#{window_size});"
        
      end
      
    end
      
  end

end