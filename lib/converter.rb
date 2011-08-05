
module Toolbox
  
  module Converter

    class Trees

      attr_reader :tree, :code, :taxa, :plain_tree, :title
      attr_writer :code

      def initialize(tree_file,title="none")
        @plain_tree = IO.readlines(tree_file).collect{ |t| t.strip }.join # a simple string
        @tree = Bio::Tree.new(Bio::Newick.new(@plain_tree).tree)          # a Bio::Tree object
        @code = Hash.new # encoding, for phylip format
        @taxa = Array.new
        @title = title
        @tree.each_node do |node|
          unless node.inspect.match(/\d/)
            @taxa.push(node.to_s.gsub(/\s/, '_'))
          end
        end

      end

      def to_nexus

        answer = Array.new

        tree_string = self.plain_tree

        answer.push("BEGIN TREES;")
        answer.push("\tTITLE #{@title}_tree;")
        answer.push("\tLINK Taxta = #{@title};")
        answer.push("\tTRANSLATE")

        @taxa.each do |taxon|
          if @taxa.index(taxon) == @taxa.nitems-1
            answer.push("\t\t#{@taxa.index(taxon)} #{taxon};")
          else
            answer.push("\t\t#{@taxa.index(taxon)} #{taxon},")
          end
          tree_string.gsub!(/#{taxon}/, "#{@taxa.index(taxon)}")
        end

        answer.push("\tTREE No_001 = #{tree_string}")
        answer.push("END;")

        return answer.join("\n")

      end

      def phylip

        answer = Array.new
        tree = self.tree
        t = self.plain_tree
        taxa = Hash.new
        counter = 1
        
        tree.each_node do |node|

          unless node.inspect.match(/[0-9]/)
            taxa[node.to_s] = "PF000#{counter}"
            t.gsub!("#{node.to_s.gsub(/\s/, '_')}", "#{taxa.fetch(node.to_s)}")
            counter += 1
          end

        end
        
        @code = taxa
        return t

      end
      
      def from_phylip(key)
        
        @keys = Hash.new
        
       IO.foreach(key) do |line|
          elements = line.strip.split(",")
          @keys[elements[1]] = elements[0].gsub(/[()]/, '')
        end
        
        answer = self.plain_tree

        @keys.each do |code,ordinary|
          answer.gsub!(code, ordinary)
        end
          
        return answer
        
      end
      
      
    end

    class Formats
      
      attr_reader :infile, :code
      attr_writer :code
      
      def initialize(infile)
        
        @infile = infile
        @code = Hash.new
        
      end
      
      def fasta_to_phylip
        
        ff = Bio::FastaFormat.open(infile)
        counter = 0
        
        ff.each do |entry|
          
          counter += 1
          if counter < 10
            zeros = "0000000"
          elsif counter >= 10
            zeros = "000000"
          elsif counter >= 100
            zeros = "00000"
          elsif counter >= 1000
            zeros = "0000"
          end
          name = "PF#{zeros}#{counter}"
          
          @code[entry.definition] = name
          
          puts ">#{name}"
          puts "#{entry.aaseq}"

        end
        
        k = File.new("#{@infile}_key.txt", "a")
        @code.each do |key,taxon|
          k.puts "#{key},#{taxon}"
        end
        k.close
        
      end
      
      def create_phylip(aHash)        
        
      end
              
    end
    
    class Phylip
    
      attr_accessor :infile, :keys
      
      def intialize(infile)
        @infile = infile
        @keys = nil
      end
      
      def convert_names        
        
      end
      
    end
    
    class Tables
      
      def initialize(file)
        f = IO.readlines(file)
        @header = f.shift.split("\t")
        @header.shift
        @rows = Hash.new
        f.each do |line|
          elements = line.strip.split("\t")
          @rows["#{elements.shift}"] = elements
        end
      end
      
      def convert
        
        new_header = []
        @rows.each_key { |k| new_header.push(k.strip) }
        new_rows = @header
        puts ",#{new_header.join(",")}"
        new_rows.each do |row|
          print "#{row.strip} - "
          @rows.each do |k,values|
            print "#{values.shift},"
            @rows[k] = values
          end
          print "\n"  
        end
      end
      
    end
    
  end

end