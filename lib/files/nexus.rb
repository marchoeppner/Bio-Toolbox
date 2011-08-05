
module Toolbox
  
  module Nexus
    

    class Writer
      
      attr_reader :title, :taxa, :characters, :tree
      
      def initialize(title,taxa,tree_file) # taxa must be a hash!
        @title = title
        @taxa = taxa
        @tree_file = tree_file
      end
      
      def to_s
        
        return "#{header}\n#{taxon_list}\n#{characters}\n#{self.tree}\n" 
      
      end
        
      def header
        return "#NEXUS"
      end
      
      def taxon_list
        
        answer = Array.new
        answer.push("\n")
        answer.push("BEGIN TAXA;")
        answer.push("\tTITLE #{@title};")
        answer.push("\tDIMENSIONS NTAX=#{@taxa.keys.nitems};")
        answer.push("\tTAXLABELS")
        answer.push("\t\t#{@taxa.keys.collect{|o| o.gsub(/\s/, '_')}.join(" ")}")
        answer.push("\t;")
        answer.push("END;")
        answer.push("\n")
        
        return answer.join("\n")
        
      end
      
      def characters
        
        answer = Array.new
        answer.push("BEGIN CHARACTERS;")
        answer.push("\tTITLE #{@title};")
        answer.push("\tDIMENSIONS NCHAR=#{@taxa.to_a[0][1].length};")
        
        symbols = ""
        @taxa.each_value { |v| symbols = symbols + "#{v}"}
        symbols = symbols.scan(/./).uniq.sort!
        
        answer.push("\tFORMAT DATATYPE = STANDARD GAP = ? MISSING = - SYMBOLS = \" #{symbols.join(" ")}\";")
        answer.push("\tCHARSTATELABELS")
        tmp = ""
        counter = 0
        @taxa.to_a[0][1].length.times do
          counter += 1
          tmp += " #{counter} group_#{counter},"
        end
        answer.push("\t\t#{tmp};")
        answer.push("\tMATRIX")
        @taxa.each { |k,v| answer.push("\t#{k.gsub(/\s/, '_')} #{v}")}
        answer.push(";")
        answer.push("\n")
        answer.push("END;")
        answer.push("\n")
        
        return answer.join("\n")
        
      end
      
      def tree
        tree = Toolbox::Converter::Trees.new(@tree_file,@title)
        return tree.to_nexus
      end
      
    end
    
  end
  
end