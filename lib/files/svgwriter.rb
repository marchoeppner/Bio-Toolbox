
module Toolbox

	module SVGWriter
	
		class Graph
			
			attr_accessor :width, :height, :elements, :output
			
			def initialize(width=744,height=1052,elements=[])
				@width = width
				@height = height
				@elements = elements
				@output = []
			end
			
			def add_rectangle(paramaters,width,height=20)
				elements << SVGWriter::Element::Rectangle.new(paramaters,width,height)	
			end
			
			def add_path(position,width)
				elements << SVGWriter::Element::Path.new(position,width)
			end
			
			def add_label(text,position,fontsize)
			  elements << SVGWriter::Element::Text.new(text,position,fontsize)
			end  
			
			def draw
			
				output << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>"
				output << "<svg"
   			output << "\sxmlns=\"http://www.w3.org/2000/svg\""
   			output << "\swidth=\"#{self.width}\""
   			output << "\sheight=\"#{self.height}\""
   			output << "\sid=\"svg#{rand(10)}\">"

  			output << "\s<g"
     		output << "\s\sid=\"layer1\">"
			
				elements.each do |element|
				
					output << element.draw
			
				end
				
				output << "\s</g>"
				output << "</svg>"
			
				return output.join("\n")
				
			end
			
		end
		
		module Element
		
		  class Text
		  
		    attr_accessor :text, :position, :fontsize
		    
		    def initialize(text,position,fontsize)
		      @text = text
		      @position = position
		      @fontsize = fontsize
		    end
		    
		    def generate_id
		      return "text#{rand(1000)}"
		    end
		    
		    def draw
		      output = []
		      output << "<text"
          output << "font-family=\"Verdana\" font-size=\"#{self.fontsize}\""
          output << "x=\"#{self.position[:x]}\""
          output << "y=\"#{self.position[:y]}\""
          output << "id=\"#{self.generate_id}\">"
          output << "#{self.text}"
          output << "</text>"
          return output.join("\n")
		    end
		    
		  end
		  
			class Rectangle
		
				attr_accessor :parameters, :height, :width
			
				def initialize(parameters,width,height=20)
					@parameters = parameters
					@height = height
					@width = width
				end
			
				def generate_id
					return "rect#{rand(1000)}"
				end
		
				def draw
					output = []
					parameters[:stroke] ? stroke = parameters[:stroke] : stroke = 1
					output << "\s\s<rect"
       		output << "\s\s\sstyle=\"fill:#{parameters[:fill]};fill-rule:evenodd;stroke:#000000;stroke-width:#{stroke}px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\""
       		output << "\s\s\sid=\"#{self.generate_id}\""
       		output << "\s\s\swidth=\"#{self.width}\""
       		output << "\s\s\sheight=\"#{self.height}\""
       		output << "\s\s\sx=\"#{self.parameters[:x]}\""
       		output << "\s\s\sy=\"#{self.parameters[:y]}\" />"
					
					return output.join("\n")
				end
	
			end
			
			class Path
				
				attr_accessor :position, :width, :angle
				def initialize(position,width,angle=0)
					
					@position = position
					@width = width
					@angle = angle
				end
			
				def draw
					output = []
					output << "\s\s<path"
					output << "\s\s\sstyle=\"fill:none;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\""
					output << "\s\s\sd=\"m #{position[:x]},#{position[:y]} #{self.width},#{self.angle}\""
					output << "\s\s\sid=\"path#{rand(100)}\" />"
					return output.join("\n")
				end
			end
		
		end

	end
	
end
