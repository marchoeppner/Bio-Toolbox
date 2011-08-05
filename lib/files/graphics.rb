module Toolbox
  
  class Table
    
    attr_accessor :data, :color_table, :header
    
    def initialize(data,header)
      @data = data # a Hash
      @header = header
      @rows = @data.to_a.nitems
      @columns = @data.to_a.first[1].nitems
      @color_table = { "0" => "#ffffff" , "1" => "#eeeeee" , "2" => "#dddddd" , "3" => "#cccccc" , "3" => "#bbbbbb" , "4" => "#aaaaaa" , "5" => "#999999" , "6" => "#888888" ,  "7" => "#777777" , 
        "8" => "666666", "9" => "#555555", "10" => "#444444" , "11" => "#333333" , "12" => "#222222" }
    end
    
    def draw_table
      
      answer = []
      answer.push("<html>")
      answer.push("<head>Some table")
      answer.push("</head>")
      answer.push("<style type='text/css'>")
      answer.push("table.sample { border-width: 1px 1px 1px 1px;border-spacing: 2px;border-style: outset outset outset outset;border-color: gray gray gray gray;border-collapse: separate;background-color: white;font-size:8px;}")
      answer.push("</style>")
      answer.push("<body>")
      answer.push(self.create_table)
      answer.push("</body>")
      answer.push("</html>")
      
      return answer.join("\n")
      
    end
    
    def get_color(value)
      
      if @color_table.has_key?("#{value}")
        return @color_table.fetch("#{value}")
      elsif value.to_i > 12
        return "#111111"
      else 
        return "#ffffff"
      end
      
    end
    
    def create_table
      
      answer = Array.new
      
      answer.push("<table class='sample'>")
      answer.push("\t<tr>")
      header = ""
      @header.each do |h|
        header = header + "<td>#{h}</td>"
      end
      answer.push("\t\t<td></td>#{header}")
      answer.push("\t</tr>")
      @data.each do |row,values|
        answer.push("\t<tr>")
        answer.push("\t\t<td>#{row}</td>")
        values.each do |value|
          answer.push("\t\t<td bgcolor='#{get_color(value)}'>#{value}</td>")
        end
        answer.push("\t</tr>")
      end
      
      return answer.join("\n")
      
    end
    
    
  end
  
end