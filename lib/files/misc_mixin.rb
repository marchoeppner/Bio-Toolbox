class String

  def phylip_name
    string = self
    answer = ""
    if string.include?("_")
      return string.capitalize[0..9]
    else
      string.split(" ").each do |s|
        answer += "_#{s[0..4]}"
      end
      return answer[0..9].gsub(/^_/, '')
    end
    
  end	

	def strip_stockholm
		return self.gsub(/^\#=G[A-Z]\s[A-Z]*(\s+)*/, '')
	end
	
  def mask(start,stop,char=".")
    return char * (stop-start+1) + self[stop+1..-1] if start == 0 # if the mask starts at the beginning, do not slice off the first nt
    return self[0..start-1] + char * (stop-start+1) + self[stop+1..-1]
  end
  
  
  def transform_date
  	
  	elements = self.split("/")
  	
  	return "#{elements[2]}-#{elements[0]}-#{elements[1]}"
  	
  end
  
  def sql_compliant
    return self.gsub(/\'/, '').gsub(/\"/, '').gsub(/\/*/, '').gsub(/\´/, '').gsub(/;/, ",").gsub(/\n/, '')
  end
  
	def to_snake_case
		return self.gsub(/\s+/, '_').downcase
	end
		
  def binary_to_ascii
    
    return self.unpack('d')
  end
  
  def to_width(width)
    answer = []
    this_seq = self
    until this_seq.length <= width
      answer << this_seq.slice!(0..width-1)
    end
    answer << this_seq
    return answer.join("\n")
  end
  
  def softmasked_seq_to_html(window_size)
    seq = self.to_s
    exon = "red"
    intron = "grey"

    answer = []

    answer << "<table width='#{window_size}px' cellspacing='1px' height='10px'>\n<tr>\n"
    answer << "<td width='#{window_size/2-1}px'></td><td width='2px' bgcolor='blue'></td><td width='#{window_size/2-1}px'></td>\n</tr>\n<tr>\n</table>"
    answer << "<table width='#{window_size}px' cellspacing='0px' height='10px'>\n<tr>\n"

    until seq.length == 0
      if seq[0..0].match(/[A-Z]/)
        this_seq = seq.slice!(/^[A-Z]*/)
        answer << "<td width='#{this_seq.length}px' bgcolor='#{exon}'></td>\n"
      else
        this_seq = seq.slice!(/^[a-z]*/)
        answer << "<td width='#{this_seq.length}px' bgcolor='#{intron}'></td>\n"
      end
    end
    answer << "\n</tr>\n</table>"

    return answer.join

  end
  
  
end
