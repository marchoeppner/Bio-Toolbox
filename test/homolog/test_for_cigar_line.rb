require 'test/unit'
require '/home/marc/scripts/tools/toolbox.rb'
include Toolbox


class HomologDBTest < Test::Unit::TestCase

  def test_for_cigar_line
		
		puts "Testing for the presence of cigar lines"
		HomologDB::DBConnection.connect(58)
		groups = HomologDB::HomologyGroup.find_all_by_snorna(true)
		groups.each do |group|
			
			hmembers = group.homology_members.select{|h| h.member.organism.active }
			assert hmembers.select{|h| h.peptide_cigar_line == ""}.empty? , "Missing cigar lines in group #{group.id}"
		
		end
  	
  	
  end
  
end
