require 'test/unit'
require '/home/marc/scripts/tools/toolbox.rb'
include Toolbox

HomologDB::DBConnection.connect(58)

class HomologDBTest < Test::Unit::TestCase

  def test_eukaryote_members
  
  	puts "Checking that eukaryote members have proteins."
  	organisms = HomologDB::Organism.find(:all, :conditions => ["database LIKE ?",'%eukaryote%'])
    organisms.each do |organism|
    
    	organism.members[0..500].each do |member|
    		
    		assert member.peptide_members.empty? == false , "#{member.stable_id} has no peptide members."
    	
    	end
    
    end
  end
end
