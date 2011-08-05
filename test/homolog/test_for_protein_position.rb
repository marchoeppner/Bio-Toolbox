require 'test/unit'
require '/home/marc/scripts/tools/toolbox.rb'
include Toolbox

HomologDB::DBConnection.connect(56)

class HomologDBTest < Test::Unit::TestCase

  def test_protein_has_position   
  
  	xrefs = HomologDB::XrefPeptideNcrna.find(:all, :conditions => ["position IS NULL"])
    assert xrefs.empty? , "#{xrefs.nitems} snorna-to-protein mappings have no position."
    puts "> All ncRNA mappings have a valid position" 
  end
end
