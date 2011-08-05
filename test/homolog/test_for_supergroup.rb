require 'test/unit'
require '/home/marc/scripts/tools/toolbox.rb'
include Toolbox

HomologDB::DBConnection.connect(56)

class HomologDBTest < Test::Unit::TestCase

  def test_must_have_supergroup_id    
    organisms = HomologDB::Organism.find_all_by_active(true)
    organisms.each do |organism|
      assert organism.supergroup_id != nil , "#{organism.name} has no supergroup_id"
    end
    puts "> All organisms have a supergroup id" 
  end
end
