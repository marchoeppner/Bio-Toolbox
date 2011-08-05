require 'test/unit'
require '/home/marc/scripts/tools/toolbox.rb'
include Toolbox

HomologDB::DBConnection.connect(58)

class HomologDBTest < Test::Unit::TestCase

  def test_homology_member_has_peptide_member
    hmembers = HomologDB::HomologyMember.find_all_by_peptide_member_id(nil).select{|hm| hm.member.organism.active}
    assert hmembers.empty? , "We have #{hmembers.nitems} homology members without peptide members :("
    puts "> All homology members have a peptide member."
  end
end
