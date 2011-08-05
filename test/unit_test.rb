#!/usr/bin/ruby

require 'test/unit'
require '/home/marc/scripts/tools/toolbox.rb'
include Toolbox

HomologDB::DBConnection.connect(56)

class HomologDBTest < Test::Unit::TestCase
  
  def test_for_supergroup    
    organisms = HomologDB::Organism.find(:all)
    organisms.each do |organism|
      assert organism.supergroup_id != nil , "#{organism.name} has no supergroup_id"
    end
    puts "> All organisms have a supergroup id"
  end
  
  def test_infernal_mapping    
    member = HomologDB::Member.find(2784127)    
    member.infernal_hits.each do |hit|
      assert_equal("#{member.nucleotide_to_fasta.naseq[hit.member_start-1..hit.member_stop-1]}","#{hit.target_seq.downcase.gsub(/u/, 't').gsub(/-/, '')}")
    end
    puts "> Infernal hits map correctly to the nucleotide sequence"
  end
  
  def test_member_has_peptide_member
    members = HomologDB::Member.find(:all).select{|m| m.peptide_members.empty?}
    assert_true members.empty? , "We have #{members.nitems} members without peptide members :("
  end
  
  def test_homology_member_has_peptide_member
    hmembers = HomologyDB::HomologyMember.find_all_by_peptide_member_id(nil)
    assert_true hmembers.empty? , "We have #{hmembers.nitems} homology members without peptide members :("
  end
  
end
  
