require 'test/unit'
require '/home/marc/scripts/tools/toolbox.rb'
include Toolbox

EukaryoteDB::DBConnection.connect(3)

class EukaryoteDBTest < Test::Unit::TestCase

  def test_produces_protein
  	genes = EukaryoteDB::Gene.find(:all, :limit => 100 )
  	genes.each do |gene|
  		gene.transcripts.each do |transcript|
  			cdna = []
  			seq = transcript.protein_seq
  			assert seq.gsub(/\*$/, '').split.select{|c| c == "*"}.nitems < 2 , "Translation product of #{gene.stable_id}, #{transcript.stable_id} contains illegal characters - #{seq}"
  		
  		end
  	end
  end
end
