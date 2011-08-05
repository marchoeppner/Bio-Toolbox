require 'test/unit'
require '/home/marc/scripts/tools/toolbox.rb'
include Toolbox

EukaryoteDB::DBConnection.connect(1)

class EukaryoteDBTest < Test::Unit::TestCase

  def test_gene_has_transcripts   
    genes = EukaryoteDB::Gene.find(:all)
    genes.each do |gene|
      assert gene.transcripts.empty? == false , "Gene #{gene.stable_id} (#{gene.description}) has no transcript"
    end
    puts "> All genes have transcripts"
  end
end