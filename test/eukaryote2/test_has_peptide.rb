require 'test/unit'
require '/home/marc/scripts/tools/toolbox.rb'
include Toolbox

EukaryoteDB::DBConnection.connect(3)

class EukaryoteDBTest < Test::Unit::TestCase
	
	def test_has_peptide
		puts "> Checking Inparanoid Clusters for peptide members."
		inparanoid_clusters = EukaryoteDB::InparanoidCluster.find(:all)
		inparanoid_clusters.each do |cluster|
			
			cluster.inparanoid_members.select{|m| m.gene_stable_id.include?("ENSG") == false}.each do |m|
				assert m.protein , "#{m.gene_stable_id}/#{m.protein_stable_id} has no peptide."
			end
		
		end

	end  
end
