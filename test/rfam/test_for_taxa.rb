require 'test/unit'
require '/home/marc/scripts/tools/toolbox.rb'
include Toolbox


class RfamDBTest < Test::Unit::TestCase

  def test_for_taxa
  	RfamDB::DBConnection.connect(10)
  	taxa = RfamDB::Taxonomy.find(:all, :conditions => ["tax_string LIKE ?", "%Eukaryota%"])
  	puts "Found #{taxa.nitems} taxa"
  	taxa.each do |taxon|
  		next unless taxon.tax_string.include?(";")
  		keyword = taxon.tax_string.split(";")[1].strip.gsub(/\./, '')
  		assert RfamDB::TaxLookup.find_by_name(keyword).nil? == false , "Unknown taxon #{taxon.tax_string}"
  	end
  	
  end
  
end
