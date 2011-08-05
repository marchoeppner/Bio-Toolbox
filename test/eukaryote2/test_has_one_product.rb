require 'test/unit'
require '/home/marc/scripts/tools/toolbox.rb'
include Toolbox

EukaryoteDB::DBConnection.connect(3)

class EukaryoteDBTest < Test::Unit::TestCase
	
	def test_has_one_product
		puts "> Checking for multiple products."
		translations = EukaryoteDB::Translation.find_by_sql("SELECT transcript_id FROM translation GROUP BY transcript_id,stable_id").collect{|t| t.transcript_id}
		puts translations[0..10].inspect
		assert translations.nitems == translations.uniq.nitems, "Difference between transcript counts encountered (#{translations.nitems}/#{translations.uniq.nitems})"

	end  
end
