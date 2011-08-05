# Marc's Bioinformatics Toolbox...

require 'rubygems'
require 'composite_primary_keys'
require 'active_record'
require 'bio'
require 'date'
#require 'ensembl'
require 'gruff'

gem 'ruby-ensembl-api', '>=0.9.6'
require 'ensembl'


require File.dirname(__FILE__) + '/lib/rfam_db.rb'
require File.dirname(__FILE__) + '/lib/pfam_db.rb'
require File.dirname(__FILE__) + '/lib/go_db.rb'
require File.dirname(__FILE__) + '/lib/expression_db.rb'
require File.dirname(__FILE__) + '/lib/mirbase_db.rb'
require File.dirname(__FILE__) + '/lib/genome_compara_db.rb'

require File.dirname(__FILE__) + '/lib/compara_db.rb'
require File.dirname(__FILE__) + '/lib/ensembl_compara_connection.rb'
require File.dirname(__FILE__) + '/lib/genome_db.rb'
require File.dirname(__FILE__) + '/lib/ncrna_db.rb'
require File.dirname(__FILE__) + '/lib/ncrna.rb'
require File.dirname(__FILE__) + '/lib/ensembl_mixin.rb'
require File.dirname(__FILE__) + '/lib/parser.rb'
require File.dirname(__FILE__) + '/lib/ensembl_compara.rb'
require File.dirname(__FILE__) + '/lib/tool_db.rb'
require File.dirname(__FILE__) + '/lib/leca_db.rb'
#require File.dirname(__FILE__) + '/lib/database_pipeline.rb'
require File.dirname(__FILE__) + '/lib/extensions.rb'
#require File.dirname(__FILE__) + '/lib/neanderthal_db.rb'
require File.dirname(__FILE__) + '/lib/homolog_db.rb'
require File.dirname(__FILE__) + '/lib/svgwriter.rb'
require File.dirname(__FILE__) + '/lib/metagenome.rb'
require File.dirname(__FILE__) + '/lib/taxonomy.rb'
require File.dirname(__FILE__) + '/lib/eukaryote_db.rb'