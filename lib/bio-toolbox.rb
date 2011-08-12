# Marc's Bioinformatics Toolbox...

require 'rubygems'
require 'composite_primary_keys'

gem 'ruby-ensembl-api', '>=0.9.6'
require 'ensembl'


# Include the source files
require File.dirname(__FILE__) + '/files/ensembl_compara_connection.rb'
require File.dirname(__FILE__) + '/files/ensembl_compara.rb'

require File.dirname(__FILE__) + '/files/bioruby_extensions.rb'
require File.dirname(__FILE__) + '/files/ensembl_mixin.rb'
