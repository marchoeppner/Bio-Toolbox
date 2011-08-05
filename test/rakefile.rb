#!/usr/bin/ruby
# HomologDB Unit tests

require '/home/marc/scripts/tools/toolbox.rb'
require "rake/testtask"
include Toolbox

HomologDB::DBConnection.connect(56)




