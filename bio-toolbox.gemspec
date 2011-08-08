Gem::Specification.new do |s|
  s.name              = "bio-toolbox"
  s.version           = "0.0.1"
  s.platform          = Gem::Platform::RUBY
  s.authors           = ["Marc P. Hoeppner"]
  s.email             = ["mphoeppner@gmail.com"]
  s.homepage          = "http://github.com/marchoeppner/Bio-Toolbox"
  s.summary           = "A toolbox for bioinformatics"
  s.description       = "See above..."

  s.required_rubygems_version = ">= 1.3.5"

  # If you have runtime dependencies, add them here
  # s.add_runtime_dependency "other", "~> 1.2"

  # If you have development dependencies, add them here
  # s.add_development_dependency "another", "= 0.9"
  s.add_development_dependency "ruby-ensembl-api", "~>1.0"
  s.add_development_dependency "composite_primary_keys", "~>3.0"
  s.add:development_dependency "mysql2", "~>0.1"

  # The list of files to be contained in the gem
  s.files         = `git ls-files`.split("\n")
  # s.executables   = `git ls-files`.split("\n").map{|f| f =~ /^bin\/(.*)/ ? $1 : nil}.compact
  # s.extensions    = `git ls-files ext/extconf.rb`.split("\n")

  s.require_path = 'lib'

  # For C extensions
  # s.extensions = "ext/extconf.rb"
end
