module Ensembl
  
  # = DESCRIPTION
  # The Ensembl::Compara module covers the compara database from
  # ensembldb.ensembl.org and includes genome alignments and orthology information.
  # For a full description of the database (and therefore the classes that
  # are available), see http://www.ensembl.org/info/software/compara/schema/index.html
  # and http://www.ensembl.org/info/software/compara/schema/schema_description.html
  #
  # = USAGE
  # To connect to the compara database, simply add 'Ensembl::Compara::DBConnection.connect'
  # to the beginning of your script. You may also want to 'include Ensembl::Compara'.
  # Note that the connection has to be spelled out fully to avoid conflicts with other 
  # database connections (Core, Vara).  
  
  module Compara
		# = DESCRIPTION
		# The Ensembl::Compara::DBConnection#connect method makes the connection
    # to the Ensembl Compara or Ensemblgenomes pan homology database. 
		# Ensembl Compara covers the original Ensembl databases (vertebrates, model systems)
		# and includes genomic alignments in addition to homology maps between protein-coding
		# genes. Ensemblgenomes includes homology maps only, but covers all of Ensembl - some of the API methods may therefore
		# not work for Ensemblgenomes data!
    # = USAGE
		#  # Connect to release 56 of the Ensembl Compara database
		# Ensembl::Compara::DBConnection.connect(56) 
		#
		# # Connect to release 56 of the Ensembl Pan Compara database
		# Ensembl::Compara::DBConnection.connect(56, { :ensembl_genomes => true })
		class DBConnection < Ensembl::DBRegistry::Base
  		self.abstract_class = true
      self.pluralize_table_names = false
      
      def self.connect(release = Ensembl::ENSEMBL_RELEASE,args = {} )
      	Ensembl::SESSION.reset
				#Use a Non-Core Compara database			
				if args[:ensembl_genomes] == true
					db_args = { :port => EG_PORT, :host => EG_HOST , :ensembl_genomes => true }
					dummy_db = DummyDBConnection.connect(db_args)
          db_name = nil
      		dummy_connection = dummy_db.connection
					# circumvent naming issues - first find all compara_dbs, then select the correct version
					if args[:fungi]
					  db_name = dummy_connection.select_values("SHOW DATABASES LIKE '%ensembl_compara_fungi%'").select{|p| p.include?("#{release}")}[0]
					elsif args[:bacteria]
					  db_name = dummy_connection.select_values("SHOW DATABASES LIKE '%ensembl_compara_bacteria%'").select{|p| p.include?("#{release}")}[0]
					elsif args[:plants]
					  db_name = dummy_connection.select_values("SHOW DATABASES LIKE '%ensembl_compara_plants%'").select{|p| p.include?("#{release}")}[0]
					elsif args[:protists]
					  db_name = dummy_connection.select_values("SHOW DATABASES LIKE '%ensembl_compara_protists%'").select{|p| p.include?("#{release}")}[0]
					elsif args[:metazoa]
					  db_name = dummy_connection.select_values("SHOW DATABASES LIKE '%ensembl_compara_metazoa%'").select{|p| p.include?("#{release}")}[0]
					else
					  db_name = dummy_connection.select_values("SHOW DATABASES LIKE '%ensembl_compara_pan_homology%'").select{|p| p.include?("#{release}")}[0]
					end
					args[:host] = EG_HOST
					args[:port] = EG_PORT
				# Use the 'normal' Ensembl Compara database					
				else
					release > 47 ? args[:port] = 5306 : args[:port] = 3306
      		dummy_db = DummyDBConnection.connect(args)
      		dummy_connection = dummy_db.connection
					db_name = dummy_connection.select_values("SHOW DATABASES LIKE '%ensembl_compara_#{release}%'")[0]
				end
				# Hack for a local installation...				
				if args[:local]
					args[:host] = 'localhost'
					#db_name = dummy_connection.select_values("SHOW DATABASES LIKE '%ensembl_compara_pan_homology%'").select{|p| p.include?("#{release}")}[0] if db_name.nil?
					db_name = "ensembl_compara_pan_homology_5_#{release}" if db_name.nil?
					args[:username] = 'tools'
					args[:password] = 'analysis'
					args[:port] = 80
				end
				
				establish_connection(
					:adapter => Ensembl::DB_ADAPTER,
					:host => args[:host] || Ensembl::DB_HOST,
					:database => db_name,
					:username => args[:username] || Ensembl::DB_USERNAME,
					:password => args[:password] || Ensembl::DB_PASSWORD,
					:port => args[:port] || 5306
				)
				self.retrieve_connection
			end   
  	      	
  	end # Compara::DBConnection

	end # Compara
  
end
