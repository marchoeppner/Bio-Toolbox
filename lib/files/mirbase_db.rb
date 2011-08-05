
Mirbase_DB_ADAPTER = 'mysql'
Mirbase_DATABASE = "mirbase"
Mirbase_DB_HOST = 'localhost'
Mirbase_DB_USERNAME = 'tools'
Mirbase_DB_PASSWORD = 'analysis'

module Toolbox 
  
  module MirbaseDB
    
    require 'composite_primary_keys'
    
    class DBConnection < ActiveRecord::Base
      self.abstract_class = true
    
      def self.connect(version="14")

        establish_connection(
                              :adapter => Mirbase_DB_ADAPTER,
                              :host => Mirbase_DB_HOST,
                              :database => "#{Mirbase_DATABASE}#{version}",
                              :username => Mirbase_DB_USERNAME,
                              :password => Mirbase_DB_PASSWORD
                              #:port => port
                            )
      end
    
    end
    
    class Mirna < DBConnection
      set_primary_key "auto_mirna"
      belongs_to :mirna_species, :foreign_key => "auto_species"
      has_one :mirna_to_prefam, :foreign_key => "auto_mirna"
      has_one :mirna_prefam, :through => :mirna_to_prefam
      has_many :mirna_pre_matures, :foreign_key => "auto_mirna"
      has_many :mirna_matures, :through => :mirna_pre_matures
      has_many :mirna_databases, :foreign_key => "auto_mirna"
      
    end
    
    class MirnaSpecies < DBConnection
      set_primary_key "auto_id"
      has_many :mirnas
    end
    
    class MirnaToPrefam < DBConnection
      set_primary_keys :auto_mirna, :auto_prefam
      belongs_to :mirna, :foreign_key => "auto_mirna"
      belongs_to :mirna_prefam, :foreign_key => "auto_prefam"
    end
    
    class MirnaPrefam < DBConnection
      set_primary_key "auto_prefam"
      has_many :mirna_to_prefams, :foreign_key => "auto_prefam"
      has_many :mirnas, :through => :mirna_to_prefams
    end
    
    class MirnaPreMature < DBConnection
    	belongs_to :mirna, :foreign_key => "auto_mirna"
    	belongs_to :mirna_mature, :foreign_key => "auto_mature"
    end
    
    class MirnaMature < DBConnection
    	set_primary_key 'auto_mature'
    	has_many :mirna_pre_matures
    	has_many :mirnas, :through => :mirna_pre_matures
    	has_many :mirna_target_links, :foreign_key => "auto_mature"
    	
    	def ebi_link
    		return self.mirna_target_links.collect{|l| l.put_db_url}.select{|d| d.include?("ebi")}[0]
			end
    end
    
    class MirnaDatabaseLink < DBConnection
    	belongs_to :mirna, :foreign_key => "auto_mirna"
    end
    
    class MirnaTargetUrl < DBConnection
    	set_primary_key "auto_db"
    	has_many :mirna_target_links, :foreign_key => "auto_db"
    
    end
    
    class MirnaTargetLink < DBConnection
    	belongs_to :mirna_mature, :foreign_key => "auto_mature"
    	belongs_to :mirna_target_url, :foreign_key => "auto_db"
    	
    	def put_db_url
    		db_string = self.mirna_target_url.url
    		return db_string.gsub(/\=\<\?\>/, "=#{self.field1}")
    		
    	end
    end
  end
  
end
