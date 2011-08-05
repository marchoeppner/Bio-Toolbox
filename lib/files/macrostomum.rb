Macrostomum_DB_ADAPTER = 'mysql'
Macrostomum_DATABASE = "macrostomum"
Macrostomum_DB_HOST = 'localhost'
Macrostomum_DB_USERNAME = 'lignano'
Macrostomum_DB_PASSWORD = 'insert heredum'

require 'active_record'

module MacrostomumDB
  
  include ActiveRecord
  
  class DBConnection < ActiveRecord::Base
    self.abstract_class = true
  	self.pluralize_table_names = false
  	
    def self.connect

      establish_connection(
                            :adapter => Macrostomum_DB_ADAPTER,
                            :host => Macrostomum_DB_HOST,
                            :database => "#{Macrostomum_DATABASE}",
                            :username => Macrostomum_DB_USERNAME,
                            :password => Macrostomum_DB_PASSWORD
                            #:port => port
                          )
    end
  
  end
  
  class Campain < DBConnection
    
  end
  
  class People < DBConnection
    
  end
  
  class Sample < DBConnection
    
  end
  
  class Site < DBConnection
    
  end
  
end