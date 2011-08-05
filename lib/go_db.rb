
GO_DB_ADAPTER = 'mysql'
GO_DB_HOST = 'mysql.ebi.ac.uk'
GO_PASSWORD = 'amigo'
GO_USER = 'go_select'
GO_DB = 'go_latest'
GO_PORT = 4085

module Toolbox
  
  module GoDB    
    
    class DBConnection < ActiveRecord::Base
      self.abstract_class = true
    
      def self.connect

        establish_connection(
                              :adapter => GO_DB_ADAPTER,
                              :host => GO_DB_HOST,
                              :database => GO_DB,
                              :username => GO_USER,
                              :password => GO_PASSWORD,
                              :port => GO_PORT
                            )
      end
    
    end
    
    class Term < DBConnection
      
      set_primary_key "id"
      has_one :term_definition
      has_many :term_descendents, :foreign_key => "descendent_id"
      
      def ancestor
        return Term.find_by_id(TermAncestor.find_by_acc(self.acc).ancestor_id)
      end
      
      def descendents
        return TermDescendent.find_all_by_descendent_id(self.id).collect{|d| d.term }
      end
      
      def root
      	if self.is_root == 1
      		return self.ancestor
      	else
      		self.ancestor.root
      	end
      end
      
    end
    
    class TermDefinition < DBConnection
      belongs_to :term, :foreign_key => "term_id"
    end
    
    class Term2Term < DBConnection
      set_table_name 'term2term'
    end
    
    class TermAncestor < DBConnection
      
    end
    
    class TermDescendent < DBConnection
     	belongs_to :term, :foreign_key => "descendent_id"
    end  
    
    class TermDefinition < DBConnection
      
    end
    
  end
  
end
