
Taxonomy_DB_ADAPTER = 'postgresql'
Taxonomy_DATABASE = "taxonomy"
Taxonomy_DB_HOST = 'localhost'
Taxonomy_DB_USERNAME = 'tools'
Taxonomy_DB_PASSWORD = 'analysis'


module Toolbox 
  module TaxonomyDB
    
    include ActiveRecord
    
    class DBConnection < ActiveRecord::Base
      self.abstract_class = true
    
      def self.connect

        establish_connection(
                              :adapter => Taxonomy_DB_ADAPTER,
                              :host => Taxonomy_DB_HOST,
                              :database => "#{Taxonomy_DATABASE}",
                              :username => Taxonomy_DB_USERNAME,
                              :password => Taxonomy_DB_PASSWORD
                              #:port => port
                            )
      end
    
    end
    
    class Name < DBConnection
    	belongs_to :node, :foreign_key => "node_id"
    	set_table_name 'name'
    	def get_name_of_rank(rank)
    		self.node.get_name_of_rank(rank)
    	end
    end
    
    class Gencode < DBConnection
    
    end
    
    class Node < DBConnection
    set_table_name 'node'
    	has_many :names
    	
    	def tax_string
    	
    		return self.get_all_parents.join(";")
    	
    	end
    	
    	def get_all_parents(parents=[])
    		
    		answer = []
    		
    		if self.get_parent and self.get_parent.scientific_name != "root"
    			parents << self.get_parent.scientific_name
    			self.get_parent.get_all_parents(parents)
    		end
    		
    		return parents
    	end
    	
    	def get_parent
    		return TaxonomyDB::Node.find(self.parent_tax_id)
    	end
    	
    	def get_children
    		return TaxonomyDB::Node.find_all_by_parent_tax_id(self.id)
    	end
    	
    	def scientific_name
    		return self.names.select{|n| n.name_class == "scientific name"}[0].name
    	end
    	
    	def child_of_node?(scientific_name,verbose=false)

    		if "#{self.scientific_name}" == "#{scientific_name}"
    			puts "\tRanks match, returning true" if verbose
					return true
					exit
				else
					puts "\tRanks don't match (#{self.scientific_name} vs #{scientific_name}) - next" if verbose
    			parent = self.get_parent
    			puts "\t\tParent: #{parent.scientific_name}" if verbose
    			if parent and parent.scientific_name != "root"
    				parent.child_of_node?(scientific_name,verbose)
    			else
    				puts "\t\tNo parent found or parent is root" if verbose
    				return false
    				break
    			end
    		end
    		
    	end  		
    	
    	def get_name_of_rank(rank)
    	
    		if self.rank == rank
    			return self.scientific_name
    		elsif self.get_parent and self.parent_tax_id > 1
    			self.get_parent.get_name_of_rank(rank)
    		end
    		
    	end
    	
    	@@answer = []
    	def get_children_by_rank(rank)
				puts "Trying #{self.scientific_name} - #{self.rank}"
    		self.get_children.each do |child|
    		
					if child.rank == rank    		
						puts child.scientific_name
						@@answer << child.scientific_name unless @@answer.include?(child.scientific_name)
    			elsif child.rank != "no rank"
    				puts "Rank is #{child.rank} - continuing"
    				child.get_children_by_rank(rank)
    			end
    			
    		end
    		
    		return @@answer

    	end

    end
    
	end
	
end
