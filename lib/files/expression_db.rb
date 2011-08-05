
Expression_DB_ADAPTER = 'postgresql'
Expression_DB_HOST = 'localhost'
Expression_PASSWORD = 'analysis'
Expression_USER = 'tools'
Expression_DB = 'expression'
#Expression_PORT = 4085

module Toolbox
  
  module ExpressionDB    
    
    class DBConnection < ActiveRecord::Base
      self.abstract_class = true
    
      def self.connect

        establish_connection(
                              :adapter => Expression_DB_ADAPTER,
                              :host => Expression_DB_HOST,
                              :database => Expression_DB,
                              :username => Expression_USER,
                              :password => Expression_PASSWORD
                              #:port => Expression_PORT
                            )
      end
    
    end
    
    class Dataset < DBConnection
      has_many :genes
      has_many :tissues
      has_many :expressions
    end
    
    class Gene < DBConnection
      belongs_to :dataset, :foreign_key => "dataset_id"
      has_many :xref_gene_tissues
      has_many :tissues, :through => :xref_gene_tissues
      
      def expression_by_tissue_group
        answer = {}
        self.xref_gene_tissues.each do |x|
          tg = x.tissue.tissue_group
          answer.has_key?(tg) ? answer[tg] << x.value : answer[tg] = [ x.value ]
        end
        answer.each do |tg,values|
          values.nitems == 2 ? answer[tg] = Math.log2((values[0]+values[1])/self.median) : answer[tg] = Math.log2(values.shift/self.median)
        end
        return answer
      end
      
      def value_by_tissue_group
        groups = {}
        self.xref_gene_tissues.each do |xref|
          tg = xref.tissue.tissue_group
          groups.has_key?(tg) ? groups[tg] << xref.value : groups[tg] = [ xref.value ]
        end       
        groups.each do |group,values|
          if values.nitems == 2
            groups[group] = (values[0]+values[1])/2
          else
            groups[group] = values[0]
          end
        end      
        return groups    
      end
      
      def tissues_by_group
        groups = {}
        self.tissues.each do |t|
          groups.has_key?(t.tissue_group) ? groups[t.tissue_group] << t : groups[t.tissue_group] = [ t ]
        end
        return groups
      end
      
      def median
        values = self.value_by_tissue_group.collect{|g,v| v}
        return values[values.nitems.to_f/2.0]
      end
      
      def sum
        sum = 0
        self.value_by_tissue_group.each do |g,v|
          sum += v
        end
        return sum
      end
      
      def mean
        sum = 0
        groups = self.expression_by_tissue_group
        groups.each {|g,v| sum += v }
        value = "#{sum/groups.keys.nitems}"
        return "#{value.slice(/[0-9]*\.[0-9][0-9][0-9]/).to_f}"
      end
          
      def entropy
        sum = self.sum
        answer = 0
        #Entropy = Sum (Pi * log (Pi)), where Pi  = Ei/total expression
        self.value_by_tissue_group.each do |t,v|
          p = v/sum
          answer += p*Math.log(p)
        end
        return answer*-1
      end
      
      def specificity
        max = nil
        sum = 0        
        groups = self.value_by_tissue_group
        groups.each do |group,val|
          sum += val
          max = val if val > max or max.nil?
        end       
        return Math.log2(max/sum)
      end
    end
    
    class XrefGeneTissue < DBConnection
      belongs_to :gene, :foreign_key => "gene_id"
      belongs_to :tissue, :foreign_key => "tissue_id"
    
      def expression
        return Math.log2(self.value/self.gene.mean)
      end
      
    end
    
    class Tissue < DBConnection
      belongs_to :tissue_group, :foreign_key => "tissue_group_id"
      has_many :xref_gene_tissues
      has_many :genes, :through => :xref_gene_tissues
    end  
    
    class TissueGroup < DBConnection
      belongs_to :dataset, :foreign_key => "dataset_id"
      has_many :tissues
      has_many :xref_gene_tissues, :through => :tissues
      
      def values
        values = []
        self.tissues.each do |t|
          t.xref_gene_tissues.each do |x|
            values << x.value
          end
        end
        return values.sort
      end
      
      def max
        self.values[-1]
      end
      
      def min
        self.values[0]
      end
      
      def percentile(value)
        values = self.values
        return values.select{|v| v >= value}.nitems.to_f/values.nitems.to_f*100
      end
        
      def median
        values = []
        self.tissues.each do |tissue|
          tissue.xref_gene_tissues.each do |x|
            values << x.value
          end
        end
        return values.sort[values.nitems/2]
      end
        
      def average   
        values = 0
        sum = 0
        self.tissues.each do |tissue|
          tissue.xref_gene_tissues.each do |x|
            values += 1
            sum += x.value
          end
        end
        return sum.to_f/values.to_f
      end
      
    end
    
  end
  
end
