#!/usr/bin/ruby

require 'bio'

lines = IO.readlines(ARGV.shift)
genomes = ARGV.shift

proteins = lines.shift.split(",").collect{|e| e.strip}
proteins.shift

lines.each do |line|
  
  elements = line.split(",").collect{|e| e.strip }
  next if elements.empty?
    
  organism = elements.shift
  
  puts "#{organism}"
  puts "#{elements.join("<->")}"
  
  proteins.each do |protein|
    
    this_acc = elements.shift.strip
    
    next if this_acc.length < 3
    
    bin = []
    IO.popen("grep \"#{this_acc}\" -A1 #{genomes}") do |f|
      while line = f.gets do
        bin << line.strip
      end
    end
    
    bin.delete_if{|x| x == "--"}
    seqs = []
    
    until bin.empty?
      name = bin.shift
      seq = bin.shift
      puts "#{name} / #{seq}" if seq.nil?
      seqs << Bio::FastaFormat.new(Bio::Sequence::AA.new(seq).to_fasta(name))
    end
    
    if seqs.empty?
      puts "#{organism} -> #{this_acc} NOT FOUND"
    else
      puts "#{organism} -> #{this_acc} FOUND (#{seqs.collect{|s| s.definition}.join(";")})"  
    end  
    
  end
end
    