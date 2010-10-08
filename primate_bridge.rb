#!/usr/bin/env ruby19
#
require 'optparse'
require 'ostruct'
require 'pp'

module Misc
  def log(msg)
    $stderr.puts "#{Time.new}: #{msg}"
  end

  def usage(msg=nil, ec=1)
    printf "\nUPS!: #{msg}\n\n" if msg
    puts DATA.read
    exit(ec)
  end
end

# Process user input arguments
#
class Arguments

  include Misc

  def initialize(arguments) 
    @arguments = arguments
    @arg_size  = arguments.size
    @o         = OpenStruct.new
    @o.indels  = FALSE
  end

  def process
    if parsed_options? && arguments_valid? 
      log "summary file : #{@o.sum_file}" if @o.sum_file
      log "snp file : #{@o.snp_file}" if @o.snp_file
      run 
    else 
      usage
    end
  end 

  private

  def run
    case @o.action
      when "to_an_snps"
        SD_Data.new(@o.sum_file, @o.snp_file).convert_to("annovar", 
                                                         @o.indels, 
                                                         @o.individual)
      when "to_en_snps"
        SD_Data.new(@o.sum_file, @o.snp_file).convert_to("ensembl",
                                                         @o.indels,
                                                         @o.individual)
    end
  end

  def parsed_options?
    # Specify options
    opts = OptionParser.new 
    opts.on('-h', '--help')      { usage }
    opts.on("-a", "--action a")  {|a| @o.action     = a }
    opts.on("-s", "--summary s") {|s| @o.sum_file   = s }
    opts.on("-f", "--snp_file f"){|f| @o.snp_file   = f }
    opts.on("-i", "--indels")    {    @o.indels     = TRUE }
    opts.on("-p", "--focus_on p"){|p| @o.individual = p}
            
    opts.parse!(@arguments) rescue return false
    true
 end

  def arguments_valid?
    case @o.action
      when /(to_en_snps|to_an_snps|to_an_indels)/
        usage "Incorrect number of parameters." unless @arg_size.to_s =~ /4|5|8|9/
        usage "Summary file doesn't exist."     if @o.individual.nil? && 
                                                   !File.exists?(@o.sum_file)
        usage "Snp file doesn't exist."         if @o.individual && 
                                                   !File.exists?(@o.snp_file)
      else
        usage "Incorrect action."
    end
    true
  end

end

# Handles data form SNP detector (sanger pipeline)
#
class SD_Data 
  include Misc

  HEADER_FIRST_COLUMN = /SNP_ID|INDEL_ID/
  # The reference used for the calls (/data/genome/Mmul2) has
  # a contig not supported in NCBI ?
  MODIFIED_CHRM       = /modified/ 

  def initialize(sum_fn, snps_fn)
    @sum_fn  = sum_fn
    @snps_fn = snps_fn
    @h_sum   = {} # We'll save summary data
    @h_snps  = {} # We'll snp data
  end

  # Converts Sanger format to the annotation tool of choice
  def convert_to(o_type, indels, individual)
    log "loading summary."
    load_summary(indels)
    log "#{@h_sum.size} snps(indels=#{indels}) loaded."

    if individual # they want us to work on a single individual
      log "Working in individual mode (#{individual})"
      col_individual = find_col_for_individual(individual)
      log "Loading Individual: #{individual} column: #{col_individual} "
      load_snps(indels, col_individual)
      log "dumping input file (For single individual) in new format."
      dump_input_file_for_individual(indels)
    else
      log "dumping input file in new format."
      dump_input_file(o_type)
    end
  end

  private
  # Add a new snp entry into the hash of snps
  # 
  def add_snp_entry(data, indels, col_individual)
    snp_entry = Struct.new(:chr, :coor, :ori, :genotype)
    # Get the columns we need 
    se    = snp_entry.new(data[1], data[2], data[3], data[col_individual])
    @h_snps[se.chr + se.coor] = se
  end

  # Hash the snp information for the specific individual
  #
  def load_snps(indels, col_individual)
    File.open(@snps_fn, "r").each_with_index do |l, i|
      data = l.split(",")
      next if data[0] =~ HEADER_FIRST_COLUMN || data.size == 1
      add_snp_entry(data, indels, col_individual)
    end
  end

  # Find what column has the genotype of the individual
  #   
  def find_col_for_individual(ind)
    col_i = nil
    File.open(@snps_fn).readlines[0].split(',').each_with_index do |i, pos|
      col_i = pos if i == ind
    end
    usage "I couldn't find individual: #{ind}" unless col_i
    return col_i
  end

  # Hash one single entry from the summary file
  #
  def add_sum_entry(d, indels)
    sum_entry = Struct.new(:id, :chr, :coor, :ori, :flank, :ref, :var, :class, :length)
    se = if indels
      sum_entry.new(d[0], d[1].gsub(/chr/i, ''), d[2], d[3], d[4], 
                    # indels have multiple '-' in the ref allele, make sure it is only 1
                    #d[9].gsub(/-+/,'-'), d[10].gsub(/-+/,'-'), 
                    d[9], d[10], d[11], d[12])
    else
      sum_entry.new(d[0], d[1].gsub(/chr/i, ''), 
                    d[2], d[3], d[4], d[7][0], d[7][2], 'subs', 0)
    end
    # Notice we can have multiple entries with the same SNP id and coor and
    # position. Since the calls are exactly the same we can just allow collisions.
    @h_sum[se.chr + se.coor] = se
  end

  # Hash all the reference values for all the snps found
  #
  def load_summary(indels)
    File.open(@sum_fn, "r").each_with_index do |l, i|
      data = l.split(",")
      next if data[0] =~ HEADER_FIRST_COLUMN
      add_sum_entry(l.split(","), indels)
    end
  end

  # Ensembl format
  # 12  1017956   1017956   T/A   +
  # 8   12601     12600     -/C   +
  # 8   12600     12602     CGT/- -
  #
  # Annovar format
  # chr2 99555161 99555161 A C comments:
  def format_line(data, o_type)
    case o_type
      when "annovar" 
        "chr" + data.chr + " " + data.coor + " " + data.coor + " " +
        data.ref  + " " + data.var + " comments: "
      when "ensembl" # chr2 99554095 99554095 G/A
        "chr" + data.chr + " " + data.coor + " " + data.coor + " " +
        data.ref  + "/" + data.var + " " + data.ori
    end
  end

  # Dump annovar input file from SNP detector data
  #
  def dump_input_file(o_type)
    @h_sum.each {|k, e| $stdout.puts format_line(e, o_type) }
  end

  def dump_input_file_for_individual(indels)
    @h_snps.each do |k, e|
      a1, a2 = indels ? process_indel_call(e) : [e.genotype[0], e.genotype[1]]
      next if skip_snp?(e, a1, a2)
      ref = @h_sum[e.chr + e.coor].ref
      $stdout.puts "#{e.chr} #{e.coor} #{e.coor} #{ref} #{a2} #{e.ori}".gsub(/-+/, '-')
    end
  end

  # Skip snps if the genotype matches the reference, or the 
  # sanger software couldn't make a call (no snps covering at that position)
  def skip_snp?(e, a1, a2)
    ref_at_loci = @h_sum[e.chr + e.coor].ref
    a1 == '.' || (ref_at_loci == a1 && ref_at_loci == a2) ? true : false
  end

  # Extract the indel genotype call from a snp file entry
  def process_indel_call(e)
    cc = e.genotype.gsub(/[\+_]$/, '') # clean genotype call
    return ['.', '.'] if cc[0] == '.'
    [ cc[0,cc.size/2], cc[cc.size/2,cc.size-1] ]
  end
end

# MAIN
#
Arguments.new(ARGV).process

__END__
Usage:
  annotate_rhesus.rb <options>    

  -h: help  
  -a: action to perform (to_en_snps|to_an_snps)
  -s: summary file from the sanger pipeline
  -i: input are indel calls
  -p: focus on a single individual
  -f: snp file from the sanger pipeline

Valid actions:

  to_en_snps : Convert SNP detector SNPs (substitutions) to ensembl format. 
  to_an_snps : Convert SNP detector SNPs (substitutions) to annovar format. 

Examples:

  # Convert to ensembl the snp calls from the sanger pipeline
  $ annotate_rhesus.rb -a to_en_snps -s summary.csv > out.txt

  # Same but with indels
  $ annotate_rhesus.rb -a to_en_snps -s indels.summary.csv -i > out.txt

  # Convert to ensembl the snp calls from the snager pipeline 
  # (for a particular individual)
  $ annotate_rhesus.rb -a to_en_snps -s summary.csv -f snps.csv -p 1-173  > out.txt

  # Same for indels
  $ annotate_rhesus.rb -a to_en_snps -s summary.csv -f snps.csv -p 1-173 -i > out.txt

