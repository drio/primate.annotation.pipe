#!/usr/bin/env ruby19
#
require 'optparse'
require 'ostruct'
require 'pp'

module Misc
  def log(msg)
    puts "#{Time.new}: #{msg}"
  end

  def usage(msg=nil, ec=1)
    printf "UPS!: #{msg}\n\n" if msg
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
  end

  def process
    if parsed_options? && arguments_valid? 
      log "Ok.. ready"
      run 
    else 
      usage
    end
  end 

  private

  def run
    case @o.action
      when "to_an_snps"
        SD_Data.new(@o.sum_file, @o.snp_file).convert_to("annovar")
      when "to_en_snps"
        SD_Data.new(@o.sum_file, @o.snp_file).convert_to("ensembl")
    end
  end

  def parsed_options?
    # Specify options
    opts = OptionParser.new 
    opts.on('-h', '--help')      { usage }
    opts.on("-a", "--action a")  {|a| @o.action     = a }
    opts.on("-s", "--summary s") {|s| @o.sum_file   = s }
    opts.on("-p", "--s_file p")  {|p| @o.snp_file   = p }
            
    opts.parse!(@arguments) rescue return false
    true
  end

  def arguments_valid?
    case @o.action
      when /(to_en_snps|to_an_snps)/
        usage "snp file not provided."          unless @o.snp_file
        usage "Incorrect number of parameters." if @arg_size != 6
        usage "SNP file doesn't exist."         if !File.exists?(@o.snp_file)
      else
        usage "Incorrect action."
    end
    true
  end

end

# Handles data form SNP detector (sanger pipeline)
#
class SD_Data 

  HEADER_FIRST_COLUMN = /SNP_ID/
  # The reference used for the calls (/data/genome/Mmul2) has
  # a contig not supported in NCBI ?
  MODIFIED_CHRM       = /modified/ 

  def initialize(sum_fn, snps_fn)
    @sum_fn  = sum_fn
    @snps_fn = snps_fn
    @h_sum   = {} # We'll save summary data
  end

  # Converts Sanger format to the annotation tool of choice
  def convert_to(o_type)
    load_summary
    dump_input_file(o_type)
    #dump_merged_sum_and_snps
  end

  private

  # Hash one single entry from the summary file
  #
  def add_sum_entry(d)
    sum_entry = Struct.new(:id, :chr, :coor, :ori, :flank, :allele)
    se        = sum_entry.new(d[0], d[1].gsub(/chr/i, ''), d[2], d[3], d[4], d[7])
    @h_sum[se.chr + se.coor] = se
  end

  # Hash all the reference values for all the snps found
  #
  def load_summary
    File.open(@sum_fn, "r").each do |l|
      data = l.split(",")
      next if data[0] =~ HEADER_FIRST_COLUMN
      add_sum_entry(data)
    end
  end

  # Given a SNP, dump the genotypes of all the samples with the
  # reference genotype at the end
  #
  def merge_line(data)
    chr, coor    = data[1], data[2]
    snps_columns = data.join(" ").chomp
    ref_column   = " (" + @h_sum[chr + coor].allele + ")"
    snps_columns + ref_column
  end

  # Iterate over all the snps in the snp file (snp detector output)
  # and print it adding the refence at the position of the SNP 
  # (last column)
  #
  def dump_merged_sum_and_snps
    puts "snp_id chrm coor orienta samples-genotype ref(alleles)"
    File.open(@snps_fn, "r").each do |l|
      data = l.split(",")
      next if data[0] =~ HEADER_FIRST_COLUMN || l =~ /^\n/
      $stdout.puts merge_line(data)
    end
  end

  # Chr	Start	End	Ref	Obs	Comments
  #
  def format_line(data, o_type)
    chr, coor    = data[1], data[2]
    ref          = @h_sum[chr + coor].allele[0]
    obs          = @h_sum[chr + coor].allele[2]
    case o_type
      when "annovar" # chr2 99555161 99555161 A C comments:
        "chr" + chr + " " + coor + " " + coor + " " + ref  + " " + obs + " comments: "
      when "ensembl" # chr2 99554095 99554095 G/A
        "chr" + chr + " " + coor + " " + coor + " " + ref  + "/" + obs
    end
  end

  # Dump annovar input file from SNP detector data
  # chr2 99560072 99560072 G A comments: 
  #
  def dump_input_file(o_type)
    skipped = 0
    File.open(@snps_fn, "r").each do |l|
      data = l.split(",")
      next if data[0] =~ HEADER_FIRST_COLUMN || l =~ /^\n/
      if data[1] =~ MODIFIED_CHRM
        skipped = skipped + 1 
        next
      end
      $stdout.puts format_line(data, o_type)
    end
    $stderr.puts "WARNING: #{skipped} calls skipped (contig names) " if skipped > 0
  end
end

# MAIN
#
Arguments.new(ARGV).process

__END__
Usage:
  annotate_rhesus.rb <options>    

  -h: help  
  -a: action to perform (xxxxx)
  -s: summary file from the sanger pipeline
  -p: snps file from the sanger pipeline
  -i: indels file from the sanger pipeline 

Valid actions:

  to_en_snps : Convert SNP detector SNPs (substitutions) to ensembl format. 
  to_an_snps : Convert SNP detector SNPs (substitutions) to annovar format. 

Examples:

  # Convert to ensembl the snp calls from the sanger pipeline
  $ annotate_rhesus.rb -a to_en_snps -s summary.csv -p snps.csv > out.txt
