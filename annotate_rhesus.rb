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
      when "a_ensembl"
        SNP_Effect_Predictor.new(@o.input_file, @o.species)
    end
  end

  def parsed_options?
    # Specify options
    opts = OptionParser.new 
    opts.on('-h', '--help')      { usage }
    opts.on("-a", "--action a")  {|a| @o.action     = a }
    opts.on("-s", "--summary s") {|s| @o.sum_file   = s }
    opts.on("-p", "--s_file p")  {|p| @o.snp_file   = p }
    opts.on("-n", "--n_file n")  {|n| @o.input_file = n }
    opts.on("-e", "--specie e")  {|e| @o.species    = e }
            
    opts.parse!(@arguments) rescue return false
    true
  end

  def arguments_valid?
    case @o.action
      when /(to_en_snps|to_an_snps)/
        usage "snp file not provided."          unless @o.snp_file
        usage "Incorrect number of parameters." if @arg_size != 6
        usage "SNP file doesn't exist."         if !File.exists?(@o.snp_file)
      when "a_ensembl"
        usage "Incorrect # of parameters." unless @arg_size == 6
        usage "Input file does not exist." unless File.exists?(@o.input_file)
      else
        usage "Incorrect action."
    end
    true
  end

end

# Encapsulates the call to the perl script from ensembl that 
# performs the necessary queries against ensembl to do the annotation
#
# snp_effect_predictor.pl -i input.ensembl.format.txt -o out.txt -s macaque
#
class SNP_Effect_Predictor
  
  include Misc

  SCRIPT_NAME = "snp_effect_predictor.pl"
  RANDOM_FILE = "/tmp/log.annotation.#{rand(1000)}.txt"

  def initialize(i_fn, specie)
    @input_fn = i_fn
    # We have to make sure we have the perl script in our path
    msg_not_found = "I can't find (#{SCRIPT_NAME}). Also, make sure it is executable."
    usage msg_not_found if script_not_in_path?

    msg_not_run = "I tried to run #{SCRIPT_NAME}, but I had a problem. "
    msg_not_run << "Details in: #{RANDOM_FILE}"
    usage msg_not_run unless running_not_fine?
  end

  def run
    puts "#{SCRIPT_NAME} -i #{i_fn} -o #{i_fn}.annotated -s #{specie}"
  end

  private
  
  def script_not_in_path?
    ENV['PATH'].split(':').each do |dir| 
      f_path = dir + '/#{SCRIPT_NAME}'
      return true if File.exists?(f_path) && File.executable?(f_path)
    end
    false
  end

  # Run the script in a shubshell and capture the stdout/err and the exit code
  #
  def running_not_fine?
    output = "" 
    # Try to run the script
    system("#{SCRIPT_NAME} &> #{RANDOM_FILE}")
    # If failed, dump the running attempt output to a file
    if $?.exitstatus != 0
      begin
        File.open(RANDOM_FILE, "w") {|f| f.puts output}
      rescue
        usage "I cannot create log file: #{RANDOM_FILE}."
      end
      false
    else
      true
    end 
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
  -n: input file for annotation
  -e: species to use when annotating

Valid actions:

  to_en_snps : Convert SNP detector SNPs (substitutions) to ensembl format. 
  to_an_snps : Convert SNP detector SNPs (substitutions) to annovar format. 
  a_ensembl  : Annotate a set of SNPs (substitions).

Examples:

  # Convert to ensembl the snp calls from the sanger pipeline
  $ annotate_rhesus.rb -a to_en_snps -s summary.csv -p snps.csv > out.txt

  # Annotate a set of snps 
  $ annotate_rhesus.rb -a a_ensembl -n ./input.txt -e macaque > input.annotated.txt
