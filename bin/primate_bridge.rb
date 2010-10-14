#!/usr/bin/env ruby19
#
require 'optparse'
require 'ostruct'
require 'pp'

$: << File.join(File.dirname(File.dirname($0)), "lib")
require 'load_libs'

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

