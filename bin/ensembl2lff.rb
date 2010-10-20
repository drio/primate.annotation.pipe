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
  end

  def process
    if parsed_options? && arguments_valid? && !STDIN.tty?
      run 
    else 
      usage
    end
  end 

  private

  def run
    e_data = Ensembl_Annotated_Data.new
    log "Iterating over input file" 
    STDIN.each {|l| e_data.add(l) unless l =~ /^Up/ }
    log "Converting to LFF"
    puts e_data.to_lff
  end

  def parsed_options?
    # Specify options
    opts = OptionParser.new 
    opts.on('-h', '--help')      { usage }
            
    opts.parse!(@arguments) rescue return false
    true
 end

  def arguments_valid?
    usage "Wrong number of arguments." unless @arg_size == 0
    true
  end
end

# MAIN
#
Arguments.new(ARGV).process

__END__
Usage:
  ensembl2lff.rb <option>

  -h: help  

Example(s):

  $ cat ensembl.annotate.txt | ensembl2lff.rb > output.lff
