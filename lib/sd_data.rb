# vim: set filetype=ruby expandtab tabstop=2 shiftwidth=2 tw=80

require 'optparse'
require 'ostruct'
require 'pp'

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
    se = snp_entry.new(data[1], data[2], data[3], data[col_individual])
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
      next if data[0] =~ HEADER_FIRST_COLUMN || data[0] =~ /^$/ 
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
      # Dump an extra line with the original data from SNP detector
      # This is useful to save some of the info lost in the translation
      # I add a hash at the beginning so it can be easily filtered out
      ref = @h_sum[e.chr + e.coor].ref
      puts "##{e.chr} #{e.coor} #{ref} #{e.genotype} #{e.ori}"
      next if skip_snp?(e, a1, a2)
      $stdout.puts "#{e.chr} #{e.coor} #{e.coor} #{ref}/#{a2} #{e.ori}".gsub(/-+/, '-')
    end
  end

  # Skip snps if the genotype matches the reference, or the 
  # sanger software couldn't make a call (no snps covering at that position)
  # If -e, the user wants then all the calls
  def skip_snp?(e, a1, a2)
    ref_at_loci = @h_sum[e.chr + e.coor].ref
    a1 == '.' || a1 == 'n' || (ref_at_loci == a1 && ref_at_loci == a2) ? true : false
  end

  # Extract the indel genotype call from a snp file entry
  # TT : reads from one stand ONLY (look Ori column) support the TT call.
  # TT+: BOTH strands support the evidence (Best case).
  # TT_: The evidence from both strands is DIFFERENT. SNPdetector pics whatever he considers is the best call.
  # nc : We had reads covering the locus, but the quality is not good so we did not call it.
  # .  : I cannot make the call because there are no traces/reads covering the locus.
  def process_indel_call(e)
    cc = e.genotype.gsub(/[\+_]$/, '') # clean genotype call
    return ['.', '.'] if cc[0] == '.'
    [ cc[0,cc.size/2], cc[cc.size/2,cc.size-1] ]
  end
end
