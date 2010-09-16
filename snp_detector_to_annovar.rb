#!/usr/bin/env ruby19
#
HEADER_FIRST_COLUMN = /SNP_ID/
MODIFIED_CHRM       = /modified/ # The reference used for the calls (/data/genome/Mmul2) has
                                 # a contig not supported in NCBI ?

@sum_fn  = "rhesus.summary.snps.txt"
@snps_fn = "rhesus.snp_detector.snps.txt"
@h_sum   = {}

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

# Given a SNP print, the genotypes of all the samples with the
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
def annovar_line(data)
  chr, coor    = data[1], data[2]
  ref          = @h_sum[chr + coor].allele[0]
  obs          = @h_sum[chr + coor].allele[2]
  "chr" + chr + " " + coor + " " + coor + " " + ref  + " " + obs + " comments: " #+ merge_line(data)
end

# Dump annovar input file from SNP detector data
#
def dump_annovar_input_file
  skipped = 0
  File.open(@snps_fn, "r").each do |l|
    data = l.split(",")
    next if data[0] =~ HEADER_FIRST_COLUMN || l =~ /^\n/
    if data[1] =~ MODIFIED_CHRM
      skipped = skipped + 1 
      next
    end
    $stdout.puts annovar_line(data)
  end
  $stderr.puts "WARNING: #{skipped} calls skipped (contig names) " if skipped > 0
end

# MAIN
#
load_summary
dump_merged_sum_and_snps
#dump_annovar_input_file
