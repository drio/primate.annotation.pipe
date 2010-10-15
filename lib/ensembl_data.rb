$: << File.join(File.dirname(File.dirname($0)), "lib")
require 'load_libs'

# Encapsulates a line of annotated data from ensembl
# Examples:
#
# 2_99515799_T/A  2:99515799      ENSMMUG00000022533 ENSMMUT00000031713      INTRONIC        -       -       -       -
#
# 2_99516432_G/A  2:99516432      ENSMMUG00000022533 ENSMMUT00000031714      SYNONYMOUS_CODING       1581    527     N -
# 2_99516432_G/A  2:99516432      ENSMMUG00000022533 ENSMMUT00000031713      SYNONYMOUS_CODING       1173    391     N -
# 
# Uploaded Variation      Location        Gene    Transcript      Consequence     
# Position in cDNA        Position in protein     Amino acid change       Corresponding Variation
#
# lff output:
#
#Annotated SNP   2_99515799_T/A  Ensembl SNP     chr2    99515799        99515799        +       .       1       .      # .       nonref allele=A;placement=intronic;gene_loc=-;aa_loc=-;aa=-;transcripts=-;      T
#
#Annotated SNP   2_99516432_G/A  Ensembl SNP     chr2    99516432        99516432        +       .       1       .      # .       nonref allele=A;placement1=synonymous_coding;placement2=synonymous_coding;
#         gene_loc1=1581;gene_loc2=1173;
#         aa_loc1=527;aa_loc2=391;
#         aa1=N;aa2=N;
#         transcripts1=-;transcripts2=-;      
# G

class Ensembl_Line
  attr_reader :up_variation

  SEP   = "\t"
  C_SEP = ";"

  def initialize(line)
    @line = line
    set_attributes
  end

  def first_part
    o = "Annotated" + SEP + "SNP" + SEP
    o << @up_variation + SEP + "Ensembl SNP" + SEP
    o << "chr" + @location.split(":")[0] + SEP
    o << @location.split(":")[1] + SEP + @location.split(":")[1] + SEP
    o << "+" + SEP + "." + SEP + "1" + SEP + "." + SEP
  end

  def non_ref
    "nonref_allele=" + @up_variation.split('/')[0][-1] + C_SEP
  end

  def second_part(i)
    o = "gene#{i}="         + @gene        + C_SEP
    o << "transcript#{i}="  + @transcript  + C_SEP
    o << "consequence#{i}=" + @consequence + C_SEP
    o << "cdna_pos#{i}="    + @cdna_pos    + C_SEP
    o << "prot_pos#{i}="    + @prot_pos    + C_SEP
    o << "a_change#{i}="    + @aa_change   + C_SEP
    o << "c_variation#{i}=" + @c_variation + C_SEP
  end

  def ref
    @up_variation.split('/')[1]
  end

private
  def set_attributes
    @columns      = @line.chomp.split("\t")
    
    @up_variation = @columns[0] # 2_99516432_G/A  
    @location     = @columns[1] # 2:99516432 
    @gene         = @columns[2] # ENSMMUG00000022533 
    @transcript   = @columns[3] # ENSMMUT00000031714    
    @consequence  = @columns[4] # SYNONYMOUS_CODING 
    @cdna_pos     = @columns[5] # 1581
    @prot_pos     = @columns[6] # 527 
    @aa_change    = @columns[7] # N
    @c_variation  = @columns[8] # -
  end
end

class Ensembl_Annotated_Data
  SEP   = "\t"
  C_SEP = ";"

  def initialize
    @h = Hash.new { |hash, key| hash[key] = [] }
  end

  def add(line)
    en_line = Ensembl_Line.new(line)
    @h[en_line.up_variation] << en_line
  end

  def to_lff
    lff_output = ""
    @h.each do |k, v|
      lff_output << v[0].first_part + SEP
      lff_output << v[0].non_ref   
      v.each_with_index {|ae, i| lff_output << ae.second_part(i)}
      lff_output << SEP +  v[0].ref + "\n"
    end
    lff_output
  end
end

