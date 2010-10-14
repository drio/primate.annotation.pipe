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
# Uploaded Variation, Location, Transcript, Consequence, Relative position in protein, Amino acid change, Corresponding Variation
#
class Ensembl_Line
  def initialize(line)
  end
end

