#!/bin/bash
#
set -e

error()
{
  echo "Ups!: $1"
  echo "`basename $0` <input_file> <path_to_ensembl_perl_libs> <species>"
  echo ""
  echo "Example:"
  echo "$0 ./input.txt /tmp/ensembl_libs macaque"
  exit 1
}

input=$1
en_mod_libs="$2"
species=$3
sd_main_path="`dirname ${BASH_SOURCE[0]}`"
sd_bin="$sd_main_path/../third-party/snp_effect_predictor.pl"

[ ".$input" == "." ]       && error "I need an input file."
[ ! -f "$input" ]          && error "Input file not found."
[ ".$en_mod_libs" == "." ] && error "path to ensembl libs not provided."
[ ".$species" == "." ]     && error "species not provided"

# Necessary Ensembl perl modules to load
#
e_libs=(
${en_mod_libs}/bioperl-1.2.3
${en_mod_libs}/ensembl/modules
${en_mod_libs}/ensembl-compara/modules
${en_mod_libs}/ensembl-variation/modules
${en_mod_libs}/ensembl-functgenomics/modules
)

# Add new libs to PERL5LIB envar
#
for m in ${e_libs[*]}
do
  echo "exporting $m"
  export PERL5LIB=$PERL5LIB:$m
done

echo "Running against ensembl ..."
$sd_bin -i $input -o $input.annotated.ensembl.txt -s $species
echo "Done."
