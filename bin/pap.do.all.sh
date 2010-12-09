#!/bin/bash
#
set -e

log()
{
  echo "`date` >> $1"
}
error()
{
  cat<<-EOF
Usage:
`basename $0` <snps.summary.csv> <indels.summary.csv>

If you want to use a pbs cluster, here is an example:
$ echo "pap.do.all.sh  snps.summary.txt indels.summary.txt" |  msub -N "pap" -q analysis -d `pwd` -e e -o o -l nodes=1:ppn=1,mem=100mb -V
EOF
  
  exit 1
}

ssf="$1" # snps summary
insf="$2" # indel summary

[ ! -f "$ssf" ]  &&  error
[ ! -f "$insf" ] &&  error

export PATH=$PATH:/stornext/snfs5/next-gen/drio-scratch/primate.annotation.pipe/bin
export PATH=/stornext/snfs5/next-gen/drio-scratch/bb/local/bin:$PATH

log "snps -> ensembl"
primate_bridge.rb -a to_en_snps -s $ssf > s.i.ensembl 2> s.pb.err
primate_bridge.rb -a to_en_snps -s $insf -i > i.i.ensembl 2> i.pb.err

log "annotate"
run_snp_effect_predictor.sh s.i.ensembl /hgsc_software/perl/perl-5.12.2/ensembl/ensembl macaque /hgsc_software/perl/perl-5.12.2/bin/perl 2> s.run.ens.err
run_snp_effect_predictor.sh i.i.ensembl /hgsc_software/perl/perl-5.12.2/ensembl/ensembl macaque /hgsc_software/perl/perl-5.12.2/bin/perl 2> i.run.ens.err

log "ensembl -> lff"
cat s.i.ensembl.annotated.ensembl.txt | ensembl2lff.rb > s.lff 2> s.lff.err
cat i.i.ensembl.annotated.ensembl.txt | ensembl2lff.rb > i.lff 2> i.lff.err
