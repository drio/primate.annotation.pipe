#!/bin/sh
#
set -e

c_dir=`pwd`
prj_name=`basename $c_dir`
db="$HOME/Dropbox/git_repo"
db_dir=$db/$prj_name
cluster_dir="$HOME/sshfs/ardmore/stornext/snfs5/next-gen/drio-scratch/synthetic.pipe/ecoli"

log()
{
  echo ">> $1"
}

if [ -d $db_dir ]
then
  log "PRJ found, pulling"
  cd $db_dir    
  git pull 
  git push github
  cd -
else
  log "PRJ not found, cloning"
  cd $db 
  git clone $c_dir 
  cd -
fi

exit # !

if [ ! -d $cluster_dir ]
then
  log "Cound't find cluster dir. ($cluster_dir)"
  exit 1
else
  log "Copying files to cluster. ($cluster_dir)"
  date > $cluster_dir/.last_transfer.txt 
  cp README bfast2novo.pl plots.R run_experiment.sh spread.sh $cluster_dir
fi
