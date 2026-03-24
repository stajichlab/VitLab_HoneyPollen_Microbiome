#!/usr/bin/bash -l
#SBATCH -N 1 -c 32 -n 1 --mem 128gb --out logs/gtkdb_classify.%a.log -J GTKDB

CPU=2
if [ $SLURM_CPUS_ON_NODE ]; then
  CPU=$SLURM_CPUS_ON_NODE
fi

GTKDB=/srv/projects/db/gtdbtk/220
MASHDB=gtkdb220.msh
module load workspace/scratch
module load gtdbtk

SAMPFILE=samples.csv
N=${SLURM_ARRAY_TASK_ID}
if [ -z $N ]; then
  N=$1
fi
if [ -z $N ]; then
  echo "cannot run without a number provided either cmdline or --array in sbatch"
  exit
fi

IFS=,
OUTFOLDER=results_bins_gtkdb
RUNDIR=classify_run
mkdir -p $OUTFOLDER $RUNDIR
tail -n +2 $SAMPFILE | sed -n ${N}p | while read STRAIN SHOTGUN
do
  PREFIX=$STRAIN
  BINFOLDER=$(realpath results/$STRAIN/bins)
  OUT=$(realpath $OUTFOLDER/$STRAIN)
  mkdir -p $OUT
  echo "$PREFIX and $BINFOLDER and $OUTFOLDER/$STRAIN"
  if [ ! -f $OUT/gtdbtk.bac120.summary.tsv ]; then
  	gtdbtk classify_wf --genome_dir $BINFOLDER --out_dir $OUT -x .fa --cpus $CPU --scratch_dir $SCRATCH --tmpdir $SCRATCH \
	  --mash_db $MASHDB 
  fi
done
