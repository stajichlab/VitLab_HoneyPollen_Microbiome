#!/usr/bin/bash -l
#SBATCH -p short -N 1 -c 48 -n 1 --mem 64gb --out logs/checkm_bins.%a.log -J checkM

CPU=2
if [ $SLURM_CPUS_ON_NODE ]; then
  CPU=$SLURM_CPUS_ON_NODE
fi

module load checkm
module load workspace/scratch
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
OUTFOLDER=results_bins_checkm
mkdir -p $OUTFOLDER
tail -n +2 $SAMPFILE | sed -n ${N}p | while read STRAIN SHOTGUN
do
  PREFIX=$STRAIN
  BINFOLDER=results/$STRAIN/bins
  mkdir -p $OUTFOLDER/$STRAIN
  if [ ! -f $OUTFOLDER/$STRAIN/lineage.ms ]; then
  	checkm lineage_wf -t $CPU -x fa $BINFOLDER $OUTFOLDER/$STRAIN
  fi
  if [ ! -f $OUTFOLDER/$STRAIN/summary_table.tsv ]; then
	  checkm qa -t $CPU $OUTFOLDER/$STRAIN/lineage.ms $OUTFOLDER/$STRAIN --tab_table -f $OUTFOLDER/$STRAIN/summary_table.tsv -o 2
  fi
done
