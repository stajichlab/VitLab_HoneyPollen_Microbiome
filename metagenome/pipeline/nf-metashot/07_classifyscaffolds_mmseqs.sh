#!/usr/bin/bash -l
#SBATCH -N 1 -c 96 -n 1 --mem 255gb --out logs/mmseqs_classify_scaffolds.%a.log -J 04_classify

CPU=2
if [ $SLURM_CPUS_ON_NODE ]; then
  CPU=$SLURM_CPUS_ON_NODE
fi

module load mmseqs2
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
DB=/srv/projects/db/ncbi/mmseqs/uniref50
#DB=/srv/projects/db/ncbi/mmseqs/swissprot
DBNAME=$(basename $DB)

IFS=,
OUTFOLDER=results_scaffold_classify_mmseqs
mkdir -p $OUTFOLDER
tail -n +2 $SAMPFILE | sed -n ${N}p | while read STRAIN SHOTGUN
do
  PREFIX=$STRAIN
  SCAFFOLDS=results/$STRAIN/scaffolds
  mkdir -p $OUTFOLDER/$STRAIN
  mmseqs touchdb $DB
  GENOME=$(ls $SCAFFOLDS/*.fa)
  echo "GENOME is $GENOME"
  OUT=$OUTFOLDER/$STRAIN/${STRAIN}_${DBNAME}
  if [ ! -s ${OUT}_tophit_aln ]; then
	mmseqs easy-taxonomy $GENOME $DB $OUT $SCRATCH --threads $CPU --lca-ranks kingdom,phylum,family  --tax-lineage 1 --db-load-mode 2
  fi
done
