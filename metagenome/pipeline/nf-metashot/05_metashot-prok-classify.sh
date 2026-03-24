#!/usr/bin/bash -l
#SBATCH -N 1 -c 48 -n 1 --mem 128gb --out logs/metashot_classify.%a.log -J GTKDB

CPU=2
if [ $SLURM_CPUS_ON_NODE ]; then
  CPU=$SLURM_CPUS_ON_NODE
fi

GTKDB=/srv/projects/db/gtdbtk/207_v2/
module load singularity
module load workspace/scratch
export NXF_SINGULARITY_CACHEDIR=/bigdata/stajichlab/shared/singularity_cache/
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
  pushd $RUNDIR
  if [ ! -f process.config ]; then ln -s ../process-classify.config process.config; ln -s ../process-classify.config ./; fi
  if [ ! -f metashot-qual.cfg ]; then ln -s ../metashot-classify.cfg ./; fi
  if [ ! -f nextflow ]; then ln -s ../nextflow ./; fi
  ./nextflow run metashot/prok-classify -c metashot-classify.cfg \
	     --genomes "$BINFOLDER/*.fa" \
  	     --gtdbtk_db /srv/projects/db/gtdbtk/207_v2 \
	     --outdir $OUT --max_cpus $CPU \
	     --scratch $SCRATCH
done
