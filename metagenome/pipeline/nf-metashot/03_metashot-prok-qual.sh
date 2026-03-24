#!/usr/bin/bash -l
#SBATCH -N 1 -c 24 -n 1 --mem 80gb --out logs/metashot_qual.%a.log

CPU=2
if [ $SLURM_CPUS_ON_NODE ]; then
  CPU=$SLURM_CPUS_ON_NODE
fi

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
OUTFOLDER=results_bins_qual
RUNDIR=qual_run
mkdir -p $OUTFOLDER $RUNDIR
tail -n +2 $SAMPFILE | sed -n ${N}p | while read STRAIN SHOTGUN
do
  PREFIX=$STRAIN
  BINFOLDER=$(realpath results/$STRAIN/bins)
  OUT=$(realpath $OUTFOLDER/$STRAIN)
  mkdir -p $OUT
  pushd $RUNDIR
  if [ ! -f process.config ]; then ln -s ../process-qual.config process.config; ln -s ../process-qual.config ./; fi
  if [ ! -f metashot-qual.cfg ]; then ln -s ../metashot-qual.cfg ./; fi
  if [ ! -f nextflow ]; then ln -s ../nextflow ./; fi
  echo "$PREFIX and $BINFOLDER and $OUTFOLDER/$STRAIN"
  if [ ! -f  $OUT/derep_info.tsv ]; then
  	./nextflow run metashot/prok-quality -c metashot-qual.cfg \
	     --genomes "$BINFOLDER/*.fa" \
	     --outdir $OUT --max_cpus $CPU \
	     --scratch $SCRATCH --resume
  fi
  popd
done
