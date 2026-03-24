#!/usr/bin/bash -l
#SBATCH -p short -N 1 -n 1 -c 8 --mem 48gb --out logs/fastp.%a.log

module load fastp

module load workspace/scratch
INPUT=input
SAMPFILE=samples.csv
WORK=working
mkdir -p $WORK
CPU=2
if [ $SLURM_CPUS_ON_NODE ]; then
  CPU=$SLURM_CPUS_ON_NODE
fi
N=${SLURM_ARRAY_TASK_ID}
if [ -z $N ]; then
  N=$1
fi
if [ -z $N ]; then
  echo "cannot run without a number provided either cmdline or --array in sbatch"
  exit
fi
IFS=,
tail -n +2 $SAMPFILE | sed -n ${N}p | while read STRAIN SHOTGUN 
do
  mkdir -p $WORK/$STRAIN
  LEFT=$(ls $INPUT/$SHOTGUN | sed -n 1p)
  RIGHT=$(ls $INPUT/$SHOTGUN | sed -n 2p)
  echo "$LEFT and $RIGHT for $INPUT/$SHOTGUN"
  if [ ! -s $WORK/$STRAIN/${STRAIN}_R1.fq.gz ]; then
  	fastp -w $CPU --detect_adapter_for_pe -j logs/$STRAIN.LIB.json -h logs/$STRAIN.LIB.html \
	      -i $LEFT -I $RIGHT -o $WORK/$STRAIN/${STRAIN}_R1.fq.gz --out2 $WORK/$STRAIN/${STRAIN}_R2.fq.gz --trim_poly_g \
	      --unpaired1 $WORK/$STRAIN/${STRAIN}_unpair1.fq.gz --unpaired2 $WORK/$STRAIN/${STRAIN}_unpair2.fq.gz --overrepresentation_analysis
  fi
done
