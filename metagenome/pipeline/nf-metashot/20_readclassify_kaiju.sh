#!/usr/bin/bash -l
#SBATCH -p short -N 1 -n 1 -c 48 --mem 256gb --out logs/kaiju.%a.log

module load kaiju

module load workspace/scratch
INPUT=input
SAMPFILE=samples.csv
WORK=working
OUT=results_kaiju
DBFOLDER=/srv/projects/db/kaiju/20240825
DBNAME=kaiju_db_nr.fmi
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
	  module load fastp
	WORK=$SCRATCH
	mkdir $WORK/$STRAIN
  	fastp -w $CPU --detect_adapter_for_pe -j logs/$STRAIN.LIB.json -h logs/$STRAIN.LIB.html \
	      -i $LEFT -I $RIGHT -o $WORK/$STRAIN/${STRAIN}_R1.fq.gz --out2 $WORK/$STRAIN/${STRAIN}_R2.fq.gz --trim_poly_g \
	      --unpaired1 $WORK/$STRAIN/${STRAIN}_unpair1.fq.gz --unpaired2 $WORK/$STRAIN/${STRAIN}_unpair2.fq.gz --overrepresentation_analysis
	module unload fastp
  fi
    if [ ! -s $OUT/${STRAIN}_kaiju.out ]; then
        kaiju -t ${DBFOLDER}/nodes.dmp -f ${DBFOLDER}/${DBNAME} -i $WORK/$STRAIN/${STRAIN}_R1.fq.gz -j $WORK/$STRAIN/${STRAIN}_R2.fq.gz \
        -o $OUT/${STRAIN}_kaiju.out -z $CPU
    fi
    if [ ! -f $OUT/${STRAIN}_kaiju_genus.tsv ]; then
        kaiju2table -t ${DBFOLDER}/nodes.dmp -n ${DBFOLDER}/names.dmp -r genus -o $OUT/${STRAIN}_kaiju_genus.tsv $OUT/${STRAIN}_kaiju.out
        kaiju2table -t ${DBFOLDER}/nodes.dmp -n ${DBFOLDER}/names.dmp -r phylum -o $OUT/${STRAIN}_kaiju_phylum.tsv $OUT/${STRAIN}_kaiju.out
        kaiju2table -t ${DBFOLDER}/nodes.dmp -n ${DBFOLDER}/names.dmp -r family -o $OUT/${STRAIN}_kaiju_family.tsv $OUT/${STRAIN}_kaiju.out
    fi
    if [ ! -s $OUT/${STRAIN}_kaiju.krona.html ]; then
        module load KronaTools
        kaiju2krona -t ${DBFOLDER}/nodes.dmp -n ${DBFOLDER}/names.dmp -o $OUT/${STRAIN}_kaiju.krona -l superkingdom,phylum,class,order,family,genus,species -u -i $OUT/${STRAIN}_kaiju.out
        ktImportText -o $OUT/${STRAIN}_kaiju.krona.html $OUT/${STRAIN}_kaiju.krona
    fi
done
