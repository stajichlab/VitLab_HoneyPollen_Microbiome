#!/usr/bin/bash -l
#SBATCH -N 1 -c 96 -n 1 --mem 255gb --out logs/mmseqs_classify_scaffolds_missing.%a.log -J mmseqs_missing

# --- CPU setup ---
CPU=2
if [ $SLURM_CPUS_ON_NODE ]; then
  CPU=$SLURM_CPUS_ON_NODE
fi

# --- Modules ---
module load mmseqs2
module load workspace/scratch

# --- Input files ---
SAMPFILE=samples_missing.csv
N=${SLURM_ARRAY_TASK_ID}
if [ -z $N ]; then
  N=$1
fi
if [ -z $N ]; then
  echo "❌ ERROR: Provide an array index or sample number (cmdline or --array in sbatch)"
  exit 1
fi

# --- Database ---
DB=/srv/projects/db/ncbi/mmseqs/uniref50
DBNAME=$(basename $DB)

# --- Output folder ---
OUTFOLDER=results_scaffold_classify_mmseqs
mkdir -p $OUTFOLDER

# --- Process each sample from samples_missing.csv ---
IFS=,
tail -n +2 $SAMPFILE | sed -n ${N}p | while read STRAIN SHOTGUN
do
  echo "🔍 Running MMseqs2 taxonomy for sample: $STRAIN"
  SCAFFOLDS=results/$STRAIN/scaffolds

  if [ ! -d "$SCAFFOLDS" ]; then
    echo "⚠️  Missing scaffolds directory: $SCAFFOLDS"
    continue
  fi

  # Find .fa or .fa.gz file
  GENOME=$(ls $SCAFFOLDS/*.fa 2>/dev/null || true)
  if [ -z "$GENOME" ]; then
    GENOME=$(ls $SCAFFOLDS/*.fa.gz 2>/dev/null || true)
  fi

  # If still empty, skip
  if [ -z "$GENOME" ]; then
    echo "⚠️  No FASTA (.fa or .fa.gz) file found in $SCAFFOLDS"
    continue
  fi

  # If gzipped, decompress to scratch
  if [[ "$GENOME" == *.gz ]]; then
    TMP_FASTA="$SCRATCH/${STRAIN}_tmp.fa"
    echo "🗜️  Decompressing $GENOME to $TMP_FASTA"
    gunzip -c "$GENOME" > "$TMP_FASTA"
    GENOME="$TMP_FASTA"
  fi

  mkdir -p $OUTFOLDER/$STRAIN
  OUT=$OUTFOLDER/$STRAIN/${STRAIN}_${DBNAME}

  echo "🧬 GENOME: $GENOME"
  echo "📂 OUTPUT: $OUT"

  mmseqs touchdb $DB

  if [ ! -s ${OUT}_tophit_aln ]; then
    mmseqs easy-taxonomy "$GENOME" "$DB" "$OUT" "$SCRATCH" \
      --threads $CPU \
      --lca-ranks kingdom,phylum,family \
      --tax-lineage 1 \
      --db-load-mode 2
    echo "✅ Finished: $STRAIN"
  else
    echo "⏩ Already exists: ${OUT}_tophit_aln"
  fi
done

