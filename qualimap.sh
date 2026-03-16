#!/bin/bash
#SBATCH --job-name=qualimap
#SBATCH --partition=normal,bigmem,long
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=48G
#SBATCH --output=log%J.out
#SBATCH --error=log%J.err

echo "loading module"

#load modules
module load Qualimap/2.2.1

#use java headless
unset DISPLAY
export JAVA_TOOL_OPTIONS="-Djava.awt.headless=true"
export _JAVA_OPTIONS="-Djava.awt.headless=true"

echo "defining paths"
#define paths
GENOME_FA_GZ=~/04.Ref/genome.fa.gz
GENOME_GTF_GZ=~/04.Ref/genome.gtf.gz
GENOME_FA=~/04.Ref/genome.fa
GENOME_GTF=~/04.Ref/genome.gtf
RAW_DATA_DIR=~/01.RawData

STAR_INDEX_DIR=~/04.Ref/star_index
STAR_RESULTS_DIR=~/results/STAR
QUALIMAP_BASE_DIR=~/results/qualimap

#qualimap for 14 samples
for i in {1..14}
do
    SAMPLE="A${i}"
    QUALIMAP_RESULTS_DIR=${QUALIMAP_BASE_DIR}/${SAMPLE}

    echo "Processing sample $SAMPLE"

    qualimap rnaseq \
      -bam "$STAR_RESULTS_DIR/${SAMPLE}_Aligned.sortedByCoord.out.bam" \
      -gtf "$GENOME_GTF" \
      -outdir "$QUALIMAP_RESULTS_DIR" \
      -a proportional \
      -p non-strand-specific \
      -pe \
      --java-mem-size=16G
done
