#!/bin/bash
#SBATCH --job-name=star
#SBATCH --partition=normal,bigmem,long
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=48G
#SBATCH --output=log%j.out
#SBATCH --error=log%j.err

set -euo pipefail

#load star module
module load STAR/2.7.10a

#directories
GENOME_FA=~/04.Ref/genome.fa
GENOME_GTF=~/04.Ref/genome.gtf
RAW_DATA_DIR=~/01.RawData
STAR_INDEX_DIR=~/04.Ref/star_index
STAR_RESULTS_DIR=~/results/STAR

#make directories for star index and alignment results
mkdir -p "$STAR_INDEX_DIR" "$STAR_RESULTS_DIR"

#build STAR index
STAR --runThreadN 6 \
    --runMode genomeGenerate \
    --genomeDir "$STAR_INDEX_DIR" \
    --genomeFastaFiles "$GENOME_FA" \
    --sjdbGTFfile "$GENOME_GTF" \
    --sjdbOverhang 149

#run STAR alignment on 14 samples
for i in {1..14}; do
    SAMPLE="A${i}"
    READ1="${RAW_DATA_DIR}/${SAMPLE}_1.fq.gz"
    READ2="${RAW_DATA_DIR}/${SAMPLE}_2.fq.gz"

    echo "Processing sample $SAMPLE"

    STAR --genomeDir "$STAR_INDEX_DIR" \
         --runThreadN 6 \
         --readFilesIn "$READ1" "$READ2" \
         --readFilesCommand zcat \
         --outFileNamePrefix "${STAR_RESULTS_DIR}/${SAMPLE}_" \
         --outSAMtype BAM SortedByCoordinate
done
