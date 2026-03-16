#!/bin/bash
#SBATCH --job-name=star
#SBATCH --partition=normal,bigmem,long
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=48G
#SBATCH --output=log%J.out
#SBATCH --error=log%J.err

#load modules
module load STAR/2.7.10a
module load Qualimap/2.2.1

#define paths
GENOME_FA_GZ=~/04.Ref/genome.fa.gz
GENOME_GTF_GZ=~/04.Ref/genome.gtf.gz
GENOME_FA=~/04.Ref/genome.fa
GENOME_GTF=~/04.Ref/genome.gtf
RAW_DATA_DIR=~/01.RawData

STAR_INDEX_DIR=~/04.Ref/star_index
STAR_RESULTS_DIR=~/results/STAR
QUALIMAP_BASE_DIR=~/results/qualimap

#create directories
mkdir -p "$STAR_INDEX_DIR"
mkdir -p "$STAR_RESULTS_DIR"
mkdir -p "$QUALIMAP_BASE_DIR"

#unzip reference genome files
gunzip -c "$GENOME_FA_GZ" > "$GENOME_FA"
gunzip -c "$GENOME_GTF_GZ" > "$GENOME_GTF"

#build star genome index
STAR --runThreadN 6 \
  --runMode genomeGenerate \
  --genomeDir "$STAR_INDEX_DIR" \
  --genomeFastaFiles "$GENOME_FA" \
  --sjdbGTFfile "$GENOME_GTF" \
  --sjdbOverhang 99

#star alignment and qualimap for 14 samples
for i in {1..14}
do
    SAMPLE="A${i}"
    READ1=${RAW_DATA_DIR}/${SAMPLE}_1.fq.gz
    READ2=${RAW_DATA_DIR}/${SAMPLE}_2.fq.gz

    QUALIMAP_RESULTS_DIR=${QUALIMAP_BASE_DIR}/${SAMPLE}

    mkdir -p "$QUALIMAP_RESULTS_DIR"

    echo "Processing sample $SAMPLE"

    #STAR alignment
    STAR --genomeDir "$STAR_INDEX_DIR" \
      --runThreadN 6 \
      --readFilesIn "$READ1" "$READ2" \
      --readFilesCommand zcat \
      --outFileNamePrefix "$STAR_RESULTS_DIR/${SAMPLE}_" \
      --outSAMtype BAM SortedByCoordinate

    #Qualimap
    qualimap rnaseq \
      -bam "$STAR_RESULTS_DIR/${SAMPLE}_Aligned.sortedByCoord.out.bam" \
      -gtf "$GENOME_GTF" \
      -outdir "$QUALIMAP_RESULTS_DIR" \
      -a proportional \
      -p non-strand-specific \
      --java-mem-size=16G
done
