#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH --partition=normal,bigmem,long
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
#SBATCH --output=log%J.out
#SBATCH --error=log%J.err

#load fastqc module
module load fastqc/0.12.1

#setup directories
RESULTS_DIR=~/fastqcresults
RAW_DATA_DIR=~/01.RawData

#create results directory
mkdir -p "$RESULTS_DIR"

#navigate to raw data directory with fastq files
cd "$RAW_DATA_DIR"

#run fastqc
fastqc -o "$RESULTS_DIR" -t 6 *.fq.gz
