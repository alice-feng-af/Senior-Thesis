CMTR2 RNA-seq Analysis Pipeline

This repository contains the complete analysis pipeline for RNA-seq data comparing Kras/Trp53/CMTR2 mutants vs Kras/Trp53 conditions in Mus musculus (GRCm38/mm10).

Overview:
14 samples (A1-A14)
A1-A7: Kras/Trp53 (KP)
A8-A14: Kras/Trp53/CMTR2 (KPC)
Mouse reference genome: GRCm38/mm10

Pipeline Steps
1. FASTQC quality control
2. STAR alignment and Qualimap
3. Feature counting
4. Differential expression analysis with DESeq2
5. Visualization
