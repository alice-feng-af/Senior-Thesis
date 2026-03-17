CMTR2 RNA-seq Analysis Pipeline

This repository contains the complete analysis pipeline for RNA-seq data comparing Kras/Trp53/CMTR2 vs Kras/Trp53 conditions in Mus musculus (GRCm38/mm10).

Overview: <br />
14 samples (A1-A14) <br />
A1-A7: Kras/Trp53 (KP) <br />
A8-A14: Kras/Trp53/CMTR2 (KPC) <br />
Mouse reference genome: GRCm38/mm10 <br />

Pipeline Steps
1. FASTQC Quality Control (1. fastqc.sh)
2. STAR Alignment (2. staralignment.sh)
3. Qualimap (3. qualimap.sh)
4. Feature Counting (4. featurecount.sh)
5. Differential Expression Analysis (5. Differential Expression Analysis.R)
