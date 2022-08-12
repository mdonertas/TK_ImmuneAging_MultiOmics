---
bibliography: references.bib
---

This repository contains the code used for the analysis of scRNAseq, kidney marrow proteomics, and plasma proteomics datasets generated as part of Morabito et al. killifish immune ageing multiomics project.

# scRNASeq

## Data preprocessing (./scripts/scRNAseq/01-dataPrep/)

N. furzeri primary genome assembly and genome annotations were downloaded from the Ensembl website (version 105) (@cunningham2021). Cell Ranger pipeline (@zheng2017) is used to create a custom reference genome for *N.furzeri* (using `cellranger mkref`) and a genome annotation file (using `cellranger mkgtf` command and only protein-coding, IG or TR genes). Count matrices were generated for each sample in each run separately using the `cellranger count` function (`./scripts/scRNAseq/01-dataPrep/cellrangercount2.sh`).

## 
