---
bibliography: references.bib
---

This repository contains the code used for the analysis of scRNAseq, kidney marrow proteomics, and plasma proteomics datasets generated as part of Morabito et al. killifish immune ageing multiomics project.

# scRNASeq

## Data preprocessing 

**./scripts/scRNAseq/01-dataPrep/**: N. furzeri primary genome assembly and genome annotations were downloaded from the Ensembl website (version 105) (@cunningham2021). Cell Ranger pipeline (@zheng2017) is used to create a custom reference genome for *N.furzeri* (using `cellranger mkref`) and a genome annotation file (using `cellranger mkgtf` command and only protein-coding, IG or TR genes). Count matrices were generated for each sample in each run separately using the `cellranger count` function (`./scripts/scRNAseq/01-dataPrep/cellrangercount2.sh`).

**./scripts/scRNAseq/02-seurat.R**: We checked the distribution of number of cells, number of features, UMI counts per cells, total read counts per cell across samples and runs. One young sample (sample B) was sequenced in only one sequencing run and had low quality (`./results/scRNAseq/qc_prenorm/allqc.png`). Thus, it is excluded from the downstream analysis. All other samples showed a comparable number of cells, features, UMI, and total read counts. All downstream steps were performed using Seurat package (@Seurat) in R. We filtered out the cells with ​​less than 200 or more than 2500 features. All samples in each run were log normalised independently using the `NormalizeData` function in the Seurat package. Next, 2000 variable features were selected using the \'vst\' method implemented in the `FindVariableFeatures` function in the Seurat package. Next, all data were integrated using canonical correlation analysis to find anchors and the first 20 dimensions for the anchor weighting process. Data were scaled, and the first 10 PCs were used to generate the tSNE, UMAP, and clustering of the cells (`./results/scRNAseq/seurat`).

## Cell-type clustering and annotations

Clustering was performed using shared nearest neighbour graph construction, using the first 10 PCs. Multiple resolutions were assessed to determine the number of clusters. Resolution of 0.3, which provided the last most stable clustering based on the clustering tree was chosen, resulting in 17 clusters (`./results/scRNAseq/seurat/cluster_resolutions.png`). Next, deferentially expressed markers for each cluster were identified. All genes that are expressed in at least 25% of the cells in a cluster were tested. Based on the resulting marker list (`./results/scRNAseq/seurat/clustermarkers.csv`), cell types are annotated (`./results/scRNAseq/seurat/clusterLabels.txt`). 

We found four clusters suggesting contamination from gonads (clusters 11, 12, 13, 14). We excluded those cells from the analysis and repeated the previous steps for better clustering and marker annotation (`./scripts/scRNAseq/03-seurat_afterQC.R`). We provide data before and after the exclusion of these cell types. Previously described steps were repeated, except we used the first 20 PCs for clustering and UMAP and tSNE plots. Cluster resolution tree graph, list of cluster markers, tSNE and UMAP plots are generated under `./results/scRNAseq/seurat_afterQC/`). New annotations for 16 clusters are given as `./results/scRNAseq/seurat_afterQC/clusterLabels.txt`.

Annotated tSNE and UMAP plots are created using `./scripts/scRNAseq/04-tSNE-UMAP.R` - tSNE plot is used in the manuscript as figure 2a.

# Helper functions/data

In order to use the same database and package versions for geneID conversions and orthology mapping across scripts, we obtained gene IDs and orthologs from Ensembl database using a small wrapper R package called [dataCollectR](https://github.com/mdonertas/dataCollectR) (in ./scripts/scRNAseq/02-seurat.R). Mapping between multiple type of gene IDs (gene external name, Ensembl Gene ID, Entrez Gene ID, Gene Biotype, and gene descriptions) and orthology mapping between human and *N. furzeri* genes were obtained on 2022.02.17 using biomaRt R Package (v2.50.2) (@durinck2009). These data are available under (`./data/processed/helperdata/`).
