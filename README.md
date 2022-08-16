---
bibliography: references.bib
---

This repository contains the code used for the analysis of scRNAseq, kidney marrow proteomics, and plasma proteomics datasets generated as part of Morabito et al. killifish immune ageing multiomics project.

# Proteomics 

## Data preprocessing

Proteomics data was analyzed using MaxQuant, version 1.6.10.43 (@Cox2008 ). Peptide fragmentation spectra were searched against the canonical sequences of the *Nothobranchius furzeri* proteome (downloaded September 2019 from UniProt). Protein names and primary gene names, corresponding to the Uniprot IDs were downloaded from Uniprot and used for annotation of the data. Methionine oxidation and protein N-terminal acetylation were set as variable modifications; cysteine carbamidomethylation was set as fixed modification. The digestion parameters were set to "specific" and "Trypsin/P". Quantification was set to "Reporter ion MS3". The isotope purity correction factors, provided by the manufacturer, were imported and included in the analysis. The minimum number of peptides and razor peptides for protein identification was 1; the minimum number of unique peptides was 0. Protein identification was performed at a peptide spectrum matches and protein false discovery rate of 0.01. The "second peptide" option was on. TMT reporter intensities were normalized using vsn (@Huber2002 ) and log2 transformed in R, version 3.4.3 (@base ).

## Exploratory data analysis and PCA

PCA done using `prcomp` function in base R, on VSN normalized, log2 transformed data with z-transformation (proteins are scaled).

## Plasma proteomics data description

We detected 474 proteins across 10 samples in plasma proteomics. One old sample (old sample 3) was an outlier (`./results/plasmaProteomics/allsamples/pca.pdf`). PC2 (27%) is mainly driven by this outlier. This PC shows enrichment for proteins involved in response to toxic substance, DNA confirmation change, regulation of catalytic activity and many other categories (we calculated enrichment based on GSEA procedure explained below). Since we do not have any technical explanations why this sample may be an outlier, it can represent biological variability thus to be more conservative we did not exclude this sample. However, we found comparable results in the downstream after exclusion of this sample (Spearman correlation coef for differential abundance=0.97, for GSEA enrichment scores=1 `./results/msfigures/proteomics_outlier_check`).

## Kidney Marrow proteomics data description

We detected 6770 proteins across 10 samples in kidney marrow proteomics. One young sample (young sample 2) was an outlier (`./results/KKMproteomics/allsamples/pca.pdf`). PC2 (20%) shows association with this outlier. This PC shows enrichment for metabolic (mostly catabolic) processes. Like in plasma proteomics, we do not have any technical explanations why this sample may be an outlier and since it can represent biological variability we wanted to be more conservative we did not exclude this sample. However, we found comparable results in the downstream after exclusion of this sample (Spearman correlation coef for differential abundance=0.99, for GSEA enrichment scores=0.97 `./results/msfigures/proteomics_outlier_check`).

## Differential expression

We calculated differences between log2 median protein abundance levels of young and old samples to determine differential expression. Statistical significance was assigned based on Wilcoxon rank sum test using `wilcox.test` function in R. DE results are given as:

-   Plasma proteomics, all samples: `./results/plasmaProteomics/allsamples/proteinDEtest.csv`

-   Plasma proteomics, after outlier removal: `./results/plasmaProteomics/outler_removed/proteinDEtest.csv`

-   Kidney marrow, all samples: `./results/KMproteomics/allsamples/proteinDEtest.csv`

-   Kidney marrow, after outlier removal: `./results/KMproteomics/outler_removed/proteinDEtest.csv`

## Gene ontology - Gene set enrichment analysis (GSEA)

Gene set enrichment analysis was performed using `gseGO` function in `clusterProfiler` package, using human gene symbols and `org.Hs.eg.db` annotation package in R. We analysed the GO Biological Processes categories with a minimum of 10 and maximum of 500 annotated genes. We used BY correction for multiple testing and considered BY-corrected p-value\<0.05 as significant. Details of Human - *N. furzeri* orthology mapping are explained under `Helper functions/data -> Gene ID & Orthology mapping`. Two outputs from `gseGO` function are used: core enrichment genes (i.e. genes that contribute most to the enrichment result.) and NES (i.e. normalized enrichment score).

### Choice of GO representatives for visualisation and summarisation

Since go enrichment results gave many significant GO categories but most of them showed overlaps, we used an in-house method to choose representative GO categories for visualisation and summarisation purposes. GO enrichment results include the whole list of GO categories, together with their representatives.

In order to choose the representatives we calculated the jaccard similarity between GO categories based on core enrichment genes. We then performed hierarchical clustering of the similarity matrix, cutting the tree at different levels (20 to 70 clusters). We calculated the median jaccard index within each cluster. We take the minimum number of clusters, k, where at least half of the clusters have median jaccard index of 0.5 or higher. This resulted in 20 clusters (=representatives) for plasma proteomics and 50 clusters (=representatives) for kidney marrow proteomics.

### Plasma proteomics - GSEA details 

438 of 474 proteins had at least one ortholog in humans. In total they map to 608 human proteins. In order to make sure 1 to many orthology does not bias the results, we repeated the analysis using only 1 random ortholog resulting in 438 proteins. We detected a spearman correlation of 0.95 (p-value \< 2.2e-16). Enrichment results are given as:

-   Plasma proteomics, all samples: `./results/plasmaProteomics/allsamples/GOenrichment.csv`

-   Plasma proteomics, after outlier removal: `./results/plasmaProteomics/outler_removed/GOenrichment.csv`

### Kidney marrow proteomics - GSEA details

5395 of 6770 proteins had at least one ortholog in humans. In total they map to 5208 human proteins. In order to make sure 1 to many orthology does not bias the results, we repeated the analysis using only 1 random ortholog resulting in 5207 proteins. We detected a spearman correlation of 0.96 (p-value \< 2.2e-16). Enrichment results are given as:

-   Kidney marrow, all samples: `./results/KMproteomics/allsamples/GOenrichment.csv`

-   Kidney marrow, after outlier removal: `./results/KMproteomics/outler_removed/GOenrichment.csv`

#### Gene overlaps across GO categories represented by 'DNA repair', 'double-strand break repair via homologous recombination', and 'response to UV'

We collected all genes associated with the significant GO categories in Kidney marrow proteomics, represented by 'DNA repair', 'double-strand break repair via homologous recombination', and 'response to UV' categories. Using only the gene set that has an absolute log2 fold change larger than 1, we calculated the jaccard index across GO categories. We plotted the jaccard similarity indices as a correlation plot, where the color and size of squares show the similarity index and the GO categories are colored by the hierarchical clustering of similarity data. We performed PCA on the same data (without scaling as jaccard index can only range between 0 and 1) to see clustering of GO categories in two dimensions.

# scRNASeq

## Data preprocessing

**./scripts/scRNAseq/01-dataPrep/**: *N. furzeri* primary genome assembly and genome annotations were downloaded from the Ensembl website (version 105) (@cunningham2021). Cell Ranger pipeline (@zheng2017) is used to create a custom reference genome for *N.furzeri* (using `cellranger mkref`) and a genome annotation file (using `cellranger mkgtf` command and only protein-coding, IG or TR genes). Count matrices were generated for each sample in each run separately using the `cellranger count` function (`./scripts/scRNAseq/01-dataPrep/cellrangercount2.sh`).

**./scripts/scRNAseq/02-seurat.R**: We checked the distribution of number of cells, number of features, UMI counts per cells, total read counts per cell across samples and runs. One young sample (sample B) was sequenced in only one sequencing run and had low quality (`./results/scRNAseq/qc_prenorm/allqc.png`). Thus, it is excluded from the downstream analysis. All other samples showed a comparable number of cells, features, UMI, and total read counts. All downstream steps were performed using Seurat package (@Seurat) in R. We filtered out the cells with ​​less than 200 or more than 2500 features. All samples in each run were log normalized independently using the `NormalizeData` function in the Seurat package. Next, 2000 variable features were selected using the 'vst' method implemented in the `FindVariableFeatures` function in the Seurat package. Next, all data were integrated using canonical correlation analysis to find anchors and the first 20 dimensions for the anchor weighting process. Data were scaled, and the first 10 PCs were used to generate the tSNE, UMAP, and clustering of the cells (`./results/scRNAseq/seurat`).

## Cell-type clustering and annotations

Clustering was performed using shared nearest neighbor graph construction, using the first 10 PCs. Multiple resolutions were assessed to determine the number of clusters. Resolution of 0.3, which provided the last most stable clustering based on the clustering tree was chosen, resulting in 17 clusters (`./results/scRNAseq/seurat/cluster_resolutions.png`). Next, deferentially expressed markers for each cluster were identified. All genes that are expressed in at least 25% of the cells in a cluster were tested. Based on the resulting marker list (`./results/scRNAseq/seurat/clustermarkers.csv`), cell types are annotated (`./results/scRNAseq/seurat/clusterLabels.txt`). 

We found four clusters suggesting contamination from gonads (clusters 11, 12, 13, 14). We excluded those cells from the analysis and repeated the previous steps for better clustering and marker annotation (`./scripts/scRNAseq/03-seurat_afterQC.R`). We provide data before and after the exclusion of these cell types. Previously described steps were repeated, except we used the first 20 PCs for clustering and UMAP and tSNE plots. Cluster resolution tree graph, list of cluster markers, tSNE and UMAP plots are generated under `./results/scRNAseq/seurat_afterQC/`). New annotations for 16 clusters are given as `./results/scRNAseq/seurat_afterQC/clusterLabels.txt`.

Annotated tSNE and UMAP plots are created using `./scripts/scRNAseq/04-tSNE-UMAP.R` - tSNE plot is used in the manuscript as figure 2a.

### Age-related changes in cell-type proportions

Since we do not have enough biological replicates in each age group, we designed a test based on re-sampling of cells of each individual and pairwise comparisons between young and old. More specifically, we rarefy each sample to the same number of cells (2,500) 1,000 times, and for each cell type, we calculate the Log2 Fold change between old and young samples using the median number of cells. Since the number of samples is small, we also calculate the pairwise Log2 difference between the number of cells between young and old sample pairs to calculate a confidence interval using the 1st and 3rd quartiles. Code for the text and EDA can be found as `./scripts/scRNAseq/05-cellType_proportions.R` and the results are under `./results/scRNAseq/celltypes/`.

## scRNAseq GO analysis

Four GO representative categories, Cellular detoxification, DNA repair, DNA replication, and Immune cell activation, were chosen for further investigation in scRNAseq data based on the results of proteomics dataset from the kidney marrow. These GO categories represent 35, 14, 55, and 15 GO categories respectively (the full list is given as `./results/scRNAseq/GO_in_scRNAseq/rep_goIDs.csv`. All human genes associated with these categories were obtained as explained under `Helper functions/data -> Gene Ontology data`. We then get the orthologs in *N. furzeri*, intersect with the scRNAseq dataset, get the genes with an absolute value of log2 fold change of 1 in kidney marrow proteomics data, and only use the non-overlapping genes. The number of genes after each step is given as `./results/scRNAseq/GO_in_scRNAseq/genes_in_categories.csv`. For each cell in the dataset, we calculated the percentage of the expressed genes in each representative category and divided this number to the percent of expressed genes among all genes detected in the scRNAseq. In this way, we aimed to normalize for the overall transcriptional profile of cells and obtained an odd's ratio showing the enrichment of genes specific to the GO categories of interest. We coloured the cells on tSNE with this odds ratio (`./results/scRNAseq/GO_in_scRNAseq/percExp_tsne.pdf`). We also calculated the correlation between the enrichment scores of i) replication and repair and ii) replication and immune cell activation. In order to correct for potential bias due to mean to the regression, we obtained p values for correlations based on a permutation test where we randomly selected N number of genes among non-overlapping scRNAseq genes, where N is the number of genes used in correlation calculation for DNA repair and Immune cell activation. We calculated correlations between DNA replication and these permutations of genes and used the number of cases where we obtained a correlation as extreme as the observed value. `'./results/scRNAseq/GO_in_scRNAseq/correlplots_percExp.pdf` shows the result where the length of each bar shows the correlation coefficient and a darker color shows the significant results based on the permutation test (corrected for multiple testing using BY procedure). The correlations are calculated using Spearman's correlation.

# Helper functions/data

## Gene ID & Orthology mapping

In order to use the same database and package versions for geneID conversions and orthology mapping across scripts, we obtained gene IDs and orthologs from Ensembl database using a small wrapper R package called [dataCollectR](https://github.com/mdonertas/dataCollectR) (in ./scripts/scRNAseq/02-seurat.R). Mapping between multiple type of gene IDs (gene external name, Ensembl Gene ID, Entrez Gene ID, Gene Biotype, and gene descriptions) and orthology mapping between human and *N. furzeri* genes were obtained on 2022.02.17 using biomaRt R Package (v2.50.2) (@durinck2009). These data are available under (`./data/processed/helperdata/`).

## Gene Ontology data

For visualization of Gene ontology categories and testing the association between repair and replication and repair and immune activation, we compiled gene ontology - gene association data that takes the GO hierarchy into account and propagates the genes included in child terms to ancestors. We obtained the Gene Ontology (GO) and gene associations on 2022.04.07 using `GO.db`, `AnnotationDbi`, and `org.Hs.eg.db` packages in R (`./scripts/helperscripts/GO2Gene.R`).

# Manuscript figures

-   Figure 1d - `./scripts/msfigures/plasmaProteomics_GOgenes.R`

-   Figure 2a & 2b - `./scripts/msfigures/scRNAseq_tSNE_celltype.R`

-   Figure 3a - `./scripts/msfigures/KMproteomics_GOresults.R`

-   Figure 3b & 3c - `./scripts/msfigures/KMproteomics_repairGOcategory_overlaps.R`

-   Figure 4 - `./scripts/scRNAseq/06-GOanalysis.R`

-   Figure S1 - `./scripts/scRNAseq/02-seurat.R`

-   Figure S2 - `./scripts/scRNAseq/03-seurat_afterQC.R`

-   Figure S3a - `./scripts/KMproteomics/02-EDA_QC.R`

-   Figure S3b & S3c - `./scripts/msfigures/KMproteomics_GOresults.R`

-   Figure SX - proteomics, effect of the outliers - `./scripts/msfigures/proteomics_outliercheck.R`
