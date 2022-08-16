source('./scripts/00-setup.R')

plasma_l2 = readRDS('./data/processed/plasmaProteomics/expression_log2.rds')

plasma_l2 = reshape2::melt(plasma_l2) %>%
  set_names(c('protein','SampleID','log2Expression')) %>%
  mutate(uppersym = toupper(protein))
plasma_l2$Age = sapply(strsplit(as.character(plasma_l2$SampleID),'_'),function(x)x[1])

plasma_wilcox = readRDS('./data/processed/plasmaProteomics/DEwilcox_allsamples.rds') %>%
  rename(uppersym = upperprotein)

objsave(plasma_l2, './data/processed/shiny/plasma_l2')
objsave(plasma_wilcox, './data/processed/shiny/plasma_wilcox')

km_l2 = readRDS('./data/processed/KMproteomics/expression_log2.rds')

km_l2 = reshape2::melt(km_l2) %>%
  set_names(c('protein','SampleID','log2Expression')) %>%
  mutate(uppersym = toupper(protein))
km_l2$Age = sapply(strsplit(as.character(km_l2$SampleID),'_'),function(x)x[1])

km_wilcox = readRDS('./data/processed/KMproteomics/DEwilcox_allsamples.rds') %>%
  rename(uppersym = upperprotein)

objsave(km_l2, './data/processed/shiny/km_l2')
objsave(km_wilcox, './data/processed/shiny/km_wilcox')

allint <-
  readRDS('./data/processed/seurat_afterQC/allintegrated.rds')

labels <- read_csv('./results/seurat_afterQC/clusterLabels',
                   col_names = c('seurat_clusters', 'label')) %>%
  mutate(seurat_clusters = as.character(seurat_clusters))

allint@meta.data = allint@meta.data %>%
  left_join(labels) %>%
  mutate(label2 = paste(label, '_cl', seurat_clusters, sep = '')) %>%
  mutate(cell = rownames(allint@meta.data))
rownames(allint@meta.data) = allint@meta.data$cell

metadata = allint@meta.data

metadata = metadata %>%
  select(cell, SampleID, Age, runID, seurat_clusters, label, label2) %>%
  unique()

objsave(metadata, './data/processed/shiny/sc_metadata')
scgenes = toupper(rownames(allint@assays$RNA))
notinsc = names(which(sapply(rownames(km_wilcox),function(x){
  !(x %in% scgenes)
})))
notinsc[2]%in% scgenes
myids = filter(myids, datagenes %in% rownames(allint@assays$RNA))
filter(myids, upperprotein %in% notinsc)$datagenes
'abat'%in%rownames(allint@assays$RNA)


#### after Ensembl team confirmed an inconsistency between web and genome version of NFU genome build v105, we used the gtf file used to create the scRNAseq count matrix to map between IDs. 

library(tidyverse)

gff = rtracklayer::readGFF('../../resources/genome_annotations/NFU/Ensembl/Nothobranchius_furzeri.Nfu_20140520.105.chr.prot_ig_tr.gtf')

head(gff)

gff2 = gff %>%
  select(gene_id, gene_name) %>%
  unique() %>%
  set_names(c('ensembl_gene_id','ensgenename'))  %>%
  mutate(datagenes = ifelse(ensgenename == '' | is.na(ensgenename), ensembl_gene_id,
                            ensgenename)) %>%
  mutate(datagenes = gsub('_','-',datagenes)) %>%
  unique()

# duplgenes = unique(gff2$datagenes[duplicated(gff2$datagenes)])

datagenes = rownames(scdat@assays$RNA)
# datagenes = gsub('-','_',datagenes)

# datagenes[which(!(sapply(strsplit(datagenes,'[.]'),function(x)paste(x[-length(x)],collapse='.')) %in% gff2$datagenes | datagenes %in% gff2$datagenes))]

# setdiff(datagenes, gff2$datagenes)


nfuids = readRDS('../TK_KidneyMarrow_ageing_multiomics/data/processed/helperdata/nfuids.rds')

nfuids = nfuids$geneIDs
nfuids$description = sapply(strsplit(nfuids$description, ' [[]Source:'),
                            function(x)
                              x[1])

gff2 = gff2 %>%
  full_join(nfuids) 

xx = apply(gff2,1,function(x){
  ensg = x['ensembl_gene_id']
  genesym = x['external_gene_name']
  genesym = ifelse(is.na(genesym) | genesym=='', ensg,genesym)
  descr = x['description']
  paste(setdiff(c(genesym,descr), c(NA,'')), collapse = ' - ')
})

gff2$displaygenes = xx

datagenesx = rownames(scdat@assays$RNA)
datagenesx2 = sapply(strsplit(datagenesx,'[.]'),function(x)paste(x[-length(x)],collapse='.'))

allgenes = gff2 %>%
  mutate(uppersym = toupper(external_gene_name)) %>%
  filter(uppersym %in% km_wilcox$uppersym | uppersym %in% plasma_wilcox$uppersym | datagenes %in% datagenesx | datagenes %in% datagenesx2) %>%
  select(external_gene_name, datagenes, displaygenes) %>%
  unique() %>%
  # filter(external_gene_name %in% 'KPNA2') %>%
  group_by(external_gene_name, datagenes) %>%
  summarise(displaygenes = paste(unique(displaygenes), collapse = ' | '))

saveRDS(allgenes, './data/geneids.rds')



