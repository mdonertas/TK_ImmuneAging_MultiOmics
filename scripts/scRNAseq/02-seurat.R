source('./scripts/00-setup.R')

# we obtain gene IDs and orthologs from Ensembl database using a small wrapper
# R package called dataCollectR. The advantage is this package will save 
# information about the ensembl version and data accessed for reproducibility 
# purposes. and the same ID data can be used across scripts. 

# humanids = dataCollectR::get_all_geneIDs('hsapiens')
# 
# humanort = dataCollectR::tidy_orthologs(organism = 'nfurzeri', target = 'hsapiens')
# 
# objsave(humanids, './data/processed/helperdata/humanids')
# objsave(humanort, './data/processed/helperdata/humanorthologs')

# nfuids = dataCollectR::get_all_geneIDs('nfurzeri')
# objsave(nfuids, './data/processed/helperdata/nfuids')
nfuids = readRDS('./data/processed/helperdata/nfuids.rds')

nfuids = nfuids$geneIDs %>%
  mutate(datagenes = ifelse(external_gene_name=='',ensembl_gene_id,external_gene_name))

humanort = readRDS('./data/processed/helperdata/humanorthologs.rds')
humanids = readRDS('./data/processed/helperdata/humanids.rds')

human = humanort$orthologs %>%
  select(2,'hsapiens_homolog_ensembl_gene') %>%
  set_names('nfu_gene','ensembl_gene_id') %>%
  left_join(humanids$geneIDs) %>%
  select(1,2) %>%
  set_names('ensembl_gene_id','human_ortholog')

nfuids$description = sapply(strsplit(nfuids$description,' [[]Source:'),function(x)x[1])

library(Seurat)

samplemeta = data.frame(SampleID = LETTERS[1:6],
                        Age = rep(c('young','old'),each=3)) %>%
  mutate(Age = factor(Age, levels = c('young', 'old')))


sampleIDs = LETTERS[1:6]
runIDs = paste('run', c(636, 640, 641), sep='')

samples = intersect(paste(rep(sampleIDs, each=3), runIDs, sep='_'),
                    list.files('./data/processed/'))

alldat = sapply(samples, function(sx){
  sampid = strsplit(sx, '_')[[1]][1]
  run = strsplit(sx, '_')[[1]][2]
  filex = paste("./data/processed/", sx,
                "/outs/filtered_feature_bc_matrix/", 
                sep='')
  mydat = CreateSeuratObject(counts = Read10X(filex), 
                             project = 'TK_HeadKidneyAgeing', min.cells = 10, 
                             min.features = 200)
  if(sampid %in% LETTERS[1:3]){
    mydat = AddMetaData(mydat, metadata = 'young', col.name = 'Age')
  } else {
    mydat = AddMetaData(mydat, metadata = 'old', col.name = 'Age')
  } 
  mydat = AddMetaData(mydat, metadata = sampid, col.name = 'SampleID')
  AddMetaData(mydat, metadata = run, col.name = 'runID')
})

# how many cells and features does each sample have?
dims = sapply(alldat,function(dat)dim(dat@assays$RNA))

dims = data.frame(t(dims)) %>%
  set_names(c('numCells', 'numFeat')) %>%
  mutate(samples = colnames(dims))

# how many UMIs in cells
numumi = sapply(alldat,function(dat)summary(colSums(dat@assays$RNA[]>0)))

numumi = data.frame(t(numumi)) %>%
  set_names(paste('UMI',c('min','Q1','median','mean','Q3','max'),sep='_')) %>%
  mutate(samples = colnames(numumi))

# total reads
totread = sapply(alldat,function(dat)summary(colSums(dat@assays$RNA[])))

totread = data.frame(t(totread)) %>%
  set_names(paste('totRead',c('min','Q1','median','mean','Q3','max'),sep='_')) %>%
  mutate(samples = colnames(totread))

sampleinfo = left_join(dims, numumi) %>%
  left_join(totread) %>%
  separate(samples, into = c('SampleID', 'runID'), remove = F) %>%
  left_join(samplemeta)

numcell_p = sampleinfo %>%
  ggplot(aes(x = SampleID, fill = Age, y = numCells)) +
  geom_bar(stat = 'identity') +
  geom_text(aes(label = numCells), vjust = 1, color = 'white', 
            size = 6/pntnorm, nudge_y = -15) +
  scale_fill_manual(values = agecol) +
  ylab('Number of Cells') +
  facet_wrap(~runID)

numfeat_p = sampleinfo %>%
  ggplot(aes(x = SampleID, fill = Age, y = numFeat)) +
  geom_bar(stat = 'identity') +
  geom_text(aes(label = numFeat), vjust = 1, color = 'white', 
            size = 6/pntnorm, nudge_y = -15) +
  scale_fill_manual(values = agecol) +
  ylab('Number of Features') +
  facet_wrap(~runID)

umicnt_p = sampleinfo %>%
  ggplot(aes(x = SampleID, y = UMI_median, color = Age)) +
  geom_point() +
  geom_segment(aes(xend = SampleID, y = UMI_Q1, yend = UMI_Q3)) +
  ylab('UMI count per cell') +
  scale_color_manual(values = agecol) +
  facet_wrap(~runID)

totread_p = sampleinfo %>%
  ggplot(aes(x = SampleID, y = totRead_median, color = Age)) +
  geom_point() +
  geom_segment(aes(xend = SampleID, y = totRead_Q1, yend = totRead_Q3)) +
  ylab('Total read count per cell') +
  scale_color_manual(values = agecol) +
  facet_wrap(~runID)

allqc = ggarrange(numcell_p, numfeat_p, umicnt_p, totread_p, ncol = 1, nrow = 4,
                  common.legend = T, legend = 'bottom', labels = 'auto', 
                  font.label = list(size = 8), align = 'hv')

plotsave(numcell_p, './results/scRNAseq/qc_prenorm/numcells', width = 16, height = 5)
plotsave(numfeat_p, './results/scRNAseq/qc_prenorm/numfeats', width = 16, height = 5)
plotsave(umicnt_p, './results/scRNAseq/qc_prenorm/umicnts', width = 16, height = 5)
plotsave(totread_p, './results/scRNAseq/qc_prenorm/totreads', width = 16, height = 5)
plotsave(allqc, './results/scRNAseq/qc_prenorm/allqc', width = 16, height = 20)

# There seems to be a problem with sample B, it has both less cells and less 
# features identified overall, and within each cell the number of UMIs and total
# counts are also low. We will exclude sample B and use the rest for analysis

# alldat = lapply(alldat, function(x){
#   x = subset(x, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
#   x
# })

# only consider the cells with number of feautres between 200 and 2500.

# allist = lapply(alldat[-4], function(x){
#   x <- NormalizeData(x, verbose = FALSE)
#   FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
# })
# 
# anchors <- FindIntegrationAnchors(object.list = allist, dims = 1:20)

# objsave(alldat, './data/processed/seurat/all_separate')
# objsave(allist, './data/processed/seurat/all_separate_normalised_noB')
# objsave(anchors, './data/processed/seurat/anchors')

alldat=readRDS('./data/processed/seurat/all_separate.rds')
allist=readRDS('./data/processed/seurat/all_separate_normalised_noB.rds')
anchors=readRDS('./data/processed/seurat/anchors.rds')

# allint <- IntegrateData(anchorset = anchors, dims = 1:20)
# 
# DefaultAssay(allint) <- "integrated"
# 
# allint <- ScaleData(allint, verbose = FALSE)
# allint <- RunPCA(allint, npcs = 30, verbose = FALSE)
# 
# allint <- RunUMAP(allint, reduction = "pca", dims = 1:10)
# allint <- FindNeighbors(allint, reduction = "pca", dims = 1:10)
# allint <- RunTSNE(allint, reduction = "pca", dims = 1:10)
# 
# for(resx in c(0.1,0.2,0.25,0.3,0.5,1,1.5)){
#   allint <- FindClusters(allint, resolution = resx)
# }
# library(clustree)
# ctree = clustree(allint@meta.data, prefix = "integrated_snn_res.")
# plotsave(ctree, './results/scRNAseq/seurat/cluster_resolutions', width = 25, height = 20)

# allint <- FindClusters(allint, resolution = 0.3)
# objsave(allint, './data/processed/seurat/allintegrated')

allint = readRDS('./data/processed/seurat/allintegrated.rds')

tsneplot = data.frame(allint@reductions$tsne@cell.embeddings) %>%
  mutate(cell = rownames(allint@reductions$tsne@cell.embeddings)) %>%
  left_join(mutate(allint@meta.data, cell = rownames(allint@meta.data))) %>%
  left_join(data.frame(cluster=Idents(allint), cell = names(Idents(allint)))) %>%
  ggplot(aes(x = tSNE_1, y = tSNE_2, color = cluster)) +
  geom_point(alpha = 0.3, size = 0.05) +
  # guides(color = 'none') +
  facet_wrap(~SampleID, ncol = 5) +
  guides(color = 'none')

plotsave(tsneplot, './results/scRNAseq/seurat/tsneplot_bysamples', width = 16, height = 7)

tsneplot2 = data.frame(allint@reductions$tsne@cell.embeddings) %>%
  mutate(cell = rownames(allint@reductions$tsne@cell.embeddings)) %>%
  left_join(mutate(allint@meta.data, cell = rownames(allint@meta.data))) %>%
  left_join(data.frame(cluster=Idents(allint), cell = names(Idents(allint)))) %>%
  ggplot(aes(x = tSNE_1, y = tSNE_2, color = cluster)) +
  geom_point(alpha = 0.3, size = 0.05) +
  # guides(color = 'none') +
  facet_wrap(~runID, ncol = 5)  +
  guides(color = 'none')

plotsave(tsneplot2, './results/scRNAseq/seurat/tsneplot_byruns', width = 16, height = 7)

umapplot = data.frame(allint@reductions$umap@cell.embeddings) %>%
  mutate(cell = rownames(allint@reductions$umap@cell.embeddings)) %>%
  left_join(mutate(allint@meta.data, cell = rownames(allint@meta.data))) %>%
  left_join(data.frame(cluster=Idents(allint), cell = names(Idents(allint)))) %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2, color = cluster)) +
  geom_point(alpha = 0.3, size = 0.05) +
  # guides(color = 'none') +
  facet_wrap(~SampleID, ncol = 5)+
  guides(color = 'none')

plotsave(umapplot, './results/scRNAseq/seurat/umapplot_bysamples', width = 16, height = 7)

umapplot2 = data.frame(allint@reductions$umap@cell.embeddings) %>%
  mutate(cell = rownames(allint@reductions$umap@cell.embeddings)) %>%
  left_join(mutate(allint@meta.data, cell = rownames(allint@meta.data))) %>%
  left_join(data.frame(cluster=Idents(allint), cell = names(Idents(allint)))) %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2, color = cluster)) +
  geom_point(alpha = 0.3, size = 0.05) +
  # guides(color = 'none') +
  facet_wrap(~runID, ncol = 5)+
  guides(color = 'none')

plotsave(umapplot2, './results/scRNAseq/seurat/umapplot_byruns', width = 16, height = 7)


umapplotly = data.frame(allint@reductions$umap@cell.embeddings) %>%
  mutate(cell = rownames(allint@reductions$umap@cell.embeddings)) %>%
  left_join(mutate(allint@meta.data, cell = rownames(allint@meta.data))) %>%
  left_join(data.frame(cluster=Idents(allint), cell = names(Idents(allint)))) %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2, color = cluster)) +
  geom_point(alpha = 0.3, size = 0.05) 
umapplotly = plotly::ggplotly(umapplotly)
htmlwidgets::saveWidget(umapplotly, "./results/scRNAseq/seurat/umapinteractive.html")

DefaultAssay(allint) <- "RNA"
markers = FindAllMarkers(allint, min.pct = 0.25, only.pos = T, logfc.threshold = 0.5)

humanids2 = rename(humanids$geneIDs, human_ortholog = ensembl_gene_id) %>% select(1,2,6) %>% set_names(c('humangene','human_ortholog','human_description'))
humanids2$human_description = sapply(strsplit(humanids2$human_description,' [[]Source:'),function(x)x[1])

markers2 = markers %>%
  rename(datagenes = gene) %>%
  left_join(nfuids) %>%
  left_join(human) %>%
  left_join(humanids2) %>%
  unique() %>%
  select(-entrezgene_id,-gene_biotype, -dataset) 

# objsave(markers, './data/processed/seurat/markers')  
# objsave(markers2, './data/processed/seurat/markers_humanorthologs')  

data.frame(allint@reductions$umap@cell.embeddings) %>%
  mutate(cell = rownames(allint@reductions$umap@cell.embeddings)) %>%
  left_join(mutate(allint@meta.data, cell = rownames(allint@meta.data))) %>%
  left_join(data.frame(cluster=Idents(allint), cell = names(Idents(allint)))) %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2, color = cluster==3)) +
  geom_point(alpha = 0.3, size = 0.05) +
  guides(color = 'none') 

markers2 %>%
  filter(cluster == 1) %>%
  arrange(-avg_log2FC) %>%
  head(10)

tablesave(markers2, './results/scRNAseq/seurat/clustermarkers')


# filter(avg_log2FC>1) %>%
arrange(markers2, -avg_log2FC) %>%
  filter(cluster == 0)  %>%
  head(15)

cl14markers = FindMarkers(allint, ident.1 = 14, ident.2 = 0, only.pos = T)
cl14markers %>%
  arrange(-avg_log2FC) %>%
  head()

allint@meta.data %>% group_by(Age, SampleID) %>% 
  summarise(n = n(), meanRead = mean(nCount_RNA), 
            meanFeat = mean(nFeature_RNA)) %>%
  arrange(-n)
