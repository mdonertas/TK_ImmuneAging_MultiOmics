source('./scripts/00-setup.R')
library(Seurat)
alldat = readRDS('./data/processed/seurat/all_separate.rds')
sort(table(unlist(sapply(alldat, function(x)rownames(x@meta.data)))),dec=T)[1:20]
# cell IDs are duplicated across samples, need to be careful.
names(alldat)

oldintegrated = readRDS('./data/processed/seurat/allintegrated.rds')
to_include = oldintegrated@meta.data %>%
  mutate(datname = paste(SampleID, runID, sep='_')) %>%
  mutate(full_cellID = rownames(oldintegrated@meta.data)) %>%
  separate(full_cellID, into = c('cellID','num'), sep = '_', remove = F) %>%
  filter(!seurat_clusters %in% 11:14) 
rm(oldintegrated)

alldat2 = sapply(names(alldat)[-4],function(nm){
  x = subset(alldat[[nm]], cells = filter(to_include, datname == nm)$cellID)
  x = subset(x, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
  x <- NormalizeData(x, verbose = FALSE)
  FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

sapply(names(alldat2), function(nm){
  length(unique(setdiff(rownames(alldat[[nm]]@meta.data),
                        rownames(alldat2[[nm]]@meta.data))))
})

# A_run636 A_run640 A_run641 C_run636 C_run640 C_run641 D_run636 D_run640 D_run641 E_run636 E_run640 E_run641 F_run636 F_run640 F_run641 
# 109      112      115       58       65       66       68       70       68       83       92       97       65       66       70
# 336 191 206 272 201

anchors <- FindIntegrationAnchors(object.list = alldat2, dims = 1:20)

# objsave(anchors, './data/processed/seurat_afterQC/anchors')


allint <- IntegrateData(anchorset = anchors, dims = 1:20)

DefaultAssay(allint) <- "integrated"

allint <- ScaleData(allint, verbose = FALSE)
allint <- RunPCA(allint, npcs = 30, verbose = FALSE)
allint <- FindNeighbors(allint, reduction = "pca", dims = 1:20)

for(resx in c(0.1,0.2,0.25,0.3,0.4,0.5,0.75,1,1.5)){
  allint <- FindClusters(allint, resolution = resx)
}

library(clustree)
ctree = clustree(allint@meta.data, prefix = "integrated_snn_res.")
ctree
plotsave(ctree, './results/scRNAseq/seurat_afterQC/cluster_resolutions', width = 25, height = 25)

allint <- FindClusters(allint, resolution = 0.3, n.iter = 100, n.start = 100)


allint <- RunUMAP(allint, reduction = "pca", dims = 1:20)
allint <- RunTSNE(allint, reduction = "pca", dims = 1:20)

# objsave(allint, './data/processed/seurat_afterQC/allintegrated')

allint = readRDS('./data/processed/seurat_afterQC/allintegrated.rds')

umapdat = data.frame(allint@reductions$umap@cell.embeddings) %>%
  mutate(cell = rownames(allint@reductions$umap@cell.embeddings)) %>%
  left_join(mutate(allint@meta.data, cell = rownames(allint@meta.data))) %>%
  left_join(data.frame(cluster=Idents(allint), cell = names(Idents(allint))))

umaptext = umapdat %>%
  group_by(cluster) %>%
  summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2))

umapplot = ggplot(umapdat, aes(x = UMAP_1, y = UMAP_2, color = cluster)) +
  geom_point(alpha = 0.5, size = 0.1) +
  geom_label(data = umaptext, aes(label = cluster), color = 'black', 
             fontface = 'bold', alpha = 0.2, label.size = NA, fill = 'gray',
             size = 6/pntnorm, label.padding = unit(1,'pt')) +
  guides(fill = 'none', color = 'none')

plotsave(umapplot, './results/scRNAseq/seurat_afterQC/umap_all', width = 8, height = 8)


tsnedat = data.frame(allint@reductions$tsne@cell.embeddings) %>%
  mutate(cell = rownames(allint@reductions$tsne@cell.embeddings)) %>%
  left_join(mutate(allint@meta.data, cell = rownames(allint@meta.data))) %>%
  left_join(data.frame(cluster=Idents(allint), cell = names(Idents(allint))))

tsnetext = tsnedat %>%
  group_by(cluster) %>%
  summarise(tSNE_1 = median(tSNE_1), tSNE_2 = median(tSNE_2))

tsneplot = ggplot(tsnedat, aes(x = tSNE_1, y = tSNE_2, color = cluster)) +
  geom_point(alpha = 0.5, size = 0.1) +
  geom_label(data = tsnetext, aes(label = cluster), color = 'black', 
             fontface = 'bold', alpha = 0.2, label.size = NA, fill = 'gray',
             size = 6/pntnorm, label.padding = unit(1,'pt')) +
  guides(fill = 'none', color = 'none') 

plotsave(tsneplot, './results/scRNAseq/seurat_afterQC/tsne_all', width = 8, height = 8)

DefaultAssay(allint) <- "RNA"
markers = FindAllMarkers(allint, min.pct = 0.25, only.pos = T, logfc.threshold = 0.5)

humanids = readRDS('./data/processed/helperdata/humanids.rds')
humanids2 = rename(humanids$geneIDs, human_ortholog = ensembl_gene_id) %>% select(1,2,6) %>% set_names(c('humangene','human_ortholog','human_description'))
humanids2$human_description = sapply(strsplit(humanids2$human_description,' [[]Source:'),function(x)x[1])

nfuids = readRDS('./data/processed/helperdata/nfuids.rds')

nfuids = nfuids$geneIDs %>%
  mutate(datagenes = ifelse(external_gene_name=='',ensembl_gene_id,external_gene_name))

nfuids$description = sapply(strsplit(nfuids$description,' [[]Source:'),function(x)x[1])

humanort = readRDS('./data/processed/helperdata/humanorthologs.rds')
human = humanort$orthologs %>%
  select(2,'hsapiens_homolog_ensembl_gene') %>%
  set_names('nfu_gene','ensembl_gene_id') %>%
  left_join(humanids$geneIDs) %>%
  select(1,2) %>%
  set_names('ensembl_gene_id','human_ortholog')

markers2 = markers %>%
  rename(datagenes = gene) %>%
  left_join(nfuids) %>%
  left_join(human) %>%
  left_join(humanids2) %>%
  unique() %>%
  select(-entrezgene_id,-gene_biotype, -dataset) 

markers_14v15 = FindMarkers(allint, ident.1 = 14, ident.2 = 15)
markers_14v15 = markers_14v15 %>%
  mutate(gene = rownames(markers_14v15))

markers_14v15 = markers_14v15 %>%
  mutate(datagenes=gene) %>%
  left_join(nfuids) %>%
  left_join(human) %>%
  left_join(humanids2) %>%
  select(-entrezgene_id,-gene_biotype,-dataset,-datagenes) %>%
  unique()

tablesave(markers_14v15, './results/scRNAseq/seurat_afterQC/DEcomp14v15')

markers_7v11 = FindMarkers(allint, ident.1 = 7, ident.2 = 11)
markers_7v11 = markers_7v11 %>%
  mutate(gene = rownames(markers_7v11))

markers_7v11 = markers_7v11 %>%
  mutate(datagenes=gene) %>%
  left_join(nfuids) %>%
  left_join(human) %>%
  left_join(humanids2) %>%
  select(-entrezgene_id,-gene_biotype,-dataset,-datagenes) %>%
  unique()

tablesave(markers_7v11, './results/scRNAseq/seurat_afterQC/DEcomp7vs11')

objsave(markers, './data/processed/seurat_afterQC/markers')  
objsave(markers2, './data/processed/seurat_afterQC/markers_humanorthologs')  

tablesave(markers2, './results/scRNAseq/seurat_afterQC/clustermarkers')


