source('./scripts/00-setup.R')

allmarkers = readRDS('./results/scRNAseq/seurat_afterQC/clustermarkers.rds') %>%
  mutate(updatagenes = toupper(datagenes))
selmarkers = read_csv('./data/raw/scRNAseq_markers_to_plot.csv', col_types = 'cc') %>%
  set_names(c('cluster','datagenes')) %>%
  mutate(updatagenes = toupper(datagenes)) 

markerdat = left_join(selmarkers, allmarkers, by = c('cluster', 'updatagenes'))

genelist = unique(markerdat$datagenes.y)

newgenenames = markerdat %>% 
  rowwise() %>%
  mutate(newname = ifelse(grepl('ENSNFUG',datagenes.y), humangene, datagenes.y)) %>%
  select(datagenes.y, humangene, newname) %>%
  unique() %>%
  filter(grepl('ENSNFUG',datagenes.y)) %>%
  filter(newname !='') %>%
  group_by(datagenes.y) %>%
  summarise(newname = paste(newname, collapse = '|'))


newgenenames$newname = paste(newgenenames$newname,'(human)')

allint <- readRDS('./data/processed/seurat_afterQC/allintegrated.rds')

labels <- read_csv('./results/scRNAseq/seurat_afterQC/clusterLabels.txt', 
                   col_names = c('seurat_clusters','label')) %>%
  mutate(seurat_clusters = as.character(seurat_clusters))

allint@meta.data = allint@meta.data %>%
  left_join(labels) %>%
  mutate(label2 = paste(label, '\n(cl',seurat_clusters,')',sep='')) %>%
  mutate(cell = rownames(allint@meta.data))  
rownames(allint@meta.data) = allint@meta.data$cell

metadat = allint@meta.data %>%
  mutate(cellID = rownames(allint@meta.data)) %>%
  select(cellID,label,label2,seurat_clusters,Age,SampleID)

mydat = reshape2::melt(as.matrix(allint@assays$RNA[genelist,])) %>%
  set_names(c('datagenes','cellID','value')) %>%
  left_join(metadat)

datsc = t(apply(as.matrix(allint@assays$RNA[genelist,]), 1 , scale))
dimnames(datsc) = dimnames(as.matrix(allint@assays$RNA[genelist,]))

mydat2 = t(sapply(unique(metadat$seurat_clusters), function(cl){
  rowMeans(datsc[genelist,unique(filter(metadat, seurat_clusters == cl)$cellID)])
}))

rownames(mydat2) = paste(unname(setNames(labels$label,labels$seurat_clusters)[rownames(mydat2)]), ' (cl.', rownames(mydat2), ')', sep='')

scdata = reshape2::melt(as.matrix(allint@assays$RNA[genelist,])) %>%
  set_names('Gene', 'cellID', 'Normalized Expression') %>%
  left_join(metadat) %>%
  group_by(Gene, label2) %>%
  summarise(meanExp = mean(`Normalized Expression`),
            percExp = mean(`Normalized Expression`>0)*100) %>%
  ungroup() %>%
  mutate(Gene = factor(Gene, labels = rev(colnames(mydat2)[hclust(dist(t(mydat2)))$order])))  

tsnedat <- data.frame(allint@reductions$tsne@cell.embeddings) %>%
  mutate(cellID = rownames(allint@reductions$tsne@cell.embeddings))

tsnelabels = left_join(metadat,tsnedat) %>%
  group_by(seurat_clusters) %>%
  summarise(tSNE_1 = median(tSNE_1),
            tSNE_2 = median(tSNE_2))

alltsne = data.frame(allint@reductions$tsne@cell.embeddings) %>%
  mutate(cellID = rownames(allint@reductions$tsne@cell.embeddings))  %>%
  inner_join((reshape2::melt(as.matrix(allint@assays$RNA[genelist,])) %>%
               set_names('Gene', 'cellID', 'Normalized Expression') %>%
               left_join(metadat))) %>%
  ggplot(aes(x = tSNE_1, y = tSNE_2, color = `Normalized Expression`)) +
  geom_point(size = 0.3, alpha = 0.5) +
  scale_color_gradient(high = '#800026', low = "#FFFFCC") +
  geom_text(data = tsnelabels, aes(label = seurat_clusters),
            color = 'black',
            size = 4/pntnorm) +
  guides(color = guide_colorbar(title= 'Normalized\nExpression')) +
  theme_void(base_size = 6) +
  facet_wrap(~Gene) +
  theme(legend.position = c(0.5,0.05), legend.justification = c(0,0), legend.direction = 'horizontal', legend.key.height = unit(6,'pt'), legend.key.width = unit(20,'pt'))

plotsave(alltsne, './results/scRNAseq/seurat_afterQC/markertSNE', width = 16, height= 12)
