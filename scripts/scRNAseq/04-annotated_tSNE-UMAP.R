# annotated tSNE plot for publication
source('./scripts/00-setup.R')

library(Seurat)

allint <- readRDS('./data/processed/seurat_afterQC/allintegrated.rds')

labels <- read_csv('./results/scRNAseq/seurat_afterQC/clusterLabels.txt', 
                   col_names = c('seurat_clusters','label')) %>%
  mutate(seurat_clusters = as.character(seurat_clusters))

allint@meta.data = allint@meta.data %>%
  left_join(labels) %>%
  mutate(label2 = paste(label, '\n(cl',seurat_clusters,')',sep='')) %>%
  mutate(cell = rownames(allint@meta.data)) 
rownames(allint@meta.data) = allint@meta.data$cell

tsnedat = data.frame(allint@reductions$tsne@cell.embeddings) %>%
  mutate(cell = rownames(allint@reductions$tsne@cell.embeddings)) %>%
  left_join(mutate(allint@meta.data, cell = rownames(allint@meta.data))) %>%
  left_join(data.frame(cluster=Idents(allint), cell = names(Idents(allint))))

tsnetext = tsnedat %>%
  group_by(cluster, label, label2) %>%
  summarise(tSNE_1 = median(tSNE_1), tSNE_2 = median(tSNE_2))

tsneplot = ggplot(tsnedat, aes(x = tSNE_1, y = tSNE_2, color = cluster)) +
  geom_point(alpha = 0.5, size = 0.1) +
  geom_label(data = tsnetext, aes(label = label), color = 'black', 
             alpha = 0.2, label.size = NA, fill = 'gray',
             size = 6/pntnorm, label.padding = unit(1,'pt')) +
  guides(fill = 'none', color = 'none') +
  xlab('tSNE 1') + ylab('tSNE 2')+ 
  geom_label(data = tsnetext, aes(label = cluster), color = 'black', 
             alpha = 0.2, label.size = NA, fill = 'gray',
             size = 4/pntnorm, label.padding = unit(1,'pt'), nudge_y = -3) + 
  theme_classic(base_size = 8)

plotsave(tsneplot, './results/scRNAseq/seurat_afterQC/tsne_annotated',width = 12, height = 8)

umapdat = data.frame(allint@reductions$umap@cell.embeddings) %>%
  mutate(cell = rownames(allint@reductions$umap@cell.embeddings)) %>%
  left_join(mutate(allint@meta.data, cell = rownames(allint@meta.data))) %>%
  left_join(data.frame(cluster=Idents(allint), cell = names(Idents(allint))))

umaptext = umapdat %>%
  group_by(cluster, label, label2) %>%
  summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2))

umapplot = ggplot(umapdat, aes(x = UMAP_1, y = UMAP_2, color = cluster)) +
  geom_point(alpha = 0.5, size = 0.1) +
  geom_label(data = umaptext, aes(label = label), color = 'black', 
             alpha = 0.2, label.size = NA, fill = 'gray',
             size = 6/pntnorm, label.padding = unit(1,'pt')) +
  guides(fill = 'none', color = 'none') +
  xlab('UMAP 1') + ylab('UMAP 2')+ 
  geom_label(data = umaptext, aes(label = cluster), color = 'black', 
             alpha = 0.2, label.size = NA, fill = 'gray',
             size = 4/pntnorm, label.padding = unit(1,'pt'), nudge_y = -0.5) + 
  theme_classic(base_size = 8)

plotsave(umapplot, './results/scRNAseq/seurat_afterQC/umap_annotated',width = 12, height = 8)
