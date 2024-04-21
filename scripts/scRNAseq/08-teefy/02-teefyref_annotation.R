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

load('data/public_scRNASeq/teefy.RData') #creates an object killi.combined.filt 

trans.anc = FindTransferAnchors(reference = killi.combined.filt , 
                                normalization.method = 'SCT',
                                query = allint, 
                                dims = 1:30, reference.reduction = "pca")

objsave(trans.anc, './data/processed/teefyref_annotation/trans.anc')

preds = TransferData(anchorset = trans.anc, 
                     refdata = killi.combined.filt$Annotation_v1, 
                     dims = 1:30)
objsave(preds, './data/processed/teefyref_annotation/preds')

pred_res = preds %>%
  mutate(cell = rownames(.)) %>%
  right_join(allint@meta.data) %>%
  rowwise() %>%
  mutate(label2 = paste(label,'_cl', seurat_clusters, sep = '')) %>%
  select(label2, predicted.id, cell) %>%
  group_by(label2, predicted.id) %>%
  summarise(n = n())

pred_mat = pred_res %>%
  spread(predicted.id, n, fill = 0) %>%
  as.data.frame()

rownames(pred_mat) = pred_mat$label2
pred_mat$label2 = NULL
pred_mat = as.matrix(pred_mat)

objsave(pred_mat, './data/processed/teefyref_annotation/pred_mat')
pred_rel_mat = pred_mat/rowSums(pred_mat) * 100
objsave(pred_rel_mat, './data/processed/teefyref_annotation/pred_rel_mat')
pheatmap::pheatmap(t(pred_rel_mat), 
                   color = RColorBrewer::brewer.pal(8,'Oranges'), 
                   cellwidth = 15, cellheight = 15,fontsize = 8, 
                   display_numbers = T, fontsize_number = 5, 
                   number_format = '%.0f', filename = './results/scRNAseq/teefyref_annotation/annotation_dist.pdf')



