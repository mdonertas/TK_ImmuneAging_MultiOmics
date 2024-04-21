source('./scripts/00-setup.R')

library(Seurat)

load('data/public_scRNASeq/teefy.RData') #creates an object killi.combined.filt 

kidney = subset(x = killi.combined.filt, subset = Tissue == "Kidney")

# choose 20% of male kidney cells
to_annotate_ids = (kidney@meta.data %>%
                     mutate(cell_id = rownames(.)) %>%
                     filter(Sex == "M") %>%
                     group_by(Annotation_v1) %>%
                     sample_n(n()/20))$cell_id


# exclude the cells we chosen
ref_ids = which(!rownames(killi.combined.filt@meta.data) %in% to_annotate_ids)

# create a data with all cells except 20% kidney male cells as reference
annot_ref = subset(killi.combined.filt, cells = ref_ids)

objsave(annot_ref, './data/processed/teefy/annot_ref')

annot_ref
#An object of class Seurat 
#38740 features across 80921 samples within 3 assays 
#Active assay: SCT (18370 features, 5000 variable features)
#3 layers present: counts, data, scale.data
#2 other assays present: RNA, integrated
#2 dimensional reductions calculated: pca, umap

# create a data with 20% male kidney cells
to_annot = subset(killi.combined.filt, cells = to_annotate_ids)

objsave(to_annot, './data/processed/teefy/to_annot')

to_annot
#An object of class Seurat 
#38740 features across 436 samples within 3 assays 
#Active assay: SCT (18370 features, 5000 variable features)
#3 layers present: counts, data, scale.data
#2 other assays present: RNA, integrated
#2 dimensional reductions calculated: pca, umap

trans.anc = FindTransferAnchors(reference = annot_ref, query = to_annot, 
                                dims = 1:30, reference.reduction = "pca")

objsave(trans.anc, './data/processed/teefy/trans.anc')
preds = TransferData(anchorset = trans.anc, refdata = annot_ref$Annotation_v1, 
                     dims = 1:30)
objsave(preds, './data/processed/teefy/preds')
pred_results = preds %>%
  mutate(cell_id = rownames(.)) %>%
  left_join(mutate(killi.combined.filt@meta.data, 
                   cell_id = rownames(killi.combined.filt@meta.data))) %>%
  select(predicted.id, cell_id, Annotation_v1)

objsave(pred_results, './data/processed/teefy/pred_results')

pred_mat = pred_results %>%
  group_by(Annotation_v1, predicted.id) %>%
  summarise(n = n()) %>%
  spread(predicted.id, n, fill = 0) %>%
  as.data.frame()

rownames(pred_mat) = pred_mat$Annotation_v1
pred_mat$Annotation_v1 = NULL
pred_mat = as.matrix(pred_mat)

objsave(pred_mat, './data/processed/teefy/pred_mat')
pred_rel_mat = pred_mat/rowSums(pred_mat) * 100
objsave(pred_rel_mat, './data/processed/teefy/pred_rel_mat')
pheatmap::pheatmap(t(pred_rel_mat), cluster_rows = F, cluster_cols = F, 
                   color = RColorBrewer::brewer.pal(8,'Oranges'), 
                   cellwidth = 15, cellheight = 15,fontsize = 8, 
                   display_numbers = T, fontsize_number = 5, 
                   number_format = '%.0f', filename = './results/scRNAseq/teefy/reannotate_20p.pdf')


