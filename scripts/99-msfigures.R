source('./scripts/00-setup.R')

# scRNAseq tsne + cell type prop change

tsne = readRDS('./results/scRNAseq/seurat_afterQC/tsne_annotated.rds')
ctype_prop = readRDS('./results/scRNAseq/celltypes/celltype_permutations.rds')

celltypes = ggarrange(tsne, ctype_prop, ncol =2 , nrow = 1, widths = c(1,1), labels = 'auto', font.label = list(size =10))

plotsave(celltypes, './results/msfigures/celltypes',width = 16,height = 6)
