source('./scripts/00-setup.R')

acc = as.matrix(read.csv('./data/processed/ml_acc.txt', row.names = 1))
rownames(acc) = c('Erythrocytes',
                  'Progenitors-1',
                  'Myelomonocytes',
                  'Lymphocytes',
                  'Polynucleated giant cells',
                  'Progenitors-2')
colnames(acc) = rownames(acc)

library(ComplexHeatmap)
basesize = 8

pdf('~/Desktop/ml_acc.pdf',width = 4, height = 3.5)
Heatmap(t(acc), cluster_rows = F, cluster_columns = F, col =RColorBrewer::brewer.pal(8,'Reds'), border = T, row_title = 'True Class', column_title = 'Predicted Class', row_names_side = 'left', column_names_side = 'top', column_names_rot = 60, cell_fun = function(j, i, x, y, width, height, fill) {  grid.text(sprintf("%.1f", acc[j, i]), x, y, gp = gpar(fontsize = basesize, col = 'gray25', font = 2))}, row_names_gp = gpar(fontsize = basesize), column_names_gp = gpar(fontsize = basesize) ,heatmap_legend_param = list(title = '%'))
dev.off()
