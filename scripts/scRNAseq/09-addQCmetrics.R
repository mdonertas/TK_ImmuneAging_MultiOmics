source('./scripts/00-setup.R')
samplemeta = data.frame(SampleID = LETTERS[1:6],
                        Age = rep(c('young','old'),each=3)) %>%
  mutate(Age = factor(Age, levels = c('young', 'old')))
alldat=readRDS('./data/processed/seurat/all_separate.rds')

dims = sapply(alldat,function(dat)dim(dat@assays$RNA))
dims = data.frame(t(dims)) %>%
  set_names(c('numFeat', 'numCells')) %>%
  mutate(samples = colnames(dims))

featnums = sapply(alldat,function(dat)summary(colSums(dat@assays$RNA[]>0)))
featnums = data.frame(t(featnums)) %>%
  set_names(paste('NumFeat',c('min','Q1','median','mean','Q3','max'),sep='_')) %>%
  mutate(samples = colnames(featnums))

totread = sapply(alldat,function(dat)summary(colSums(dat@assays$RNA[])))

totread = data.frame(t(totread)) %>%
  set_names(paste('totRead',c('min','Q1','median','mean','Q3','max'),sep='_')) %>%
  mutate(samples = colnames(totread))

sampleinfo = left_join(dims, featnums) %>%
  left_join(totread) %>%
  separate(samples, into = c('SampleID', 'runID'), remove = F) %>%
  left_join(samplemeta)


allint = readRDS('./data/processed/seurat/allintegrated.rds')

sampleinfo2 = as.data.frame(t(sapply(LETTERS[c(1,3:6)], function(samp){
  mymat = allint@assays$RNA[,rownames(filter(allint@meta.data, SampleID == samp))]
  numcell = ncol(mymat)
  numfeat = sum(rowMeans(mymat>0)>0)
  numfeat_perCell = setNames(summary(colSums(mymat>0)),paste('NumFeat',c('min','Q1','median','mean','Q3','max'),sep='_'))
  numread_perCell = setNames(summary(colSums(mymat)), paste('totRead',c('min','Q1','median','mean','Q3','max'),sep='_'))
  c(numCells = numcell, numFeat = numfeat, numfeat_perCell, numread_perCell)
}))) %>%
  mutate(SampleID = rownames(.)) %>%
  left_join(samplemeta) %>%
  mutate(runID = 'integrated')


allint_afterQC = readRDS('./data/processed/seurat_afterQC/allintegrated.rds')

sampleinfo3 = as.data.frame(t(sapply(LETTERS[c(1,3:6)], function(samp){
  mymat = allint_afterQC@assays$RNA[,rownames(filter(allint_afterQC@meta.data, SampleID == samp))]
  numcell = ncol(mymat)
  numfeat = sum(rowMeans(mymat>0)>0)
  numfeat_perCell = setNames(summary(colSums(mymat>0)),paste('NumFeat',c('min','Q1','median','mean','Q3','max'),sep='_'))
  numread_perCell = setNames(summary(colSums(mymat)), paste('totRead',c('min','Q1','median','mean','Q3','max'),sep='_'))
  c(numCells = numcell, numFeat = numfeat, numfeat_perCell, numread_perCell)
}))) %>%
  mutate(SampleID = rownames(.)) %>%
  left_join(samplemeta) %>%
  mutate(runID = 'integrated_afterQC')


sampleinfo = add_row(sampleinfo, sampleinfo2) %>%
  add_row(sampleinfo3) %>%
  left_join(read_tsv('./data/processed/seurat/mapping_rate.tsv')) %>%
  mutate(runID = factor(runID, levels = c('run636','run640','run641','integrated','integrated_afterQC'))) 

numcell_p = sampleinfo %>%
  ggplot(aes(x = SampleID, fill = Age, y = numCells)) +
  geom_bar(stat = 'identity') +
  geom_text(aes(label = numCells), vjust = 1, color = 'white', 
            size = 6/pntnorm, nudge_y = -15) +
  scale_fill_manual(values = agecol) +
  ylab('Number of Cells') +
  facet_wrap(~runID, ncol = 5)

numfeat_p = sampleinfo %>%
  ggplot(aes(x = SampleID, fill = Age, y = numFeat)) +
  geom_bar(stat = 'identity') +
  geom_text(aes(label = numFeat), vjust = 1, color = 'white', 
            size = 6/pntnorm, nudge_y = -15) +
  scale_fill_manual(values = agecol) +
  ylab('Number of Features') +
  facet_wrap(~runID, ncol = 5)

maprate_p = sampleinfo %>%
  ggplot(aes(x = SampleID, fill = Age, y = mapped_to_genome)) +
  geom_bar(stat = 'identity') +
  geom_text(aes(label = mapped_to_genome), vjust = 1, color = 'white', 
            size = 6/pntnorm, nudge_y = -10) +
  scale_fill_manual(values = agecol) +
  ylab('% Mapped to Genome') +
  facet_wrap(~runID, ncol = 5)

numfeat_cell_p = sampleinfo %>%
  ggplot(aes(x = SampleID, y = NumFeat_median, color = Age)) +
  geom_point() +
  geom_segment(aes(xend = SampleID, y = NumFeat_Q1, yend =NumFeat_Q3)) +
  ylab('Feature count per cell') +
  scale_color_manual(values = agecol) +
  facet_wrap(~runID, ncol = 5)

totread_p = sampleinfo %>%
  ggplot(aes(x = SampleID, y = totRead_median, color = Age)) +
  geom_point() +
  geom_segment(aes(xend = SampleID, y = totRead_Q1, yend = totRead_Q3)) +
  ylab('Total read count per cell') +
  scale_color_manual(values = agecol) +
  facet_wrap(~runID, ncol = 5) 

allqc = ggarrange(numcell_p, numfeat_p, maprate_p, numfeat_cell_p, totread_p, ncol = 1, nrow = 5,
                  common.legend = T, legend = 'bottom', labels = 'auto', 
                  font.label = list(size = 8), align = 'hv')

totread_p2 = sampleinfo %>%
  mutate(totRead_median = ifelse(runID %in% c('integrated','integrated_afterQC'), NA, totRead_median),
         totRead_Q1 = ifelse(runID %in% c('integrated','integrated_afterQC'), NA, totRead_Q1),
         totRead_Q3 = ifelse(runID %in% c('integrated','integrated_afterQC'), NA, totRead_Q3)) %>%
  ggplot(aes(x = SampleID, y = totRead_median, color = Age)) +
  geom_point() +
  geom_segment(aes(xend = SampleID, y = totRead_Q1, yend = totRead_Q3)) +
  ylab('Total read count per cell') +
  scale_color_manual(values = agecol) +
  facet_wrap(~runID, ncol = 5) 

allqc2 = ggarrange(numcell_p, numfeat_p, maprate_p, numfeat_cell_p, totread_p2, ncol = 1, nrow = 5,
                  common.legend = T, legend = 'bottom', labels = 'auto', 
                  font.label = list(size = 8), align = 'hv')


allqc2

plotsave(numcell_p, './results/scRNAseq/qc_prenorm/numcells', width = 20, height = 5)
plotsave(maprate_p, './results/scRNAseq/qc_prenorm/maprate_p', width = 20, height = 5)
plotsave(numfeat_p, './results/scRNAseq/qc_prenorm/numfeats', width = 20, height = 5)
plotsave(totread_p, './results/scRNAseq/qc_prenorm/totreads', width = 20, height = 5)
plotsave(allqc, './results/scRNAseq/qc_prenorm/allqc', width = 20, height = 20)
plotsave(allqc2, './results/scRNAseq/qc_prenorm/allqc2', width = 20, height = 20)

tablesave(sampleinfo,'./results/scRNAseq/qc_prenorm/sampleinfo')









