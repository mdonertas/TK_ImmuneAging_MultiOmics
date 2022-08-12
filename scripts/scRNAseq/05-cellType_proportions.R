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

celltypenums = allint@meta.data %>%
  group_by(SampleID, Age, seurat_clusters, label, label2) %>%
  summarise(n = length(unique(cell)))

celltypenums = celltypenums %>%
  ungroup() %>% group_by(SampleID) %>% summarise( tot = sum(n)) %>%
  right_join(celltypenums) %>%
  mutate(prop = n/tot)

celltypenums = celltypenums %>%
  mutate(Age = factor(Age, levels = c('young', 'old'))) %>%
  mutate(seurat_clusters = as.numeric(as.character(seurat_clusters))) %>%
  mutate(label2 = fct_reorder(label2, seurat_clusters))

allcelltype_props = ggplot(celltypenums, aes(x = reorder(paste(label, ' (',seurat_clusters,')',sep=''), prop), y = prop*100)) +
  geom_point(aes( color = Age )) +
  scale_color_manual(values = agecol) +
  coord_flip() +
  ylab('% abundance') +
  xlab(NULL)

plotsave(allcelltype_props, './results/scRNAseq/celltypes/allcelltype_props',width = 12, height = 8)

allcelltype_props = celltypenums %>%
  group_by(SampleID,Age,label) %>%
  summarise(prop = sum(prop)) %>%
  ggplot(aes(x = reorder(label, prop), y = prop*100)) +
  geom_point(aes( color = Age )) +
  scale_color_manual(values = agecol) +
  coord_flip() +
  ylab('% abundance') +
  xlab(NULL)

plotsave(allcelltype_props, './results/scRNAseq/celltypes/allcelltype_props2',width = 12, height = 8)

celltypenums %>%
  group_by(SampleID,Age,label) %>%
  summarise(prop = sum(prop)) %>%
  group_by(label) %>%
  summarise(prop = mean(prop)) %>%
  arrange(-prop)

# label                   prop
# <chr>                  <dbl>
# 1 Neutrophils           0.399 
# 2 Monocytes             0.149 
# 3 RBCs                  0.120 
# 4 RBC Progenitors       0.104 
# 5 Progenitors           0.0569
# 6 B Cells               0.0511
# 7 Macrophages           0.0453
# 8 T Cells               0.0315
# 9 Macrophages/Monocytes 0.0156
# 10 HSCs                  0.0139
# 11 NK cells              0.0139

celltypeprops = ggplot(celltypenums, aes(x = SampleID, y = prop, fill = label2)) +
  geom_bar(stat = 'identity', color = 'gray50') +
  guides(fill = guide_legend(title=NULL)) +
  ylab('Cell type proportions') 

plotsave(celltypeprops, './results/scRNAseq/celltypes/celltypeprops_bySample',width = 12, height = 10)

celltypeprops = ggplot(celltypenums, aes(x = SampleID, y = prop, fill = label)) +
  geom_bar(stat = 'identity', color = 'gray50') +
  guides(fill = guide_legend(title=NULL)) +
  ylab('Cell type proportions') 

plotsave(celltypeprops, './results/scRNAseq/celltypes/celltypeprops_bySample_summarised',width = 12, height = 10)


myeloid_bar = celltypenums %>%
  filter(seurat_clusters %in% c(0, 1, 5, 12)) %>%
  ggplot(aes(x = SampleID, y = prop, fill = label2)) +
  geom_bar(stat = 'identity', position = 'stack') +
  scale_fill_brewer(type = 'qual', palette = 5) +
  guides(fill = guide_legend(title = NULL)) +
  ylab('Proportion of cell type') +
  annotate('text', label = 'young', x = 1.5, y = 0.77, vjust = 0, size = 6/pntnorm, fontface = 'bold') +
  annotate('text', label = 'old', x = 4, y = 0.77, vjust = 0, size = 6/pntnorm, fontface = 'bold') +
  annotate('segment', x = 1, xend = 2, y = 0.75, yend = 0.75) +
  annotate('segment', x = c(1,2,3,5), xend = c(1,2,3,5), y = 0.73, yend = 0.75) +
  annotate('segment', x = 3, xend = 5, y = 0.75, yend = 0.75) +
  theme(legend.box.background = element_rect(color=NA))

plotsave(myeloid_bar, './results/scRNAseq/celltypes/myeloid_bar',width = 8, height = 8)

myeloid_box = celltypenums %>%
  filter(seurat_clusters %in% c(0, 1, 5, 12)) %>%
  group_by(SampleID, Age, tot) %>%
  summarise(n = sum(n)) %>%
  mutate(prop = n / tot) %>%
  ggplot(aes(x = Age, y = prop)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, size = 1) +
  stat_compare_means(size = 6/pntnorm) +
  ylab('Proportion of myeloid cells') +
  xlab(NULL)

plotsave(myeloid_box, './results/scRNAseq/celltypes/myeloid_box',width = 8, height = 8)


lymph_bar = celltypenums %>%
  filter(seurat_clusters %in% c(4, 13, 8)) %>%
  ggplot(aes(x = SampleID, y = prop, fill = label2)) +
  geom_bar(stat = 'identity', position = 'stack') +
  scale_fill_brewer(type = 'qual', palette = 5) +
  guides(fill = guide_legend(title = NULL)) +
  ylab('Proportion of cell type') +
  annotate('text', label = 'young', x = 1.5, y = 0.16, vjust = 0, size = 6/pntnorm, fontface = 'bold') +
  annotate('text', label = 'old', x = 4, y = 0.16, vjust = 0, size = 6/pntnorm, fontface = 'bold') +
  annotate('segment', x = 1, xend = 2, y = 0.15, yend = 0.15) +
  annotate('segment', x = c(1,2,3,5), xend = c(1,2,3,5), y = 0.14, yend = 0.15) +
  annotate('segment', x = 3, xend = 5, y = 0.15, yend = 0.15) +
  theme(legend.box.background = element_rect(color=NA))

plotsave(lymph_bar, './results/scRNAseq/celltypes/lymph_bar',width = 8, height = 8)

lymph_box = celltypenums %>%
  filter(seurat_clusters %in% c(4, 13, 8)) %>%
  group_by(SampleID, Age, tot) %>%
  summarise(n = sum(n)) %>%
  mutate(prop = n / tot) %>%
  ggplot(aes(x = Age, y = prop)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, size = 1) +
  stat_compare_means(size = 6/pntnorm) +
  ylab('Proportion of lymphocytes') +
  xlab(NULL)

plotsave(lymph_box, './results/scRNAseq/celltypes/lymph_box',width = 8, height = 8)

progenitor_bar = celltypenums %>%
  filter(seurat_clusters %in% c(3, 7, 9, 10, 11, 14, 15)) %>%
  ggplot(aes(x = SampleID, y = prop, fill = label2)) +
  geom_bar(stat = 'identity', position = 'stack') +
  scale_fill_brewer(type = 'qual', palette = 4) +
  guides(fill = guide_legend(title = NULL)) +
  ylab('Proportion of cell type') +
  annotate('text', label = 'young', x = 1.5, y = 0.28, vjust = 0, size = 6/pntnorm, fontface = 'bold') +
  annotate('text', label = 'old', x = 4, y = 0.28, vjust = 0, size = 6/pntnorm, fontface = 'bold') +
  annotate('segment', x = 1, xend = 2, y = 0.27, yend = 0.27) +
  annotate('segment', x = c(1,2,3,5), xend = c(1,2,3,5), y = 0.25, yend = 0.27) +
  annotate('segment', x = 3, xend = 5, y = 0.27, yend = 0.27) +
  theme(legend.box.background = element_rect(color=NA))

plotsave(progenitor_bar, './results/scRNAseq/celltypes/progenitor_bar',width = 8, height = 8)

tablesave(celltypenums,'./results/scRNAseq/celltypes/celltypenums')

progenitor_bar2 = celltypenums %>%
  filter(seurat_clusters %in% c(3, 7, 9, 10, 11, 14, 15)) %>%
  ggplot(aes(x = SampleID, y = prop, fill = label)) +
  geom_bar(stat = 'identity', position = 'stack') +
  scale_fill_brewer(type = 'qual', palette = 4) +
  guides(fill = guide_legend(title = NULL)) +
  ylab('Proportion of cell type') +
  annotate('text', label = 'young', x = 1.5, y = 0.28, vjust = 0, size = 6/pntnorm, fontface = 'bold') +
  annotate('text', label = 'old', x = 4, y = 0.28, vjust = 0, size = 6/pntnorm, fontface = 'bold') +
  annotate('segment', x = 1, xend = 2, y = 0.27, yend = 0.27) +
  annotate('segment', x = c(1,2,3,5), xend = c(1,2,3,5), y = 0.25, yend = 0.27) +
  annotate('segment', x = 3, xend = 5, y = 0.27, yend = 0.27) +
  theme(legend.box.background = element_rect(color=NA))

plotsave(progenitor_bar2, './results/scRNAseq/celltypes/progenitor_bar2',width = 8, height = 8)

progenitor_box = celltypenums %>%
  filter(seurat_clusters %in% c(3, 7, 9, 10, 11, 14, 15)) %>%
  group_by(SampleID, Age, tot) %>%
  summarise(n = sum(n)) %>%
  mutate(prop = n / tot) %>%
  ggplot(aes(x = Age, y = prop)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, size = 1) +
  stat_compare_means(size = 6/pntnorm) +
  ylab('Proportion of progenitor cells') +
  xlab(NULL)

plotsave(progenitor_box, './results/scRNAseq/celltypes/progenitor_box',width = 8, height = 8)

###
# new R session

source('./scripts/00-setup.R')

celltype = readRDS('./results/scRNAseq/celltypes/celltypenums.rds')

cellmat = celltype %>%
  select(SampleID,seurat_clusters,n) %>%
  spread(SampleID,n) %>%
  as.data.frame()

rownames(cellmat) = cellmat$seurat_clusters
cellmat$seurat_clusters = NULL
cellmat = as.matrix(cellmat)

allrares = lapply(1:1000,function(i){
  raredcell = apply(cellmat,2,function(x){
    table(sample(unlist(apply(cbind(0:15,x),1,function(x){
      rep(x[1],x[2])
    })),2500))
  })
  x=raredcell[1,]
  changes = data.frame(t(apply(raredcell,1,function(x){
    change = log2((median(x[3:5]))/median(x[1:2]))
    oldifs = sapply(x[1:2],function(y){
      sapply(x[3:5],function(o){
        log2(o/y)
      })
    })
    Qs=quantile(oldifs,probs = c(0.25,0.75))
    c(change,Qs)
  }))) %>%
    set_names(c('change', 'Q1','Q3')) 
  changes %>%
    mutate(seurat_clusters = as.numeric(rownames(changes))) %>%
    filter(seurat_clusters %in% 0:13)
})

allrares = reshape2::melt(allrares, id.vars = c('change','Q1','Q3','seurat_clusters')) %>%
  rename(perm = L1)

ctypechanges = allrares %>%
  left_join(unique(select(celltype,seurat_clusters, label))) %>%
  mutate(label2 = paste(label,seurat_clusters,sep="_")) %>%
  mutate(label2 = fct_reorder(label2, change)) %>%
  ggplot(aes(y = label2)) +
  geom_vline(xintercept = 0, color = 'darkred', linetype = 'dotted', size = 0.1) +
  geom_errorbar(data = . %>%
                  group_by(label2) %>%
                  summarise(Q1 = min(Q1), Q3 = max(Q3)), aes(xmin=Q1, xmax = Q3), width = 0.3, size = 0.1, color = 'gray35') +
  geom_point(size = 0.1, aes(x = change)) +
  xlab('Log2 Fold Change') +
  ylab(NULL) 

ctype_change_p = ctypechanges +
  annotate('text', x = -0.05 , y = 15, label = 'young', hjust = 1, size = 6/pntnorm)+ 
  annotate('text', x = 0.05 , y = 15, label = 'old', hjust = 0, size = 6/pntnorm)+ 
  theme_classic(base_size = 8)+
  coord_cartesian(clip = 'off')

plotsave(ctype_change_p, './results/scRNAseq/celltypes/celltype_permutations', width = 12, height = 8)
