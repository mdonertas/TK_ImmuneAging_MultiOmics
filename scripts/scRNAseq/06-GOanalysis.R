source('./scripts/00-setup.R')

#### ID conversion data ####

humanids = readRDS('./data/processed/helperdata/humanids.rds')
humanids2 = rename(humanids$geneIDs, human_ortholog = ensembl_gene_id) %>%
  select(1, 2, 6) %>%
  set_names(c('humangene', 'human_ortholog', 'human_description'))
humanids2$human_description = sapply(strsplit(humanids2$human_description,
                                              ' [[]Source:'), function(x)
                                                x[1])

nfuids = readRDS('./data/processed/helperdata/nfuids.rds')

nfuids = nfuids$geneIDs %>%
  mutate(datagenes = ifelse(external_gene_name == '', ensembl_gene_id,
                            external_gene_name))

nfuids$description = sapply(strsplit(nfuids$description, ' [[]Source:'),
                            function(x)
                              x[1])

humanort = readRDS('./data/processed/helperdata/humanorthologs.rds')

human = humanort$orthologs %>%
  dplyr::select(2, 'hsapiens_homolog_ensembl_gene') %>%
  set_names('nfu_gene', 'ensembl_gene_id') %>%
  left_join(humanids$geneIDs) %>% 
  dplyr::select(1, 2) %>%
  set_names('ensembl_gene_id', 'human_ortholog')

human = left_join(human, humanids2)  %>% 
  select(1, 2) %>%
  unique()

myids = left_join(nfuids, human) %>%
  select(1,datagenes,human_ortholog) %>%
  mutate(upperprotein = toupper(external_gene_name))

rm(human,humanids,humanids2,humanort,nfuids)


#### Read & prep scRNAseq data ####
library(Seurat)

allint <-
  readRDS('./data/processed/seurat_afterQC/allintegrated.rds')

labels <- read_csv('./results/scRNAseq/seurat_afterQC/clusterLabels.txt',
                   col_names = c('seurat_clusters', 'label')) %>%
  mutate(seurat_clusters = as.character(seurat_clusters))

allint@meta.data = allint@meta.data %>%
  left_join(labels) %>%
  mutate(label2 = paste(label, '_cl', seurat_clusters, sep = '')) %>%
  mutate(cell = rownames(allint@meta.data))
rownames(allint@meta.data) = allint@meta.data$cell
rm(labels)

DefaultAssay(allint) <- "RNA"

#### Read & prep GO data ####
gogenes = readRDS('./data/processed/GO2Gene/GO2GeneBP_20220407.rds')

goterm = readRDS('./data/processed/GO2Gene/goterms_20220407.rds')


## this analysis is performed after the proteomics analysis and uses the results of KM proteomics data.
gores = readRDS('./data/processed/KMproteomics/DE_GOenrich_allsamples.rds')
tsnegos = read_csv('./data/processed/gocats_for_tsne.csv', col_names = c('group','representative_description'))

#### KM proteomics results ####

wiltest = readRDS('./data/processed/KMproteomics/DEwilcox_allsamples.rds')

#### GO in scRNAseq ####

#### 1. get all go ids for representatives. ####

goids = gores %>%
  right_join(tsnegos) %>%
  group_by(group) %>%
  summarise(goids = list(unique(ID)))

goids %>%
  group_by(group) %>%
  summarise(goterms = paste(goterm[unlist(goids)], collapse = ', '),
            goids = paste(unlist(goids), collapse = ', ')) %>% 
  tablesave('./results/scRNAseq/GO_in_scRNAseq/rep_goIDs')

goids

#### 2. get all genes under go categories -> ENSEMBL human ####

goids$fullgenelist = sapply(1:nrow(goids), function(i){
  gogenes = unique(unlist(gogenes[goids$goids[[i]]]))
})

#### 3. get NFU orthologs for all genes in the go categories ####

goids$nfugenes = sapply(goids$fullgenelist, function(x){
  setdiff(filter(myids,human_ortholog %in% x)$datagenes,c(NA,''))
})

#### 4. get all NFU genes in scRNAseq ####

goids$allgenes_in_scRNAseq = sapply(goids$nfugenes,function(x){
  intersect(rownames(allint@assays$RNA),x)
})

#### 5.1. get non-overlapping gene set for these genes ####

goids$no_genes_in_scRNAseq = sapply(1:nrow(goids),function(i){
  setdiff(goids$allgenes_in_scRNAseq[[i]],unique(unlist(goids$allgenes_in_scRNAseq[setdiff(1:nrow(goids),i)])))
})

#### 5.2.1. get genes in KM proteomics with |log2FC|>=1 ####

l2fc_1genes = filter(myids, upperprotein %in% filter(wiltest,abs(l2FC)>=1)$upperprotein)$datagenes

goids$l2fc_1genes = sapply(1:nrow(goids), function(i){
  intersect(goids$allgenes_in_scRNAseq[[i]], l2fc_1genes)
})

#### 5.2.2. and non-overlapping ####

goids$no_l2fc_1genes = sapply(1:nrow(goids),function(i){
  setdiff(goids$l2fc_1genes[[i]],unique(unlist(goids$l2fc_1genes[setdiff(1:nrow(goids),i)])))
})

genenums = data.frame(representatives = goids$group, apply(goids[,-1], 2, function(x) sapply(x, length)))

tablesave(genenums, './results/scRNAseq/GO_in_scRNAseq/genes_in_categories')

#### 6.1. calculate % expressed genes in go categories using 5.1. and 5.2. ####

percExpressed = sapply(goids$no_genes_in_scRNAseq, function(x){
  mat = as.matrix((allint@assays$RNA[x,]>0)+1-1)
  colMeans(mat)
})
percExpressed = data.frame(percExpressed, totExpressed = colMeans(as.matrix((allint@assays$RNA[]>0)+1-1)))

colnames(percExpressed) = paste(c(goids$group,'totExpressed'),'_allnonoverlap_perc',sep="")

percExpressed2 = sapply(goids$no_l2fc_1genes, function(x){
  mat = as.matrix((allint@assays$RNA[x,]>0)+1-1)
  colMeans(mat)
})
colnames(percExpressed2) = paste(goids$group,'_l2FC1_perc',sep="")

#### 6.2. calculate mean scaled expression for 5.1. and 5.2. ####

scExp = sapply(goids$no_genes_in_scRNAseq, function(x){
  mat = as.matrix(allint@assays$RNA[x,])
  colMeans(t(scale(t(mat))))
})
colnames(scExp) = paste(goids$group,'_allnonoverlap_scExp',sep="")

scExp2 = sapply(goids$no_l2fc_1genes, function(x){
  mat = as.matrix(allint@assays$RNA[x,])
  colMeans(t(scale(t(mat))))
})
colnames(scExp2) = goids$group

colnames(scExp2) = paste(goids$group,'_l2FC1_scExp',sep="")

go_exp_dat = cbind(percExpressed,percExpressed2,scExp,scExp2)

go_exp_dat = go_exp_dat %>%
  mutate(cell = rownames(go_exp_dat)) %>%
  gather(key = 'type', value = 'value', -cell, -totExpressed_allnonoverlap_perc) %>%
  separate(type, into = c('GO', 'genes', 'expType'), remove = F, sep = '_') %>%
  left_join(select(allint@meta.data, cell, label, seurat_clusters, label2, Age, SampleID)) 

#### 7.1. tSNE for 4 categories using % expressed genes ####

tsnedat = data.frame(allint@reductions$tsne@cell.embeddings) %>%
  mutate(cell = rownames(allint@reductions$tsne@cell.embeddings))

tsnelabels = left_join(allint@meta.data,tsnedat) %>% 
  group_by(seurat_clusters) %>%
  summarise(tSNE_1 = median(tSNE_1),
            tSNE_2 = median(tSNE_2))

percExp_tsne = tsnedat %>%
  left_join(go_exp_dat) %>%
  filter(expType == 'perc' & genes == 'allnonoverlap') %>%
  mutate(xx = value/totExpressed_allnonoverlap_perc) %>%
  ggplot(aes(x = tSNE_1, y = tSNE_2, color = xx)) +
  facet_wrap(~GO) +
  geom_point(size = 0.01, alpha = 0.1) +
  geom_text(data = tsnelabels, aes(label = seurat_clusters), color = 'black', 
            alpha = 0.5,
            size = 4/pntnorm)  + 
  scale_color_gradient2(low = '#08306B', high = '#800026', mid = "gray99", midpoint = 1) +
  theme_classic(base_size = 6) + 
  guides(color = guide_colorbar(title = 'OR')) +
  theme(legend.key.width = unit(4,'pt'), legend.key.height = unit(8,'pt'), legend.background = element_rect(fill=NA, color = NA)) 

labels = allint@meta.data %>% select(seurat_clusters, label) %>% unique() %>% mutate(seurat_clusters = as.numeric(seurat_clusters)) %>% arrange(seurat_clusters) %>%
  mutate(label = paste(seurat_clusters, '-',label,sep=''))

labels = paste(paste(labels$label[1:5],collapse=', '), 
               paste(labels$label[6:10],collapse=', '), 
               paste(labels$label[11:14],collapse=', '),
               paste(labels$label[15:16],collapse=', '),  sep = '\n')


tsneplots2 = annotate_figure(percExp_tsne, bottom = text_grob(labels, size = 6))


plotsave(tsneplots2, './results/scRNAseq/GO_in_scRNAseq/percExp_tsne', width = 8, height = 9)


go_ages_summary = tsnedat %>%
  left_join(go_exp_dat) %>%
  filter(expType == 'perc' & genes == 'allnonoverlap') %>%
  mutate(xx = value/totExpressed_allnonoverlap_perc) 

go_ages_summary = sapply(unique(go_ages_summary$GO),function(go){
  sapply(unique(go_ages_summary$label2), function(celltype){
    dat = filter(go_ages_summary, GO == go & label2 == celltype)
    wilcox.test(xx~Age, data = dat)$p.value
  })
})

go_ages_summary = reshape2::melt(go_ages_summary) %>%
  set_names(c('label2','GO','p')) %>%
  mutate(padj = p.adjust(p, method = 'fdr')) %>%
  mutate(label = ifelse(padj<=0.05,'*',''))

go_ages_scRNAdat = tsnedat %>%
  left_join(go_exp_dat) %>%
  filter(expType == 'perc' & genes == 'allnonoverlap') %>%
  mutate(xx = value/totExpressed_allnonoverlap_perc) %>%
  # mutate(GO = paste(GO,Age,sep='_')) %>%
  group_by(GO, label2,Age,seurat_clusters) %>%
  summarise(meanOR = mean(xx)) %>%
  mutate(Age = factor(Age, levels = c('young','old'))) %>%
  left_join(go_ages_summary) %>%
  ungroup() %>%
  mutate(label2 = fct_reorder(factor(label2), -as.numeric(seurat_clusters))) 
go_ages_scRNAdat$GO = gsub(' ', '\n',go_ages_scRNAdat$GO)
go_ages_scRNA = ggplot(go_ages_scRNAdat, aes(x = Age, y = label2)) +
  facet_wrap(~GO, ncol = 4) +
  geom_tile(aes(fill = meanOR)) +
  geom_text(aes(label = label), size = 6/pntnorm) +
  scale_fill_gradient2(low = '#08306B', high = '#800026', mid = "gray99", midpoint = 1) +
  theme_minimal(base_size = 6) +
  xlab(NULL) + ylab(NULL) +
  guides(fill = guide_colorbar(title = 'mean\nOR')) +
  theme(legend.key.width = unit(4,'pt'), legend.key.height = unit(8,'pt'), legend.background = element_rect(fill=NA, color = NA)) 

plotsave(go_ages_scRNA, './results/scRNAseq/GO_in_scRNAseq/go_ages_scRNA_percExp', width = 8, height = 8)

#### 7.2. tSNE for 4 categories using mean scaled expression ####

scExp_tsne = tsnedat %>%
  left_join(go_exp_dat) %>%
  filter(expType == 'scExp' & genes == 'allnonoverlap') %>%
  ggplot(aes(x = tSNE_1, y = tSNE_2, color = value)) +
  facet_wrap(~GO) +
  geom_point(size = 0.01, alpha = 0.1) +
  geom_text(data = tsnelabels, aes(label = seurat_clusters), color = 'black', 
            alpha = 0.5,
            size = 4/pntnorm)  + 
  scale_color_gradient(high = '#800026', low = "#FFFFCC") +
  theme_classic(base_size = 6) + 
  guides(color = guide_colorbar(title = 'Mean\nExp')) +
  theme(legend.key.width = unit(4,'pt'), legend.key.height = unit(8,'pt'), legend.background = element_rect(fill=NA, color = NA)) 

labels = allint@meta.data %>% select(seurat_clusters, label) %>% unique() %>% mutate(seurat_clusters = as.numeric(seurat_clusters)) %>% arrange(seurat_clusters) %>%
  mutate(label = paste(seurat_clusters, '-',label,sep=''))

labels = paste(paste(labels$label[1:5],collapse=', '), 
               paste(labels$label[6:10],collapse=', '), 
               paste(labels$label[11:14],collapse=', '),
               paste(labels$label[15:16],collapse=', '),  sep = '\n')


tsneplots_scExp = annotate_figure(scExp_tsne, bottom = text_grob(labels, size = 6))

plotsave(tsneplots_scExp, './results/scRNAseq/GO_in_scRNAseq/scExp_tsne', width = 8, height = 9)


go_ages_summary = tsnedat %>%
  left_join(go_exp_dat) %>%
  filter(expType == 'scExp' & genes == 'allnonoverlap') %>%
  mutate(xx = value) 

go_ages_summary = sapply(unique(go_ages_summary$GO),function(go){
  sapply(unique(go_ages_summary$label2), function(celltype){
    dat = filter(go_ages_summary, GO == go & label2 == celltype)
    wilcox.test(xx~Age, data = dat)$p.value
  })
})

go_ages_summary = reshape2::melt(go_ages_summary) %>%
  set_names(c('label2','GO','p')) %>%
  mutate(padj = p.adjust(p, method = 'fdr')) %>%
  mutate(label = ifelse(padj<=0.05,'*',''))

go_ages_scRNAdat = tsnedat %>%
  left_join(go_exp_dat) %>%
  filter(expType == 'scExp' & genes == 'allnonoverlap') %>%
  mutate(xx = value) %>%
  # mutate(GO = paste(GO,Age,sep='_')) %>%
  group_by(GO, label2,Age,seurat_clusters) %>%
  summarise(meanOR = mean(xx)) %>%
  mutate(Age = factor(Age, levels = c('young','old'))) %>%
  left_join(go_ages_summary) %>%
  ungroup() %>%
  mutate(label2 = fct_reorder(factor(label2), -as.numeric(seurat_clusters))) 
go_ages_scRNAdat$GO = gsub(' ', '\n',go_ages_scRNAdat$GO)
go_ages_scExp = ggplot(go_ages_scRNAdat, aes(x = Age, y = label2)) +
  facet_wrap(~GO, ncol = 4) +
  geom_tile(aes(fill = meanOR)) +
  geom_text(aes(label = label), size = 6/pntnorm) +
  scale_fill_gradient(high = '#800026', low = "#FFFFCC") +
  theme_minimal(base_size = 6) +
  xlab(NULL) + ylab(NULL) +
  guides(fill = guide_colorbar(title = 'mean\nExp')) +
  theme(legend.key.width = unit(4,'pt'), legend.key.height = unit(8,'pt'), legend.background = element_rect(fill=NA, color = NA)) 

plotsave(go_ages_scExp, './results/scRNAseq/GO_in_scRNAseq/go_ages_scRNA_scExp', width = 8, height = 8)

#### 8.1.1 correlation between repair and replication using % expression ####

percExp_correlation_p = go_exp_dat %>%
  filter(expType == 'perc' & genes == 'allnonoverlap') %>%
  mutate(OR = value/totExpressed_allnonoverlap_perc) %>%
  select(-type,-value) %>%
  spread(GO,OR) %>%
  ggplot(aes(y = `DNA repair`, x = `DNA replication`)) +
  geom_point(size = 0.3 , alpha = 0.3) +
  geom_smooth(method = 'lm') +
  facet_wrap(~label2) +
  stat_cor(method='spearman', size = 6/pntnorm, cor.coef.name = 'rho')

#### 8.1.2. calcualte correlations using permutations - % expression ####

allgenes = intersect(unique(filter(myids, 
                                   human_ortholog %in% unique(unlist(gogenes)))$datagenes),
                     rownames(allint@assays$RNA))

allbg = setdiff(allgenes, goids$allgenes_in_scRNAseq[[3]])


# percExpressed_perm = sapply(1:1000,function(i){
#   print(i)
#   mat = as.matrix((allint@assays$RNA[sample(allbg,length(goids$no_genes_in_scRNAseq[[2]])),]>0)+1-1)
#   colMeans(mat)
# })
# 
# objsave(percExpressed_perm,'./data/processed/GO_in_scRNAseq/repair_percExp_perm')

percExpressed_perm = readRDS('./data/processed/GO_in_scRNAseq/repair_percExp_perm.rds')

cell_list = allint@meta.data %>%
  group_by(label2) %>%
  summarise(cells = list(unique(cell)))

cors = sapply(1:nrow(cell_list), function(i){
  cells=as.character(cell_list$cells[[i]])
  replication = percExpressed[cells,'DNA replication_allnonoverlap_perc']/percExpressed[cells,'totExpressed_allnonoverlap_perc']
  apply(percExpressed_perm, 2, function(x){
    x2 = x[cells]/percExpressed[cells,'totExpressed_allnonoverlap_perc']
    cor(replication,x2, method = 's')
  })
})
names(cors) = cell_list$label2

realcor = sapply(1:nrow(cell_list), function(i){
  cells=as.character(cell_list$cells[[i]])
  replication = percExpressed[cells,'DNA replication_allnonoverlap_perc']/percExpressed[cells,'totExpressed_allnonoverlap_perc']
  repair = percExpressed[cells,'DNA repair_allnonoverlap_perc']/percExpressed[cells,'totExpressed_allnonoverlap_perc']
  cor(replication, repair, method = 's')
})
names(realcor) = cell_list$label2

pvals = sapply(1:ncol(cors),function(i){
  ifelse(realcor[i]>0,mean(cors[,i]>=realcor[i]),mean(cors[,i]<=realcor[i]))
})

names(pvals) = cell_list$label2

repli_repair_corplot = data.frame(permmean =colMeans(cors), label2= cell_list$label2, real = realcor, p = pvals) %>%
  mutate(padj = p.adjust(p, method = 'BY')) %>%
  arrange(padj) %>%
  mutate(padj<=0.05) %>%
  ggplot(aes(y = reorder(label2,-real) , x = real, fill = padj<=0.05)) +
  geom_bar(stat = 'identity') +
  xlab('Spearman Correlation Coefficient') +
  ylab(NULL) +
  ggtitle('DNA Replication - DNA Repair') +
  theme_classic(base_size = 6) +
  scale_fill_brewer(palette = 3) +
  guides(fill = guide_legend(title = 'BY-corrected\np<=0.05')) +
  theme(legend.key.width = unit(4,'pt'), legend.key.height = unit(4,'pt'), legend.background = element_rect(fill=NA, color = NA),
        legend.position = c(0.005,0.005), legend.justification = c(0,0)) +
  xlim(-0.5,0.1)


## immune activation
## 
# percExpressed_perm_immune = sapply(1:1000,function(i){
#   print(i)
#   mat = as.matrix((allint@assays$RNA[sample(allbg,length(goids$no_genes_in_scRNAseq[[4]])),]>0)+1-1)
#   colMeans(mat)
# })
# objsave(percExpressed_perm_immune,'./data/processed/GO_in_scRNAseq/immune_percExp_perm')
percExpressed_perm_immune = readRDS('./data/processed/GO_in_scRNAseq/immune_percExp_perm.rds')


cors_immune = sapply(1:nrow(cell_list), function(i){
  cells=as.character(cell_list$cells[[i]])
  replication = percExpressed[cells,'DNA replication_allnonoverlap_perc']/percExpressed[cells,'totExpressed_allnonoverlap_perc']
  apply(percExpressed_perm_immune, 2, function(x){
    x2 = x[cells]/percExpressed[cells,'totExpressed_allnonoverlap_perc']
    cor(replication,x2, method = 's')
  })
})
names(cors_immune) = cell_list$label2

realcor_immune = sapply(1:nrow(cell_list), function(i){
  cells=as.character(cell_list$cells[[i]])
  replication = percExpressed[cells,'DNA replication_allnonoverlap_perc']/percExpressed[cells,'totExpressed_allnonoverlap_perc']
  repair = percExpressed[cells,'Immune cell activation_allnonoverlap_perc']/percExpressed[cells,'totExpressed_allnonoverlap_perc']
  cor(replication, repair, method = 's')
})
names(realcor_immune) = cell_list$label2

pvals_immune = sapply(1:ncol(cors_immune),function(i){
  ifelse(realcor_immune[i]>0,mean(cors_immune[,i]>=realcor_immune[i]),mean(cors_immune[,i]<=realcor_immune[i]))
})

names(pvals_immune) = cell_list$label2

repli_immune_corplot = data.frame(permmean =colMeans(cors_immune), label2= cell_list$label2, real = realcor_immune, p = pvals_immune) %>%
  mutate(padj = p.adjust(p, method = 'BY')) %>%
  arrange(padj) %>%
  mutate(padj<=0.05) %>%
  ggplot(aes(y = reorder(label2,-real) , x = real, fill = padj<=0.05)) +
  geom_bar(stat = 'identity') +
  xlab('Spearman Correlation Coefficient') +
  ylab(NULL) +
  ggtitle('DNA Replication - Immune Cell Activation') +
  theme_classic(base_size = 6) +
  scale_fill_brewer(palette = 3) +
  guides(fill = guide_legend(title = 'BY-corrected\np<=0.05')) +
  theme(legend.key.width = unit(4,'pt'), legend.key.height = unit(4,'pt'), legend.background = element_rect(fill=NA, color = NA),
        legend.position = c(0.005,0.005), legend.justification = c(0,0))+
  xlim(-0.5,0.1)

plotsave(repli_immune_corplot, './results/scRNAseq/GO_in_scRNAseq/repli_immune_corplot_percExp',width = 8, height = 4.5)
plotsave(repli_repair_corplot, './results/scRNAseq/GO_in_scRNAseq/repli_repair_corplot_percExp',width = 8, height = 4.5)


correlplots = ggarrange(repli_repair_corplot,repli_immune_corplot+guides(fill='none'), ncol = 1, nrow = 2)

plotsave(correlplots, './results/scRNAseq/GO_in_scRNAseq/correlplots_percExp',width = 8, height = 9)

# 10.2.1. correlation between repair and replication using mean scaled expression 

scExp_correlation_p = go_exp_dat %>%
  filter(expType == 'scExp' & genes == 'allnonoverlap') %>%
  # mutate(OR = value/totExpressed_allnonoverlap_perc) %>%
  select(-type) %>%
  spread(GO,value) %>%
  ggplot(aes(y = `DNA repair`, x = `DNA replication`)) +
  geom_point(size = 0.3 , alpha = 0.3) +
  geom_smooth(method = 'lm') +
  facet_wrap(~label2) +
  stat_cor(method='spearman', size = 6/pntnorm, cor.coef.name = 'rho')

# 10.2.2. calcualte correlations using permutations - mean scaled expression 

# scExpressed_perm = sapply(1:1000,function(i){
#   print(i)
#   mat = as.matrix(allint@assays$RNA[sample(allbg,length(goids$no_genes_in_scRNAseq[[2]])),])
#   colMeans(t(scale(t(mat))))
# })
# 
# objsave(scExpressed_perm,'./data/processed/GO_in_scRNAseq/repair_scExp_perm')

scExpressed_perm = readRDS('./data/processed/GO_in_scRNAseq/repair_scExp_perm.rds')

sccors = sapply(1:nrow(cell_list), function(i){
  cells=as.character(cell_list$cells[[i]])
  replication = scExp[cells,'DNA replication_allnonoverlap_scExp']
  apply(scExpressed_perm, 2, function(x){
    cor(replication,x[cells], method = 's')
  })
})
names(sccors) = cell_list$label2

screalcor = sapply(1:nrow(cell_list), function(i){
  cells=as.character(cell_list$cells[[i]])
  replication = scExp[cells,'DNA replication_allnonoverlap_scExp']
  repair = scExp[cells,'DNA repair_allnonoverlap_scExp']
  cor(replication, repair, method = 's')
})
names(screalcor) = cell_list$label2

scpvals = sapply(1:ncol(sccors),function(i){
  ifelse(screalcor[i]>0,mean(sccors[,i]>=screalcor[i]),mean(sccors[,i]<=screalcor[i]))
})

names(scpvals) = cell_list$label2

sc_repli_repair_corplot = data.frame(permmean =colMeans(sccors), label2= cell_list$label2, real = screalcor, p = scpvals) %>%
  mutate(padj = p.adjust(p, method = 'BY')) %>%
  arrange(padj) %>%
  mutate(padj<=0.05) %>%
  ggplot(aes(y = reorder(label2,-real) , x = real, fill = padj<=0.05)) +
  geom_bar(stat = 'identity') +
  xlab('Spearman Correlation Coefficient') +
  ylab(NULL) +
  ggtitle('DNA Replication - DNA Repair') +
  theme_classic(base_size = 6) +
  scale_fill_brewer(palette = 3) +
  guides(fill = guide_legend(title = 'BY-corrected\np<=0.05')) +
  theme(legend.key.width = unit(4,'pt'), legend.key.height = unit(4,'pt'), legend.background = element_rect(fill=NA, color = NA),
        legend.position = c(0.995,0.995), legend.justification = c(1,1)) +
  xlim(0,1)
sc_repli_repair_corplot
### immune 
# scExpressed_perm_immune = sapply(1:1000,function(i){
#   print(i)
#   mat = as.matrix(allint@assays$RNA[sample(allbg,length(goids$no_genes_in_scRNAseq[[4]])),])
#   colMeans(t(scale(t(mat))))
# })
# objsave(scExpressed_perm_immune,'./data/processed/GO_in_scRNAseq/immune_scExp_perm')

scExpressed_perm_immune = readRDS('./data/processed/GO_in_scRNAseq/immune_scExp_perm.rds')

sccors_immune = sapply(1:nrow(cell_list), function(i){
  cells=as.character(cell_list$cells[[i]])
  replication = scExp[cells,'DNA replication_allnonoverlap_scExp']
  apply(scExpressed_perm_immune, 2, function(x){
    cor(replication,x[cells], method = 's')
  })
})
names(sccors_immune) = cell_list$label2

screalcor_immune = sapply(1:nrow(cell_list), function(i){
  cells=as.character(cell_list$cells[[i]])
  replication = scExp[cells,'DNA replication_allnonoverlap_scExp']
  repair = scExp[cells,'Immune cell activation_allnonoverlap_scExp']
  cor(replication, repair, method = 's')
})
names(screalcor_immune) = cell_list$label2

scpvals_immune = sapply(1:ncol(sccors_immune),function(i){
  ifelse(screalcor_immune[i]>0,mean(sccors_immune[,i]>=screalcor_immune[i]),mean(sccors_immune[,i]<=screalcor_immune[i]))
})

names(scpvals_immune) = cell_list$label2

sc_repli_immune_corplot = data.frame(permmean =colMeans(sccors_immune), label2= cell_list$label2, real = screalcor_immune, p = scpvals_immune) %>%
  mutate(padj = p.adjust(p, method = 'BY')) %>%
  arrange(padj) %>%
  mutate(padj<=0.05) %>%
  ggplot(aes(y = reorder(label2,-real) , x = real, fill = padj<=0.05)) +
  geom_bar(stat = 'identity') +
  xlab('Spearman Correlation Coefficient') +
  ylab(NULL) +
  ggtitle('DNA Replication - Immune cell activation') +
  theme_classic(base_size = 6) +
  scale_fill_brewer(palette = 3) +
  guides(fill = guide_legend(title = 'BY-corrected\np<=0.05')) +
  theme(legend.key.width = unit(4,'pt'), legend.key.height = unit(4,'pt'), legend.background = element_rect(fill=NA, color = NA),
        legend.position = c(0.995,0.995), legend.justification = c(1,1)) +
  xlim(0,1)

plotsave(sc_repli_immune_corplot, './results/scRNAseq/GO_in_scRNAseq/repli_immune_corplot_scExp',width = 8, height = 4.5)
plotsave(sc_repli_repair_corplot, './results/scRNAseq/GO_in_scRNAseq/repli_repair_corplot_scExp',width = 8, height = 4.5)
correlplots_sc = ggarrange(sc_repli_repair_corplot,sc_repli_immune_corplot+guides(fill='none'), ncol = 1, nrow = 2)

plotsave(correlplots_sc, './results/scRNAseq/GO_in_scRNAseq/correlplots_scExp',width = 8, height = 9)

