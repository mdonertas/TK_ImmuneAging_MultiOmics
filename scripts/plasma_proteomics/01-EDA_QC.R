source('./scripts/00-setup.R')

dat = readxl::read_xlsx('../../raw_data/multiomics/TK_KidneyMarrowAndPlasma_proteomics_metabolomics/plasma_proteomics.xlsx',
                        col_names = c('Protein','Age',paste('s',1:5, sep = '')), skip = 1) %>%
  na.omit() 

dat = dat %>%
  gather('sampleID','level',-Protein,-Age) %>%
  mutate(fullID = paste(Age,sampleID,sep='_')) %>%
  group_by(Protein,fullID) %>%
  summarise(level = median(level)) %>%
  spread(fullID, level, fill=NA)

mat = as.matrix(dat[,-1]) 
rownames(mat) = dat$Protein
mat_l2 = log2(mat)

dim(mat_l2)
# 474 proteins, 10 samples


allsamples_boxp = reshape2::melt(mat_l2) %>%
  ggplot(aes(x = Var2, y = value)) +
  geom_boxplot(outlier.size = 0.2)  +
  xlab(NULL) + ylab('Log2 Expression')

plotsave(allsamples_boxp, './results/plasmaProteomics/allsamples/exp_dist_box', width = 12, height = 6)


# correlations
pheatmap::pheatmap(cor(mat_l2), cellwidth = 10, cellheight = 10, color = RColorBrewer::brewer.pal(8,'Oranges')[-1], filename = './results/plasmaProteomics/allsamples/correlationplot.pdf')
pheatmap::pheatmap(cor(mat_l2), cellwidth = 10, cellheight = 10, color = RColorBrewer::brewer.pal(8,'Oranges')[-1], filename = './results/plasmaProteomics/allsamples/correlationplot.png')
#old sample 3 is distant to all but still have a correlation of 0.86

pcx = prcomp(t(mat_l2), scale. = T)

imp = summary(pcx)$imp[2,]*100
imptext = paste('PC', 1:10,' (',round(imp,2),'%)' , sep='')

pcaplot = data.frame(pcx$x ) %>%
  mutate(fullID = rownames(pcx$x)) %>%
  separate(fullID,into=c('Age','sampleID'),remove = F) %>%
  ggplot(aes(x = PC1, y = PC2, color = Age) ) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'gray35') +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'gray35') +
  geom_point(size = 3) +
  # geom_text_repel(aes(label = fullID), color = 'black', size = 7/pntnorm) +
  scale_color_manual(values = agecol) +
  xlab(imptext[1]) +
  ylab(imptext[2]) +
  # theme_rvis(base_size = 6, legend.pos = c(0,1), just.par = c(0,1)) +
  # coord_fixed(ratio = imp[2]/imp[1]) +
  theme_minimal() 
pcaplot

plotsave(pcaplot, './results/plasmaProteomics/allsamples/pca', width = 8, height = 8)

#old s3 is an outlier clearly, why? check the proteins in second dim

humanids = readRDS('./data/processed/helperdata/humanids.rds')
humanids2 = rename(humanids$geneIDs, human_ortholog = ensembl_gene_id) %>% select(1,2,6) %>% set_names(c('humangene','human_ortholog','human_description'))
humanids2$human_description = sapply(strsplit(humanids2$human_description,' [[]Source:'),function(x)x[1])

nfuids = readRDS('./data/processed/helperdata/nfuids.rds')

nfuids = nfuids$geneIDs %>%
  mutate(datagenes = ifelse(external_gene_name=='',ensembl_gene_id,external_gene_name))

nfuids$description = sapply(strsplit(nfuids$description,' [[]Source:'),function(x)x[1])

humanort = readRDS('./data/processed/helperdata/humanorthologs.rds')
human = humanort$orthologs %>%
  dplyr::select(2,'hsapiens_homolog_ensembl_gene') %>%
  set_names('nfu_gene','human_ortholog') %>%
  left_join(humanids2) %>%
  dplyr::select(1,3) %>%
  set_names('ensembl_gene_id','human_ortholog')

myids = left_join(nfuids,human)

myids2 = myids %>%
  filter(external_gene_name %in% rownames(pcx$rotation)) %>%
  filter(!(human_ortholog =='' | is.na(human_ortholog)) ) %>%
  select(1,human_ortholog) %>%
  unique() 

humanpc2 = data.frame(pc2 = pcx$rotation[,2], external_gene_name = rownames(pcx$rotation)) %>%
  left_join(myids2) %>%
  dplyr::select(1,3) %>%
  na.omit() %>%
  unique()

exclude = which(table(humanpc2$human_ortholog)>1)

humanpc2 = humanpc2 %>%
  filter(!human_ortholog%in%exclude)

humanpc2 = sort(setNames(humanpc2$pc2, humanpc2$human_ortholog), dec=T)

gorespc2 = clusterProfiler::gseGO(humanpc2, OrgDb = 'org.Hs.eg.db',keyType = 'SYMBOL', pvalueCutoff = 1, minGSSize = 10, maxGSSize = 500)

gores = gorespc2@result %>%
  mutate(BY_corrected_p = p.adjust(pvalue, method = 'BY'))  %>%
  mutate(FDR_corrected_p = p.adjust(pvalue, method = 'fdr'))  %>%
  dplyr::select(1:6,FDR_corrected_p, BY_corrected_p)

gores %>%
  tablesave('./results/plasmaProteomics/allsamples/outlierPC2_GOres')

objsave(mat, './data/processed/plasmaProteomics/expression')
objsave(mat_l2, './data/processed/plasmaProteomics/expression_log2')

#### repeat the same after excluding the outlier. 

mat = mat[,!colnames(mat)%in%c('old_s3')]
mat_l2 = mat_l2[,!colnames(mat_l2)%in%c('old_s3')]

outremoved_boxp = reshape2::melt(mat_l2) %>%
  ggplot(aes(x = Var2, y = value)) +
  geom_boxplot(outlier.size = 0.2)  +
  xlab(NULL) + ylab('Log2 Expression')

plotsave(outremoved_boxp, './results/plasmaProteomics/outlier_removed/exp_dist_box', width = 12, height = 6)


# correlations
pheatmap::pheatmap(cor(mat_l2), cellwidth = 10, cellheight = 10, color = RColorBrewer::brewer.pal(8,'Oranges')[-1], filename = './results/plasmaProteomics/outlier_removed/correlationplot.pdf')
pheatmap::pheatmap(cor(mat_l2), cellwidth = 10, cellheight = 10, color = RColorBrewer::brewer.pal(8,'Oranges')[-1], filename = './results/plasmaProteomics/outlier_removed/correlationplot.png')
#old sample 3 is distant to all but still have a correlation of 0.86

pcx = prcomp(t(mat_l2), scale. = T)

imp = summary(pcx)$imp[2,]*100
imptext = paste('PC', 1:10,' (',round(imp,2),'%)' , sep='')

pcaplot = data.frame(pcx$x ) %>%
  mutate(fullID = rownames(pcx$x)) %>%
  separate(fullID,into=c('Age','sampleID'),remove = F) %>%
  ggplot(aes(x = PC1, y = PC2, color = Age)) +
  geom_point(size = 3) +
  # geom_text_repel(aes(label = fullID), color = 'black', size = 7/pntnorm) +
  scale_color_manual(values = agecol) +
  xlab(imptext[1]) +
  ylab(imptext[2]) +
  theme_rvis(base_size = 6, legend.pos = c(0,0), just.par = c(0,0)) +
  coord_fixed(ratio = imp[2]/imp[1])
pcaplot

plotsave(pcaplot, './results/plasmaProteomics/outlier_removed/pca', width = 8, height = 8)


pcaplot = data.frame(pcx$x ) %>%
  mutate(fullID = rownames(pcx$x)) %>%
  separate(fullID,into=c('Age','sampleID'),remove = F) %>%
  ggplot(aes(x = PC1, y = PC2, color = Age)) +
  geom_point(size = 3) +
  # geom_text_repel(aes(label = fullID), color = 'black', size = 7/pntnorm) +
  scale_color_manual(values = agecol) +
  xlab(imptext[1]) +
  ylab(imptext[2]) +
  theme_rvis(base_size = 6, legend.pos = c(0,0), just.par = c(0,0)) 
pcaplot

plotsave(pcaplot, './results/plasmaProteomics/outlier_removed/pca_noratio', width = 8, height = 8)

objsave(mat, './data/processed/plasmaProteomics/expression_outlierremoved')
objsave(mat_l2, './data/processed/plasmaProteomics/expression_log2_outlierremoved')
