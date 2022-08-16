source('./scripts/00-setup.R')

gores = readRDS('./data/processed/KMproteomics/DE_GOenrich_allsamples.rds')

wiltest = readRDS('./data/processed/KMproteomics/DEwilcox_allsamples.rds')

goreps = c('DNA repair', 'double-strand break repair via homologous recombination', 'response to UV')

repcats = gores %>%
  filter(representative_description %in% goreps)

gogenes = readRDS('./data/processed/GO2Gene/GO2GeneBP_20220407.rds')

gogenes = reshape2::melt(gogenes[repcats$ID]) %>%
  set_names(c('human_ensembl','ID'))

humanids = readRDS('./data/processed/helperdata/humanids.rds')
humanids2 = rename(humanids$geneIDs, human_ensembl = ensembl_gene_id) %>% select(1,2,6) %>% set_names(c('humangene','human_ensembl','human_description'))
humanids2$human_description = sapply(strsplit(humanids2$human_description,' [[]Source:'),function(x)x[1])

nfuids = readRDS('./data/processed/helperdata/nfuids.rds')

nfuids = nfuids$geneIDs %>%
  mutate(datagenes = ifelse(external_gene_name=='',ensembl_gene_id,external_gene_name))

nfuids$description = sapply(strsplit(nfuids$description,' [[]Source:'),function(x)x[1])

humanort = readRDS('./data/processed/helperdata/humanorthologs.rds')

human = humanort$orthologs %>%
  dplyr::select(2,'hsapiens_homolog_ensembl_gene') %>%
  set_names('nfu_gene','ensembl_gene_id') %>%
  left_join(humanids$geneIDs) %>%
  dplyr::select(1,2) %>%
  set_names('ensembl_gene_id','human_ensembl')

human = left_join(human,humanids2)  %>%
  select(1,2,3) %>%
  unique()

myids = left_join(nfuids,human)

rm(human,humanids,nfuids,humanort,humanids2)


# collect all genes in GO categories represented by the three categories of interest

allgenesinGO = select(myids,external_gene_name,human_ensembl) %>%
  inner_join(gogenes) %>%
  left_join(unique(select(repcats,ID,Description,representative_description))) %>%
  mutate(upperprotein = toupper(external_gene_name)) %>%
  left_join(wiltest) %>%
  filter(!is.na(l2FC))

# go-gene association summarisation

gogenemat = allgenesinGO %>% filter(abs(l2FC)>=1) %>% select(upperprotein, Description) %>% unique() %>% mutate(val = 1) %>% spread(Description, val, fill = 0) %>% as.data.frame()
rownames(gogenemat) = gogenemat$upperprotein
gogenemat$upperprotein = NULL
gogenemat = as.matrix(gogenemat)
dnagooverlaps = apply(gogenemat,2,function(x){
  apply(gogenemat,2,function(y){
    sum(x & y) / sum( x | y)
  })
})

library(corrplot)
order.hc = corrMatOrder(dnagooverlaps, order = 'hclust',hclust.method = 'average')
clustersx = cutree(hclust(dist(dnagooverlaps)),5)

pdf('./results/msfigures/kmproteomics_repairgocat_jaccard.pdf', useDingbats = T, width = 5, height = 5)
corrplot(dnagooverlaps, method = 'square',
         diag = T, addCoef.col = 'gray25', 
         order = 'hclust', addrect = 5, 
         hclust.method = 'average',
         type = 'lower', cl.cex = 0.5,
         bg = NA, tl.cex = 0.7, number.cex = 0.3, tl.pos = 'l',
         tl.col = ggthemes_data$gdocs$colors$value[clustersx[colnames(dnagooverlaps)[order.hc]]], 
         number.digits = 1, col.lim = c(0,1), 
         col = colorRampPalette(c('blue','white','#b82e2e'))(200))
dev.off()

pcax = prcomp(t(dnagooverlaps))
axesnms = paste(colnames(pcax$x),' (', round(summary(pcax)$imp[2,]*100,2),'%)', sep='')
prx = pcax$x
pcaplot = as.data.frame(prx) %>%
  mutate(gocats = colnames(dnagooverlaps)) %>%
  mutate(cols = ggthemes_data$gdocs$colors$value[1:5][clustersx[gocats]]) %>%
  ggplot(aes(x = PC1, y = PC2, color = cols)) +
  geom_point(size = 3) +
  geom_text_repel(size = 6/pntnorm, aes(label = gocats)) +
  scale_color_identity() +
  theme_classic(base_size = 8) +
  ggforce::geom_mark_circle(aes(fill = cols), col = NA, alpha = 0.05) +
  scale_fill_identity()  +
  xlab(axesnms[1]) + ylab(axesnms[2]) 
pcaplot
plotsave(pcaplot, './results/msfigures/kmproteomics_repairgocat_overlap', width = 8, height = 8)
