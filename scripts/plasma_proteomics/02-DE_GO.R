source('./scripts/00-setup.R')

humanids = readRDS('./data/processed/helperdata/humanids.rds')
humanids2 = rename(humanids$geneIDs, human_ortholog = ensembl_gene_id) %>% select(1, 2, 6) %>% set_names(c('humangene', 'human_ortholog', 'human_description'))
humanids2$human_description = sapply(strsplit(humanids2$human_description, ' [[]Source:'), function(x)
  x[1])

nfuids = readRDS('./data/processed/helperdata/nfuids.rds')

nfuids = nfuids$geneIDs %>%
  mutate(datagenes = ifelse(external_gene_name == '', ensembl_gene_id, external_gene_name))

nfuids$description = sapply(strsplit(nfuids$description, ' [[]Source:'), function(x)
  x[1])

humanort = readRDS('./data/processed/helperdata/humanorthologs.rds')

human = humanort$orthologs %>%
  dplyr::select(2, 'hsapiens_homolog_ensembl_gene') %>%
  set_names('nfu_gene', 'ensembl_gene_id') %>%
  left_join(humanids$geneIDs) %>%
  dplyr::select(1, 2) %>%
  set_names('ensembl_gene_id', 'human_ortholog')

human = left_join(human, humanids2)  %>%
  select(1, 3) %>%
  unique()

myids = left_join(nfuids, human)

rm(human, humanids, nfuids, humanort, humanids2)

mat_l2 = readRDS('./data/processed/plasmaProteomics/expression_log2.rds')

wiltest = data.frame(t(apply(mat_l2, 1, function(x) {
  y = x[grepl('young', names(x))]
  o = x[grepl('old', names(x))]
  l2fc = median(o) - median(y)
  difs = c(sapply(o, function(xx) {
    sapply(y, function(z)
      xx - z)
  }))
  quants = quantile(difs, c(0.25, 0.75))
  wi = wilcox.test(y, o)
  c(wi$stat, wi$p.val, l2fc, quants)
}))) %>%
  set_names(c('w', 'p', 'l2FC', 'Q1', 'Q3')) %>%
  mutate(padj = p.adjust(p, method = 'fdr'))

wiltest$protein = rownames(wiltest)
wiltest$upperprotein = toupper(wiltest$protein)

wiltest %>%
  select(protein, 3, 4, 5, 2, 6) %>%
  tablesave('./results/plasmaProteomics/allsamples/proteinDEtest')

objsave(wiltest,
        './data/processed/plasmaProteomics/DEwilcox_allsamples')

myids2 = myids %>%
  mutate(upperprotein = toupper(external_gene_name)) %>%
  select(upperprotein, humangene) %>%
  na.omit() %>%
  filter(upperprotein != '' & humangene != '')

wiltest2 = wiltest %>%
  left_join(myids2)

wiltesthuman = wiltest2 %>%
  select(l2FC, humangene) %>%
  unique() %>%
  na.omit() %>%
  group_by(humangene) %>%
  summarise(l2fc = median(l2FC))

genelist = sort(setNames(wiltesthuman$l2fc, wiltesthuman$humangene), dec =T)

goenrich = clusterProfiler::gseGO(
  genelist,
  OrgDb = 'org.Hs.eg.db',
  keyType = 'SYMBOL',
  ont = 'BP',
  pvalueCutoff = 1,
  minGSSize = 10,
  maxGSSize = 500
)

gores = goenrich@result %>%
  mutate(BY_corrected_p = p.adjust(pvalue, method = 'BY'))  %>%
  mutate(FDR_corrected_p = p.adjust(pvalue, method = 'fdr'))
detach('package:org.Hs.eg.db')
detach('package:AnnotationDbi')

genelist = filter(gores, BY_corrected_p <= 0.05) %>% select(ID, core_enrichment)

genelist = setNames(strsplit(genelist$core_enrichment, '/'), genelist$ID)

incmat = sapply(genelist, function(x) {
  sapply(genelist, function(y) {
    length(unique(intersect(x, y))) / length(unique(union(x, y)))
  })
})

clmeds = (sapply(20:70, function(k) {
  cls = cutree(hclust(as.dist(1 - incmat)), k)
  mean(sapply(unique(cls), function(i) {
    mycl = names(which(cls == i))
    if (length(mycl) > 1) {
      mymat = incmat[mycl, mycl]
      median(mymat[upper.tri(mymat)])
    } else{
      1
    }
  }) >= 0.5)
}))

setSize = setNames(gores$setSize, gores$ID)
k = which(clmeds >= 0.5)[1] + 19
cls = cutree(hclust(as.dist(1 - incmat)), k)

reps = sapply(1:k, function(i) {
  genes = names(which(cls == i))
  if (length(genes) > 1) {
    selectgo = names(which(rowMeans(incmat[genes, genes]) == max(rowMeans(incmat[genes, genes]))))
    setNames(rep(names(which.max(setSize[selectgo])), length(genes)), genes)
  } else{
    setNames(genes, genes)
  }
})

reps = unlist(reps)

idnames = setNames(gores$Description, gores$ID)

gores = gores %>%
  # filter(ID %in% reps) %>%
  # filter(BY_corrected_p <= 0.05) %>%
  mutate(representative = reps[as.character(ID)]) %>%
  mutate(representative_description = idnames[as.character(representative)])

goresplot = gores %>%
  filter(BY_corrected_p <= 0.05) %>%
  filter(representative == ID) %>%
  ggplot(aes(
    x = NES,
    y = reorder(Description, NES),
    fill = NES > 0
  )) +
  geom_vline(xintercept = 0, color = 'gray') +
  geom_bar(stat = 'identity') +
  ylab(NULL) +
  xlab("Normalized enrichment score") +
  scale_fill_manual(values = unname(agecol)) +
  guides(fill = 'none') +
  annotate(
    'text',
    x = 0.1,
    y = 21.5,
    label = 'Enriched in the old',
    hjust = 0,
    vjust = 1.4,
    size = 6 / pntnorm
  ) +
  annotate(
    'text',
    x = -0.1,
    y = 21.5,
    label = 'Enriched in the young',
    hjust = 1,
    vjust = 1.4,
    size = 6 / pntnorm
  ) +
  coord_cartesian(clip = 'off')


plotsave(
  goresplot,
  './results/plasmaProteomics/allsamples/GOreps',
  width = 10,
  height = 8
)


objsave(gores,
        './data/processed/plasmaProteomics/DE_GOenrich_allsamples')

gores %>%
  select(-7, -8) %>%
  select(1:6, 10:13, 7:9) %>%
  tablesave('./results/plasmaProteomics/allsamples/GOenrichment')

repgenelist = sapply(unique(gores$representative_description), function(repdesc) {
  genes = filter(gores, representative_description == repdesc)$core_enrichment
  unique(unlist(strsplit(genes, '/')))
})
repgenelist = repgenelist[!sapply(repgenelist,is.null)]
repgenedat = reshape2::melt(repgenelist) %>%
  set_names('humangene', 'GO representative') %>%
  left_join(myids2) %>%
  left_join(wiltest2) %>%
  select(2, 6, 7, 8, 9, 10) %>%
  na.omit() %>%
  unique()
lim = max(c(abs(repgenedat$Q1), abs(repgenedat$Q3)))

gogeneplots = lapply(unique(repgenedat$`GO representative`), function(gorep) {
  repgenedat %>%
    filter(`GO representative` == gorep) %>%
    filter(abs(l2FC) >= 1) %>%
    mutate(type = c('young', 'old')[1 + (l2FC > 0)]) %>%
    mutate(protein = fct_reorder(protein, l2FC)) %>%
    ggplot(aes(x = l2FC, y = protein, fill = type)) +
    geom_bar(stat = 'identity') +
    geom_errorbar(
      aes(xmin = Q1, xmax = Q3),
      width = 0.25,
      color = 'gray35',
      size = 0.25
    ) +
    scale_fill_manual(values = agecol) +
    # xlim(-lim,lim) +
    ggtitle(gorep) +
    ylab(NULL) +
    xlab('Log2 Fold Change') +
    guides(fill = 'none') +
    geom_text(
      data = . %>% filter(padj <= 0.1) %>% mutate(labx = ifelse(abs(Q3) > abs(Q1), Q3 +
                                                                  0.1, Q1 - 0.1)),
      aes(x = labx),
      label = '*',
      vjust = 0.8,
      size = 8 / pntnorm,
      color = 'gray35'
    )
})
names(gogeneplots) = unique(repgenedat$`GO representative`)

sapply(names(gogeneplots), function(gorep) {
  fname = paste('./results/plasmaProteomics/allsamples/GOgenes/',
                gsub(' ', '_', gorep),
                sep = '')
  plotsave(gogeneplots[[gorep]],
           fname,
           width = 8,
           height = max(4, length(unique(
             gogeneplots[[gorep]]$data$protein
           )) / 3))
})


drivergene_plot = repgenedat %>%
  # filter(`GO representative` == gorep) %>%
  select(-1) %>%
  unique() %>%
  filter(abs(l2FC) >= 1) %>%
  mutate(type = c('young', 'old')[1 + (l2FC > 0)]) %>%
  mutate(protein = fct_reorder(protein, l2FC)) %>%
  ggplot(aes(x = l2FC, y = protein, fill = type)) +
  geom_bar(stat = 'identity') +
  # geom_errorbar(aes(xmin = Q1, xmax = Q3), width = 0.25, color = 'gray35', size = 0.25) +
  scale_fill_manual(values = agecol) +
  ggtitle('Driver genes') +
  ylab(NULL) +
  xlab('Log2 Fold Change') +
  guides(fill = 'none') +
  # xlim(-4.5,4.5) +
  # geom_text(aes(label = protein), x= 0, fill = 'white', size = 6/pntnorm) +
  geom_text(
    data = . %>% filter(padj <= 0.1) %>% mutate(labx = ifelse(abs(Q3) > abs(Q1), l2FC +
                                                                0.1, l2FC - 0.1)),
    aes(x = labx),
    label = '*',
    vjust = 0.8,
    size = 8 / pntnorm,
    color = 'gray35'
  )

plotsave(
  drivergene_plot,
  './results/plasmaProteomics/allsamples/drivergenes',
  width = 8,
  height = 10
)



####
#check if the use of human orthologs may bias go enrichment

source('./scripts/00-setup.R')

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
  set_names('nfu_gene','ensembl_gene_id') %>%
  left_join(humanids$geneIDs) %>%
  dplyr::select(1,2) %>%
  set_names('ensembl_gene_id','human_ortholog')

human = left_join(human,humanids2)  %>%
  select(1,3) %>%
  unique()

myids = left_join(nfuids,human)

rm(human,humanids,nfuids,humanort,humanids2)

wiltest = readRDS('./data/processed/plasmaProteomics/DEwilcox_allsamples.rds')

myids2 = myids %>%
  mutate(upperprotein = toupper(external_gene_name)) %>%
  select(upperprotein, humangene) %>%
  na.omit() %>%
  filter(upperprotein!='' & humangene!='')

wiltest2 = wiltest %>%
  left_join(myids2) 

wiltesthuman = wiltest2 %>%
  group_by(protein,l2FC) %>%
  summarise(humangene = sample(humangene,1)) %>%
  ungroup() %>%
  dplyr::select(l2FC, humangene) %>%
  unique() %>%
  na.omit() %>%
  group_by(humangene) %>%
  summarise(l2fc = median(l2FC))

genelist = sort(setNames(wiltesthuman$l2fc, wiltesthuman$humangene),dec=T)

goenrich = clusterProfiler::gseGO(genelist, OrgDb = 'org.Hs.eg.db', keyType = 'SYMBOL', ont = 'BP', pvalueCutoff = 1, minGSSize = 10, maxGSSize = 500)

gores = goenrich@result %>%
  mutate(BY_corrected_p = p.adjust(pvalue, method = 'BY'))  %>%
  mutate(FDR_corrected_p = p.adjust(pvalue, method = 'fdr'))  

allortgores = readRDS('./data/processed/plasmaProteomics/DE_GOenrich_allsamples.rds') %>%
  # filter(BY_corrected_p<=0.05) %>%
  dplyr::select(ID,NES,BY_corrected_p) %>%
  dplyr::rename(allortNES = NES, allortp=BY_corrected_p)

signifgos = inner_join(gores,allortgores) %>%
  filter(BY_corrected_p<=0.05 | allortp<=0.05)

cor.test(signifgos$NES,signifgos$allortNES,method='spearman')
# Spearman's rank correlation rho
# 
# data:  signifgos$NES and signifgos$allortNES
# S = 2548.5, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.9513583 
