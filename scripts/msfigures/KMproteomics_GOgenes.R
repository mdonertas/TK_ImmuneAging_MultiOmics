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
gores = readRDS('./data/processed/KMproteomics/DE_GOenrich_allsamples.rds')

wiltest = readRDS('./data/processed/KMproteomics/DEwilcox_allsamples.rds')

myids2 = myids %>%
  mutate(upperprotein = toupper(external_gene_name)) %>%
  select(upperprotein, humangene) %>%
  na.omit() %>%
  filter(upperprotein!='' & humangene!='')

wiltest2 = wiltest %>%
  left_join(myids2) 

repgenelist = sapply(unique(gores$representative_description), function(repdesc){
  genes = filter(gores, representative_description == repdesc)$core_enrichment
  unique(unlist(strsplit(genes, '/')))
})
repgenelist = repgenelist[!sapply(repgenelist,is.null)]
repgenedat = reshape2::melt(repgenelist) %>%
  set_names(c('humangene','GO representative')) %>%
  left_join(myids2) %>%
  left_join(wiltest2) %>%
  select(2,6,7,8,9,10) %>%
  na.omit() %>%
  unique()

repmat = repgenedat %>%
  filter(abs(l2FC)>=1.5) %>%
  select(1,2,protein) %>%
  spread(protein,l2FC,fill=0) %>% 
  as.data.frame()

rownames(repmat) = repmat$`GO representative`
repmat$`GO representative` = NULL
repmat = as.matrix(repmat)
repmat = t((repmat!=0)+1-1)

toadd = setdiff(unique(repgenedat$`GO representative`), colnames(repmat))
toadd = matrix(data = 0, nrow = nrow(repmat), ncol = length(toadd), dimnames = list(rownames(repmat),toadd))
repmat = cbind(repmat,toadd)

goresx = filter(gores, Description %in% colnames(repmat))
goresx = setNames(goresx$NES, goresx$Description)[colnames(repmat)]
generesx = repgenedat %>%
  filter(abs(l2FC)>=1.5) %>%
  select(protein, l2FC) %>% unique()
generesx = setNames(generesx$l2FC, generesx$protein)[rownames(repmat)]
generessignif = repgenedat %>%
  # filter((l2FC)<=-0.5) %>%
  filter(abs(l2FC)>=1.5) %>%
  select(protein, padj) %>% unique()
generessignif = 1+(setNames(generessignif$padj<=0.1, generessignif$protein)[rownames(repmat)])


library(ComplexHeatmap)
nescols = agecol[(generesx>0)+1]
gocols = agecol[(goresx>0)+1]
basesize = 6
column_ha = HeatmapAnnotation(`Normalized\nEnrichment\nScore` = anno_barplot(abs(goresx), height = unit(1,'cm'), gp = gpar(fill = gocols)), annotation_name_gp = gpar(fontsize = basesize+1))
row_ha = rowAnnotation(Log2FC = anno_barplot(abs(generesx), width = unit(1,'cm'), gp = gpar(fill = nescols)), annotation_name_gp = gpar(fontsize = basesize+1))
inchconv = 0.393701
pdf(file="./results/msfigures/KMproteomics_GOgenesl2fc15.pdf", width = 16*inchconv, height = 17*inchconv, useDingbats = F)
Heatmap(repmat, top_annotation = column_ha, right_annotation = row_ha, col = c('gray95','gray25'), 
        rect_gp = gpar(col = "white", lwd = 1), 
        column_names_max_height = max_text_width(
          colnames(repmat), 
          gp = gpar(fontsize = basesize)),
        row_names_gp = gpar(fontsize = basesize, fontface = generessignif), column_names_gp = gpar(fontsize = basesize), 
        # row_km = 3, 
        # column_km = 4, 
        row_title = NULL, 
        row_gap = unit(basesize, "pt"), column_gap = unit(basesize, "pt"), border = 'gray35',
        width = ncol(repmat)*unit(basesize, "pt"),
        row_order = order(-generesx),
        column_order = order(goresx),
        column_title_gp = gpar(fontsize = basesize + 2),
        height = nrow(repmat)*unit(basesize, "pt"), show_heatmap_legend = F, column_title = 'Kidney Marrow Proteomics\nCore Enrichment Genes (|Log2FC| >= 1.5)')
dev.off()


####

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

allgenesinGO = select(myids,external_gene_name,human_ensembl) %>%
  inner_join(gogenes) %>%
  left_join(unique(select(repcats,ID,Description,representative_description))) %>%
  mutate(upperprotein = toupper(external_gene_name)) %>%
  left_join(wiltest) %>%
  filter(!is.na(l2FC))

repmat = allgenesinGO %>%
  filter(abs(l2FC)>=1) %>%
  select(representative_description,l2FC,protein) %>%
  unique() %>%
  spread(protein,l2FC,fill=0) %>% 
  as.data.frame()

rownames(repmat) = repmat$representative_description
repmat$representative_description = NULL
repmat = as.matrix(repmat)
repmat = t((repmat!=0)+1-1)

toadd = setdiff(unique(allgenesinGO$representative_description), colnames(repmat))
toadd = matrix(data = 0, nrow = nrow(repmat), ncol = length(toadd), dimnames = list(rownames(repmat),toadd))
repmat = cbind(repmat,toadd)

goresx = filter(gores, representative_description %in% colnames(repmat))
goresx = setNames(goresx$NES, goresx$representative_description)[colnames(repmat)]
generesx = allgenesinGO %>%
  filter(abs(l2FC)>=1) %>%
  select(protein, l2FC) %>% unique()
generesx = setNames(generesx$l2FC, generesx$protein)[rownames(repmat)]
generessignif = allgenesinGO %>%
  # filter((l2FC)<=-0.5) %>%
  filter(abs(l2FC)>=1) %>%
  select(protein, padj) %>% unique()
generessignif = 1+(setNames(generessignif$padj<=0.1, generessignif$protein)[rownames(repmat)])


library(ComplexHeatmap)
nescols = agecol[(generesx>0)+1]
gocols = agecol[(goresx>0)+1]
basesize = 6
column_ha = HeatmapAnnotation(`Normalized\nEnrichment\nScore` = anno_barplot(abs(goresx), height = unit(1,'cm'), gp = gpar(fill = gocols)), annotation_name_gp = gpar(fontsize = basesize+1))
row_ha = rowAnnotation(Log2FC = anno_barplot(abs(generesx), width = unit(1,'cm'), gp = gpar(fill = nescols)), annotation_name_gp = gpar(fontsize = basesize+1))
inchconv = 0.393701
geneorder = (data.frame(repmat) %>% mutate(protein = rownames(repmat)) %>% 
               mutate(n = rowSums(repmat)) %>% 
               set_names(c('gr1','gr2','gr3','protein','n')) %>% 
               arrange(-gr2,-gr3,-gr1))$protein

pdf(file="./results/msfigures/KMproteomics_GOrepairgenes_l2fc1.pdf", width = 10*inchconv, height = 25*inchconv, useDingbats = F)
Heatmap(repmat, top_annotation = column_ha, right_annotation = row_ha, col = c('gray95','gray25'), 
        rect_gp = gpar(col = "white", lwd = 1), 
        column_names_max_height = max_text_width(
          colnames(repmat), 
          gp = gpar(fontsize = basesize)),
        row_names_gp = gpar(fontsize = basesize, fontface = generessignif), column_names_gp = gpar(fontsize = basesize), 
        # row_km = 3, 
        # column_km = 4, 
        row_title = NULL, 
        row_gap = unit(basesize, "pt"), column_gap = unit(basesize, "pt"), border = 'gray35',
        width = ncol(repmat)*unit(basesize+1, "pt"),
        # row_order = geneorder,
        row_order = order(-generesx),
        column_order = order(-goresx),
        column_title_gp = gpar(fontsize = basesize + 2),
        height = nrow(repmat)*unit(basesize, "pt"), show_heatmap_legend = F, column_title = 'Kidney Marrow Proteomics\nGenes in DNA Repair-related\nGO categories (|Log2FC| >= 1)')
dev.off()
