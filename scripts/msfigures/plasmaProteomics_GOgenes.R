source('./scripts/00-setup.R')

gores = readRDS('./data/processed/plasmaProteomics/DE_GOenrich_allsamples.rds')

wiltest = readRDS('./data/processed/plasmaProteomics/DEwilcox_allsamples.rds')

# goreps = c('DNA repair', 'double-strand break repair via homologous recombination', 'response to UV')

repcats = gores %>%
  filter(BY_corrected_p <= 0.05) 

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

row_ha= rowAnnotation(`Normalized\nEnrichment\nScore` = anno_barplot(abs(goresx), height = unit(1,'cm'), gp = gpar(fill = gocols)), annotation_name_gp = gpar(fontsize = basesize+1))

column_ha = HeatmapAnnotation(Log2FC = anno_barplot(abs(generesx), width = unit(1,'cm'), gp = gpar(fill = nescols)), annotation_name_gp = gpar(fontsize = basesize+1))

inchconv = 0.393701

pdf(file="./results/msfigures/plasmaProteomics_GOgenesl2fc1.pdf", width = 25*inchconv, height = 10*inchconv, useDingbats = F)
Heatmap(t(repmat), top_annotation = column_ha, right_annotation = row_ha, 
        col = c('gray95','gray25'), 
        rect_gp = gpar(col = "white", lwd = 1), 
        column_names_max_height = max_text_width(
          rownames(repmat), 
          gp = gpar(fontsize = basesize)),
        row_names_gp = gpar(fontsize = basesize, 
                            fontface = colnames(repmat) %in% 
                              c("regulation of response to external stimulus", 
                                "regulation of coagulation", 
                                "response to stimulus", 
                                "regulation of response to stimulus", 
                                "defense response", 
                                "leukocyte mediated immunity") +1), 
        column_names_gp = gpar(fontsize = basesize, 
                               fontface = rownames(repmat) %in% c('C8G','F5',
                                                                  'C8B','FGB',
                                                                  'FGG') +1), 
        # row_km = 3, 
        # column_km = 4, 
        row_title = NULL, 
        row_gap = unit(basesize, "pt"), 
        column_gap = unit(basesize, "pt"), border = 'gray35',
        height = ncol(repmat)*unit(basesize+1, "pt"),
        column_order = order(generesx),
        # row_order = order(-generesx),
        row_order = order(-goresx),
        column_title_gp = gpar(fontsize = basesize + 2),
        width = nrow(repmat)*unit(basesize, "pt"), show_heatmap_legend = F, 
        column_title = 'Plasma Proteomics\nGenes in significant GO categories (representative)\n(|Log2FC| >= 1)')
dev.off()
