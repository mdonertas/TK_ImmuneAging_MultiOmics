source('./scripts/00-setup.R')

gores = readRDS('./data/processed/KMproteomics/DE_GOenrich_allsamples.rds')

goresplot = gores %>%
  filter(BY_corrected_p<=0.05) %>%
  filter(representative == ID) %>%
  mutate(NESres = factor(c('Enriched in young','Enriched in old')[1+(NES>0)],levels = c('Enriched in young','Enriched in old'))) %>%
  ggplot(aes(y = abs(NES), x = reorder(Description,NES), fill = NESres)) +
  geom_bar(stat = 'identity') +
  ylab("Normalized enrichment score") + xlab(NULL) +
  scale_fill_manual(values = unname(agecol)) +
  theme_pubr(base_size = 6) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))  
goresplot

plotsave(goresplot,'./results/msfigures/KMproteomics_GOresults',width = 16, height = 10)

