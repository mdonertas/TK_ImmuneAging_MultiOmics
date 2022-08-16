source('./scripts/00-setup.R')

km_all = readRDS('./data/processed/KMproteomics/DEwilcox_allsamples.rds')
km_outrem = readRDS('./data/processed/KMproteomics/DEwilcox_outlier_removed.rds')

go_km_all = readRDS('./data/processed/KMproteomics/DE_GOenrich_allsamples.rds')
go_km_outrem = readRDS('./data/processed/KMproteomics/DE_GOenrich_outlier_removed.rds')

pl_all = readRDS('./data/processed/plasmaProteomics/DEwilcox_allsamples.rds')
pl_outrem = readRDS('./data/processed/plasmaProteomics/DEwilcox_outlier_removed.rds')

go_pl_all = readRDS('./data/processed/plasmaProteomics/DE_GOenrich_allsamples.rds')
go_pl_outrem = readRDS('./data/processed/plasmaProteomics/DE_GOenrich_outlier_removed.rds')


p_c = km_all %>%
  select(upperprotein, l2FC) %>%
  rename(all = l2FC) %>%
  left_join(select(km_outrem, upperprotein, l2FC)) %>%
  rename(outrem = l2FC) %>%
  ggplot(aes(x = all, y = outrem)) +
  geom_point(size = 0.2, alpha = 0.2) +
  # geom_smooth(method = 'lm', color = 'darkred', se = F) +
  stat_cor(size = 8/pntnorm, cor.coef.name = 'rho') +
  xlab('with all samples') + ylab('after outlier removal') +
  ggtitle('KM Proteomics - Differential Abundance')

p_d = go_km_all %>%
  select(ID, NES, BY_corrected_p) %>%
  rename(all = NES) %>%
  left_join(select(go_km_outrem, ID, NES)) %>%
  rename(outrem = NES) %>%
  ggplot(aes(x = all, y = outrem)) +
  geom_point(size = 0.2, alpha = 0.2) +
  # geom_smooth(method = 'lm', color = 'darkred', se = F) +
  stat_cor(size = 8/pntnorm, cor.coef.name = 'rho') +
  xlab('with all samples') + ylab('after outlier removal') +
  ggtitle('KM Proteomics - GSEA') 

p_a = pl_all %>%
  select(upperprotein, l2FC) %>%
  rename(all = l2FC) %>%
  left_join(select(pl_outrem, upperprotein, l2FC)) %>%
  rename(outrem = l2FC) %>%
  ggplot(aes(x = all, y = outrem)) +
  geom_point(size = 0.2, alpha = 0.2) +
  # geom_smooth(method = 'lm', color = 'darkred', se = F) +
  stat_cor(size = 8/pntnorm, cor.coef.name = 'rho') +
  xlab('with all samples') + ylab('after outlier removal') +
  ggtitle('Plasma Proteomics - Differential Abundance')

p_b = go_pl_all %>%
  select(ID, NES, BY_corrected_p) %>%
  rename(all = NES) %>%
  left_join(select(go_pl_outrem, ID, NES)) %>%
  rename(outrem = NES) %>%
  ggplot(aes(x = all, y = outrem)) +
  geom_point(size = 0.2, alpha = 0.2) +
  # geom_smooth(method = 'lm', color = 'darkred', se = F) +
  stat_cor(size = 8/pntnorm, cor.coef.name = 'rho') +
  xlab('with all samples') + ylab('after outlier removal') +
  ggtitle('Plasma Proteomics - GSEA') 

outlierplot = ggarrange(p_a, p_b, p_c, p_d, ncol=2, nrow=2, labels = 'auto', font.label = list(size = 8))

plotsave(outlierplot, './results/msfigures/proteomics_outlier_check', width = 16, height = 12)


