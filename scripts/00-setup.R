## compilation of common libraries, functions, visualisation settings across scripts. each script starts by running this setup script.

# ### Folder Organisation ####
# sapply(c('docs', 'data/processed', 'data/raw', 'scripts', 'results'),
#        function(folder){
#          system(paste('mkdir -p',folder))
#        }) # run only at the beginning of the project

#### Libraries used across all scripts ####
easypackages::libraries(
  'tidyverse'
)

#### Visualization helpers ####
easypackages::libraries(
  'ggrepel',
  'ggpubr',
  'ggthemes',
  'Rvislib'
  # 'gghalves',
  # 'ggridges'
)
theme_set(theme_pubr(base_size = 8))
pntnorm <- (1/0.352777778)
geom_point2 <- function(...)geom_point(shape=21,...)

## Colors ##
agecol = setNames(c('cornflowerblue','indianred2'),c('young','old'))

#### Functions ####

plotsave <- function(ggobj, prefix, width, height, ...){
  path = strsplit(prefix,'/')[[1]]
  path = paste(path[-length(path)], collapse = '/')
  if(!file.exists(path)){
    system(paste('mkdir -p',path))
  }
  saveRDS(object = ggobj, file = paste(prefix,'.rds',sep=''))
  ggsave(file = paste(prefix, '.pdf', sep = ''), plot = ggobj, units = 'cm', 
         width = width, height = height, useDingbats = F, limitsize = F)
  ggsave(file = paste(prefix, '.png', sep = ''), plot = ggobj, units = 'cm', 
         width = width, height = height, limitsize = F)
}

tablesave <- function(tib, prefix, ...){
  path = strsplit(prefix,'/')[[1]]
  path = paste(path[-length(path)], collapse = '/')
  if(!file.exists(path)){
    system(paste('mkdir -p',path))
  }
  readr::write_csv(tib, path = paste(prefix,'.csv', sep = ''), append = F)
  readr::write_tsv(tib, path = paste(prefix,'.tsv', sep = ''), append = F)
  saveRDS(object = tib, file = paste(prefix,'.rds',sep=''))
}

objsave <- function(obj, prefix){
  path = strsplit(prefix,'/')[[1]]
  path = paste(path[-length(path)], collapse = '/')
  if(!file.exists(path)){
    system(paste('mkdir -p',path))
  }
  saveRDS(object = obj, file = paste(prefix,'.rds', sep = ''))
}


### Folder Organisation ####
# sapply(c('scRNAseq','KM_proteomics','plasma_proteomics','helperscripts','imageStream'),
#        function(folder){
#          system(paste('mkdir -p ./scripts/',folder, sep=''))
#        }) # run only at the beginning of the project
