# GO2Gene
library(GO.db)
library(AnnotationDbi)
library(org.Hs.eg.db)
term=select(GO.db,columns=c("GOID","TERM"),keytype="GOID",keys=keys(GO.db))
terms=as.character(term[,2])
names(terms)=as.character(term[,1])
# save(terms,file="./data/processed/GO2Gene/goterms_20220407.RData")
saveRDS(terms,'./data/processed/GO2Gene/goterms_20220407.rds')
x<-c(as.list(GOBPOFFSPRING))
xx<-sapply(1:length(x),function(i)x[[i]]=c(names(x)[i],x[[i]]))
xxx<-sapply(1:length(xx),function(i)xx[[i]][complete.cases(xx[[i]])])
length(x)
length(xx)
length(xxx)
names(xx)<-names(x)
names(xxx)<-names(x)
xxx[[1]]
go2gene=select(org.Hs.eg.db,keys=names(xxx),columns = c("GO","ENSEMBL"),keytype = "GO")
GO2Gene=list()
for(i in 1:length(xxx)){
  print(i/length(xxx))
  myg=go2gene[go2gene[,1]%in%xxx[[i]],4]
  myg=unique(myg[complete.cases(myg)])
  GO2Gene[[i]]=myg
}
names(GO2Gene)=names(xxx)
head(GO2Gene)
length(GO2Gene)

saveRDS(GO2Gene,file="./data/processed/GO2Gene/GO2GeneBP_20220407.rds")

rm(list=ls())

x<-c(as.list(GOMFOFFSPRING))
xx<-sapply(1:length(x),function(i)x[[i]]=c(names(x)[i],x[[i]]))
xxx<-sapply(1:length(xx),function(i)xx[[i]][complete.cases(xx[[i]])])
length(x)
length(xx)
length(xxx)
names(xx)<-names(x)
names(xxx)<-names(x)
xxx[[1]]
go2gene=select(org.Hs.eg.db,keys=names(xxx),columns = c("GO","ENSEMBL"),keytype = "GO")
GO2Gene=list()
for(i in 1:length(xxx)){
  print(i/length(xxx))
  myg=go2gene[go2gene[,1]%in%xxx[[i]],4]
  myg=unique(myg[complete.cases(myg)])
  GO2Gene[[i]]=myg
}
names(GO2Gene)=names(xxx)
head(GO2Gene)
length(GO2Gene)

saveRDS(GO2Gene,file="./data/processed/GO2Gene/GO2GeneMF_20220407.rds")
rm(list=ls())

x<-c(as.list(GOCCOFFSPRING))
xx<-sapply(1:length(x),function(i)x[[i]]=c(names(x)[i],x[[i]]))
xxx<-sapply(1:length(xx),function(i)xx[[i]][complete.cases(xx[[i]])])
length(x)
length(xx)
length(xxx)
names(xx)<-names(x)
names(xxx)<-names(x)
xxx[[1]]
go2gene=select(org.Hs.eg.db,keys=names(xxx),columns = c("GO","ENSEMBL"),keytype = "GO")
GO2Gene=list()
for(i in 1:length(xxx)){
  print(i/length(xxx))
  myg=go2gene[go2gene[,1]%in%xxx[[i]],4]
  myg=unique(myg[complete.cases(myg)])
  GO2Gene[[i]]=myg
}
names(GO2Gene)=names(xxx)
head(GO2Gene)
length(GO2Gene)

saveRDS(GO2Gene,file="./data/processed/GO2Gene/GO2GeneCC_20220407.rds")
