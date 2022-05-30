######################
# This prepares SNP data for genomic selection imputing missing values
# Authors: Matteo Dell'Acqua
# Date: May 30th, 2022
########################

options(stringsAsFactors = F)

library(rrBLUP)

wd<-"C:/Users/admin/OneDrive - Scuola Superiore Sant'Anna/projects/EtNAM.GS/analysis/2.impute.snps"
setwd(wd)

#load in the datasets
load("../../data/diversity.panel.data.gp.Rdata")
load("../../data/etnam.data.gp.Rdata")


#start by fixing genoyptic data
#set a function doing the conversion
convert<-function(x){
  tmp<-as.character(x)
  tmp[tmp=="N"]<-NA
  tmp[tmp  %in% c("R","Y","S","W","K","M") ]<-0
  uniqus<-sort(unique(tmp))
  uniqus<-uniqus[!uniqus %in% c("0", NA)]
  if(length(uniqus)>1){
    tmp[tmp==uniqus[1]]<- 1
    tmp[tmp==uniqus[2]]<- -1
  }    
  tmp[tmp==uniqus[1]]<-1
  #overwrite
  return(tmp)
}

###########
#convert snps in -1,0,1 and NA

###do some data preparation
#diversity panel data
geno[geno=="N"]<-NA
geno[which(geno %in% c("R","Y","S","W","K","M"))]<-0

geno2<-apply(geno, 1, convert)
geno2<-data.frame(geno2)
rownames(geno2)<-colnames(geno)
geno2<-apply(geno2, 2, as.numeric)

#diversity panel data subsetted to NAM data
genored<-geno[rownames(geno) %in% rownames(genonam),]
genored<-as.matrix(genored[order(rownames(genored)),])

genored[genored=="N"]<-NA
genored[which(genored %in% c("R","Y","S","W","K","M"))]<-0

genored2<-apply(genored, 1, convert)
genored2<-data.frame(genored2)
rownames(genored2)<-colnames(genored)
genored2<-apply(genored2, 2, as.numeric)

#NAM data
genonam[genonam=="N"]<-NA
genonam[which(genonam %in% c("R","Y","S","W","K","M"))]<-0

genonam2<-apply(genonam, 1, convert)
genonam2<-data.frame(genonam2)
rownames(genonam2)<-colnames(genonam)
genonam2<-apply(genonam2, 2, as.numeric)

###perform imputation with rrBLUP
#diversity panel data
imputed<-A.mat(geno2, max.missing=0.5, impute.method="mean", return.imputed=T)
failrate<-apply(imputed$imputed,2, function(x) length(which(is.na(x))))
todrop<-which(failrate>0)
if(length(todrop)>1){
  imputed.p<-imputed$imputed[,-which(failrate>0)]
} else {
  imputed.p<-imputed$imputed
}

#reduced dp panel data
imputedred<-A.mat(genored2, max.missing=0.5, impute.method="mean", return.imputed=T)
failrate<-apply(imputedred$imputed,2, function(x) length(which(is.na(x))))
todrop<-which(failrate>0)
if(length(todrop)>1){
  imputed.p.sub<-imputedred$imputed[,-which(failrate>0)]
} else {
  imputed.p.sub<-imputedred$imputed
}

#NAM data
#only polymorphic markers will be retained!
imputednam<-A.mat(genonam2, max.missing=0.8, impute.method="mean", return.imputed=T)
failrate<-apply(imputednam$imputed,2, function(x) length(which(is.na(x))))
todrop<-which(failrate>0)
if(length(todrop)>1){
  imputed.nam<-imputednam$imputed[,-which(failrate>0)]
} else {
  imputed.nam<-imputednam$imputed
}

#fix object names
save(geno, genonam, imputed.p, imputed.p.sub, imputed.nam, file="genotypic.data.nam.panel.Rdata")
