######################
# This script runs GWAS with GAPIT3
# Authors: Matteo Dell'Acqua
# Date: May 30th, 2022
########################

wd<-"C:/Users/admin/OneDrive - Scuola Superiore Sant'Anna/projects/EtNAM.GS/analysis/11.GWAS/results"
setwd(wd)

#install.packages("devtools")
#devtools::install_github("jiabowang/GAPIT3",force=TRUE)
library(GAPIT3)

#get pheno
load("../../1.calculate.BLUPs/BLUP.values.EtNAM.Rdata")
myY<-blups.met.nam
myYfarm<-blups.farm.nam

#get geno
myG <- read.table("../../../data/clean.full.NAM.geno.data.15k.refseq.positions.q10.hmp", head = FALSE)

#run GWAS for metric traits
for (i in 2:ncol(myY)){
  print(colnames(myY)[i])
  myGAPIT<-GAPIT(
    Y=myY[,c(1,i)], #fist column is ID
    G=myG,
    PCA.total=3,
    model=c("FarmCPU", "Blink"),
    Multiple_analysis=TRUE)
}

#run GWAS for farmer traits
for (i in 2:ncol(myYfarm)){
  print(colnames(myYfarm)[i])
  myGAPIT<-GAPIT(
    Y=myYfarm[,c(1,i)], #fist column is ID
    G=myG,
    PCA.total=3,
    model=c("FarmCPU", "Blink"),
    Multiple_analysis=TRUE)
}
