######################
# This script performs genomic selection
# Authors: Matteo Dell'Acqua, Jesse Poland
# Date: May 30th, 2022
########################

options(stringsAsFactors = F)

library(rrBLUP)
library(ggplot2)

wd<-"C:/Users/admin/OneDrive - Scuola Superiore Sant'Anna/projects/EtNAM.GS/analysis/3.make.predictions"
setwd(wd)

#get in EtNAM BLUP data
load("../1.calculate.BLUPs/BLUP.values.EtNAM.Rdata")

#get in diversity panel BLUP data used in the 3D breeding paper experiment
load("C:/Users/admin/OneDrive - Scuola Superiore Sant'Anna/projects/tricot.genomics/analyses/0a.make.blups/BLUPs.h2.Rdata")

#get in stability data
load("../1.calculate.BLUPs/genoytpe.stability.analysis.Rdata")
stability<-cbind(rownames(stability), stability)

## fix phenotypes
phenolist<-list(blups.farm, blups.farm.nam, blups.met, blups.met.nam, stability)

for (i in 1:length(phenolist)){
  rownames(phenolist[[i]])<-phenolist[[i]][,1]
  if(length(grep("^ID_", rownames(phenolist[[i]])))>0){
    rownames(phenolist[[i]])<-sub("^ID_","",rownames(phenolist[[i]]))
  }
  phenolist[[i]]<-phenolist[[i]][order(rownames(phenolist[[i]])),]
  phenolist[[i]]<-phenolist[[i]][,-1]
  apply(phenolist[[i]], 2, as.numeric)
}

blups.farm<-phenolist[[1]]
blups.farm.nam<-phenolist[[2]]
blups.met<-phenolist[[3]]
blups.met.nam<-phenolist[[4]]
stability<-phenolist[[5]]

#put together phenotypes
pheno.p<-cbind(blups.met, blups.farm)
pheno.nam<-cbind(blups.met.nam, blups.farm.nam)

## get in genotypes
load("../2.impute.snps/genotypic.data.nam.panel.Rdata")

#make sure that the marker sets are equal for nam and panel
imputed.p<-imputed.p[,colnames(imputed.p) %in% colnames(imputed.nam)]
imputed.nam<-imputed.nam[,colnames(imputed.nam) %in% colnames(imputed.p)]
dim(imputed.p)
dim(imputed.nam)
rownames(imputed.p)<-colnames(geno)
rownames(imputed.nam)<-colnames(genonam)

#reorder the dataset
imputed.p<-imputed.p[,order(colnames(imputed.p))]
imputed.nam<-imputed.nam[,order(colnames(imputed.nam))]
imputed.p<-imputed.p[order(rownames(imputed.p)),]
imputed.nam<-imputed.nam[order(rownames(imputed.nam)),]

stopifnot(all(colnames(imputed.nam) == colnames(imputed.p)))

imputed.p[1:4,1:4]
imputed.nam[1:4,1:4]

#make sure to keep only samples that have been phenotyped AND genoyped
imputed.p<-imputed.p[which(rownames(imputed.p) %in% rownames(pheno.p)),]
pheno.p<-pheno.p[which(rownames(pheno.p) %in% rownames(imputed.p)),]
stopifnot(all(rownames(imputed.p) == rownames(pheno.p)))

imputed.nam<-imputed.nam[which(rownames(imputed.nam) %in% rownames(pheno.nam)),]
pheno.nam<-pheno.nam[which(rownames(pheno.nam) %in% rownames(imputed.nam)),]
stopifnot(all(rownames(imputed.nam) == rownames(pheno.nam)))

###########
# test different GS scenarios

###
## conduct GS using diversity panel to predict the NAM

## subset phenotypes to relevant ones
traits<-c("DB", "DF", "DM", "PH", "NET", "SPL", "SPS", "BM", "GY", "EARLINESS", "OVERALL", "SPIKE", "TILLER")

#set number of sample runs
nsubs<-100

#set training popualtion size
fraction<-0.8

#run the prediction on different subsets of the panel
tsize<-floor(nrow(imputed.p)*1)
vsize<-floor(nrow(imputed.nam)*fraction)
  
predresults<-list()
for(zz in 1:nsubs){ ##zz=1
  print(zz)
  #extract the training set
  training<-sample(1:nrow(imputed.p), tsize)             
  valid<-sample(1:nrow(imputed.nam), vsize)             
  
  #identify training and validation populations
  t.pheno<-pheno.p[training,traits]
  t.geno<-imputed.p[training,]
  
  v.pheno<-pheno.nam[valid,]
  v.geno<-imputed.nam[valid,]
  
  #run the loop
  traitout<-list()
  for (p in 1:ncol(t.pheno)){ # p=1
    curp<-colnames(t.pheno)[p]
    #print(curp)
    trainpheno<-t.pheno[,curp]
 
    model<-mixed.solve(trainpheno, Z=t.geno, K=NULL, SE= F, return.Hinv = F )
    effects<-model$u
    
    genoeffect <- v.geno %*% effects
    GEBV<-genoeffect[,1] + rep(model$beta, ncol(genoeffect))
    
    save(model, GEBV, file=paste0("./output/GEBV.", curp, ".panel.on.nam.Rdata"))
    
    accuracy<-cor(GEBV, v.pheno, use="complete")
    accuracy<-data.frame(t(accuracy))
    traitout[[p]]<-accuracy
    names(traitout[[p]])<-curp
  }
  
  traitoutdf<-do.call(cbind, traitout)
  colnames(traitoutdf)<-colnames(t.pheno)
  predresults[[zz]]<-traitoutdf
}#for zz

save(predresults, file=paste0("./output/panel.to.nam.gs.", fraction, ".subset.Rdata"))

###
## conduct GS using diversity panel to predict STABILITY of NAM phenotypes

## subset phenotypes to relevant ones
traits<-c("DB", "DF", "DM", "PH", "NET", "SPL", "SPS", "BM", "GY", "EARLINESS", "OVERALL", "SPIKE", "TILLER")

#set number of sample runs
nsubs<-100

#set training popualtion size
fraction<-0.8

#run the prediction on different subsets of the panel
tsize<-floor(nrow(imputed.p)*1)
vsize<-floor(nrow(imputed.nam)*fraction)

predresults<-list()
for(zz in 1:nsubs){ ##zz=1
  print(zz)
  #extract the training set
  training<-sample(1:nrow(imputed.p), tsize)             
  valid<-sample(1:nrow(imputed.nam), vsize)             
  
  #identify training and validation populations
  t.pheno<-pheno.p[training,traits]
  t.geno<-imputed.p[training,]
  
  v.pheno<-stability[valid,]
  v.geno<-imputed.nam[valid,]
  
  #run the loop
  traitout<-list()
  for (p in 1:ncol(t.pheno)){ # p=1
    curp<-colnames(t.pheno)[p]
    #print(curp)
    trainpheno<-t.pheno[,curp]
    
    model<-mixed.solve(trainpheno, Z=t.geno, K=NULL, SE= F, return.Hinv = F )
    effects<-model$u
    
    genoeffect <- v.geno %*% effects
    GEBV<-genoeffect[,1] + rep(model$beta, ncol(genoeffect))
    
    save(model, GEBV, file=paste0("./output/GEBV.", curp, ".panel.on.stability.nam.Rdata"))
    
    accuracy<-cor(GEBV, v.pheno, use="complete")
    accuracy<-data.frame(t(accuracy))
    traitout[[p]]<-accuracy
    names(traitout[[p]])<-curp
  }
  
  traitoutdf<-do.call(cbind, traitout)
  colnames(traitoutdf)<-colnames(t.pheno)
  predresults[[zz]]<-traitoutdf
}#for zz

save(predresults, file=paste0("./output/panel.to.stability.nam.gs.", fraction, ".subset.Rdata"))

###
## conduct GS using NAM to predict NAM

#set number of sample runs
nsubs<-100

#set training popualtion size
fraction<-0.8

#run the prediction on different subsets of the panel
tsize<-floor(nrow(imputed.nam)*fraction)
vsize<-floor(nrow(imputed.nam)*(1-fraction))

predresults<-list()
for(zz in 1:nsubs){ ##zz=1
  print(zz)
  #extract the training set
  training<-sample(1:nrow(imputed.nam), tsize)             
  valid<-sample(1:nrow(imputed.nam), vsize)             
  
  #identify training and validation populations
  t.pheno<-pheno.nam[training,]
  t.geno<-imputed.nam[training,]
  
  v.pheno<-pheno.nam[valid,]
  v.geno<-imputed.nam[valid,]
  
  #run the loop
  traitout<-list()
  for (p in 1:ncol(t.pheno)){ # p=1
    curp<-colnames(t.pheno)[p]
    #print(curp)
    trainpheno<-t.pheno[,curp]
    
    model<-mixed.solve(trainpheno, Z=t.geno, K=NULL, SE= F, return.Hinv = F )
    effects<-model$u
    
    genoeffect <- v.geno %*% effects
    GEBV<-genoeffect[,1] + rep(model$beta, ncol(genoeffect))
    
    save(model, GEBV, file=paste0("./output/GEBV.", curp, ".nam.on.nam.Rdata"))
    
    accuracy<-cor(GEBV, v.pheno, use="complete")
    accuracy<-data.frame(t(accuracy))
    traitout[[p]]<-accuracy
    names(traitout[[p]])<-curp
  }
  
  traitoutdf<-do.call(cbind, traitout)
  colnames(traitoutdf)<-colnames(t.pheno)
  predresults[[zz]]<-traitoutdf
}#for zz

save(predresults, file=paste0("./output/nam.to.nam.gs.", fraction, ".subset.Rdata"))

###
## conduct GS using NAM to predict NAM across locations

kul.t<-pheno.nam[,grep("kul", colnames(pheno.nam))]
ger.t<-pheno.nam[,grep("ger", colnames(pheno.nam))]
ade.t<-pheno.nam[,grep("ade", colnames(pheno.nam))]

kul.v<-pheno.nam[,-grep("kul", colnames(pheno.nam))]
ger.v<-pheno.nam[,-grep("ger", colnames(pheno.nam))]
ade.v<-pheno.nam[,-grep("ade", colnames(pheno.nam))]

tlist<-list(kul.t, ger.t, ade.t)
vlist<-list(kul.v, ger.v, ade.v)
lapply(tlist, dim)
lapply(vlist, dim)

names(tlist)<-c("kulumsa", "geregera", "adet")

#run the prediction on different subsets of the panel
predresults<-list()
for(zz in 1:length(tlist)){ ##zz=1
  print(names(tlist)[[zz]])

  #identify training and validation populations
  t.pheno<-tlist[[zz]]
  t.geno<-imputed.nam
  
  v.pheno<-vlist[[zz]]
  v.geno<-imputed.nam
  
  #run the loop
  traitout<-list()
  for (p in 1:ncol(t.pheno)){ # p=1
    curp<-colnames(t.pheno)[p]
    #print(curp)
    trainpheno<-t.pheno[,curp]
    
    model<-mixed.solve(trainpheno, Z=t.geno, K=NULL, SE= F, return.Hinv = F )
    effects<-model$u
    
    genoeffect <- v.geno %*% effects
    GEBV<-genoeffect[,1] + rep(model$beta, ncol(genoeffect))
    
    save(model, GEBV, file=paste0("./output/GEBV.", curp , ".", names(tlist)[[zz]], "nam.on.nam.Rdata"))
    
    accuracy<-cor(GEBV, v.pheno, use="complete")
    accuracy<-data.frame(t(accuracy))
    traitout[[p]]<-accuracy
    names(traitout[[p]])<-curp
  }
  
  traitoutdf<-do.call(cbind, traitout)
  colnames(traitoutdf)<-colnames(t.pheno)
  predresults[[zz]]<-traitoutdf
}#for zz

save(predresults, file=paste0("./output/nam.to.nam.gs.byloc.subset.Rdata"))


## conduct gs to predict panel from NAM

#set training popualtion size
fraction<-0.8

#set number of sample runs
nsubs<-20

#run the prediction on different subsets of the panel
tsize<-floor(nrow(imputed.nam)*fraction)
vsize<-floor(nrow(imputed.p)*1)

predresults<-list()
for(zz in 1:nsubs){ ##zz=1
  print(zz)
  #extract the training set
  training<-sample(1:nrow(imputed.nam), tsize)             
  valid<-sample(1:nrow(imputed.p), vsize)             
  
  #identify training and validation populations
  t.pheno<-pheno.nam[training,]
  t.geno<-imputed.nam[training,]
  
  v.pheno<-pheno.p[valid,]
  v.geno<-imputed.p[valid,]
  
  #run the loop
  traitout<-list()
  for (p in 1:ncol(t.pheno)){ # p=1
    curp<-colnames(t.pheno)[p]
    #print(curp)
    trainpheno<-t.pheno[,curp]
    
    model<-mixed.solve(trainpheno, Z=t.geno, K=NULL, SE= F, return.Hinv = F )
    effects<-model$u
    
    genoeffect <- v.geno %*% effects
    GEBV<-genoeffect[,1] + rep(model$beta, ncol(genoeffect))
    
    save(model, GEBV, file=paste0("./output/GEBV.", curp, ".nam.on.panel.Rdata"))
    
    accuracy<-cor(GEBV, v.pheno, use="complete")
    accuracy<-data.frame(t(accuracy))
    traitout[[p]]<-accuracy
    names(traitout[[p]])<-curp
  }
  
  traitoutdf<-do.call(cbind, traitout)
  colnames(traitoutdf)<-colnames(t.pheno)
  predresults[[zz]]<-traitoutdf
}#for zz

save(predresults, file=paste0("./output/nam.to.panel.gs.", fraction, ".subset.Rdata"))
