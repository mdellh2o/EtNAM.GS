######################
# This script runs QTL mappung with r/QTL2
# Authors: Matteo Dell'Acqua
# Date: May 30th, 2022
########################

#set the higher level that will be used in this script
maindir<-"C:/Users/admin/OneDrive - Scuola Superiore Sant'Anna/projects/EtNAM.GS/analysis/9.QTL.mapping"
setwd(maindir)

options(stringsAsFactors=F)
options(scipen=999)
library(plyr)
library(qtl2)


filen<-list.files("rqtl.files/")
#get yamls
filen<-filen[grep("yaml", filen)]


#run genome reconstruction and QTL mapping for all families
for (i in 1:length(filen)){ #i=1
  
  tmpfam<-sub("\\.yaml", "", filen[i])  
  print(tmpfam)
  cross<-read_cross2(paste0("rqtl.files/",filen[i]))
  print(cross)  
  
  #insert pseudomarkers
  map <- insert_pseudomarkers(cross$gmap, step=1)
  
  #get phenos
  pheno<- cross$pheno
  
  #get genotype probs
  gp<-calc_genoprob(cross)
  
  #get kinship
  #kinship <- calc_kinship(gp)
  kinship_loco <- calc_kinship(gp, "loco", cores=4)
  
  #out <- scan1(gp, pheno, kinship)
  out <- scan1(gp, pheno,  kinship_loco)
  
  #save output
  save(out, kinship_loco, gp, map, cross, file=paste0("./rqtl.files/", tmpfam, ".QTL.mapping.data.Rdata"))
  
}#for i

#make permutations
for (i in 1:length(filen)){ #i=1
  tmpfam<-sub("\\.yaml", "", filen[i])  
  print(tmpfam)
  load(paste0("./rqtl.files/", tmpfam, ".QTL.mapping.data.Rdata"))
  perm <- scan1perm(gp, cross$pheno, n_perm=1000) #############
  save(perm, file=paste0("./rqtl.files/", tmpfam, ".permutations.Rdata") )
}


#extract results from QTL
for (i in 1:length(filen)){ #i=1
  
  tmpfam<-sub("\\.yaml", "", filen[i])  
  print(tmpfam)
  
  if(!dir.exists(paste0("./rqtl.files/", tmpfam, ".results"))){
    dir.create(paste0("./rqtl.files/", tmpfam, ".results"))
  }
  
  #load QTL mapping and permutations
  load(paste0("./rqtl.files/", tmpfam, ".QTL.mapping.data.Rdata"))
  load(paste0("./rqtl.files/", tmpfam, ".permutations.Rdata"))
  
  #extract thresholds for each trait
  thr<-apply(perm, 2, function(x) quantile(x, .9))
  
  #save and plot QTL peaks conditional on threshold
  pks<-find_peaks(out, map, threshold=thr,  peakdrop=1, prob=0.9)
  save(pks,file=paste0("./rqtl.files/", tmpfam, ".results/", tmpfam, ".QTL.peaks.Rdata") )
  
  #get phenos
  pheno<- cross$pheno
  
  for(j in 1:ncol(pheno)){ #j=1
 #   ymx <- maxlod(out[,j]) # overall maximum LOD score
    
    png(paste0("./rqtl.files/", tmpfam, ".results/", colnames(pheno)[j], ".png"), width=1200, height = 600)
      par(mar=c(5.1, 4.1, 1.1, 1.1))
      plot(out, map, lodcolumn=j, col="slateblue", main=colnames(pheno)[j], las=2,xlab="")
      abline(h=thr[j], col="red", lty=2)
    dev.off()  

  }#for j
  
}# for i


#perform mapping again but using DH as a covariate

for (i in 1:length(filen)){ #i=1
  tmpfam<-sub("\\.yaml", "", filen[i])  
  print(tmpfam)
  
  #load QTL mapping and permutations
  load(paste0("./rqtl.files/", tmpfam, ".QTL.mapping.data.Rdata"))
  load(paste0("./rqtl.files/", tmpfam, ".permutations.Rdata"))
  
  #get phenos
  pheno<- cross$pheno
  
  #define covariate trait
  covtrait<-"DB"
  
  #extract thresholds for each trait
  thr<-apply(perm, 2, function(x) quantile(x, .9))

  #add covariate to QTL scan
  outcov <- scan1(gp, pheno,  addcovar=pheno[,covtrait], kinship_loco)
  
  #save and plot QTL peaks conditional on threshold
  pks<-find_peaks(outcov, map, threshold=thr,  peakdrop=1, prob=0.9)
  save(pks,file=paste0("./rqtl.files/", tmpfam, ".results/", colnames(pheno)[j], ".QTL.", covtrait, "covariate.peaks.Rdata") )
  
  for(j in 1:ncol(pheno)){ #j=1
    #   ymx <- maxlod(out[,j]) # overall maximum LOD score
    
    png(paste0("./rqtl.files/", tmpfam, ".results/", colnames(pheno)[j], ".DB.covariate.png"), width=1200, height = 600)
    par(mar=c(5.1, 4.1, 1.1, 1.1))
    plot(out, map, lodcolumn=j, col="slateblue", main=colnames(pheno)[j], las=2,xlab="")
    abline(h=thr[j], col="red", lty=2)
    dev.off()  
    
  }#for j
  
}# for i




