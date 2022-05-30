######################
# This script takes the output from joinmap and parses it, including to produce consensus maps
# Authors: Matteo Dell'Acqua
# Date: May 30th, 2022
########################
options(stringsAsFactors = F)

wd<-"C:/Users/admin/OneDrive - Scuola Superiore Sant'Anna/projects/EtNAM.GS/analysis/9.QTL.mapping"
setwd(wd)

library(tidyverse)
library(plyr)
library(LPmerge)

#get names of chromsomes
nameschr<-sort(paste0(rep(1:7,2), c("A", "B")))

#read in the phenotypes
load("../1.calculate.BLUPs/BLUP.values.EtNAM.Rdata")

#read in dictionary for markers
markername<-read_delim("../../data/mapmaker/marker.info.txt")
markername$new.marker.name<-sub("^\\*", "", markername$new.marker.name)
markername

#read in maps for each family
filen<-list.files("../../data/maps/")
filen<-filen[grep("^N", filen)]

#create a list to store map lengths
maplengths<-list()

for (i in 1:length(filen)){ #i = 1
  tmp<-filen[i]
  curfam<-sub("\\.map.csv", "", tmp)  
  print(curfam)
  names(maplengths)[i]<-curfam
  
  #load the corresponding sample names
  tmpnames<-read_delim(paste0("../../data/mapmaker/NAM.family.", curfam ,".ids.correspondance.txt"))
  
  #load the file
  df<-read.csv(paste0("../../data/maps/", tmp))
  head(df)  
  
  #make sure everything is as it should be
  df[1:5,1:10]
  
  lgs<-unique(df$LG)
  lgdict<-data.frame(lg=lgs, newlg=1:length(unique(lgs)))
  lgdict
  df$LG<-mapvalues(df$LG, from=lgdict$lg, to=lgdict$newlg)
  df$Pos<-as.numeric(df$Pos)
  head(df)
  print(range(df$Pos, na.rm=T))
  
  #fix sample names
  colnames(df)<-sub("^X", "", colnames(df))
  colnames(df) <- mapvalues(colnames(df), from=tmpnames$idx, to=tmpnames$ID)
  
  #order the file by chr, lg, position
  df<-df[order(df$CHR, df$LG, df$Pos),]
  
  #split by chr to fix the LG number
  bychr<-split(df, df$CHR)
  
  #for each chromosome sort LG within CHR
  outordered<-list()
  for(cr in 1:length(bychr)){ #cr=1
    curchr<-names(bychr)[cr]
    tmpchr<-bychr[[cr]]

    #get mean physical position for each LG
    lgs<-split(tmpchr, tmpchr$LG)
    
    #prepare output for lg physical pos 
    mpos<-c()
    for(j in 1:length(lgs)){#j=1
      tmplg<-lgs[[j]]
      markerpos<-as.numeric(sub("^.*\\.m\\.", "", tmplg$Locus))
      mpos[j]<-mean(markerpos)
    }#for lg
    names(mpos)<-names(lgs)
    mpos<-sort(mpos)
    #create new names and map them back to the original file
    mapper<-data.frame(oldLG=names(mpos), newLG=paste0(curchr, "_LG_", 1:length(mpos)))
    
    #fix the LG names
    tmpchr$LG<-mapvalues(tmpchr$LG, from=mapper$oldLG, to=mapper$newLG, warn_missing = F)
    #put everything backe where it belongs
    outordered[[cr]]<-tmpchr
  }#
  
  #assemble the df with the new LG names
  dfoldlg<-df
  df<-do.call(rbind, outordered)
  
  #sort the file
  df<-df[order(df$CHR, df$LG, df$Pos),]
  head(df)
  
  #fix locus names
  df$Locus <- mapvalues(df$Locus, from=markername$new.marker.name, to=markername$marker, warn_missing = F)
  df[1:5,1:20]
  
  #fix duplicated markers if any
  dups<-which(duplicated(df$Locus))
  if(length(dups)>0){
    #remove duplicates for now
    df<-df[-dups,]
  }

  #####slice the file to produce different outputs useful for r/QTL
  
  #genetic map
  map<-df[,c("Locus", "LG", "Pos")]
  names(map)<-c("marker", "chr", "pos")
  head(map)
  maplengths[[i]]<-map
  write.csv(map, file = paste0("./rqtl.files/", curfam, "_gmap.csv"), row.names = F, quote=F)
  
  #get out genotypeing data
  #genotypic file
  geno<-df[,grep("^[0-9]", colnames(df))]
  geno<-data.frame(t(geno))
  colnames(geno)<-df$Locus
  geno[1:5,1:5]
  class(geno[,1])
  
  geno[geno=="a"]<-"A"
  geno[geno=="b"]<-"B"
  geno[geno=="h"]<-"H"
  
  geno<-data.frame(id = rownames(geno), geno)
  geno[1:5,1:5]
  write.csv(geno, file = paste0("./rqtl.files/", curfam, "_geno.csv"), row.names = F, quote=F)
  
  #phenotypic file
  dfphemet<-blups.met.nam[which(blups.met.nam$ID %in% geno$id),]
  dfphefarm<-blups.farm.nam[which(blups.farm.nam$ID %in% geno$id),]
  
  head(dfphemet)
  head(dfphefarm)
  
  dfphe<-merge(dfphemet, dfphefarm, by="ID")
  head(dfphe)
  
  write.csv(dfphe, file = paste0("./rqtl.files/", curfam, "_pheno.csv"), row.names = F, quote=F)
  
  #write the corresponding yaml file
  fileConn<-file(paste0("./rqtl.files/", curfam, ".yaml"))
    l1<-paste0("# Data from EtNAM family ", curfam)  
    l2<-"crosstype: riself"
    l3<-paste0("geno: ", curfam,"_geno.csv")
    l4<-paste0("pheno: ", curfam,"_pheno.csv")
    l5<-paste0("gmap: ", curfam,"_gmap.csv")
    l6<-"alleles:"
    l7<-"- A"
    l8<-"- B"
    l9<-"genotypes:"
    l10<-"   A:  1"
    l11<-"   B:  2"
    l12<-"na.strings:"
    l13<-"-  '-'"
    writeLines(c(l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,l13), fileConn)
  close(fileConn)

}

#summarize map lengths
names(maplengths)
lgs<-lapply(maplengths, function(x) unique(x[,2]))
lgnumber<-unlist(lapply(lgs, length))
cMlength<-unlist(lapply(maplengths, function(x) max(x[,3])))
nmarks<-unlist(lapply(maplengths, nrow))
  
lgnumber
cMlength
nmarks



#create an input for MergeMap
#read in maps for each family
filen<-list.files("../../analysis/9.QTL.mapping/rqtl.files/")
filen<-filen[grep("_gmap.csv",filen)]

for (i in 1:length(filen)){ #i = 1
  tmp<-filen[i]
  curfam<-sub("\\_gmap.csv", "", tmp)  
  tmpmap<-read.csv(paste0("./rqtl.files/", tmp))
  
  #split by chromosome, deriving chr name from LG naming
  tmpmap$chr<-sub("_LG_.*$", "", tmpmap$chr)
  bylg<-split(tmpmap, tmpmap$chr)
  length(bylg)  

  #open a file and write in append mode
  tmpname<-paste0("./mergemap.files/", curfam, ".mergemap.txt")

  #initialize file and keep writing it
  tmplg<-bylg[[1]]
  tmplg<-as.matrix(tmplg[,c(1,3)])
  head(tmplg)
  
  cat(paste0("group lg1 \n"),file=tmpname)
  cat(";BEGINOFGROUP \n", file=tmpname, append = T)
  write.table(tmplg, tmpname, sep="\t",append=TRUE, col.names=F, row.names = F, quote=F)
  cat(";ENDOFGROUP \n", file=tmpname, append = T)
  
  
  for(j in 2:length(bylg)){ #j =2
    
    tmplg<-bylg[[j]]
    tmplg<-as.matrix(tmplg[,c(1,3)])
    head(tmplg)
    
    cat(paste0("\ngroup lg", j, "\n"),file=tmpname, append = T)
    cat(";BEGINOFGROUP \n", file=tmpname, append = T)
    write.table(tmplg, tmpname, sep="\t",append=TRUE, col.names=F, row.names = F, quote=F)
    cat(";ENDOFGROUP \n", file=tmpname, append = T)
    
  } 
      
}

#create an input for LPmerge
#read in maps for each family
filen<-list.files("../../analysis/9.QTL.mapping/rqtl.files/")
filen<-filen[grep("_gmap.csv",filen)]
filen

#prepare list of list of LGs
listout<-list()

for (i in 1:length(filen)){ #i = 1  
  tmp<-filen[i]
  curfam<-sub("\\_gmap.csv", "", tmp)  
  tmpmap<-read.csv(paste0("./rqtl.files/", tmp))
  
  #split by chromosome, deriving chr name from LG naming
  tmpmap$chr<-sub("_LG_.*$", "", tmpmap$chr)
  bylg<-split(tmpmap, tmpmap$chr)
  length(bylg)  
  
  #create an internal list
  chrlist<-list()
  
  for(j in 1:length(bylg)){ #j =1
    
    tmplg<-bylg[[j]]
    tmplg<-tmplg[,c(1,3)]
    head(tmplg)
    
    chrlist[[j]]<-tmplg
  } 
  names(chrlist)<-nameschr
  listout[[i]]<-chrlist
  
  names(listout)[i]<-curfam
  
}

#extract elements from lists and run consensus
for (i in 1:14){ # i=1
  print(paste("Working on chr", nameschr[i]))
  tmplist<-lapply(listout, function(x) x[[i]])
  tmpout<-LPmerge(tmplist, max.interval = 4, weights = NULL)
  save(tmpout, file=paste0("./consensus/consensus.map.int.4.chr.", nameschr[i],".Rdata"))
}

