######################
# This script parses QTL mapping results
# Authors: Matteo Dell'Acqua
# Date: May 30th, 2022
########################

wd<-"C:/Users/admin/OneDrive - Scuola Superiore Sant'Anna/projects/EtNAM.GS/analysis/10.check.map.QTL"
setwd(wd)

library(tidyverse)
library(qtl2)
library(plyr)
library(qtlTools)
library(lsmeans)

#load physical position of markers
markerdictionary<-read.delim("../../data/mapmaker/deprecated/marker.info.txt")
markerdictionary<-markerdictionary[order(markerdictionary$chr.phy, markerdictionary$pos.phy),]

#create a cumulative genetic position
bychr<-split(markerdictionary, markerdictionary$chr.phy)
starter<-0
outlist<-list()
for (i in 1:length(bychr)){ #i=1
  tmpchr<-bychr[[i]]
  tmpchr$cumulative.pos<-tmpchr$pos.phy + starter
  tmpchr$chr.start<-min(tmpchr$cumulative.pos)
  tmpchr$chr.stop<-max(tmpchr$cumulative.pos)
  
  outlist[[i]]<-tmpchr
  
  #update the starter
  starter<-max(tmpchr$cumulative.pos)
}
#build back dictionary
markerdictionary<-do.call(rbind, outlist)
head(markerdictionary)
tail(markerdictionary)


#set family names
families<-c("N1", "N3", "N5", "N8", "N10", "N16", "N19", "N32", "N36", "N45", "N46", "N51" )


#### Parse QTL mapping results
outlist<-list()

snp2keep<-list() #this whill hold ONLY those SNPs having a map position in any of the EtNAM families

for (i in families){ #i="N1"
  tmpfam<-i
  
  #load QTL mapping and permutations
  load(paste0("../9.QTL.mapping/rqtl.files/", tmpfam, ".QTL.mapping.data.Rdata"))
  load(paste0("../9.QTL.mapping/rqtl.files/", tmpfam, ".permutations.Rdata"))
  
  #get all markers having a map position
  tmpmap<-cross$gmap
  tmpsnps<-lapply(tmpmap, function(x) names(x))
  names(tmpsnps)<-NULL
  snp2keep[[i]]<-tmpsnps
  
  
  #extract thresholds for each trait
  thr<-apply(perm, 2, function(x) quantile(x, .9))
  
  #save and plot QTL peaks conditional on threshold
  pks<-find_peaks(out, map, threshold=thr,  peakdrop=1, prob=0.9)
  
  #extract markers relative to each peak
  pks$marker<-NA
  for(j in 1:nrow(pks)){#j=1
    tmpLG<-pks[j,3]
    topLG<-map[[tmpLG]]
    
    diff<-topLG-pks[j,4]
    hit<-which(diff==min(abs(diff)))
    
    #if more than one, take the first
    if(length(hit)>1){
      hit<-hit[1]
    }
    
    pks$marker[j]<-names(hit)
  }#for j
  pks
  
  #now get chr and pos on physical map for each marker
  pks$chr.p<-mapvalues(pks$marker, from=markerdictionary$marker, to=markerdictionary$chr.phy, warn_missing =F)
  pks$pos.p<-mapvalues(pks$marker, from=markerdictionary$marker, to=markerdictionary$pos.phy, warn_missing =F)
  pks$cumpos<-mapvalues(pks$marker, from=markerdictionary$marker, to=markerdictionary$cumulative.pos, warn_missing =F)
  
  pks$family<-tmpfam
  
  outlist[[i]]<-pks
}

pklist<-do.call(rbind, outlist)
lapply(snp2keep, length)

save(pklist, file="parsed.QTL.peaks.Rdata")

#get unique markers having map position
snps<-unlist(snp2keep)
length(snps)
snps<-snps[!duplicated(snps)]
length(snps)

#add a flag to the markerdictionary reporting markers included in EtNAM maps
markerdictionary$EtNAM.map<-0
hit<-which(markerdictionary$marker %in% snps)
markerdictionary$EtNAM.map[hit]<-1
sum(markerdictionary$EtNAM.map)









#make a plotter function to show co-location of QTL across populations

#start by selecting the QTL that you want to see on the plot
#bring everything in Mb
pklist$pos.p<-as.numeric(pklist$pos.p)
pklist$cumpos<-as.numeric(pklist$cumpos)

if(max(pklist$pos.p>1e6)==T){
  pklist$pos.p<-pklist$pos.p*1e-6
  pklist$cumpos<-pklist$cumpos*1e-6
}#if max

#plot for traits of relevance
qtls<-pklist[grep("OA$|DB$|DF$|DH$|DM$|GY$|TGW$|SPL$|NSPKPS$|SPS$|PH$|NTPP$", pklist$lodcolumn),]
#qtls<-pklist[grep("adet$", pklist$lodcolumn),]

qtls

#extract those chromosomes that have at least one QTL hit
chroms<-sort(unique(qtls$chr.p))
chroms

#create a chromosome index
chridx<-data.frame(chr=chroms, idx=1:length(chroms))

#now get the marker database and move everything to Mb
pmrk<-markerdictionary

if(max(pmrk$pos.phy>1e6)==T){
  pmrk$pos.phy<-pmrk$pos.phy*1e-6
  pmrk$cumulative.pos<-pmrk$cumulative.pos*1e-6
  pmrk$chr.start<-pmrk$chr.start*1e-6
  pmrk$chr.stop<-pmrk$chr.stop*1e-6
}#if max

#get maximum length of chr
maxchr<-max(pmrk$pos.phy)

#intialize plot

pdf("qtl.chromosome.map.pdf", height=15, width=10)

    subTMP<-pmrk[which(pmrk$chr.phy == "1A"),]
    subTMPonEtNAM<-subTMP[which(subTMP$EtNAM.map == 1), ]
    plot(x = rep(1, nrow(subTMPonEtNAM)), y = subTMPonEtNAM$pos.phy, 
         type="n",
         xlim = c(0,(length(chroms)+1)), ylim=c(0, maxchr*1.1),
         ylab = "Mbp position",
         xlab = "Chromosome",
         xaxt="n")
    
    #fix the x axis lables
    xtick<-c(1:length(chroms))
    axis(side=1, at=xtick, labels = chroms)
    
    #add chromosome bars
    for (j in 1:length(chroms)){ #j=1
      tmpchr<-chroms[j]
      subTMP<-pmrk[which(pmrk$chr.phy == tmpchr),]
      
      chrstart<-min(subTMP$pos.phy)
      chrstop<-max(subTMP$pos.phy)
      
      segments(y0=chrstart, x0=j, y1=chrstop, x1=j, col="gray")
      
    }
    
    #add all marker ticks
    for (j in 1:length(chroms)){ #j=
      tmpchr<-chroms[j]
      subTMP<-pmrk[which(pmrk$chr.phy == tmpchr),]
      subTMPonEtNAM<-subTMP[which(subTMP$EtNAM.map == 1), ]
      
      points(x = rep(j, nrow(subTMPonEtNAM)), y = subTMPonEtNAM$pos.phy, 
             pch="-", col="gray20", cex=1.5)
    }
    
    #now add those QTL that have been mapped
    qtls
    
    #bulk all OAs into a single phenotype
    oas<-grep("OA",qtls$lodcolumn)
    qtls$lodcolumn[oas]<-"Farmer's appreciation"
    
    #bulk all flowering times into a single phenotype
    phenol<-grep("DB|DF|DH|DM",qtls$lodcolumn)
    qtls$lodcolumn[phenol]<-"Phenology"
    
    #bulk all yield components into a single phenotype
    yld<-grep("GY|TGW|SPL|NSPKPS|SPS|PH|NTPP|BM",qtls$lodcolumn)
    qtls$lodcolumn[yld]<-"Yield components"
    
    #loop through pheonotype
    tmptraits<-unique(qtls$lodcolumn)
    tmptraits
    
    #set parametrs depeding on number of traits
    cols<-c("#648FFF",
            "#FFB000",
            "#DC267F")
    
    params<-data.frame(traits=tmptraits, 
                       offset=c(-0.2, 0, +0.2), 
                       colors=alpha(cols, 0.6),
                       pch= c(18,18,18))
    
    for (t in tmptraits){ #t=tmptraits[1]
      tmp<-qtls[which(qtls$lodcolumn == t),]
      
      #loop through chromosomes
      tmpchrs<-unique(tmp$chr.p)
      
      #add information
      for (z in tmpchrs){ #z=tmpchrs[1]
        
        #get index for plotting
        tmpchridx<-chridx[which(chridx$chr == z), 2]
        
        #extract QTL for the right trait in the right chromosome
        tmpchrqtl<-tmp[which(tmp$chr.p == z),]
          
        #loop through qlts in chromosomes
        for (q in 1:nrow(tmpchrqtl)){
          points(x = tmpchridx + params[which(params$traits==t),2], 
                 y= tmpchrqtl[q,10], 
                 col=params[which(params$traits==t),3], 
                 cex=2.5,
                 pch=params[which(params$traits==t),4])
        }
      }
    }
    
    legend("topleft", legend=params$traits,
           col=params$colors,pch=params$pch, cex=1.5)
    
dev.off()





#genomerange<-c(0, max(markerdictionary$cumulative.pos))
#chrstarts<-unique(markerdictionary$chr.start)



head(markerdictionary)

plot(x=chrstarts, y=rep(1, length(chrstarts)), xlim=genomerange)


ggplot(markerdictionary, aes(x=cumulative.pos, y=1))+geom_point()

