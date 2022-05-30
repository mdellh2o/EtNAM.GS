######################
# This script parses GWAS outputs
# Authors: Matteo Dell'Acqua
# Date: May 30th, 2022
########################

wd<-"C:/Users/admin/OneDrive - Scuola Superiore Sant'Anna/projects/EtNAM.GS/analysis/11.GWAS"
setwd(wd)

filen<-list.files("results/")
filen<-filen[grep("GWAS.Results", filen)]
filen<-filen[grep("Blink", filen)]
filen


#combine all
tmp<-read.csv(paste0("./results/",filen[1]))

tmpname<-sub("^GAPIT.Blink.", "", filen[1])
tmpname<-sub(".GWAS.Results.csv$", "", tmpname)
tmpname


out<-tmp[order(tmp[,1]),]
out<-out[,1:4]
colnames(out)[ncol(out)]<-tmpname

for (i in 2:length(filen)){ #i=1
  tmp<-read.csv(paste0("./results/",filen[i]))
  
  tmpname<-sub("^GAPIT.Blink.", "", filen[i])
  tmpname<-sub(".GWAS.Results.csv$", "", tmpname)
  tmpname
  
  
  tmp<-tmp[order(tmp[,1]),]
  out<-cbind(out, tmp$P.value)
  colnames(out)[ncol(out)]<-tmpname
  
}

out<-out[order(out[,2], out[,3]),]
head(out)

write.csv(out, file="supplementary.table.gwas.pvalues.csv", quote=F, row.names=F)

#get threshold
thr<-0.05/nrow(out)
-log10(thr)

#make a parsing table
out1<-out
class(out1)

for (i in 4:ncol(out1)){ #i=1
  hit<-which(out1[,i]<=thr)
  not<-which(out1[,i]>thr)
  if(length(hit)>0){
    out1[hit,i]<-1
    out1[not,i]<-0
  } else {
    out1[,i]<-0
  }
}

head(out1)


mtas<-apply(out1[,4:ncol(out1)],1,sum)
mtas
range(mtas)

tokeep<-which(mtas>0)
out2<-out1[tokeep,]
tail(out2)

write.csv(out2, file="MTA.hits.table.csv", quote=F, row.names=F)
