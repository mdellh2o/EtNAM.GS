######################
# This scripts takes in phenotypic data and creates BLUPs and H2 estimates, by trait
# Authors: Matteo Dell'Acqua, Kesse Poland
# Date: May 30th, 2022
########################

options(stringsAsFactors = F)

library(asreml)
library(ggplot2)
library(tidyr)
library(metan)
library(ggcorrplot)

wd<-"C:/Users/admin/OneDrive - Scuola Superiore Sant'Anna/projects/EtNAM.GS/analysis/1.calculate.BLUPs"
setwd(wd)

#get in EtNAM phenotypic data
load("../0.check.phenos/phenotypic.data.EtNAM.Rdata")

## state functions to derive BLUPs and variances
getBLUPs = function(m){
  res = data.frame(coef(m)$rand)
  res = res[grepl("ID",row.names(res)) & !grepl(":",row.names(res)) ,  , drop=FALSE]
  ##res = res[grepl("GY:ID_",row.names(res)), , drop=FALSE]
  return(res)
}

getVariance = function(asr, comp){
  var = summary(asr)$varcomp
  idx = which(rownames(var)==comp)
  v = var$component[idx]
  print(paste("variance component", v))
  return(v)
}

######################
## fix trait dataframes

## metric first
## check data classes and assign factors
head(metN)
sapply(metN, class)

f =c("LOCATION", "REP", "COL", "ROW", "ID")
head(metN[f])

metN[f] = lapply(metN[f], as.factor)
sapply(metN, class)

metN[!colnames(metN) %in% f] = lapply(metN[!colnames(metN) %in% f], as.numeric)
sapply(metN, class)

head(metN)

## farmer traits now

## collapse farmer data to key-value pairs
head(farmN)
f.names = colnames(farmN)[-c(1:5)] ## get names of farmer traits

## transform data to format for mixed model ##
farm2 = gather(farmN, key="farmer", value = 'OA', f.names)
head(farm2)

farm2 = separate(farm2, col='farmer', into=c('GENDER', 'FARMER'), sep="_")
head(farm2)

sapply(farm2, class) ## check if factors

n = colnames(farm2)[1:7]
farm2[n] = lapply(farm2[n], as.factor) ## change to factors
sapply(farm2, class)

traits = colnames(farm2)[8]
farm2[traits] = lapply(farm2[traits], as.numeric) ## change phenotypes to numeric
sapply(farm2, class)

head(farm2)

##############
## make BLUPs for metric data
##############

## get list of traits and number of trials
traits = colnames(metN)[-c(1:5)] 
nLoc = 3
nRep = 2

## set up dataframe for holding BLUPs and h2
blups.met.nam = data.frame(ID = unique(metN$ID), row.names = unique(metN$ID))
head(blups.met.nam)
h2.met.nam = data.frame()

## loop through traits with model for combined year, within year, within location, within location & year
for(t in traits){  ## t = traits[1]
  ## model for across year analysis  
  asr= asreml(fixed = as.formula(paste(t,'~ 1')), random= ~ id(ID) + 
                  id(LOCATION) + 
                  id(ID):id(LOCATION) , data=metN, maxit=60)#, G.param = iv)

  print(summary(asr)$varcomp)
  
  b = predict(asr, "ID")
  b = b$pvals
  b.tmp = data.frame(b$ID, b$predicted.value)
  names(b.tmp) = c("ID", t)
  head(b.tmp)
  
  blups.met.nam = merge(blups.met.nam, b.tmp, by="ID")
  head(blups.met.nam)
  
  Vg = getVariance(asr, 'ID')
  Vl = getVariance(asr, 'LOCATION')
  Ve = getVariance(asr, 'units!R')
  
  h2 = Vg/(Vg + Vl/nLoc + Ve/nRep)  ## NOT SURE, NEED TO CHECK THIS CALCULATION
  print(paste("H2 for", t, ":", h2))
  
  h2.met.nam = rbind(h2.met.nam, data.frame(trait=t, location="ALL", h2=h2))
  
  for(l in levels(metN$LOCATION)){  ## l='adet'
    
    writeLines("\n\n"); print(paste("RUNNING TRAIT:", t, 'for', l))
    
    asr=asreml(fixed = as.formula(paste(t,'~ 1')), random= ~ id(ID) + 
                  id(REP), data=metN, subset= LOCATION==l, maxit=60) ##, G.param = iv)
    print(summary(asr)$varcomp)
    
    name = paste(t, l, sep=".")

    b = predict(asr, "ID")
    b = b$pvals
    b.tmp = data.frame(b$ID, b$predicted.value)
    names(b.tmp) = c("ID", name)
    head(b.tmp)
    
    blups.met.nam = merge(blups.met.nam, b.tmp, by="ID")
    head(blups.met.nam)
    
    Vg = getVariance(asr, 'ID')
    Ve = getVariance(asr, 'units!R')
    
    h2 = Vg/(Vg + Ve/nRep)  ## NOT SURE, NEED TO CHECK THIS CALCULATION
    print(paste("H2 for", t, "in", l, ":", h2))
    
    h2.met.nam = rbind(h2.met.nam, data.frame(trait=t, location=l, h2=h2))
    
  } # by location   
  
} #across locations
head(blups.met.nam)
head(h2.met.nam)

##############
## make BLUPs for farmer data
##############

## get list of traits and number of trials
traits = colnames(farm2)[-c(1:7)] 
nLoc = 3
nRep = 2

## set up dataframe for holding BLUPs and h2
blups.farm.nam = data.frame(ID = unique(metN$ID), row.names = unique(metN$ID))
head(blups.farm.nam)
h2.farm.nam = data.frame()

## loop through all traits
for(t in traits){  ## t = traits[1]
  ## model for across year analysis  
  asr=asreml(fixed = as.formula(paste(t,'~ 1')), random= ~ id(ID) + 
                id(LOCATION) +
                id(GENDER) +
                id(FARMER) +
                id(ID):id(LOCATION) +
                id(ID):id(GENDER) +
                id(ID):id(FARMER) +
                id(GENDER):id(LOCATION) +
                id(ID):id(FARMER):id(GENDER):id(LOCATION), data=farm2, maxit=60) ##G.param=iv)

  print(summary(asr)$varcomp)
  
  Vg = getVariance(asr, 'ID')
  Vgl = getVariance(asr, 'ID:LOCATION')
  Vgm = getVariance(asr, 'ID:GENDER')
  Vgf = getVariance(asr, 'ID:FARMER')
  Ve = getVariance(asr, 'units!R')
  
  h2 = Vg/(Vg + Vgl/nLoc + Vgm/2 + Vgf/10 + Ve/(nRep*nLoc*10*2))  ## DON'T KNOW WHAT DO DO WITH FARMER AND GROUP VARIANCE ! .... NEED TO CHECK THIS CALCULATION
  print(paste("H2 for", t, ":", h2))
  
  pred = predict(asr, 'ID')
  a = pred$avsed
  
  h2.farm.nam = rbind(h2.farm.nam, data.frame(trait=t, gender='both', location="ALL", h2=h2, avsed=a))
  
  b = predict(asr, "ID")
  b = b$pvals
  b.tmp = data.frame(b$ID, b$predicted.value)
  names(b.tmp) = c("ID", t)
  head(b.tmp)
  
  blups.farm.nam = merge(blups.farm.nam, b.tmp, by="ID")
  head(blups.farm.nam)
  
  for(l in levels(farm2$LOCATION)){  ## l="adet"
    
    asr=asreml(fixed = as.formula(paste(t,'~ 1')), random= ~ id(ID) + 
                 id(GENDER) +  
                 id(FARMER) +
                 id(ID):id(GENDER) + 
                 id(ID):id(FARMER) + 
                 id(ID):id(FARMER):id(GENDER), data=farm2, subset= LOCATION==l, maxit=60) ##, G.param = iv)
    
    print(summary(asr)$varcomp)
    
    Vg = getVariance(asr, 'ID')
    Vgm = getVariance(asr, 'ID:GENDER')
    Vgf = getVariance(asr, 'ID:FARMER')
    Ve = getVariance(asr, 'units!R')
    
    h2 = Vg/(Vg + Vgm/2 + Vgf/10 + Ve/(nRep*2)) ## NOT SURE, NEED TO CHECK THIS CALCULATION !!!
    print(paste("H2 for", t, 'in', l, ":", h2))
    
    pred = predict(asr, 'ID')
    a = pred$avsed
    
    h2.farm.nam = rbind(h2.farm.nam, data.frame(trait=t, gender='both', location=l, h2=h2, avsed=a))
    
    ## blups 
    name = paste(t, l, sep=".")
    
    b = predict(asr, "ID")
    b = b$pvals
    b.tmp = data.frame(b$ID, b$predicted.value)
    names(b.tmp) = c("ID", name)
    head(b.tmp)
    
    blups.farm.nam = merge(blups.farm.nam, b.tmp, by="ID")
    head(blups.farm.nam)
    
  }   
  
  for(g in levels(farm2$GENDER)){  ## l='hagreselam', g='M'
    writeLines("\n\n"); print(paste("RUNNING TRAIT:", t, 'for', g))
    
    asr=asreml(fixed = as.formula(paste(t,'~ 1')), random= ~ id(ID) + 
                 id(LOCATION) +
                 id(FARMER) +
                 id(ID):id(LOCATION) +
                 id(ID):id(FARMER) +
                 id(ID):id(FARMER):id(LOCATION) , data=farm2, subset=GENDER==g, maxit=60) ##, G.param = iv)
    
    print(summary(asr)$varcomp)
    
    Vg = getVariance(asr, 'ID')
    Vgl = getVariance(asr, 'ID:LOCATION')
    Vgf = getVariance(asr, 'ID:FARMER')
    Ve = getVariance(asr, 'units!R')
    
    h2 = Vg/(Vg + Vgl/nLoc + Vgf/5 + Ve/(nRep*nLoc))
    print(paste("H2 for", t, 'for', g, ":", h2))
    
    pred = predict(asr, 'ID')
    a = pred$avsed
    
    h2.farm.nam = rbind(h2.farm.nam, data.frame(trait=t, gender=g, location='ALL', h2=h2, avsed=a))
    
    ## blups
    name = paste(t, g, sep=".")
    
    b = predict(asr, "ID")
    b = b$pvals
    b.tmp = data.frame(b$ID, b$predicted.value)
    names(b.tmp) = c("ID", name)
    head(b.tmp)
    
    blups.farm.nam = merge(blups.farm.nam, b.tmp, by="ID")
    head(blups.farm.nam)
    
  }     
  
  for(l in levels(farm2$LOCATION)){  ## l='adet', g='M'
    for(g in levels(farm2$GENDER)){
      
      writeLines("\n\n"); print(paste("RUNNING TRAIT:", t, 'for', g, 'in', l))
      
      asr=asreml(fixed = as.formula(paste(t,'~ 1')), random= ~ id(ID) + 
                    id(FARMER) + 
                    id(ID):id(FARMER), data=farm2, subset= GENDER==g & LOCATION==l, maxit=60) ##, G.param = iv)
      
      print(summary(asr)$varcomp)
      
      Vg = getVariance(asr, 'ID')
      Vgf = getVariance(asr, 'ID:FARMER')
      Ve = getVariance(asr, 'units!R')
      
      h2 = Vg/(Vg + Vgf/5 + Ve/(nRep))
      print(paste("H2 for", t, 'for', g, 'in', l,  ":", h2))
      
      pred = predict(asr, 'ID')
      a = pred$avsed
      
      h2.farm.nam = rbind(h2.farm.nam, data.frame(trait=t, gender=g, location=l, h2=h2, avsed=a))
      
      
      ## blups
      name = paste(t, l, g, sep=".")
      
      b = predict(asr, "ID")
      b = b$pvals
      b.tmp = data.frame(b$ID, b$predicted.value)
      names(b.tmp) = c("ID", name)
      head(b.tmp)
      
      blups.farm.nam = merge(blups.farm.nam, b.tmp, by="ID")
      head(blups.farm.nam)
      
    }
  }  
}  

########
#check correlations
cr<-cor(blups.met.nam[,2:ncol(blups.met.nam)], blups.farm.nam[,2:ncol(blups.farm.nam)], use = "complete.obs")
ggcorrplot(cr)


##################
# save all relevant results
save(blups.farm.nam, h2.farm.nam, blups.met.nam, h2.met.nam, file="BLUP.values.EtNAM.Rdata")

#############
## derive Shukla's stability parameter (Shukla, G.K. 1972)
#############
#get farmers'OA
bf<-blups.farm.nam[,1:2]
bm<-blups.met.nam[,c("ID", "DF", "SPS", "NSPKPS", "PH", "TGW", "GY", "SPL")]
rownames(bf)<-bf[,1]
rownames(bm)<-bm[,1]
phenos<-cbind(bf,bm)
phenos<-phenos[,-c(1,3)]

head(phenos)

#correlate OA with all the rest
cor(phenos[,1], phenos[,2:ncol(phenos)])

#impute phenotypes with mean value
forsh<-metN
#impute missing data
for(i in 6:ncol(forsh)){
  forsh[is.na(forsh[,i]),i]<-mean(forsh[,i], na.rm=T)
}#for i

ge_plot(forsh, env=LOCATION, gen=ID, resp = GY, type=1)

#Calculate multi-trait stability as in Olivoto et al 2019
#This index works for stability only
MTSI <- forsh %>%
        waasb(LOCATION, ID, REP,
        resp = c(DF, SPS, NSPKPS, PH, TGW, GY, SPL)) %>%
        mtsi(verbose = FALSE, index = "waasb")
#get MTSI score on which the selection is based
MTSIidx<-MTSI$MTSI

#This index works for stability AND performance
MTSI_perf <- forsh %>%
        waasb(LOCATION, ID, REP,
        resp = c(DF, SPS, NSPKPS, PH, TGW, GY, SPL),
        mresp=c(0,100,0,100,100,100,0),  #use values coming from correlations
        wresp = c(50, 70, 50, 60, 70, 80, 50)) %>%
        mtsi()

#get MTSI score on which the selection is based
MTSIperf<-MTSI_perf$MTSI
colnames(MTSIperf)[2]<-"MTSIperf"

#now get stability on GY alone
allstat<-ge_stats(forsh, env=LOCATION, gen=ID, rep=REP, resp = c(GY), verbose = TRUE)
moddat<-get_model_data(allstat)

#combine everything in a stability dataframe
stability<-merge(moddat, MTSIidx, by.x = "gen", by.y="Genotype", all=T)
stability<-merge(stability, MTSIperf, by.x = "gen", by.y="Genotype", all=T)
rownames(stability)<-stability[,1]
stability<-stability[,-c(1:2)]

#remove all stability indexes with lots of NAs
nas<-apply(stability, 2, function(x) length(which(is.na(x)==T)))
nas

stability<-stability[,-which(nas > 500)]

####
save(stability, file="genoytpe.stability.analysis.Rdata")
