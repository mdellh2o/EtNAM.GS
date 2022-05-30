######################
# This scripts takes genomic selection outputs, parses them and produces most relevant plots
# Authors: Matteo Dell'Acqua
# Date: May 30th, 2022
########################

options(stringsAsFactors = F)

wd<-"C:/Users/admin/OneDrive - Scuola Superiore Sant'Anna/projects/EtNAM.GS/analysis/5.plot.GS"
setwd(wd)

#setting up a plotting function to be used in GS plotting
plotGP<-function(gpout, idx, predictors, predicted, cols=NA, gpname="", ymin=NULL, ymax=NULL){
  gpout<-gpout
  idx<-idx
  gp<-gpout[[idx]]
  
  gp<-gp[which(names(gp) %in% predictors)]
  
  for (i in 1:length(gp)){
    tmp<-gp[[i]]
    tmp<-tmp[rownames(tmp) %in% predicted, ]
    tmpout<-data.frame(Accuracy=apply(tmp, 1, mean), SE=apply(tmp, 1, function(x) sd(x)/sqrt(length(x))))
    tmpoutdf<-data.frame(predictor=names(gp)[i], predicted=rownames(tmpout), tmpout)
    gp[[i]]<-tmpoutdf
  }#for i
  
  #assemble a new dataset to plot
  out<-do.call(rbind, gp)
  rownames(out)<-NULL
  out
  #get mean and SE
  
  #write input tables to file
  write.csv(out, file=paste(names(gpout)[idx], paste(predictors, collapse="."), "csv", sep="."), quote=F, row.names=F)
  
  if(is.na(gpname)==T){
    gpname<-names(gpout)[idx]
  }
  
  if(length(which(is.na(cols)))>0){
    cols<-gray.colors(length(predictors))
  }
  
  #get to plotting
  ggplot(out, aes(x=predicted, y=Accuracy, fill=predictor)) + 
    geom_bar(stat="identity", position = position_dodge(), colour="black") + 
    geom_errorbar(aes(ymin = Accuracy-SE, ymax = Accuracy+SE), width=0.4,
                  position=position_dodge(.9)) +
    scale_fill_manual(values=cols)+
    labs(fill = "Predictor", x="Predictedm") +
    theme(panel.background = element_rect(fill = 'white', colour = 'black'),axis.text.x = element_text(angle = 45, hjust = 1)) +
    coord_cartesian(ylim = c(ymin, ymax))
}### end of function

#include function to reduce legend size
addSmallLegend <- function(myPlot, pointSize = 0.9, textSize = 9, spaceLegend = 0.8) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_text(size = textSize), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
}

library(ggplot2)
library(plyr)
library(reshape2)
library(ggplot2)
library(cowplot)
library(tidyr)
library(ggridges)
library(patchwork)
library(viridis)
library(hrbrthemes)
library(ggcorrplot)
library(sjPlot)
library(vcd)

#####################
# set colors by locations
cola<-"#697ed5"
colg<-"#bc7d39"
colk<-"#b94663"

# set colors by gender
colmen<-"#6fac5d"
colwom<-"#9350a1"

# set colors by traits
colea<-"#b5002d"
coloa<-"#E6E6E6"
colsp<-"#9a5cd2"
coltl<-"#ffc65f"

#load all needed files
#get phenotypes
load("../1.calculate.BLUPs/BLUP.values.EtNAM.Rdata")
#lickert scale plots by farmer
load("../0.check.phenos/phenotypic.data.EtNAM.Rdata")
## get in genotypes
load("../2.impute.snps/genotypic.data.nam.panel.Rdata")
#get files from previous GS iteration
load("prediction.scenarios.by.trait.Rdata")


#######################
# PHENOTYPES PLOT

met<-blups.met.nam
met[,1]<-sub("^ID_", "", met[,1])
rownames(met)<-met[,1]

farm<-blups.farm.nam
farm[,1]<-sub("^ID_", "", farm[,1])
rownames(farm)<-farm[,1]

# get to long format 
oa<-gather(farm, factor, value, ID:OA.kulumsa.M, factor_key=TRUE)
id<-subset(oa, factor == "ID")
oa<-subset(oa, factor != "ID")
oa$location<-"Combined"
oa[grep("adet", oa$factor), "location"]<-"Adet"
oa[grep("geregera", oa$factor), "location"]<-"Geregera"
oa[grep("kulumsa", oa$factor), "location"]<-"Kulumsa"
oa$gender<-"Combined"
oa[grep(".M$", oa$factor), "gender"]<-"Men"
oa[grep(".F$", oa$factor), "gender"]<-"Women"
oa$factor<-"OA"
oa<-oa[,c("factor", "location", "gender", "value")]
oa<-cbind(id$value, oa)
oa[,1]<-as.factor(oa[,1])
oa[,2]<-as.factor(oa[,2])
oa[,3]<-as.factor(oa[,3])
oa[,4]<-as.factor(oa[,4])
oa$value<-as.numeric(oa$value)
str(oa) 

#add interaction of gender and location
oa$genderloc<-interaction(oa$gender, oa$location)

#add info about type of materials
oa$type<-"MV"
oa[grep("^[0-9]", oa[,1]), "type"]<-"EtNAM RIL"
oa[grep("^RF$", oa[,1]), "type"]<-"Recurrent parental"
oa[grep("^N[0-9]", oa[,1]), "type"]<-"Landrace parental"

#work with metric traits
m<-gather(met, factor, value, ID:GY.kulumsa, factor_key=TRUE)
id<-subset(m, factor == "ID")
m<-subset(m, factor != "ID")
m$location<-"Combined"
m[grep("adet", m$factor), "location"]<-"Adet"
m[grep("geregera", m$factor), "location"]<-"Geregera"
m[grep("kulumsa", m$factor), "location"]<-"Kulumsa"
m$factor<-sub("\\..*$", "", m$factor)
m<-data.frame(m[,c("factor", "location", "value")])
m<-cbind(id$value, m)
m[,1]<-as.factor(m[,1])
m[,2]<-as.factor(m[,2])
m[,3]<-as.factor(m[,3])
m[,4]<-as.numeric(m[,4])
str(m)

#add info about type of materials
m$type<-"MV"
m[grep("^[0-9]", m[,1]), "type"]<-"EtNAM RIL"
m[grep("^RF$", m[,1]), "type"]<-"Recurrent parental"
m[grep("^N[0-9]", m[,1]), "type"]<-"Landrace parental"

#########
#reshape datasets to make the easily manageble
men<-subset(oa, gender == "Men")
women<-subset(oa, gender == "Women")
comb<-subset(oa, gender == "Combined")
stopifnot(all(men[,1]==women[,1]))
stopifnot(all(men[,1]==comb[,1]))
toplot<-data.frame(ID= men$`id$value`, type=men$type, location=men$location, Men = men$value, Women = women$value, Combined = comb$value)
#add metric traits
bytrait<-split(m, m$factor)

for (t in 1:length(bytrait)){
  tmp<-bytrait[[t]]
  #check that everything is as expected
  stopifnot(all(tmp[,1]==toplot[,1]))
  stopifnot(all(tmp[,3]==toplot[,3]))
  
  #add to the plotting dataframe
  toplot<-cbind(toplot, tmp$value)
  colnames(toplot)[ncol(toplot)]<-names(bytrait)[t]
  
} #for t in bytrait

toplot[,1:3]<-apply(toplot[,1:3], 2, as.factor)
toplot[,4:ncol(toplot)]<-apply(toplot[,4:ncol(toplot)], 2, as.numeric)
str(toplot)

#reorder factor levels
toplot$location <- factor(toplot$location, levels = c("Adet", "Geregera", "Kulumsa", "Combined"))
combpheno<-subset(toplot, location=="Combined")

##################
#make plots
##

#################  FIGURE 1 DEPRECATED #####################
#average farmer scores across replicas
farmNl<-split(farmN, farmN[,1])
byfarm<-list()

#set factor levels
f<-c("W1", "W2", "W3", "W4", "W5", "M1", "M2", "M3", "M4", "M5")
locnames<-c("Adet", "Geregera", "Kulumsa")

for (i in 1:length(farmNl)){ #i=1
  tmp<-farmNl[[i]]
  r1<-tmp[which(tmp[,2]==1),]
  r2<-tmp[which(tmp[,2]==2),]
  stopifnot(all(r1[,"ID"]==r2[,"ID"]))
  tmpout<-data.frame(W1=rowMeans(cbind(r1[,6], r2[,6]), na.rm=T),
                     W2=rowMeans(cbind(r1[,7], r2[,7]), na.rm=T),
                     W3=rowMeans(cbind(r1[,8], r2[,8]), na.rm=T),
                     W4=rowMeans(cbind(r1[,9], r2[,9]), na.rm=T),
                     W5=rowMeans(cbind(r1[,10], r2[,10]), na.rm=T),
                     M1=rowMeans(cbind(r1[,11], r2[,11]), na.rm=T),
                     M2=rowMeans(cbind(r1[,12], r2[,12]), na.rm=T),
                     M3=rowMeans(cbind(r1[,13], r2[,13]), na.rm=T),
                     M4=rowMeans(cbind(r1[,14], r2[,14]), na.rm=T),
                     M5=rowMeans(cbind(r1[,15], r2[,15]), na.rm=T))
  rownames(tmpout)<-r1[,"ID"]
  tmpout<-data.frame(apply(tmpout, 2, round))
  head(tmpout)
  tmpout[,f] <- lapply(tmpout[,f], factor)
  head(tmpout)
  byfarm[[i]]<-tmpout
}
names(byfarm)<-locnames

#make a table with all farmer scores
allPVS<-do.call(cbind, byfarm)
xq<-sapply(allPVS, table)
xqdf<-do.call(cbind, xq)

write.table(xqdf, file="./figures/summary.lickert.txt", sep="\t", quote=F, row.names=T)

#make a complex plot
assoc(xqdf[1:5,], shade = TRUE, las=1, labeling= labeling_border(rot_labels = c(90,0,0,0)),  main = NULL, xlab = NULL, ylab = NULL)

lkad<-plot_likert(byfarm[["Adet"]], catcount=5, expand.grid=TRUE,
            reverse.scale=TRUE, show.n=F, values="hide", geom.colors = "RdBu",
            legend.labels = c("Very poor","Poor", "Average","Good", "Very good"), 
            legend.title="Scores", title="Adet", show.legend=F)
lkad<-lkad+theme_bw()+labs(tag="a")

lkge<-plot_likert(byfarm[["Geregera"]], catcount=5, expand.grid=TRUE,
                  reverse.scale=TRUE, show.n=F, values="hide", geom.colors = "RdBu",
                  legend.labels = c("Very poor","Poor", "Average","Good", "Very good"), 
                  legend.title="Scores", title="Geregera", show.legend=F)
lkge<-lkge+theme_bw()

lkku<-plot_likert(byfarm[["Kulumsa"]], catcount=5, expand.grid=TRUE,
                  reverse.scale=TRUE, show.n=F, values="hide", geom.colors = "RdBu",
                  legend.labels = c("Very poor","Poor", "Average","Good", "Very good"), 
                  legend.title="Scores", title="Kulumsa", show.legend=T)
lkku<-lkku+theme_bw()+theme(legend.position = "bottom", legend.title = element_text( size=8), legend.text=element_text(size=8))

lickertplot<-lkad / lkge / lkku
lickertplot

#XY plot showing covariance of men and women
scat<-ggplot(subset(toplot, location != "Combined"), aes(x=Men, y=Women, color=location)) + 
      geom_point(alpha=0.3) +
      xlim(1,5) +
      ylim(1,5) +
      #geom_density_2d() +
      #facet_wrap(.~location) +
      scale_color_manual(values=c(cola,colg, colk)) +
      geom_smooth(method=lm, se=FALSE, fullrange=TRUE) +
      coord_fixed() +
      theme_bw()+
      theme(legend.position = "none")+
      labs(fill="Gender")

scat<-scat+labs(tag="b")

#make a ridge chart
rcharoa<-ggplot(subset(oa, gender != "Combined" & location != "Combined"), aes(x = value, y = genderloc, fill=gender)) +
  geom_density_ridges()+
  scale_fill_manual(values=c("Men"=colmen, "Women"=colwom))+
  theme_bw() +
  scale_y_discrete(labels= c("", "Adet", "", "Geregera", "", "Kulumsa")) +
  theme(axis.ticks.y = element_blank())+
  theme(legend.position="top", legend.title = element_text( size=8), legend.text=element_text(size=8)) +
  labs(x ="OA", y = "", fill="Gender", tag="c")
rcharoa 

rchargy<-ggplot(subset(m, location != "Combined" & factor =="GY"), aes(x = value, y = location, fill=location)) +
  geom_density_ridges() +
  scale_fill_manual(values=c("Geregera"=colg, "Adet"=cola, "Kulumsa"=colk))+
  theme_bw() +
  theme(legend.position="top", legend.title = element_text(size=8), legend.text=element_text(size=8)) +
  labs(x ="GY", y = "", fill="Location")
rchargy 

ridgeplot<-rcharoa | rchargy
ridgeplot

#assemble FIGURE 1 plot
fig1<- lickertplot | scat / ridgeplot 
fig1

ggsave(fig1, file="./figures/fig1.deprecated.pdf", width = 11, height = 10)

################# FIGURE 1 ENDS  #####################

##check gender differences in locations
tmpad<-subset(toplot, location=="Adet")
t.test(tmpad[,"Men"], tmpad[,"Women"], paired = T)

tmpge<-subset(toplot, location=="Geregera")
t.test(tmpge[,"Men"], tmpge[,"Women"], paired = T)

tmpku<-subset(toplot, location=="Kulumsa")
t.test(tmpku[,"Men"], tmpku[,"Women"], paired = T)

################


################# FIGURE 2 DEPRECATED #####################
#########

#make  plots about type of materials
bxpltoa<-ggplot(subset(toplot, type != "Recurrent parental"), aes(y=Combined, x=type, fill=type)) + 
  geom_boxplot(width=1, size=0.2) +
  scale_fill_viridis(discrete=TRUE)+
  coord_flip()+
 #ylim(1,5)+
  facet_wrap(.~location) +
 #theme_ipsum() +
  theme_bw() +
  theme(legend.position="none")+
  labs(fill="Type", y="OA", x="", tag="a")

bxpltoa

#now yield
bxpltgy<-ggplot(subset(toplot, type != "Recurrent parental"), aes(y=GY, x=type, fill=type)) + 
  geom_boxplot(width=1, size=0.2) +
  scale_fill_viridis(discrete=TRUE)+
  coord_flip()+
  facet_wrap(.~location) +
  # theme_ipsum() +
  theme_bw() +
  theme(legend.position="none")+
  labs(fill="Type", y="GY", x="")
bxpltgy

types<-bxpltoa / bxpltgy
types

#make PCA based on molecular diversity
dim(imputed.nam)
pca<-prcomp(imputed.nam, scale=T, center=T)
pcsc<-pca$x[,1:3]
rownames(pcsc)<-colnames(genonam)
#get variaces
eigs<-pca$sdev^2
pcvars<-(eigs/sum(eigs))[1:3]

pctoplot<-merge(combpheno, pcsc, by.x="ID",by.y="row.names")
head(pctoplot)

#get subsets based on quantiles for plotting
subm<-pctoplot[which(pctoplot$Men> quantile(pctoplot$Men, 0.95)),]
dim(subm)

subf<-pctoplot[which(pctoplot$Women> quantile(pctoplot$Women, 0.95)),]
dim(subf)

#plot it!
pcplot<-ggplot(pctoplot, aes(x=PC1, y=PC2))+
  geom_point(alpha=0.9,size=1.5) + theme_bw() + theme(aspect.ratio = 1, legend.title = element_blank())
pcplot<-pcplot + geom_point(data=pctoplot, aes(x=PC1, y=PC2), alpha=0.9,size=1.5, col = "gray")
pcplot<-pcplot + geom_point(data=subm, aes(x=PC1, y=PC2), col = colmen, alpha=0.3, size=2)
pcplot<-pcplot + geom_point(data=subf, aes(x=PC1, y=PC2), col = colwom, alpha=0.3, size=2)
pcplot<-pcplot  + labs(tag= "b", title="Molecular diversity", x=paste0("PC1 (",round(pcvars[1]*100,2), "%)" ) , y=paste0("PC2 (",round(pcvars[2]*100,2), "%)" ))
pcplot 

#now pca of phneotypes
pcaph<-prcomp(combpheno[,7:ncol(combpheno)], scale=T, center=T)
pcscph<-pcaph$x[,1:3]
#get variaces
eigs<-pcaph$sdev^2
pcphvars<-(eigs/sum(eigs))[1:3]

rownames(pcscph)<-combpheno[,"ID"]
pcphtoplot<-merge(combpheno, pcscph, by.x="ID",by.y="row.names")
head(pcphtoplot)

#get subsets based on quantiles for plotting
submph<-pcphtoplot[which(pcphtoplot$Men> quantile(pcphtoplot$Men, 0.95)),]
dim(submph)

subfph<-pcphtoplot[which(pcphtoplot$Women> quantile(pcphtoplot$Women, 0.95)),]
dim(subfph)

pcphplot<-ggplot(pcphtoplot, aes(x=PC1, y=PC2, col=location))+
  geom_point(alpha=0.9,size=1.5) + theme_bw() + theme(aspect.ratio = 1, legend.title = element_blank(), legend.position = "bottom", legend.text=element_text(size=8)) + 
  scale_color_manual(values = c("gray", colmen, colwom),
                     limits = c("1", "2", "3"), 
                     labels = c("Wheat genotype", "Men choice", "Women choice"))
pcphplot<-pcphplot + geom_point(data=pcphtoplot, aes(x=PC1, y=PC2), alpha=0.9,size=1.5, col = "gray")
pcphplot<-pcphplot + geom_point(data=submph, aes(x=PC1, y=PC2), col = colmen, alpha=0.3, size=2)
pcphplot<-pcphplot + geom_point(data=subfph, aes(x=PC1, y=PC2), col = colwom, alpha=0.3, size=2)
pcphplot<-pcphplot + labs(title="Phenotypic diversity", x=paste0("PC1 (",round(pcphvars[1]*100,2), "%)" ) , y=paste0("PC2 (",round(pcphvars[2]*100,2), "%)" ))

pcphplot 

pcplots<-pcplot / pcphplot
####make FIGURE 2

#fig2<-corplot / (types + PCAplot) + plot_layout(ncol = 1, heights = c(1,1,2))

fig2<-types | pcplots #+  plot_layout(heights = c(1,2))
#Üfig2

ggsave(fig2, file="./figures/fig2.deprecated.pdf", width = 10, height = 10)

################# FIGURE 2 ENDS  #####################


################ UPDATED FIGURE 1 ######################
rcharoa <-ggplot(subset(oa, gender != "Combined" & location != "Combined"), aes(x = value, y = genderloc, fill=gender)) +
  geom_density_ridges()+
  scale_fill_manual(values=c("Men"=colmen, "Women"=colwom))+
  theme_bw() +
  scale_y_discrete(labels= c("", "Adet", "", "Geregera", "", "Kulumsa")) +
  theme(axis.ticks.y = element_blank())+
 # theme(legend.title = element_text( size=8), legend.text=element_text(size=8)) +
  labs(x ="OA", y = "", fill="Gender", tag="A")
rcharoa

rchargy <- ggplot(subset(m, location != "Combined" & factor =="GY"), aes(x = value, y = location, fill=location)) +
  geom_density_ridges() +
  scale_fill_manual(values=c("Geregera"=colg, "Adet"=cola, "Kulumsa"=colk))+
  theme_bw() +
 # theme(legend.title = element_text(size=8), legend.text=element_text(size=8)) +
  labs(x ="GY", y = "", fill="Location", tag="B")
rchargy 


scat<- scat<-ggplot(subset(toplot, location != "Combined"), aes(x=Men, y=Women, color=location)) + 
  geom_point(alpha=0.3) +
  xlim(1,5) +
  ylim(1,5) +
  #geom_density_2d() +
  #facet_wrap(.~location) +
  scale_color_manual(values=c(cola,colg, colk)) +
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE) +
  coord_fixed() +
  theme_bw()+
  theme(legend.position = "none")+
  labs(tag="C", fill="Gender")

pcplot <- ggplot(pctoplot, aes(x=PC1, y=PC2))+
    geom_point(alpha=0.9,size=1.5) + theme_bw() + theme(legend.title = element_blank())
    pcplot<-pcplot + geom_point(data=pctoplot, aes(x=PC1, y=PC2), alpha=0.9,size=1.5, col = "gray90")
    pcplot<-pcplot + geom_point(data=subm, aes(x=PC1, y=PC2), col = colmen, alpha=0.4, size=3)
    pcplot<-pcplot + geom_point(data=subf, aes(x=PC1, y=PC2), col = colwom, alpha=0.4, size=3)
    pcplot<-pcplot  + labs(tag= "D", title="Molecular diversity", x=paste0("PC1 (",round(pcvars[1]*100,2), "%)" ) , y=paste0("PC2 (",round(pcvars[2]*100,2), "%)" ))
pcplot 

pcphplot <- ggplot(pcphtoplot, aes(x=PC1, y=PC2, col=location))+
    geom_point(alpha=0.9,size=1.5) + theme_bw() + theme(legend.title = element_blank()) + 
    scale_color_manual(values = c("gray", colmen, colwom),
                         limits = c("1", "2", "3"), 
                         labels = c("Wheat genotype", "Men choice", "Women choice"))
    pcphplot<-pcphplot + geom_point(data=pcphtoplot, aes(x=PC1, y=PC2), alpha=0.9,size=1.5, col = "gray90")
    pcphplot<-pcphplot + geom_point(data=submph, aes(x=PC1, y=PC2), col = colmen, alpha=0.4, size=3)
    pcphplot<-pcphplot + geom_point(data=subfph, aes(x=PC1, y=PC2), col = colwom, alpha=0.4, size=3)
    pcphplot<-pcphplot + labs(tag="E", title="Phenotypic diversity", x=paste0("PC1 (",round(pcphvars[1]*100,2), "%)" ) , y=paste0("PC2 (",round(pcphvars[2]*100,2), "%)" ))
pcphplot

f1new<- (rcharoa | rchargy| scat) / (pcplot | pcphplot)
f1new<- f1new + plot_layout(guides = "collect") 
f1new

ggsave(f1new, file="./figures/fig1.pdf", width = 12, height = 9)
################ END OF UPDATED FIGURE 1 ######################


################# FIGURE 3 #####################

######################################
# GP PLOTS
#
names(listpred) # check the order of prediction files!!!

########
#panel over nam plot
##
idx<-4
names(listpred[[idx]])
rownames(listpred[[idx]][[1]])

 
# GY and OA by location
predictors<-c("OVERALL", "GY")
predicted<-c("GY","GY.geregera", "GY.kulumsa", "GY.adet")
p2<-plotGP(gpout=listpred, idx=idx, predictors=predictors, predicted=predicted)
#p2<-addSmallLegend(p2) + geom_abline(intercept=0, slope=0, linetype="dashed", color="gray80")
p2<- p2 + labs(title = "DP to EtNAM yield", fill="Predictors", x ="Predicted", tag="a") #+ theme(aspect.ratio=0.5)
#fix lables
p2<- p2 + ylim(-0.3,0.3)+ theme(legend.position = "none")  +
  scale_x_discrete(labels=c(expression('GY'[EtNAM]),expression('GY'[NAM]*'Ad'),expression('GY'[NAM]*'Ge'),expression('GY'[NAM]*'Ku'))) +
  scale_fill_manual(labels = c(expression('GY'[DP]), expression('OA'[DP])), values=gray.colors(length(predictors))) + theme(legend.text.align = 0)
p2<- p2 +  geom_rect(aes(xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf),
            fill = "pink", alpha = 0.03)
p2 

# GY and OA by location
predictors<-c("OVERALL", "GY")
predicted<-c("OA","OA.geregera", "OA.kulumsa", "OA.adet")
po<-plotGP(gpout=listpred, idx=idx, predictors=predictors, predicted=predicted)
#po<-addSmallLegend(p2) + geom_abline(intercept=0, slope=0, linetype="dashed", color="gray80")
po<- po + labs(title = "DP to EtNAM appreciation", fill="Predictors", x ="Predicted") #+ theme(aspect.ratio=0.5)
#fix lables
po<- po + ylim(-0.3,0.3) +
  scale_x_discrete(labels=c(expression('OA'[NAM]),expression('OA'[NAM]*'Ad'),expression('OA'[NAM]*'Ge'),expression('OA'[NAM]*'Ku'))) + 
  scale_fill_manual(labels = c(expression('GY'[DP]), expression('OA'[DP])), values=gray.colors(length(predictors))) + theme(legend.text.align = 0)
po<- po +  geom_rect(aes(xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf),
                      fill = "pink", alpha = 0.03)
po

pred<-p2 | po 
pred


#what do farmer predict?
predictors<-c("OVERALL", "EARLINESS", "SPIKE", "TILLER")
predicted<-c("DF", "SPS", "TGW", "GY", "BM", "NSPKPS", "PH", "SPL")
p3<-plotGP(gpout=listpred, idx=idx, predictors=predictors, predicted=predicted)
p3<-addSmallLegend(p3) + geom_abline(intercept=0, slope=0, linetype="dashed", color="gray80")
p3<- p3 + labs(title = "PVS to metric", fill="Predictors", x ="Predicted", tag="b") #+ theme(aspect.ratio=0.5)
#fix lables
p3<- p3 + ylim(-0.3,0.3) + 
  scale_x_discrete(labels=c(expression('BM'[NAM]),expression('DF'[NAM]),
                            expression('GY'[NAM]),expression('NSPKPS'[NAM]),
                            expression('PH'[NAM]),expression('SPL'[NAM]),
                            expression('SPS'[NAM]),expression('TGW'[NAM]))) + 
  scale_fill_manual(labels = c(expression('EA'[DP]), expression('OA'[DP]), expression('SP'[DP]), expression('TL'[DP])), 
                      values=c(colea, coloa, colsp, coltl)) + theme(legend.text.align = 0)
p3

#Men and women prediction of EtNAM over DP
idx<-3
names(listpred[[idx]])
rownames(listpred[[idx]][[1]])


# gender predictions on combined data
predictors<-c("OA.F", "OA.M", "OA", "GY")
predicted<-c("GY","OVERALL")
pgend<-plotGP(gpout=listpred, idx=idx, predictors=predictors, predicted=predicted)
#p2<-addSmallLegend(p2) + geom_abline(intercept=0, slope=0, linetype="dashed", color="gray80")
pgend<- pgend + labs(title = "EtNAM to DP", fill="Predictors", x ="Predicted", tag="c") #+ theme(aspect.ratio=0.5)
#fix lables
pgend<- pgend +
  scale_x_discrete(labels=c(expression('GY'[DP]),expression('OA'[DP]))) +
  scale_fill_manual(labels = c(expression('GY'[NAM]), expression('OA'[NAM]), 
                               expression('OA'[NAM]*'W'), expression('OA'[NAM]*'M')), 
                    values=c(gray.colors(2),colwom,colmen)) + theme(legend.text.align = 0)
pgend 


#assemble Figure 3
fig3<- (pred / p3)| pgend
fig3

ggsave(fig3, file="./figures/fig3.deprecated.pdf", width = 16, height = 10)

#measure of stability compared with GY -> do OA prefer stable genotypes over high yielding?

################# FIGURE 3 ENDS #####################


################# UPDATED FIGURE 3 ####################
########
#panel over nam plot
##
idx<-4
names(listpred[[idx]])
rownames(listpred[[idx]][[1]])

# GY and OA by location
predictors<-c("OVERALL", "GY")
predicted<-c("GY","GY.geregera", "GY.kulumsa", "GY.adet")
p2<-plotGP(gpout=listpred, idx=idx, predictors=predictors, predicted=predicted)
#p2<-addSmallLegend(p2) + geom_abline(intercept=0, slope=0, linetype="dashed", color="gray80")
p2<- p2 + labs(title = "DP to EtNAM: grain yield", fill="Predictors", x ="Predicted", tag="A") #+ theme(aspect.ratio=0.5)
#fix lables
p2<- p2 + ylim(-0.3,0.3)+ #theme(legend.position = "none")  +
  scale_x_discrete(labels=c(expression('GY'[EtNAM]),expression('GY'[EtNAM]*'Ad'),expression('GY'[EtNAM]*'Ge'),expression('GY'[EtNAM]*'Ku'))) +
  scale_fill_manual(labels = c(expression('GY'[DP]), expression('OA'[DP])), values=gray.colors(length(predictors))) + theme(legend.text.align = 0)
p2<- p2 +  geom_rect(aes(xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf),
                     fill = "pink", alpha = 0.03)
p2 

# GY and OA by location
predictors<-c("OVERALL", "GY")
predicted<-c("OA","OA.geregera", "OA.kulumsa", "OA.adet")
po<-plotGP(gpout=listpred, idx=idx, predictors=predictors, predicted=predicted)
#po<-addSmallLegend(p2) + geom_abline(intercept=0, slope=0, linetype="dashed", color="gray80")
po<- po + labs(title = "DP to EtNAM: farmer appreciation", fill="Predictors", x ="Predicted", tag="B") #+ theme(aspect.ratio=0.5)
#fix lables
po<- po + ylim(-0.3,0.3) +
  scale_x_discrete(labels=c(expression('OA'[EtNAM]),expression('OA'[EtNAM]*'Ad'),expression('OA'[EtNAM]*'Ge'),expression('OA'[EtNAM]*'Ku'))) + 
  scale_fill_manual(labels = c(expression('GY'[DP]), expression('OA'[DP])), values=gray.colors(length(predictors))) + theme(legend.text.align = 0)
po<- po +  geom_rect(aes(xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf),
                     fill = "pink", alpha = 0.03)
po

#Men and women prediction of EtNAM over DP
idx<-3
names(listpred[[idx]])
rownames(listpred[[idx]][[1]])

# gender predictions on combined data
predictors<-c("OA.F", "OA.M", "OA", "GY")
predicted<-c("GY","OVERALL")
pgend<-plotGP(gpout=listpred, idx=idx, predictors=predictors, predicted=predicted)
#p2<-addSmallLegend(p2) + geom_abline(intercept=0, slope=0, linetype="dashed", color="gray80")
pgend<- pgend + labs(title = "EtNAM to DP", fill="Predictors", x ="Predicted", tag="C") #+ theme(aspect.ratio=0.5)
#fix lables
pgend<- pgend +
  scale_x_discrete(labels=c(expression('GY'[DP]),expression('OA'[DP]))) +
  scale_fill_manual(labels = c(expression('GY'[EtNAM]), expression('OA'[EtNAM]), 
                               expression('OA'[EtNAM]*'W'), expression('OA'[EtNAM]*'M')), 
                    values=c(gray.colors(2),colwom,colmen)) + theme(legend.text.align = 0)
pgend 

fig3new<-p2 / po / pgend
fig3new

ggsave(fig3new, file="./figures/fig.3.pdf", height = 10, width =6)


################# END OF UPDATED FIGURE 3 #####################






################ SUPPLEMENTAL PLOTS #################

################################
#lickert panel
lkad<-plot_likert(byfarm[["Adet"]], catcount=5, expand.grid=TRUE,
                  reverse.scale=TRUE, show.n=F, values="hide", geom.colors = "RdBu",
                  legend.labels = c("Very poor","Poor", "Average","Good", "Very good"), 
                  legend.title="Scores", title="Adet", show.legend=F)
lkad<-lkad+theme_bw()

lkge<-plot_likert(byfarm[["Geregera"]], catcount=5, expand.grid=TRUE,
                  reverse.scale=TRUE, show.n=F, values="hide", geom.colors = "RdBu",
                  legend.labels = c("Very poor","Poor", "Average","Good", "Very good"), 
                  legend.title="Scores", title="Geregera", show.legend=F)
lkge<-lkge+theme_bw()

lkku<-plot_likert(byfarm[["Kulumsa"]], catcount=5, expand.grid=TRUE,
                  reverse.scale=TRUE, show.n=F, values="hide", geom.colors = "RdBu",
                  legend.labels = c("Very poor","Poor", "Average","Good", "Very good"), 
                  legend.title="Scores", title="Kulumsa", show.legend=T)
lkku<-lkku+theme_bw()+theme(legend.position = "bottom", legend.title = element_text( size=8), legend.text=element_text(size=8))

lickertplot<-lkad / lkge / lkku
lickertplot

ggsave(lickertplot, file="./figures/supplemental.lickerplot.pdf", width = 10, height = 10)


##################################
#score distributions by class
bxpltoa<-ggplot(subset(toplot, type != "Recurrent parental"), aes(y=Combined, x=type, fill=type)) + 
  geom_boxplot(width=1, size=0.2) +
  scale_fill_viridis(discrete=TRUE)+
  coord_flip()+
  #ylim(1,5)+
  facet_wrap(.~location) +
  #theme_ipsum() +
  theme_bw() +
  theme(legend.position="none")+
  labs(fill="Type", y="OA", x="", tag="A")

bxpltoa

#now yield
bxpltgy<-ggplot(subset(toplot, type != "Recurrent parental"), aes(y=GY, x=type, fill=type)) + 
  geom_boxplot(width=1, size=0.2) +
  scale_fill_viridis(discrete=TRUE)+
  coord_flip()+
  facet_wrap(.~location) +
  # theme_ipsum() +
  theme_bw() +
  theme(legend.position="none")+
  labs(fill="Type", y="GY", x="", tag = "B")
bxpltgy

bytype<-bxpltoa | bxpltgy
bytype

ggsave(bytype, file = "./figures/trait.distribution.by.type.pdf")


###################################
#correlation across locations
forcor<-cbind(farm[,2:ncol(farm)], met[,2:ncol(met)])
forcor<-forcor[,grep("ad|ge|ku", colnames(forcor))]
#remove combined OAs
forcor<-forcor[,-c(1:3, grep("DB|DH|DM", colnames(forcor)))]
crs<-cor(forcor)
colnames(crs)
crspval<- cor_pmat(forcor)

fixedlab<-c(expression('OA'[Adet]*'W'),
            expression('OA'[Adet]*'M'),
            expression('OA'[Geregera]*'W'),
            expression('OA'[Geregera]*'M'),
            expression('OA'[Kulumsa]*'W'),
            expression('OA'[Kulumsa]*'M'),
            expression('BM'[Adet]),
            expression('BM'[Geregera]),
            expression('BM'[Kulumsa]),
            expression('DF'[Adet]),
            expression('DF'[Geregera]),
            expression('DF'[Kulumsa]),
            expression('GY'[Adet]),
            expression('GY'[Geregera]),
            expression('GY'[Kulumsa]),
            expression('NSPKS'[Adet]),
            expression('NSPKS'[Geregera]),
            expression('NSPKS'[Kulumsa]),
            expression('NTPP'[Adet]),
            expression('NTPP'[Geregera]),
            expression('NTPP'[Kulumsa]),
            expression('PH'[Adet]),
            expression('PH'[Geregera]),
            expression('PH'[Kulumsa]),
            expression('SPL'[Adet]),
            expression('SPL'[Geregera]),
            expression('SPL'[Kulumsa]),
            expression('SPS'[Adet]),
            expression('SPS'[Geregera]),
            expression('SPS'[Kulumsa]),
            expression('TGW'[Adet]),
            expression('TGW'[Geregera]),
            expression('TGW'[Kulumsa]))

corplot<-ggcorrplot(crs,  p.mat =crspval)
corplot<-corplot + 
  scale_x_discrete(labels=fixedlab) + scale_y_discrete(labels=fixedlab)
corplot

ggsave(corplot, file="./figures/correlation.plot.pdf", height=10, width=10)

###############################
#heritabilities plot
head(h2.farm.nam)
head(h2.met.nam)
h2s<-rbind.fill(h2.farm.nam, h2.met.nam)
h2s<-h2s[,-5]
write.table(h2s, file="./figures/heritability.table.txt", row.names=F, quote=F, sep="\t")

#reorder factor levels
h2.farm.nam$location <- factor(h2.farm.nam$location, levels = c("adet", "geregera", "kulumsa", "ALL"))

h2farmplt<-ggplot(h2.farm.nam, aes(y=h2, x=gender, fill=location))+
          geom_bar(position="dodge", stat="identity") +
          scale_fill_manual(name = "Location",values=c("geregera"=colg, "adet"=cola, "kulumsa"=colk, "ALL"="gray"), labels = c("Adet", "Geregera", "Kulumsa", "Combined"))+
          theme_bw()+
          ylim(0,1)+
          scale_x_discrete(labels=c("Combined", "Women", "Men"))+
          labs(y=expression(italic("H"^"2")), x="", tag="A", title="Overal appreciation")
h2farmplt

#reorder factor levels
h2.met.nam$location <- factor(h2.met.nam$location, levels = c("adet", "geregera", "kulumsa", "ALL"))

h2metplt<-ggplot(subset(h2.met.nam, !trait %in% c("DB", "DH", "DM") ), aes(y=h2, x=trait, fill=location))+
  geom_bar(position="dodge", stat="identity") +
  scale_fill_manual(name = "Location",values=c("geregera"=colg, "adet"=cola, "kulumsa"=colk, "ALL"="gray"), labels = c("Adet", "Geregera", "Kulumsa", "Combined"))+
  theme_bw()+
  ylim(0,1)+
  theme(axis.text.x=element_text(angle=90))+ 
  labs(y=expression(italic("H"^"2")), x="", tag="B", title="Agronomic traits")
h2metplt

h2plot<-h2farmplt / h2metplt
h2plot

ggsave(h2plot, file="./figures/h2plots.pdf")


###############################
##test shuklas stabikilty and OA

load("../1.calculate.BLUPs/genoytpe.stability.analysis.Rdata")

#subset the stability dataset
stability<-stability[,c("Y", "CV", "Var", "Shukla", "Wi_g", "MTSI", "MTSIperf")]

teststab<-merge(combpheno, stability, by.x="ID", by.y="row.names")
cr<-cor(teststab[,c(4:6, 13:ncol(teststab))], method="spearman", use = "pairwise.complete")
crspval<- cor_pmat(teststab[,c(4:6, 13:ncol(teststab))])

corstab<-ggcorrplot(cr[,1:3], p.mat=crspval[,1:3], lab = TRUE)+
      scale_y_discrete(labels=c(expression('OA'[Men]), expression('OA'[Women]), expression('OA'))) +
      scale_x_discrete(labels=c(expression('OA'[Men]), 
                                expression('OA'[Women]), 
                                expression('OA'),
                                expression('GY'[Model]),
                                expression('CV'[GY]),
                                expression(sigma[g][GY]),
                                expression('S'[GY]),
                                expression('WI'[GY]),
                                expression('MTSI'),
                                expression('MTSI'[P])))

corstab

ggsave(corstab, file="./figures/stability.correlation.pdf", width=8, height=4)

quantstab<-quantile(teststab[,"MTSI"], 0.75)
beststab<-teststab[which(teststab[,"MTSI"]> quantstab), "ID"]

quantoa<-quantile(teststab[,"Combined"], 0.75)
bestoa<-teststab[which(teststab[,"Combined"]>quantoa), "ID"]

topgen<-beststab[which(beststab %in% bestoa)]
topgen

plot(teststab[,c("Combined", "MTSIperf")])
cor(teststab[,c("Combined", "MTSIperf")])

#ggplot(subset(toplot, location != "Combined"), aes(y=GY, x=location))+ geom_boxplot(col="gray20") +
#  geom_point(data=subset(toplot, ID %in% topgen & location != "Combined"), aes(y=GY, x=location, size=Combined), col="goldenrod")+
#  theme_bw()

##########
#predicting Stability
idx<-5
names(listpred[[idx]])
rownames(listpred[[idx]][[1]])

#select traits to plot
predictors<-c("OVERALL", "GY")
predicted<-c("Y", "CV", "Var", "Shukla", "Wi_g", "MTSI", "MTSIperf")
pstab<-plotGP(gpout=listpred, idx=idx, predictors=predictors, predicted=predicted)
pstab

pstab<-addSmallLegend(pstab) + geom_abline(intercept=0, slope=0, linetype="dashed", color="gray80")
pstab<- pstab + labs(title = "Panel to EtNAM Stability", fill="Predictors", x ="Predicted") #+ theme(aspect.ratio=0.5)

#fix lables
pstab<- pstab + 
  scale_x_discrete(labels=c(expression('CV'[GY]), expression('MTSI'),
                            expression('MTSI'[P]), expression('S'[GY]),
                            expression('CV'[GY]),
                            expression('WI'[GY]),
                            expression('GY'[Model]))) + 
  scale_fill_manual(labels = c(expression('GY'[DP]), expression('OA'[DP])), values=gray.colors(length(predictors))) + 
  theme(legend.text.align = 0)
pstab

ggsave(pstab, file="./figures/predicting.stability.pdf")


###############################
#Men and women prediction by location 
idx<-1
names(listpred[[idx]])
rownames(listpred[[idx]][[1]])

#ADET
#select traits to plot
predictors<-c("OA.adet.F", "OA.adet.M", "OA.adet", "GY.adet")
predicted<-c("GY.geregera", "GY.kulumsa")
padet<-plotGP(gpout=listpred, idx=idx, predictors=predictors, predicted=predicted)
padet<-addSmallLegend(padet) + geom_abline(intercept=0, slope=0, linetype="dashed", color="gray80")
padet<- padet + ylim(-0.2,0.4) + labs(title = "EtNAM Adet over EtNAM", fill="Predictors", x ="Predicted", tag="A") #+ theme(aspect.ratio=0.5)
#fix lables
padet<- padet + 
  scale_x_discrete(labels=c(expression('GY'[Ge]), expression('GY'[Ku]))) + 
  scale_fill_manual(labels = c(expression('GY'[Ad]), expression('OA'[Ad]), expression('OA'[Ad]*"W"), expression('OA'[Ad]*"M")), values=gray.colors(length(predictors))) + theme(legend.text.align = 0)
padet

#GEREGERA
#select traits to plot
predictors<-c("OA.geregera.F", "OA.geregera.M", "OA.geregera", "GY.geregera")
predicted<-c("GY.adet", "GY.kulumsa")
pger<-plotGP(gpout=listpred, idx=idx, predictors=predictors, predicted=predicted)
pger<-addSmallLegend(pger) + geom_abline(intercept=0, slope=0, linetype="dashed", color="gray80")
pger<- pger + ylim(-0.2,0.4)+ labs(title = "EtNAM Geregera over EtNAM", fill="Predictors", x ="Predicted", tag="B") #+ theme(aspect.ratio=0.5)
#fix lables
pger<- pger + 
  scale_x_discrete(labels=c(expression('GY'[Ad]), expression('GY'[Ku]))) + 
  scale_fill_manual(labels = c(expression('GY'[Ge]), expression('OA'[Ge]), expression('OA'[Ge]*"W"), expression('OA'[Ge]*"M")), values=gray.colors(length(predictors))) + theme(legend.text.align = 0)
pger

#KULUMSA
#select traits to plot
predictors<-c("OA.kulumsa.F", "OA.kulumsa.M", "OA.kulumsa", "GY.kulumsa")
predicted<-c("GY.adet", "GY.geregera")
pkul<-plotGP(gpout=listpred, idx=idx, predictors=predictors, predicted=predicted)
pkul<-addSmallLegend(pkul) + geom_abline(intercept=0, slope=0, linetype="dashed", color="gray80")
pkul<- pkul + ylim(-0.2,0.4)+ labs(title = "EtNAM Kulumsa over EtNAM", fill="Predictors", x ="Predicted", tag="C") #+ theme(aspect.ratio=0.5)
#fix lables
pkul<- pkul + 
  scale_x_discrete(labels=c(expression('GY'[Ad]), expression('GY'[Ge]))) + 
  scale_fill_manual(labels = c(expression('GY'[Ku]), expression('OA'[Ku]), expression('OA'[Ku]*"W"), expression('OA'[Ku]*"M")), values=gray.colors(length(predictors))) + theme(legend.text.align = 0)
pkul

nambyloc<-padet + pger + pkul  +  plot_layout(ncol = 2, nrow=2)

ggsave(nambyloc, file="./figures/nambyloc.bygender.pdf", height=8, width=8)

#########################################
## what do farmers predict exactly?

predictors<-c("OVERALL", "EARLINESS", "SPIKE", "TILLER")
predicted<-c("DF", "SPS", "TGW", "GY", "BM", "NSPKPS", "PH", "SPL")
ppred1<-plotGP(gpout=listpred, idx=idx, predictors=predictors, predicted=predicted)
ppred1<-addSmallLegend(ppred1) + geom_abline(intercept=0, slope=0, linetype="dashed", color="gray80")
ppred1<- ppred1 + labs(title = "PVS to metric", fill="Predictors", x ="Predicted") #+ theme(aspect.ratio=0.5)
#fix lables
ppred1<- ppred1 + ylim(-0.3,0.3) + 
  scale_x_discrete(labels=c(expression('BM'[EtNAM]),expression('DF'[EtNAM]),
                            expression('GY'[EtNAM]),expression('NSPKPS'[EtNAM]),
                            expression('PH'[EtNAM]),expression('SPL'[EtNAM]),
                            expression('SPS'[EtNAM]),expression('TGW'[EtNAM]))) + 
  scale_fill_manual(labels = c(expression('EA'[DP]), expression('OA'[DP]), expression('SP'[DP]), expression('TL'[DP])), 
                    values=c(colea, coloa, colsp, coltl)) + theme(legend.text.align = 0)
ppred1

ggsave(ppred1, file="./figures/fig.prediction.PVS.trait.pdf", width = 10, height =8)

####################################
## make supplementary tables
for (idx in 1:length(listpred)){
  gp<-listpred[[idx]]
  tmpname<-names(listpred)[[idx]]
  
  for (i in 1:length(gp)){
    tmp<-gp[[i]]
    tmpout<-data.frame(Accuracy=apply(tmp, 1, mean), SE=apply(tmp, 1, function(x) sd(x)/sqrt(length(x))))
    tmpoutdf<-data.frame(predictor=names(gp)[i], predicted=rownames(tmpout), tmpout)
    gp[[i]]<-tmpoutdf
  }#for i
  
  #assemble a new dataset to plot
  out<-do.call(rbind, gp)
  rownames(out)<-NULL
  out
  
  #write input tables to file
  write.csv(out, file=paste("full.prediction.stats", tmpname, "csv", sep="."), quote=F, row.names=F)
}




######################################
#ADDITIONAL PLOTS
#######################################
idx<-3
names(listpred[[idx]])
rownames(listpred[[idx]][[1]])

# GY and OA
#select traits to plot
predictors<-c("OA", "GY")
predicted<-c("OVERALL", "GY")
p1<-plotGP(gpout=listpred, idx=idx, predictors=predictors, predicted=predicted)
p1<-addSmallLegend(p1) + geom_abline(intercept=0, slope=0, linetype="dashed", color="gray80")
p1<- p1 + labs(title = "Panel to EtNAM", fill="Predictors", x ="Predicted") #+ theme(aspect.ratio=0.5)
#fix lables
p1<- p1 + 
  scale_x_discrete(labels=c(expression('GY'[NAM]), expression('OA'[NAM]))) + 
  scale_fill_manual(labels = c(expression('GY'[DP]), expression('OA'[DP])), values=gray.colors(length(predictors))) + theme(legend.text.align = 0)
p1

ggsave(p1, file="./figures/combined.GS.OA.GY.pdf")


########
#nam over nam GS plot
##
idx<-1

names(listpred[[idx]])
rownames(listpred[[idx]][[1]])

# OA over traits
#select traits to plot
predictors<-c("OA.M", "OA.F", "GY")
predicted<-c("DF", "GY", "SPS","TGW")
p4<-plotGP(gpout=listpred, idx=idx, predictors=predictors, predicted=predicted)
p4<-addSmallLegend(p4) + geom_abline(intercept=0, slope=0, linetype="dashed", color="gray80")
p4<- p4 + labs(title = "EtNAM to EtNAM", fill="Predictors", x ="Predicted") #+ theme(aspect.ratio=0.5)
#fix lables
p4<- p4 + 
  scale_x_discrete(labels=c(expression('DF'[NAM]),
                            expression('GY'[NAM]),
                            expression('SPS'[NAM]),
                            expression('TGW'[NAM]))) + 
  scale_fill_manual(labels = c(expression('GY'[NAM]), expression('OA'[NAM]*'F'),expression('OA'[NAM]*'M')), values=gray.colors(length(predictors))) + theme(legend.text.align = 0)

p4

ggsave(p1, file="./figures/combined.nam.on.nam.pdf")


