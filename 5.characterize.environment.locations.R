######################
# This scripts extracts climatic variation data for the locations under study
# Authors: Kaue de Sousa, Matteo Dell'Acqua
# Date: May 30th, 2022
########################

wd <- "C:/Users/admin/OneDrive - Scuola Superiore Sant'Anna/projects/EtNAM.GS/analysis/8.climatrends"
setwd(wd)

#install.packages("climatrends")
#install.packages("nasapower")
#install.packages("chirps")

library(climatrends)
library(nasapower)
library(chirps)
library(patchwork)

info<-read.delim("../../data/experiments.date.location.txt")
str(info)

inf1<-info[1:2,]
inf2<-info[3,]

inl<-list(inf1, inf2)
outr<-list()
outt<-list()

for (i in 1:2){ #i=1
  tmp<-inl[[i]]
  latlon<-tmp[,2:3]
  
  outr[[i]] <- rainfall(latlon, 
                   day.one = tmp$planting,
                   last.day = tmp$harvesting,
                   #span= 60,
                   timeseries = T,
                   intervals = 7)
  
  outt[[i]] <- temperature(latlon, 
                      day.one = tmp$planting,
                      last.day = tmp$harvesting,
                      #span= 60,
                      timeseries = T,
                      intervals = 7)
}


#fix names before putting everything back together
outr[[1]]$id[which(outr[[1]]$id == 1)]<-"Adet"
outr[[1]]$id[which(outr[[1]]$id == 2)]<-"Geregera"
outr[[2]]$id[which(outr[[2]]$id == 1)]<-"Kulumsa"

rain<-do.call(rbind, outr)
rain$week<-format(rain$date, "%U")

outt[[1]]$id[which(outt[[1]]$id == 1)]<-"Adet"
outt[[1]]$id[which(outt[[1]]$id == 2)]<-"Geregera"
outt[[2]]$id[which(outt[[2]]$id == 1)]<-"Kulumsa"

temp<-do.call(rbind, outt)
temp$week<-format(temp$date, "%U")

ggplot(subset(rain, index==c("Rtotal")), aes(y=value, x= week))+
  geom_point()+ geom_smooth(method = "loess", formula = value ~ week, se=FALSE)+
  facet_wrap(.~id)+
  theme_bw()


ggplot(subset(temp, index==c("maxDT")), aes(y=value, x= week))+
  geom_point()+ geom_smooth(method = "loess", formula = value ~ week, se=FALSE)+
  facet_wrap(.~id)+
  theme_bw()


df<-temp[which(temp$index %in% c("maxDT", "minDT")),]
dfwide<-pivot_wider(df, names_from = index, values_from = value)
dfwide$week<-as.numeric(dfwide$week)
dfwide<-dfwide[order(dfwide$week),]


tmpplot<-ggplot(dfwide) +
          geom_segment(aes(x=week, xend=week, y=minDT, yend=maxDT), color="grey") +
          geom_point( aes(x=week, y=minDT, colour="Min"), size=2.5)+
          geom_point( aes(x=week, y=maxDT, colour="Max"), size=2.5)+
          #coord_flip()+
          scale_color_manual(values = c("Min" = rgb(0.2,0.7,0.1,0.5),'Max' = rgb(0.7,0.2,0.1,0.5))) +
          theme_bw() +
          scale_x_continuous(breaks=seq(24, 50, 3)) +
          labs(x="Week", y="°C", colour="Daily Temperature")+
          facet_wrap(.~id)
tmpplot

rain

df<-rain[which(rain$index %in% c("Rtotal", "SDII")),]
dfwide<-pivot_wider(df, names_from = index, values_from = value)
dfwide$week<-as.numeric(dfwide$week)
dfwide<-dfwide[order(dfwide$week),]
head(dfwide)

rnplot<-ggplot(dfwide, aes(x=week, y = SDII, fill=SDII)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(breaks=seq(24, 50, 3)) +
  scale_fill_continuous(low="lightblue", high="navyblue", na.value = NA) +
  theme_bw() +
  labs(x="Week", y="mm/days")+
  facet_wrap(.~id)

rnplot


combplot<-tmpplot/rnplot + plot_layout(guides = "collect")

combplot      
ggsave(combplot, file="climate.plot.pdf")

