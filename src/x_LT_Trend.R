## Title:      FCE LTER Regional Long-term trend analysis
## Created by: Paul Julian (pauljulianphd@gmail.com)
## Created on: 01/03/2023

## BAD ## https://www.tidyverse.org/articles/2017/12/workflow-vs-script/
#Clears Everything...start fresh.
rm(list=ls(all=T));cat("\014");dev.off()

# Libraries
# devtools::install_github("SwampThingPaul/AnalystHelper")
library(AnalystHelper);
library(plyr)
library(reshape2)
library(openxlsx)

library(mblm)

#GIS Libraries
library(sp)
library(rgdal)
library(PROJ)
library(rgeos)
library(tmap)
library(raster)

library(fields)
library(spdep)

# Paths
wd="C:/Julian_LaCie/_GitHub/EVER_FKNMS_WQTrend"
paths=paste0(wd,c("/Plots/FCE_Retreat","/Export/","/Data/","/GIS"))
# Folder.Maker(paths);#One and done. Creates folders in working directory.
plot.path=paths[1]
export.path=paths[2]
data.path=paths[3]
GIS.path=paths[4]


# -------------------------------------------------------------------------
## Dataset at https://doi.org/10.6073/pasta/b2c8e104ecac44a44bcc02ffef3dc6f3
## Currently Embargoed, manuscript in prep
## preprint available; https://www.researchsquare.com/article/rs-1753636/v1

dat=read.csv(paste0(export.path,"20230103_GM_allparam.csv"))
site.locs=read.csv(paste0(export.path,"Sites_TableS1.csv"))

ann.trend=ddply(subset(dat,N.val>3),c("STATION","variable"),summarise,
                est=as.numeric(cor.test(WY,GM,method="kendall")$estimate),
                pval=cor.test(WY,GM,method="kendall")$p.value,
                sen.slope=as.numeric(zyp::zyp.sen(GM~WY)$coefficients[2]),
                N.WY=N.obs(WY))
subset(ann.trend,pval<0.05)
ann.trend$stat.sig=with(ann.trend,ifelse(pval<0.05,"sig","not-sig"))
ann.trend$stat.sig=with(ann.trend,ifelse(is.na(stat.sig)==T,"Empty",as.character(stat.sig)))


AGM.all=ddply(dat,c("STATION","variable"),summarise,
              mean.GM=mean(GM,na.rm=T),
              N.val=N.obs(GM),
              SE.val=SE(GM),
              var.val=var(GM,na.rm=T),# sample variance s^2; use this one (variance within site)
              pop.var=var.val*((N.val-1)/N.val), # population variance sigma^2 
              sd.val=sd(GM,na.rm=T),
              CV.val=cv.per(GM)*100)

vars=c("TN","DIN","TP","SRP","Chla","TOC")
col.vars=c("STATION", "variable", "est", "pval", "sen.slope", "N.WY")
