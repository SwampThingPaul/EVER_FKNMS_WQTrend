## Title:      South Florida WQ Trend re-analysis
## Created by: Paul Julian (pauljulianphd@gmail.com)
## Created on: 05/11/2022

## BAD ## https://www.tidyverse.org/articles/2017/12/workflow-vs-script/
#Clears Everything...start fresh.
rm(list=ls(all=T));cat("\014");dev.off()

# Libraries
# devtools::install_github("SwampThingPaul/AnalystHelper")
library(AnalystHelper);
library(plyr)
library(reshape2)
library(openxlsx)

# #GIS Libraries
# library(sp)
# library(rgdal)
# library(PROJ)
# library(rgeos)
# library(tmap)
# library(raster)
# 
# library(fields)
# library(spdep)
# 
# # GAM
# library(ggplot2)
# library(viridis)
# library(mgcv)
# library(gratia)

#Paths
wd="C:/Julian_LaCie/_GitHub/EVER_FKNMS_WQTrend"

paths=paste0(wd,c("/Plots/","/Export/","/Data/","/GIS","/src/"))
# Folder.Maker(paths);#One and done. Creates folders in working directory.
plot.path=paths[1]
export.path=paths[2]
data.path=paths[3]
GIS.path=paste0(wd,"/GIS")

gen.GIS="C:/Julian_LaCie/_GISData"
db.path=paste(gen.GIS,"/SFER_GIS_Geodatabase.gdb",sep=""); 

#Helper variables 
N.mw=14.0067
P.mw=30.973762
C.mw=12.0107

# Map Porjects
wgs84=CRS("+init=epsg:4326")
utm17=CRS("+init=epsg:26917")


# -------------------------------------------------------------------------
# need updated SFWMD WQ data if we want to extend analysis to 2022
# dates=date.fun(c("1995-05-01","2022-04-31"))
dates=date.fun(c("1995-05-01","2019-04-30"))

# Hydrologic Analysis -----------------------------------------------------
RF.sites=data.frame(SITE=c(rep("S9_R",2),rep("S12D_R",2),rep("3A-S_R",2),rep("3A-SW_R",2),rep('S331_R',2),"S331W"),
                    DBKEY=c("13027","UJ621","06055","LS269","05865","HC941","05882","JA344","05967","P6930","16261"),
                    REGION=c(rep("WCA3",8),rep("TS",3)))

RF.dat=data.frame()
pb=txtProgressBar(max=length(RF.sites$DBKEY),style=2)
for(i in 1:length(RF.sites$DBKEY)){
  tmp=DBHYDRO_daily(dates[1],dates[2],RF.sites$DBKEY[i])
  tmp$DBKEY=as.character(RF.sites$DBKEY[i])
  RF.dat=rbind(RF.dat,tmp)
  setTxtProgressBar(pb, i)
}

RF.dat=merge(RF.dat,RF.sites,"DBKEY")
RF.dat$Date=date.fun(RF.dat$Date)
RF.dat$WY=WY(RF.dat$Date)
range(RF.dat$Data.Value,na.rm=T)

RF.dat.da=ddply(RF.dat,c('Date',"REGION"),summarise,mean.RF=mean(Data.Value,na.rm=T),N.val=N.obs(Data.Value))
RF.dat.da$WY=WY(RF.dat.da$Date)

subset(RF.dat.da,REGION=="TS"&WY==2005)
subset(RF.dat,REGION=="TS"&WY==2005)

RF.dat.WY=ddply(RF.dat.da,c("WY","REGION"),summarise,TRF.cm=sum(in.to.cm(mean.RF),na.rm=T))

plot(TRF.cm~WY,subset(RF.dat.WY,REGION=="TS"))
plot(TRF.cm~WY,subset(RF.dat.WY,REGION=="WCA3"))
