## 
## Everglades/South Florida Nitrogen Synthesis
##   Modified Code from EvergladesRegionalTrend_v2.R 
##   in .../UF/Eveglades_Nitrogen/... folder
## 
##
## Code was compiled by Paul Julian
## contact info: pjulian@ufl.edu/pjulian.sccf.org

## BAD ## https://www.tidyverse.org/articles/2017/12/workflow-vs-script/
#Clears Everything...start fresh.
rm(list=ls(all=T));cat("\014");dev.off()

#Libraries
library(AnalystHelper);#devtools::install_github("SwampThingPaul/AnalystHelper")
library(plyr)
library(reshape)
library(openxlsx)
library(classInt)
library(zoo)

#GIS Libraries
library(rgdal)
library(rgeos)
library(tmap)
library(raster)

library(fields)
library(spdep)

# GAM
library(ggplot2)
library(viridis)
library(mgcv)
library(gratia)

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

utm17=wkt(CRS(SRS_string="EPSG:26917"))
wgs84=wkt(CRS(SRS_string="EPSG:4326"))

# custom functions
read.lter=function(data.package,PASTA,DOI,na.string.val=c("-9999","-9999.00","-9999.000")){
  prefix="http://pasta.lternet.edu/package/data/eml/"
  
  infile1=paste0(prefix,data.package,PASTA,DOI)
  dt1=read.csv(infile1,na.strings=na.string.val)
  return(dt1)
}
notidy_glance_gam<-function(model,...){
  data.frame(
    df=sum(model$edf),
    df.residual=stats::df.residual(model),
    logLik=as.numeric(stats::logLik(model)),
    AIC = stats::AIC(model),
    BIC = stats::BIC(model),
    adj.r.squared=summary(model)$r.sq,
    deviance=summary(model)$dev.expl,
    nobs = stats::nobs(model),
    method=as.character(summary(model)$method),
    sp.crit=as.numeric(summary(model)$sp.criterion),
    scale.est=summary(model)$scale
  )
}

notidy_tidy_gam<-function(model,dig.num=2,...){
  ptab <- data.frame(summary(model)$p.table)
  ptab$term<-rownames(ptab)
  rownames(ptab)=NULL
  ptab$Component="A. parametric coefficients"
  ptab<-ptab[,c(6,5,1:4)]
  colnames(ptab) <- c("Component","Term", "Estimate", "Std.Error", "t.value", "p.value")
  ptab$p.value=with(ptab,ifelse(p.value<0.01,"<0.01",round(p.value,2)))
  ptab[,3:5]=format(round(ptab[,3:5],dig.num),nsmall=dig.num)
  ptab
  
  stab= data.frame(summary(model)$s.table)
  stab$term<-rownames(stab)
  rownames(stab)=NULL
  stab$Component="B. smooth terms"
  stab<-stab[,c(6,5,1:4)]
  colnames(stab) <- c("Component","Term", "edf", "Ref. df", "F.value", "p.value")
  stab$p.value=with(stab,ifelse(p.value<0.01,"<0.01",round(p.value,2)))
  stab[,3:5]=format(round(stab[,3:5],dig.num),nsmall=dig.num)
  stab
  
  ptab.cnames = c("Component","Term", "Estimate", "Std Error", "t-value", "p-value")
  stab.cnames = c("Component","Term", "edf", "Ref. df", "F-value", "p-value")
  
  colnames(ptab) = c("A", "B", "C", "D")
  if (ncol(stab) != 0) {
    colnames(stab) = colnames(ptab)
  }
  tab = rbind(ptab, stab)
  colnames(tab) = ptab.cnames
  
  tab2 = rbind(c(ptab.cnames), tab[1:nrow(ptab), ])
  if (nrow(stab) > 0) {
    tab2 = rbind(tab2, c(stab.cnames), tab[(nrow(ptab) + 1):nrow(tab), ])
  }
  
  tab2
}

notidy_as_flextable_gam<-function(x,data_t=NULL,data_g=NULL,dig.num=2,r2dig=2,...){
  # needs flextable
  # magrittr
  if(sum(class(x)%in%c("gam"))==1&is.null(data_t)&is.null(data_g)){
    data_t <- notidy_tidy_gam(x)
    data_g <- notidy_glance_gam(x)
  }
  
  std_border=officer::fp_border(color = "black", style = "solid", width = 2)
  data.frame(data_t)%>%
    flextable()%>%
    delete_part(part="header")%>%
    hline(i=which(data_t=="Component"),border=std_border)%>%
    hline(i=which(data_t=="Component")[2]-1,border=std_border)%>%
    bold(i=which(data_t=="Component"))%>%
    align(j=1,part="all")%>%
    hline_top(border=std_border)%>%
    hline_bottom(border=std_border)%>%
    merge_v(j=1)%>%valign(j=1,valign="top")%>%fix_border_issues()%>%
    autofit(part = c("header", "body"))%>%
    add_footer_lines(values = c(
      sprintf("Adjusted R-squared: %s, Deviance explained %s", formatC(data_g$adj.r.squared,digits = r2dig,format="f"), formatC(data_g$deviance,digits = r2dig,format="f")),
      paste0(data_g$method,": ",format(round(data_g$sp.crit,dig.num),dig.num),", Scale est.: ",format(round(data_g$scale.est,dig.num),dig.num),", N: ",data_g$nobs)
    ))
}

tmap_mode("view")
# GIS ---------------------------------------------------------------------
ogrListLayers(paste(gen.GIS,"/SFER_GIS_Geodatabase.gdb",sep=""))
ogrListLayers(paste(gen.GIS,"/AHED_release/AHED_20171102.gdb",sep=""))

structures=spTransform(readOGR(paste(gen.GIS,"/AHED_release/AHED_20171102.gdb",sep=""),"STRUCTURE"),utm17)
wmd.mon=spTransform(readOGR(paste0(gen.GIS,"/SFWMD_Monitoring_20200221"),"Environmental_Monitoring_Stations"),utm17)
lter=spTransform(readOGR(paste0(gen.GIS,"/LTER"),"ltersites_utm"),utm17)

canal=spTransform(readOGR(paste(gen.GIS,"/SFER_GIS_Geodatabase.gdb",sep=""),"SFWMD_Canals"),utm17)
ENP.shore=spTransform(readOGR(paste(gen.GIS,"/SFER_GIS_Geodatabase.gdb",sep=""),"Shoreline_ENPClip"),utm17)
ENP=spTransform(readOGR(paste(gen.GIS,"/SFER_GIS_Geodatabase.gdb",sep=""),"ENP"),utm17)
wca=spTransform(readOGR(paste(gen.GIS,"/SFER_GIS_Geodatabase.gdb",sep=""),"WCAs"),utm17)
bcnp=spTransform(readOGR(paste(gen.GIS,"/SFER_GIS_Geodatabase.gdb",sep=""),"BCNP"),utm17)

regions=spTransform(readOGR(paste0(GIS.path,"/Segments"),"SFL_NNC"),utm17)
regions=merge(regions,data.frame(ESTUARY_SE=c(paste0("ENRE",9:15),paste0("ENRE",5:7),paste0("ENRF",1:6),"ENRH2",paste0("ENRG",1:7)),
                                 Region=c(rep("Coastal_Mangroves",7),rep("Shelf",3),rep("FLBay",7),rep("Keys",7))),"ESTUARY_SE")

shore=spTransform(readOGR(paste(gen.GIS,"/SFER_GIS_Geodatabase.gdb",sep=""),"SFWMD_Shoreline"),utm17)
sloughs=spTransform(readOGR(paste(gen.GIS,"/LTER",sep=""),"sloughs_utm2"),utm17)

sloughs.clp=gIntersection(ENP.shore,sloughs)
sloughs.clp=SpatialPolygonsDataFrame(sloughs.clp,data.frame(ID=1))
plot(sloughs.clp)

##
ddply(regions@data,c("Region","ESTUARY"),summarise,N.val=N.obs(ESTUARY_SE))
subset(regions@data,is.na(Region)==F)

unique(wmd.mon$ACTIVITY_S)
wq.mon=subset(wmd.mon,ACTIVITY_S=="Surface Water Grab")

ENP.sites=c(paste0("P",33:37),"EP","TSB","CR2","RG1","G-3273","NE1","SRS1C","SRS1B","SRS2","NP201",paste0("S12",LETTERS[1:4]),"S333","S332DX","S332D","S332","S197","S177","S18C")
TTI.sites=c("TTI51B",paste0("TTI",52:66))
# Include S-177/S-199? inflow from aerojet (c-111) and/or S18C for coastal basins inputs?
# ENP_FLB=subset(wmd.mon,substr(SITE,1,3)=="SWS"|substr(SITE,1,4)=="FLAB"|SITE%in%ENP.sites)
ENP_FLB=subset(wmd.mon,ACTIVITY_S=="Surface Water Grab"&substr(SITE,1,4)=="FLAB"|
                 ACTIVITY_S=="Surface Water Grab"&STATION%in%c(ENP.sites,TTI.sites))

# combine ENP and regions
region.mask=bind(subset(regions,is.na(Region)==F),ENP)
region.mask$dis=1

# Dissolved
region.mask=gUnaryUnion(region.mask,id=region.mask@data$dis)
region.mask=SpatialPolygonsDataFrame(region.mask,data.frame(ID=1))
# remove empty spaces (could just join shoreline?)
region.mask=gBuffer(gBuffer(region.mask,width=4000,id=1),width=-4100,id=1)
row.names(region.mask)
region.mask=SpatialPolygonsDataFrame(region.mask,data.frame(ID=1))

regions2=gUnaryUnion(regions,id=regions@data$Region)
regions2=SpatialPolygonsDataFrame(regions2,data.frame(Region=c("Coastal_Mangroves", "FLBay", "Keys", "Shelf")),FALSE)
plot(regions2)

# Data --------------------------------------------------------------------
WYs=seq(1996,2019,1)

#from Caccia and Boyer & Rookery Bay SERC report
mdls=data.frame(param=c("TN","TP","TOC","DOC","NOX","NO2","NH4","SRP","CHLA"),
                mdl=c(0.05,0.0012,0.04,0.06,0.0024,0.0003,0.0057,0.0022,0.1))
# mdls$mdl.uM=with(mdls,ifelse(param%in%c("TN","NOX","NO2","NH4"),(mdl*1000)/N.mw,
#                              ifelse(param%in%c("TP","SRP"),(mdl*1000)/P.mw,
#                                     ifelse(param%in%c("TOC","DOC"),(mdl*1000)/C.mw,mdl))))


# FIU SERC ----------------------------------------------------------------
# http://serc.fiu.edu/wqmnetwork/
serc=read.xlsx(paste0(data.path,"FIU_SERC/WQFloridaKeys&Shelf (ppm) UPDATED 6-6-2020.xlsx"),sheet=3)
serc$DATE=date.fun(convertToDate(serc$DATE))
serc=subset(serc,DATE>date.fun("1995-01-01"))
# serc=rename(serc,c("TN-S"="TN"))

serc.sites=ddply(serc,c("STATION"),summarise,LATDEC=mean(LATDEC,na.rm=T),LONDEC=mean(LONDEC,na.rm=T),n.val=N.obs(SITE),min.date=min(DATE,na.rm=T),max.date=max(DATE,na.rm=T))
serc.sites=subset(serc.sites,is.na(LONDEC)==F)
serc.sites.shp=spTransform(SpatialPointsDataFrame(serc.sites[,c("LONDEC","LATDEC")],data=serc.sites,proj4string=CRS(SRS_string="EPSG:4326")),utm17)

tm_shape(serc.sites.shp)+tm_dots()+
  tm_shape(ENP_FLB)+tm_dots(col="red",alpha=0.5)+
  tm_shape(lter)+tm_dots(col="yellow",alpha=0.5)

serc.sites.region=spatialEco::point.in.poly(serc.sites.shp,regions)@data[,c("STATION","ESTUARY","SEGMENT_NA","Region")]
# plot(serc.sites.shp)
# plot(subset(serc.sites.shp,STATION%in%subset(serc.sites.region,is.na(Region)==F)$STATION),add=T,pch=21,bg="red")

serc.sites.region=subset(serc.sites.region,is.na(Region)==F)
serc.sites.region$source="SERC"
serc.sites.region=rbind(serc.sites.region,
                        data.frame(STATION=c(c(316,508,509,506,253,254,504),c(311,303,507,505),500,502,501,503),
                                   ESTUARY=c(rep("Florida Keys",13),"Florida Bay","Florida Bay"),
                                   SEGMENT_NA=c(rep("Lower Keys",7),rep("Back Bay",4),"Upper Keys","Middle Keys","East Central Florida Bay","Southern Florida Bay"),
                                   Region=c(rep("Keys",13),"FLBay","FLBay"),
                                   source="SERC"))
subset(serc.sites.region,is.na(Region)==T)

# Surface data only
names(serc)
serc2=serc[,c("DATE","STATION","TN-S","NOX-S","NH4-S","TP-S","SRP-S","CHLA-S")]
serc2=rename(serc2,c("TN-S"="TN","DIN-S"="DIN","NOX-S"="NOx","NH4-S"="NH4","TP-S"="TP","SRP-S"="SRP","CHLA-S"="Chla"))
serc2[,c("TN","NOx","NH4","TP","SRP","Chla")]=sapply(serc2[,c("TN","NOx","NH4","TP","SRP","Chla")],as.numeric)
serc2$TN=with(serc2,ifelse(TN<0.05,0.05/2,TN))
serc2$NOx=with(serc2,ifelse(NOx<0.0024,0.0024/2,NOx))
serc2$NH4=with(serc2,ifelse(NH4<0.0057,0.0057/2,NH4))
serc2$DIN=with(serc2,ifelse(is.na(NH4)==T|is.na(NOx)==T,NA,NH4+NOx))
serc2$TP=with(serc2,ifelse(TP<0.0012,0.0012/2,TP))
serc2$SRP=with(serc2,ifelse(SRP<0.0022,0.0022/2,SRP))
serc2=subset(serc2,is.na(TN)==F)
serc2$WY=WY(serc2$DATE)
serc2$season=FL.Hydroseason(serc2$DATE)

chk=ddply(serc2,c("STATION"),summarise,n.val=N.obs(TN),min.date=min(DATE,na.rm=T),max.date=max(DATE,na.rm=T))
chk[chk$max.date>date.fun("2011-07-30"),"STATION"]

head(serc2)
serc2=ddply(serc2,c("STATION","DATE","WY","season"),summarise,
            TN=mean(TN,na.rm=T),DIN=mean(DIN,na.rm=T),
            TP=mean(TP,na.rm=T),SRP=mean(SRP,na.rm=T),
            Chla=mean(Chla,na.rm=T))
# tm_shape(subset(serc.sites.shp,STATION%in%chk[chk$max.date>date.fun("2011-07-30"),"STATION"]))+tm_dots()

# Sanity Check
# plot(TN~DATE,subset(serc2,STATION==200))
# plot(DIN~DATE,subset(serc2,STATION==200))
# plot(TP~DATE,subset(serc2,STATION==200))
# plot(SRP~DATE,subset(serc2,STATION==200))
# plot(Chla~DATE,subset(serc2,STATION==200))

serc2=subset(serc2,WY%in%WYs)

# LTER --------------------------------------------------------------------
fce.boyer=read.lter("knb-lter-fce/","1055/9/","acc2be4a77ab1e55b740efdef27648bd",na.string.val=c("-9999","-9999.00","-9999.000","-9999.0"))
fce.boyer=rename(fce.boyer,c("LO.9999EC"="LONDEC"))
fce.boyer.sites=ddply(fce.boyer,c("SITE","LATDEC","LONDEC"),summarise,n.val=N.obs(SITE))
fce.boyer.sites=subset(fce.boyer.sites,n.val>300)
fce.boyer.sites=subset(fce.boyer.sites,is.na(LATDEC)==F|is.na(LONDEC)==F)
fce.boyer.sites=spTransform(SpatialPointsDataFrame(fce.boyer.sites[,c("LONDEC","LATDEC")],data=fce.boyer.sites,proj4string=CRS(SRS_string="EPSG:4326")),utm17)
# only TS/PH9,10,11

tm_shape(serc.sites.shp)+tm_dots()+
  tm_shape(ENP_FLB)+tm_dots(col="red",alpha=0.5)+
  tm_shape(lter)+tm_dots(col="yellow",alpha=0.5)+
  tm_shape(fce.boyer.sites)+tm_dots(col="green",alpha=0.5)

# tm_shape(regions)+tm_polygons()+
#   tm_shape(fce.boyer.sites)+tm_dots(col="green",alpha=0.5)

## Boyer dataset cleanup 
fce.boyer=subset(fce.boyer,SITE%in%c("TS/Ph9","TS/Ph10","TS/Ph11"))
fce.boyer$SITENAME=as.character(fce.boyer$SITE)
fce.boyer$N.N=fce.boyer$NOX
fce.boyer$DOC=NA
fce.boyer$DATE=fce.boyer$Date
fce.boyer$SALINITY=fce.boyer$SAL_S

vars=c("SITENAME", "DATE", "TIME", "SALINITY", "TN", "TP", "TOC","NH4", "N.N", "NO2", "SRP", "DOC", "NO3")
fce.boyer=fce.boyer[,vars]
fce.boyer$clean=rowSums(is.na(fce.boyer[,4:13]))
fce.boyer=subset(fce.boyer,clean!=10)[,vars]
##

fce.grahl=read.lter("knb-lter-fce/","1073/13/","15e6ff92f875a272ba6d98db3f738026")
fce.losada=read.lter("knb-lter-fce/","1075/9/","a3a0e59ec99adb9b15483fda8b504d4b")
fce.rondeau=read.lter("knb-lter-fce/","1077/3/","832e5493edf5a5d79d8a60535e26a012")
fce.rubio=read.lter("knb-lter-fce/","1080/9/","8d000635ce154d208975d15eb7e4f0db")

names(fce.grahl)
names(fce.losada)
names(fce.rondeau)
names(fce.rubio)
fce.grahl=rename(fce.grahl,c("NandN"="N.N","Salinity"="SALINITY","Date"="DATE","Time"="TIME"));

fce.wq=rbind(fce.grahl,fce.losada,fce.rondeau,fce.rubio,fce.boyer)
fce.wq$DATE=date.fun(fce.wq$DATE)
fce.wq$TN=(fce.wq$TN/1000)*N.mw
fce.wq$N.N=(fce.wq$N.N/1000)*N.mw
fce.wq$NH4=(fce.wq$NH4/1000)*N.mw
fce.wq$TP=(fce.wq$TP/1000)*P.mw
fce.wq$SRP=(fce.wq$SRP/1000)*P.mw
fce.wq$DOC=(fce.wq$DOC/1000)*C.mw
range(fce.wq$DOC,na.rm=T)

fce.wq$TN=with(fce.wq,ifelse(TN<0.05,0.05/2,TN))
fce.wq$N.N=with(fce.wq,ifelse(N.N<0.0024,0.0024/2,N.N))
fce.wq$NH4=with(fce.wq,ifelse(NH4<0.0057,0.0057/2,NH4))
fce.wq$DIN=with(fce.wq,ifelse(is.na(NH4)==T|is.na(N.N)==T,NA,NH4+N.N))
fce.wq$TP=with(fce.wq,ifelse(TP<0.0012,0.0012/2,TP))
fce.wq$SRP=with(fce.wq,ifelse(SRP<0.0022,0.0022/2,SRP))

# Sanity Check
# plot(TN~DATE,subset(fce.wq,SITENAME=="SRS2"))
# plot(TN~DATE,subset(fce.wq,SITENAME=="TS/Ph11"))
# plot(DIN~DATE,subset(fce.wq,SITENAME=="SRS2"))

unique(fce.wq$SITENAME)
fce.wq$SITENAME=with(fce.wq,ifelse(substr(SITENAME,1,5)=="TS/Ph",paste0("TS/PH",substr(SITENAME,6,7)),as.character(SITENAME)))

fce.wq$WY=WY(fce.wq$DATE)
fce.wq=subset(fce.wq,WY%in%WYs)
fce.wq$season=FL.Hydroseason(fce.wq$DATE)
fce.wq$Chla=NA
fce.wq$STATION=fce.wq$SITENAME

# write.csv(fce.wq,paste0(export.path,"fce_dat.csv"),row.names = F)

LTER.sites.region=spatialEco::point.in.poly(lter,regions)@data[,c("SITE","ESTUARY","SEGMENT_NA","Region")]
LTER.sites.region=subset(LTER.sites.region,is.na(Region)==F)
LTER.sites.region$source="LTER"
LTER.sites.region=rbind(LTER.sites.region,
                        data.frame(SITE=c(paste("SRS",c("1a","1b","1c","1d",2:3),sep="-"),paste("TS/Ph",c("1a","1b","2","2a","2b",4,3,5,"6a","6b","7a","7b",8:11),sep="-")),
                                   ESTUARY=c(rep("SRS",6),rep("TS",13),rep("FLBay",3)),
                                   SEGMENT_NA=c(rep("SRS_FW",5),rep("SRS_est",1),rep("TS_FW",6),rep("TS_est",7),rep("FLBay",3)),
                                   Region=c(rep("ENP",22)),
                                   source="LTER"))
LTER.sites.region=rename(LTER.sites.region,c("SITE"="STATION"))
LTER.sites.region
spl=strsplit(as.character(LTER.sites.region$STATION),"-")
LTER.sites.region$STATION=paste0(toupper(sapply(spl,"[",1)),(sapply(spl,"[",2)))

# SFWMD -------------------------------------------------------------------
dates=date.fun(c(paste(min(WYs)-1,"05-01",sep="-"),paste(max(WYs),"05-01",sep="-")))
# ENP_FLB$SITE
# ENP_FLB$STATION

ENP_FLB$STATION

# params=data.frame(Test.Number=c(21,20,18,80,25,23,61,179),param=c("TKN","NH4","NOx","TN","TP","SRP","Chla","Chla"))
# wmd.dat=data.frame()
# for(i in 1:length(ENP_FLB$SITE)){
#   tmp=DBHYDRO_WQ(dates[1],dates[2],ENP_FLB$STATION[i],params$Test.Number)
#   wmd.dat=rbind(tmp,wmd.dat)
#   print(i)
# }
# wmd.dat=merge(wmd.dat,params,"Test.Number")
# wmd.dat=subset(wmd.dat,Collection.Method%in%c("G","GP"))
# write.csv(wmd.dat,paste0(export.path,"wmd_dat.csv"),row.names = F)
wmd.dat=read.csv(paste0(export.path,"wmd_dat.csv"))
wmd.dat$Date.EST=date.fun(wmd.dat$Date.EST)

wmd.dat.xtab=reshape2::dcast(wmd.dat,Station.ID+Date.EST~param,value.var="HalfMDL",mean)
wmd.dat.xtab$TN=with(wmd.dat.xtab,TN_Combine(NOx,TKN,TN))
wmd.dat.xtab$DIN=with(wmd.dat.xtab,NH4+NOx)

subset(wmd.dat.xtab,SRP==0)
wmd.dat.xtab$SRP[wmd.dat.xtab$SRP==0]<-NA

wmd.dat.xtab$WY=WY(wmd.dat.xtab$Date.EST)
wmd.dat.xtab=subset(wmd.dat.xtab,WY%in%WYs)
wmd.dat.xtab$season=FL.Hydroseason(wmd.dat.xtab$Date.EST)
wmd.dat.xtab$STATION=wmd.dat.xtab$Station.ID
wmd.dat.xtab$DATE=wmd.dat.xtab$Date.EST
#wmd.station.xwalk=data.frame(Station.ID=ENP_FLB$STATION)
#wmd.station.xwalk$STATION=with(wmd.station.xwalk,ifelse(Station.ID=="S332","S332D",as.character(Station.ID)))
#wmd.dat.xtab=merge(wmd.dat.xtab,wmd.station.xwalk,"Station.ID",all.x=T)

# For analysis purposes, S332 monitoring was terminated and moved to S332D (i.e. inflow to Taylor slough shifted)
wmd.dat.xtab$STATION=with(wmd.dat.xtab,ifelse(STATION=="S332DX","S332D",as.character(STATION)))
# wmd.dat.xtab$STATION=with(wmd.dat.xtab,ifelse(STATION=="S332","S332D",as.character(STATION)))

head(subset(wmd.dat.xtab,Station.ID=="S332"))
head(subset(wmd.dat.xtab,Station.ID=="S332D"))
head(subset(wmd.dat.xtab,Station.ID=="S332DX"))

# remove obvious outliers
wmd.dat.xtab$TN=with(wmd.dat.xtab,ifelse(TN>10,NA,TN))
wmd.dat.xtab$DIN=with(wmd.dat.xtab,ifelse(DIN>10,NA,TN))
wmd.dat.xtab$Chla=with(wmd.dat.xtab,ifelse(Chla>100,NA,TN))

plot(TN~Date.EST,wmd.dat.xtab,type="n")
with(subset(wmd.dat.xtab,Station.ID=="S332"),points(Date.EST,TN,pch=21,bg="red"))
with(subset(wmd.dat.xtab,Station.ID=="S332DX"),points(Date.EST,TN,pch=21,bg="green"))
with(subset(wmd.dat.xtab,Station.ID=="S332D"),points(Date.EST,TN,pch=21,bg="blue"))

wmd.sites.region=spatialEco::point.in.poly(ENP_FLB,regions)@data[,c("STATION","ESTUARY","SEGMENT_NA","Region")]
wmd.sites.region=subset(wmd.sites.region,is.na(Region)==F)
wmd.sites.region$source="WMD"
wmd.sites.region=rbind(wmd.sites.region,
                       data.frame(STATION=c(paste0("S12",LETTERS[1:4]),"S333","SRS1B","SRS1C","NP201","NE1","SRS2",paste0("P",33:36),"G-3273","RG1","CR2","S332D","TSB","S177","S197","S18C","EP","P37"),
                                  ESTUARY=c(rep("SRS",17),rep("TS",7)),
                                  SEGMENT_NA=c(rep("SRS_FW",17),rep("TS_FW",7)),
                                  Region=rep("ENP",24),
                                  source="WMD"))

# Combine datasets
names(serc2)
names(fce.wq)
names(wmd.dat.xtab)

vars=c("STATION","DATE","WY","season","TN","DIN","TP","SRP","Chla")
dat.all=rbind(serc2[,vars],fce.wq[,vars],wmd.dat.xtab[,vars])
# write.csv(dat.all,paste0(export.path,"SERC_FCE_WMD_data.csv"),row.names = F)

all.sites.region=rbind(serc.sites.region,LTER.sites.region,wmd.sites.region)

# test=merge(dat.all,all.sites.region,"STATION",all.x=T)
# nrow(test)
# nrow(dat.all)
# test$dup=as.numeric(duplicated(test))
# test2=ddply(test,c("STATION","Region"),summarise,N.val=N.obs(STATION))
# test2=merge(test2,all.sites.region,c("STATION","Region"),all.y=T)
# subset(test2,is.na(N.val))
# 
# test3=ddply(test,c("STATION","Region","WY","dup","source"),summarise,N.val=N.obs(STATION))
# subset(test3,dup==1)
# subset(dat.all,STATION=="S332D"&WY==2003)

# dat.all2= dat.all
dat.all2=merge(dat.all,all.sites.region,"STATION",all.x=T)
head(dat.all2)




# TN ----------------------------------------------------------------------
TN.samp.size=cast(dat.all2,STATION+WY+ESTUARY+Region+source~season,value="TN",fun.aggregate = function(x)N.obs(x))
TN.samp.size$TSamp=rowSums(TN.samp.size[,c("A_Wet","B_Dry")],na.rm=T)
#for GM calc.
TN.samp.size$sea.screen=with(TN.samp.size, ifelse(A_Wet>0&B_Dry>0&TSamp>=4,1,0))

# number of consecutive years https://stackoverflow.com/a/23786245/5213091
TN.consec.WY=ddply(TN.samp.size,c("STATION","ESTUARY","Region","source"),summarise,
                   max.consec=max(rle(cumsum(c(1,diff(WY)>1)))$lengths),
                   min.consec=min(rle(cumsum(c(1,diff(WY)>1)))$lengths))
TN.consec.WY$ceonsec.diff=with(TN.consec.WY,max.consec-min.consec)
TN.consec.WY$consec.screen=with(TN.consec.WY,ifelse(max.consec>=10,1,0))
TN.consec.WY

vars_join1=c("STATION","WY","ESTUARY","Region","source","sea.screen")
vars_by1=c("STATION","WY","ESTUARY","Region","source")
vars_join2=c("STATION","ESTUARY","Region","source","consec.screen")
vars_by2=c("STATION","ESTUARY","Region","source")
dat.all.TN=merge(dat.all2,TN.samp.size[,vars_join1],vars_by1)
dat.all.TN=merge(dat.all.TN,TN.consec.WY[,vars_join2],vars_by2)

# DIN ---------------------------------------------------------------------
DIN.samp.size=cast(dat.all2,STATION+WY+ESTUARY+Region+source~season,value="DIN",fun.aggregate = function(x)N.obs(x))
DIN.samp.size$TSamp=rowSums(DIN.samp.size[,c("A_Wet","B_Dry")],na.rm=T)
#for GM calc.
DIN.samp.size$sea.screen=with(DIN.samp.size, ifelse(A_Wet>0&B_Dry>0&TSamp>=4,1,0))

# number of consecutive years https://stackoverflow.com/a/23786245/5213091
DIN.consec.WY=ddply(DIN.samp.size,c("STATION","ESTUARY","Region","source"),summarise,
                    max.consec=max(rle(cumsum(c(1,diff(WY)>1)))$lengths),
                    min.consec=min(rle(cumsum(c(1,diff(WY)>1)))$lengths))
DIN.consec.WY$ceonsec.diff=with(DIN.consec.WY,max.consec-min.consec)
DIN.consec.WY$consec.screen=with(DIN.consec.WY,ifelse(max.consec>=10,1,0))
DIN.consec.WY

dat.all.DIN=merge(dat.all2,DIN.samp.size[,vars_join1],vars_by1)
dat.all.DIN=merge(dat.all.DIN,DIN.consec.WY[,vars_join2],vars_by2)

# TP ----------------------------------------------------------------------
TP.samp.size=cast(dat.all2,STATION+WY+ESTUARY+Region+source~season,value="TP",fun.aggregate = function(x)N.obs(x))
TP.samp.size$TSamp=rowSums(TP.samp.size[,c("A_Wet","B_Dry")],na.rm=T)
#for GM calc.
TP.samp.size$sea.screen=with(TP.samp.size, ifelse(A_Wet>0&B_Dry>0&TSamp>=4,1,0))

# number of consecutive years https://stackoverflow.com/a/23786245/5213091
TP.consec.WY=ddply(TP.samp.size,c("STATION","ESTUARY","Region","source"),summarise,
                   max.consec=max(rle(cumsum(c(1,diff(WY)>1)))$lengths),
                   min.consec=min(rle(cumsum(c(1,diff(WY)>1)))$lengths))
TP.consec.WY$ceonsec.diff=with(TP.consec.WY,max.consec-min.consec)
TP.consec.WY$consec.screen=with(TP.consec.WY,ifelse(max.consec>=10,1,0))
TP.consec.WY

dat.all.TP=merge(dat.all2,TP.samp.size[,vars_join1],vars_by1)
dat.all.TP=merge(dat.all.TP,TP.consec.WY[,vars_join2],vars_by2)

# SRP ---------------------------------------------------------------------
SRP.samp.size=cast(dat.all2,STATION+WY+ESTUARY+Region+source~season,value="SRP",fun.aggregate = function(x)N.obs(x))
SRP.samp.size$TSamp=rowSums(SRP.samp.size[,c("A_Wet","B_Dry")],na.rm=T)
#for GM calc.
SRP.samp.size$sea.screen=with(SRP.samp.size, ifelse(A_Wet>0&B_Dry>0&TSamp>=4,1,0))

# number of consecutive years https://stackoverflow.com/a/23786245/5213091
SRP.consec.WY=ddply(SRP.samp.size,c("STATION","ESTUARY","Region","source"),summarise,
                    max.consec=max(rle(cumsum(c(1,diff(WY)>1)))$lengths),
                    min.consec=min(rle(cumsum(c(1,diff(WY)>1)))$lengths))
SRP.consec.WY$ceonsec.diff=with(SRP.consec.WY,max.consec-min.consec)
SRP.consec.WY$consec.screen=with(SRP.consec.WY,ifelse(max.consec>=10,1,0))
SRP.consec.WY

dat.all.SRP=merge(dat.all2,SRP.samp.size[,vars_join1],vars_by1)
dat.all.SRP=merge(dat.all.SRP,SRP.consec.WY[,vars_join2],vars_by2)

# Chla --------------------------------------------------------------------
Chla.samp.size=cast(dat.all2,STATION+WY+ESTUARY+Region+source~season,value="Chla",fun.aggregate = function(x)N.obs(x))
Chla.samp.size$TSamp=rowSums(Chla.samp.size[,c("A_Wet","B_Dry")],na.rm=T)
#for GM calc.
Chla.samp.size$sea.screen=with(Chla.samp.size, ifelse(A_Wet>0&B_Dry>0&TSamp>=4,1,0))

# number of consecutive years https://stackoverflow.com/a/23786245/5213091
Chla.consec.WY=ddply(Chla.samp.size,c("STATION","ESTUARY","Region","source"),summarise,
                     max.consec=max(rle(cumsum(c(1,diff(WY)>1)))$lengths),
                     min.consec=min(rle(cumsum(c(1,diff(WY)>1)))$lengths))
Chla.consec.WY$ceonsec.diff=with(Chla.consec.WY,max.consec-min.consec)
Chla.consec.WY$consec.screen=with(Chla.consec.WY,ifelse(max.consec>=10,1,0))
Chla.consec.WY

dat.all.Chla=merge(dat.all2,Chla.samp.size[,vars_join1],vars_by1)
dat.all.Chla=merge(dat.all.Chla,Chla.consec.WY[,vars_join2],vars_by2)

### 

dat.all.TN.GM=ddply(subset(dat.all.TN,sea.screen==1&consec.screen==1),c("STATION","WY"),summarise,TN.GM=exp(mean(log(TN),na.rm=T)),N.val=N.obs(TN))
dat.all.DIN.GM=ddply(subset(dat.all.DIN,sea.screen==1&consec.screen==1),c("STATION","WY"),summarise,DIN.GM=exp(mean(log(DIN),na.rm=T)),N.val=N.obs(DIN))
dat.all.TP.GM=ddply(subset(dat.all.TP,sea.screen==1&consec.screen==1),c("STATION","WY"),summarise,TP.GM=exp(mean(log(TP*1000),na.rm=T)),N.val=N.obs(TP))
dat.all.SRP.GM=ddply(subset(dat.all.SRP,sea.screen==1&consec.screen==1),c("STATION","WY"),summarise,SRP.GM=exp(mean(log(SRP*1000),na.rm=T)),N.val=N.obs(SRP))
dat.all.Chla.GM=ddply(subset(dat.all.Chla,sea.screen==1&consec.screen==1),c("STATION","WY"),summarise,Chla.GM=exp(mean(log(Chla),na.rm=T)),N.val=N.obs(Chla))

plot(TP.GM~WY,dat.all.TP.GM)
plot(SRP.GM~WY,dat.all.SRP.GM)
plot(TN.GM~WY,dat.all.TN.GM)
plot(DIN.GM~WY,dat.all.DIN.GM)
plot(Chla.GM~WY,dat.all.Chla.GM)

# Trend -------------------------------------------------------------------
TN.ann.trend=ddply(dat.all.TN.GM,"STATION",summarise,
                   est=as.numeric(cor.test(WY,TN.GM,method="kendall")$estimate),
                   pval=cor.test(WY,TN.GM,method="kendall")$p.value,
                   sen.slope=as.numeric(zyp::zyp.sen(TN.GM~WY)$coefficients[2]),
                   N.WY=N.obs(WY))
subset(TN.ann.trend,pval<0.05)
TN.AGM.all=ddply(dat.all.TN.GM,"STATION",summarise,mean.TN.GM=mean(TN.GM,na.rm=T),N.val=N.obs(TN.GM),SE.val=SE(TN.GM))

DIN.ann.trend=ddply(dat.all.DIN.GM,"STATION",summarise,
                    est=as.numeric(cor.test(WY,DIN.GM,method="kendall")$estimate),
                    pval=cor.test(WY,DIN.GM,method="kendall")$p.value,
                    sen.slope=as.numeric(zyp::zyp.sen(DIN.GM~WY)$coefficients[2]),
                    N.WY=N.obs(WY))
# DIN.ann.trend$sen.slope
subset(DIN.ann.trend,pval<0.05)
DIN.AGM.all=ddply(dat.all.DIN.GM,"STATION",summarise,mean.DIN.GM=mean(DIN.GM,na.rm=T),N.val=N.obs(DIN.GM),SE.val=SE(DIN.GM))

TP.ann.trend=ddply(dat.all.TP.GM,"STATION",summarise,
                   est=as.numeric(cor.test(WY,TP.GM,method="kendall")$estimate),
                   pval=cor.test(WY,TP.GM,method="kendall")$p.value,
                   sen.slope=as.numeric(zyp::zyp.sen(TP.GM~WY)$coefficients[2]),
                   N.WY=N.obs(WY))
subset(TP.ann.trend,pval<0.05)
TP.AGM.all=ddply(dat.all.TP.GM,"STATION",summarise,mean.TP.GM=mean(TP.GM,na.rm=T),N.val=N.obs(TP.GM),SE.val=SE(TP.GM))

SRP.ann.trend=ddply(dat.all.SRP.GM,"STATION",summarise,
                    est=as.numeric(cor.test(WY,SRP.GM,method="kendall")$estimate),
                    pval=cor.test(WY,SRP.GM,method="kendall")$p.value,
                    sen.slope=as.numeric(zyp::zyp.sen(SRP.GM~WY)$coefficients[2]),
                    N.WY=N.obs(WY))
subset(SRP.ann.trend,pval<0.05)
subset(SRP.ann.trend,is.na(pval))
SRP.AGM.all=ddply(dat.all.SRP.GM,"STATION",summarise,mean.SRP.GM=mean(SRP.GM,na.rm=T),N.val=N.obs(SRP.GM),SE.val=SE(SRP.GM))

Chla.ann.trend=ddply(dat.all.Chla.GM,"STATION",summarise,
                     est=as.numeric(cor.test(WY,Chla.GM,method="kendall")$estimate),
                     pval=cor.test(WY,Chla.GM,method="kendall")$p.value,
                     sen.slope=as.numeric(zyp::zyp.sen(Chla.GM~WY)$coefficients[2]),
                     N.WY=N.obs(WY))
subset(Chla.ann.trend,pval<0.05)
Chla.AGM.all=ddply(dat.all.Chla.GM,"STATION",summarise,mean.Chla.GM=mean(Chla.GM,na.rm=T),N.val=N.obs(Chla.GM),SE.val=SE(Chla.GM))

# Spatial -----------------------------------------------------------------

serc.shp=cbind(data.frame(STATION=serc.sites.shp$STATION),coordinates(serc.sites.shp))
colnames(serc.shp)=c("STATION","UTMX","UTMY")

spl=strsplit(as.character(lter$SITE),"-")
lter.shp=cbind(data.frame(STATION=paste0(toupper(sapply(spl,"[",1)),(sapply(spl,"[",2)))),coordinates(lter))
tmp=subset(lter.shp,STATION=="TS/PH2a")
tmp$STATION="TS/PH2"
lter.shp=rbind(lter.shp,tmp)
colnames(lter.shp)=c("STATION","UTMX","UTMY")

wmd.shp=cbind(data.frame(STATION=ENP_FLB$STATION),coordinates(ENP_FLB))
# wmd.shp$STATION=with(wmd.shp,ifelse(STATION=="S332DX","S332D",as.character(STATION)))
colnames(wmd.shp)=c("STATION","UTMX","UTMY")

sites.shp=rbind(serc.shp,lter.shp,wmd.shp)

sites.shp.TN.trend=merge(sites.shp,TN.ann.trend,"STATION",all.y=T)
# sites.shp.TN.trend[is.na(sites.shp.TN.trend$UTMY),]
sites.shp.TN.trend=SpatialPointsDataFrame(sites.shp.TN.trend[,c("UTMX","UTMY")],data=sites.shp.TN.trend,proj4string=CRS(SRS_string="EPSG:26917"))
sites.shp.TN.trend$stat.sig=with(sites.shp.TN.trend@data,ifelse(pval<0.05,"sig","not-sig"))
sites.shp.TN.trend$stat.sig=with(sites.shp.TN.trend@data,ifelse(is.na(stat.sig)==T,"Empty",as.character(stat.sig)))
sites.shp.TN.trend$stat.sig=as.factor(sites.shp.TN.trend$stat.sig)
levels(sites.shp.TN.trend$stat.sig)

sites.shp.DIN.trend=merge(sites.shp,DIN.ann.trend,"STATION",all.y=T)
# sites.shp.DIN.trend[is.na(sites.shp.DIN.trend$UTMY),]
sites.shp.DIN.trend=SpatialPointsDataFrame(sites.shp.DIN.trend[,c("UTMX","UTMY")],data=sites.shp.DIN.trend,proj4string=CRS(SRS_string="EPSG:26917"))
sites.shp.DIN.trend$stat.sig=with(sites.shp.DIN.trend@data,ifelse(pval<0.05,"sig","not-sig"))
sites.shp.DIN.trend$stat.sig=with(sites.shp.DIN.trend@data,ifelse(is.na(stat.sig)==T,"Empty",as.character(stat.sig)))
sites.shp.DIN.trend$stat.sig=as.factor(sites.shp.DIN.trend$stat.sig)
levels(sites.shp.DIN.trend$stat.sig)

sites.shp.TP.trend=merge(sites.shp,TP.ann.trend,"STATION",all.y=T)
# sites.shp.TP.trend[is.na(sites.shp.TP.trend$UTMY),]
sites.shp.TP.trend=SpatialPointsDataFrame(sites.shp.TP.trend[,c("UTMX","UTMY")],data=sites.shp.TP.trend,proj4string=CRS(SRS_string="EPSG:26917"))
sites.shp.TP.trend$stat.sig=with(sites.shp.TP.trend@data,ifelse(pval<0.05,"sig","not-sig"))
sites.shp.TP.trend$stat.sig=with(sites.shp.TP.trend@data,ifelse(is.na(stat.sig)==T,"Empty",as.character(stat.sig)))
sites.shp.TP.trend$stat.sig=as.factor(sites.shp.TP.trend$stat.sig)
levels(sites.shp.TP.trend$stat.sig)

# STATION 286 SRP all non-detect therefore no trend could be calculated. 
sites.shp.SRP.trend=merge(sites.shp,subset(SRP.ann.trend,is.na(pval)==F),"STATION",all.y=T)
# sites.shp.SRP.trend[is.na(sites.shp.SRP.trend$UTMY),]
sites.shp.SRP.trend=SpatialPointsDataFrame(sites.shp.SRP.trend[,c("UTMX","UTMY")],data=sites.shp.SRP.trend,proj4string=CRS(SRS_string="EPSG:26917"))
sites.shp.SRP.trend$stat.sig=with(sites.shp.SRP.trend@data,ifelse(pval<0.05,"sig","not-sig"))
sites.shp.SRP.trend$stat.sig=with(sites.shp.SRP.trend@data,ifelse(is.na(stat.sig)==T,"Empty",as.character(stat.sig)))
sites.shp.SRP.trend$stat.sig=as.factor(sites.shp.SRP.trend$stat.sig)
levels(sites.shp.SRP.trend$stat.sig)

sites.shp.Chla.trend=merge(sites.shp,Chla.ann.trend,"STATION",all.y=T)
# sites.shp.Chla.trend[is.na(sites.shp.Chla.trend$UTMY),]
sites.shp.Chla.trend=SpatialPointsDataFrame(sites.shp.Chla.trend[,c("UTMX","UTMY")],data=sites.shp.Chla.trend,proj4string=CRS(SRS_string="EPSG:26917"))
sites.shp.Chla.trend$stat.sig=with(sites.shp.Chla.trend@data,ifelse(pval<0.05,"sig","not-sig"))
sites.shp.Chla.trend$stat.sig=with(sites.shp.Chla.trend@data,ifelse(is.na(stat.sig)==T,"Empty",as.character(stat.sig)))
sites.shp.Chla.trend$stat.sig=as.factor(sites.shp.Chla.trend$stat.sig)
levels(sites.shp.Chla.trend$stat.sig)

sites.shp2=SpatialPointsDataFrame(sites.shp[,c("UTMX","UTMY")],data=sites.shp,proj4string = CRS(SRS_string="EPSG:26917"))


# SamplingMap -------------------------------------------------------------
serc.shp2=SpatialPointsDataFrame(serc.shp[,c("UTMX","UTMY")],data=serc.shp,proj4string =CRS(SRS_string="EPSG:26917") )
wmd.shp2=SpatialPointsDataFrame(wmd.shp[,c("UTMX","UTMY")],data=wmd.shp,proj4string =CRS(SRS_string="EPSG:26917") )
lter.shp2=SpatialPointsDataFrame(lter.shp[,c("UTMX","UTMY")],data=lter.shp,proj4string =CRS(SRS_string="EPSG:26917") )

cols=wesanderson::wes_palette("Zissou1",4,"continuous")
# png(filename=paste0(plot.path,"SamplingMap.png"),width=6.5,height=4,units="in",res=200,type="windows",bg="white")
par(family="serif",oma=c(0.25,0.25,0.25,0.25),mar=c(0.1,0.1,0.1,0.1),xpd=F)
layout(matrix(c(1:2),1,2,byrow=T),widths = c(1,0.3))
bbox.lims=bbox(region.mask)

plot(shore,col="cornsilk",border="grey",bg="lightblue",ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.01)
plot(wca,add=T,col="grey90",border=NA)
plot(bcnp,add=T,col="grey90",border=NA)
plot(canal,add=T,col="lightblue",lwd=1)
plot(gSimplify(subset(regions2,Region=="Coastal_Mangroves"),500),add=T,col=adjustcolor(cols[1],0.5),lwd=0.5)
plot(gSimplify(subset(regions2,Region=="FLBay"),500),add=T,col=adjustcolor(cols[2],0.5),lwd=0.5)
plot(gSimplify(subset(regions2,Region=="Shelf"),500),add=T,col=adjustcolor(cols[3],0.5),lwd=0.5)
plot(gSimplify(subset(regions2,Region=="Keys"),500),add=T,col=adjustcolor(cols[4],0.5),lwd=0.5)
plot(ENP,add=T,bg=NA,lwd=1)
plot(serc.shp2,add=T,lwd=0.1,pch=21,bg="white",cex=0.75)
plot(wmd.shp2,add=T,lwd=0.1,pch=22,bg="grey",cex=0.75)
plot(lter.shp2,add=T,lwd=0.1,pch=24,bg="black",cex=0.75)
mapmisc::scaleBar(utm17,"bottomright",bty="n",cex=1,seg.len=4)
box(lwd=1)
plot(0:1,0:1,ann=F,axes=F,type="n")
legend(0.5,0.8,legend=c(paste0("SERC (",length(serc.shp2$STATION),")"),
                        paste0("SFWMD (",length(wmd.shp2$STATION),")"),
                        paste0("FCE LTER (",length(lter.shp2$STATION),")")),
       pch=c(21,22,24),lty=c(NA),lwd=c(0.1),
       col=c("black"),pt.bg=c("white","grey","black"),
       pt.cex=1.25,ncol=1,cex=0.8,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=1,
       title.adj = 0,title="Monitoring\n(Number of Sites)")
legend(0.5,0.4,legend=c("ENP","Mangrove Fringe","Florida Bay","W. Florida Shelf","Keys"),
       pch=c(22),lty=c(NA),lwd=c(0.1),
       col=c("black"),pt.bg=c("white",adjustcolor(cols,0.5)),
       pt.cex=1.25,ncol=1,cex=0.8,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=1,
       title.adj = 0,title="Regions")
dev.off()

# Spline ------------------------------------------------------------------
#thin plate spline https://rspatial.org/raster/analysis/4-interpolation.html

region.buf.r=raster(region.mask)
res(region.buf.r)=1000

## Spatial Trend
# TN
m=Tps(coordinates(sites.shp.TN.trend),sites.shp.TN.trend$sen.slope)
tps=interpolate(region.buf.r,m)
tps.TN.trend=mask(tps,region.mask)

tm_shape(tps.TN.trend)+tm_raster(title="Sen Slope (mg N L\u207B\u00B9 Yr\u207B\u00B9)",alpha=0.5,palette = "viridis")+
  tm_shape(sites.shp.TN.trend)+tm_dots(col="stat.sig",palette=c("green","red","grey"),title="Statistical Significant (\u03C1<0.05)")

# DIN
m=Tps(coordinates(sites.shp.DIN.trend),sites.shp.DIN.trend$sen.slope)
tps=interpolate(region.buf.r,m)
tps.DIN.trend=mask(tps,region.mask)

tm_shape(tps.DIN.trend)+tm_raster(title="Sen Slope (mg N L\u207B\u00B9 Yr\u207B\u00B9)",alpha=0.5,palette = "viridis")+
  tm_shape(sites.shp.DIN.trend)+tm_dots(col="stat.sig",palette=c("green","red","grey"),title="Statistical Significant (\u03C1<0.05)")

# TP
m=Tps(coordinates(sites.shp.TP.trend),sites.shp.TP.trend$sen.slope)
tps=interpolate(region.buf.r,m)
tps.TP.trend=mask(tps,region.mask)

tm_shape(tps.TP.trend)+tm_raster(title="Sen Slope (ug P L\u207B\u00B9 Yr\u207B\u00B9)",alpha=0.5,palette = "viridis")+
  tm_shape(sites.shp.TP.trend)+tm_dots(col="stat.sig",palette=c("green","red","grey"),title="Statistical Significant (\u03C1<0.05)")

# SRP
m=Tps(coordinates(sites.shp.SRP.trend),sites.shp.SRP.trend$sen.slope)
tps=interpolate(region.buf.r,m)
tps.SRP.trend=mask(tps,region.mask)

tm_shape(tps.SRP.trend)+tm_raster(title="Sen Slope (ug P L\u207B\u00B9 Yr\u207B\u00B9)",alpha=0.5,palette = "viridis")+
  tm_shape(sites.shp.SRP.trend)+tm_dots(col="stat.sig",palette=c("green","red","grey"),title="Statistical Significant (\u03C1<0.05)")

# Chla
m=Tps(coordinates(sites.shp.Chla.trend),sites.shp.Chla.trend$sen.slope)
tps=interpolate(region.buf.r,m)
tps.Chla.trend=mask(tps,region.mask)

tm_shape(tps.Chla.trend)+tm_raster(title="Sen Slope (ug L\u207B\u00B9 Yr\u207B\u00B9)",alpha=0.5,palette = "viridis")+
  tm_shape(sites.shp.Chla.trend)+tm_dots(col="stat.sig",palette=c("green","red","grey"),title="Statistical Significant (\u03C1<0.05)")

## Spatial AGM
# TN
sites.shp.TN.GM=merge(sites.shp,TN.AGM.all,"STATION",all.y=T)
sites.shp.TN.GM=SpatialPointsDataFrame(sites.shp.TN.GM[,c("UTMX","UTMY")],data=sites.shp.TN.GM,proj4string=CRS(SRS_string ="EPSG:26917"))
sites.shp.TN.GM=subset(sites.shp.TN.GM,is.na(mean.TN.GM)==F)

m.GM=Tps(coordinates(sites.shp.TN.GM),sites.shp.TN.GM$mean.TN.GM)
tps.GM=interpolate(region.buf.r,m.GM)
tps.TN.GM=mask(tps.GM,region.mask)

tm_shape(tps.TN.GM)+tm_raster(title="Average Annual GM TN (mg N L\u207B\u00B9)",alpha=0.5,palette = "viridis")+
  tm_shape(sites.shp.TN.GM)+tm_dots(col="white",alpha=0.5)

m.GM.SE=Tps(coordinates(sites.shp.TN.GM),sites.shp.TN.GM$SE.val)
tps.GM.SE=interpolate(region.buf.r,m.GM.SE)
tps.TN.GM.SE=mask(tps.GM.SE,region.mask)

tm_shape(tps.TN.GM.SE)+tm_raster(title="SE of average GM TN (mg N L\u207B\u00B9)",alpha=0.5,palette = "viridis")+
  tm_shape(sites.shp.TN.GM)+tm_dots(col="white",alpha=0.5)

# DIN
sites.shp.DIN.GM=merge(sites.shp,DIN.AGM.all,"STATION",all.y=T)
sites.shp.DIN.GM=SpatialPointsDataFrame(sites.shp.DIN.GM[,c("UTMX","UTMY")],data=sites.shp.DIN.GM,proj4string=CRS(SRS_string ="EPSG:26917"))
sites.shp.DIN.GM=subset(sites.shp.DIN.GM,is.na(mean.DIN.GM)==F)

m.GM=Tps(coordinates(sites.shp.DIN.GM),sites.shp.DIN.GM$mean.DIN.GM)
tps.GM=interpolate(region.buf.r,m.GM)
tps.DIN.GM=mask(tps.GM,region.mask)

tm_shape(tps.DIN.GM)+tm_raster(title="Average Annual GM DIN (mg N L\u207B\u00B9)",alpha=0.5,palette = "viridis")+
  tm_shape(sites.shp.DIN.GM)+tm_dots(col="white",alpha=0.5)

m.GM.SE=Tps(coordinates(sites.shp.DIN.GM),sites.shp.DIN.GM$SE.val)
tps.GM.SE=interpolate(region.buf.r,m.GM.SE)
tps.DIN.GM.SE=mask(tps.GM.SE,region.mask)

tm_shape(tps.DIN.GM.SE)+tm_raster(title="SE of average GM DIN (mg N L\u207B\u00B9)",alpha=0.5,palette = "viridis")+
  tm_shape(sites.shp.DIN.GM)+tm_dots(col="white",alpha=0.5)

# TP
sites.shp.TP.GM=merge(sites.shp,TP.AGM.all,"STATION",all.y=T)
sites.shp.TP.GM=SpatialPointsDataFrame(sites.shp.TP.GM[,c("UTMX","UTMY")],data=sites.shp.TP.GM,proj4string=CRS(SRS_string ="EPSG:26917"))
sites.shp.TP.GM=subset(sites.shp.TP.GM,is.na(mean.TP.GM)==F)

m.GM=Tps(coordinates(sites.shp.TP.GM),sites.shp.TP.GM$mean.TP.GM,method="REML")
tps.GM=interpolate(region.buf.r,m.GM)
tps.TP.GM=mask(tps.GM,region.mask)

tm_shape(tps.TP.GM)+tm_raster(title="Average Annual GM TP (ug P L\u207B\u00B9)",alpha=0.5,palette = "viridis")+
  tm_shape(sites.shp.TP.GM)+tm_dots(col="white",alpha=0.5)

m.GM.SE=Tps(coordinates(sites.shp.TP.GM),sites.shp.TP.GM$SE.val)
tps.GM.SE=interpolate(region.buf.r,m.GM.SE)
tps.TP.GM.SE=mask(tps.GM.SE,region.mask)

tm_shape(tps.TP.GM.SE)+tm_raster(title="SE of average GM TP (ug P L\u207B\u00B9)",alpha=0.5,palette = "viridis")+
  tm_shape(sites.shp.TP.GM)+tm_dots(col="white",alpha=0.5)

# SRP
sites.shp.SRP.GM=merge(sites.shp,SRP.AGM.all,"STATION",all.y=T)
sites.shp.SRP.GM=SpatialPointsDataFrame(sites.shp.SRP.GM[,c("UTMX","UTMY")],data=sites.shp.SRP.GM,proj4string=CRS(SRS_string ="EPSG:26917"))
sites.shp.SRP.GM=subset(sites.shp.SRP.GM,is.na(mean.SRP.GM)==F)

m.GM=Tps(coordinates(sites.shp.SRP.GM),sites.shp.SRP.GM$mean.SRP.GM)
tps.GM=interpolate(region.buf.r,m.GM)
tps.SRP.GM=mask(tps.GM,region.mask)

tm_shape(tps.SRP.GM)+tm_raster(title="Average Annual GM SRP (ug P L\u207B\u00B9)",alpha=0.5,palette = "viridis")+
  tm_shape(sites.shp.SRP.GM)+tm_dots(col="white",alpha=0.5)

m.GM.SE=Tps(coordinates(sites.shp.SRP.GM),sites.shp.SRP.GM$SE.val)
tps.GM.SE=interpolate(region.buf.r,m.GM.SE)
tps.SRP.GM.SE=mask(tps.GM.SE,region.mask)

tm_shape(tps.SRP.GM.SE)+tm_raster(title="SE of average GM SRP (ug P L\u207B\u00B9)",alpha=0.5,palette = "viridis")+
  tm_shape(sites.shp.TP.GM)+tm_dots(col="white",alpha=0.5)

# Chla
sites.shp.Chla.GM=merge(sites.shp,Chla.AGM.all,"STATION",all.y=T)
sites.shp.Chla.GM=SpatialPointsDataFrame(sites.shp.Chla.GM[,c("UTMX","UTMY")],data=sites.shp.Chla.GM,proj4string=CRS(SRS_string ="EPSG:26917"))
sites.shp.Chla.GM=subset(sites.shp.Chla.GM,is.na(mean.Chla.GM)==F)

m.GM=Tps(coordinates(sites.shp.Chla.GM),sites.shp.Chla.GM$mean.Chla.GM)
tps.GM=interpolate(region.buf.r,m.GM)
tps.Chla.GM=mask(tps.GM,region.mask)

tm_shape(tps.Chla.GM)+tm_raster(title="Average Annual GM Chla (ug L\u207B\u00B9)",alpha=0.5,palette = "viridis")+
  tm_shape(sites.shp.Chla.GM)+tm_dots(col="white",alpha=0.5)

m.GM.SE=Tps(coordinates(sites.shp.Chla.GM),sites.shp.Chla.GM$SE.val)
tps.GM.SE=interpolate(region.buf.r,m.GM.SE)
tps.Chla.GM.SE=mask(tps.GM.SE,region.mask)

tm_shape(tps.Chla.GM.SE)+tm_raster(title="SE of average GM Chla (ug N L\u207B\u00B9)",alpha=0.5,palette = "viridis")+
  tm_shape(sites.shp.TN.GM)+tm_dots(col="white",alpha=0.5)


# Climate -----------------------------------------------------------------
## AMO
vars=c('year',month.abb)
row.count=length(seq(1856,2021,1))
noaa.amo.path="https://psl.noaa.gov/data/correlation/amon.us.long.data"

# AMO.dat=read.table("https://psl.noaa.gov/data/correlation/amon.us.long.data",header=F,skip=1,col.names=vars,nrows=row.count,na.string="-99.990")
AMO.dat=read.table("https://psl.noaa.gov/data/correlation/amon.sm.long.data",header=F,skip=1,col.names=vars,nrows=row.count,na.string="-99.990")
AMO.dat.melt=melt(AMO.dat,id.vars="year")
AMO.dat.melt=merge(AMO.dat.melt,data.frame(variable=month.abb,month=1:12))
AMO.dat.melt$Date.mon=with(AMO.dat.melt,date.fun(paste(year,month,"01",sep="-")))
AMO.dat.melt=AMO.dat.melt[order(AMO.dat.melt$Date.mon),c("Date.mon","value")]
AMO.dat.melt$warm=with(AMO.dat.melt,ifelse(value>0,value,0))
AMO.dat.melt$dry=with(AMO.dat.melt,ifelse(value<0,value,0))
AMO.dat.melt$ma=with(AMO.dat.melt,c(rep(NA,120),zoo::rollapply(value,width=121,FUN=function(x)mean(x,na.rm=T))))
head(AMO.dat.melt)
tail(AMO.dat.melt)

layout(matrix(c(1:2),2,1,byrow=T))
ylim.val=c(-0.4,0.4);by.y=0.5;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=date.fun(c("1870-01-01","2016-12-01"));xmaj=seq(xlim.val[1],xlim.val[2],"20 years");xmin=seq(xlim.val[1],xlim.val[2],"1 years")

plot(value~Date.mon,AMO.dat.melt,xlim=xlim.val,ylim=ylim.val,type="n")
#with(AMO.dat.melt,lines(Date.mon,ma,col="red"))
with(subset(AMO.dat.melt,is.na(value)==F),shaded.range(Date.mon,rep(0,length(Date.mon)),ifelse(value>0,value,0),"indianred1",lty=1))
with(subset(AMO.dat.melt,is.na(value)==F),shaded.range(Date.mon,ifelse(value<0,value,0),rep(0,length(Date.mon)),"dodgerblue1",lty=1))
abline(h=0)

# PDO
# https://www.ncdc.noaa.gov/teleconnections/pdo/
# pdo=read.csv("https://www.ncdc.noaa.gov/teleconnections/pdo/data.csv",skip=1)
nrow.val=length(seq(1854,2021,1))
head.val=c("yr",month.abb)
pdo=read.table(paste0(data.path,"NOAA/PDO/data.txt"),skip=2,header=F,col.names=head.val,nrows = nrow.val-1)
head(pdo);tail(pdo)
pdo$yr=as.numeric(pdo$yr)
pdo=melt(pdo,id.var="yr")
pdo$month.num=with(pdo,as.numeric(match(variable,month.abb)))
pdo$monCY.date=with(pdo,date.fun(paste(yr,month.num,1,sep="-")))
pdo$WY=WY(pdo$monCY.date)
pdo$dec.WY=decimal.WY(pdo$monCY.date)
pdo=pdo[order(pdo$dec.WY),]

pdo.WY.dat=ddply(pdo,"WY",summarise,mean.PDO=mean(as.numeric(value),na.rm=T),sd.PDO=sd(as.numeric(value),na.rm=T),N.val=N.obs(value))
pdo.WY.dat$UCI=with(pdo.WY.dat,mean.PDO+qnorm(0.975)*sd.PDO/sqrt(N.val))
pdo.WY.dat$LCI=with(pdo.WY.dat,mean.PDO-qnorm(0.975)*sd.PDO/sqrt(N.val))

plot(mean.PDO~WY,pdo.WY.dat)
with(pdo.WY.dat,lines(WY,UCI))
with(pdo.WY.dat,lines(WY,LCI))

# png(filename=paste0(plot.path,"ClimateIndex.png"),width=6.5,height=4,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,3,0.5,1.5),oma=c(2,1,0.25,0.25),xpd=F);
layout(matrix(c(1:2),2,1,byrow=T))

xlim.val=date.fun(c("1870-01-01","2021-12-01"));xmaj=seq(xlim.val[1],xlim.val[2],"20 years");xmin=seq(xlim.val[1],xlim.val[2],"1 years")
ylim.val=c(-0.4,0.4);by.y=0.2;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(value~Date.mon,AMO.dat.melt,xlim=xlim.val,ylim=ylim.val,type="n",axes=F,ann=F)
abline(h=ymaj,v=xmaj,lty=1,col=adjustcolor("grey",0.5))
#with(AMO.dat.melt,lines(Date.mon,ma,col="red"))
with(subset(AMO.dat.melt,is.na(value)==F),shaded.range(Date.mon,rep(0,length(Date.mon)),ifelse(value>0,value,0),"indianred1",lty=1))
with(subset(AMO.dat.melt,is.na(value)==F),shaded.range(Date.mon,ifelse(value<0,value,0),rep(0,length(Date.mon)),"dodgerblue1",lty=1))
abline(h=0)
axis_fun(1,xmaj,xmin,NA)
axis_fun(2,ymaj,ymin,format(ymaj,nsmall=1));box(lwd=1)
abline(v=date.fun(c("1995-05-01","2019-05-01")),lty=2)
text(date.fun(date.fun("1995-05-01")+lubridate::ddays(4383)),ylim.val[2],"Study Period",cex=0.75,font=3)
mtext(side=2,line=2.5,"AMO Index")

ylim.val=c(-4,4);by.y=2;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(value~monCY.date,pdo,xlim=xlim.val,ylim=ylim.val,type="n",axes=F,ann=F)
abline(h=ymaj,v=xmaj,lty=1,col=adjustcolor("grey",0.5))
with(subset(pdo,is.na(value)==F),shaded.range(monCY.date,rep(0,length(monCY.date)),ifelse(value>0,value,0),"indianred1",lty=1))
with(subset(pdo,is.na(value)==F),shaded.range(monCY.date,ifelse(value<0,value,0),rep(0,length(monCY.date)),"dodgerblue1",lty=1))
abline(h=0)
axis_fun(1,xmaj,xmin,format(xmaj,"%Y"),line=-0.5)
axis_fun(2,ymaj,ymin,format(ymaj,nsmall=1));box(lwd=1)
abline(v=date.fun(c("1995-05-01","2019-05-01")),lty=2)
mtext(side=2,line=2.5,"PDO Index")
mtext(side=1,line=1.5,"Year")
dev.off()

# png(filename=paste0(plot.path,"ClimateIndex_v2.png"),width=6.5,height=4,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,3,0.5,1.5),oma=c(2,1,0.25,0.25),xpd=F);
layout(matrix(c(1:2),2,1,byrow=T))

xlim.val=date.fun(c("1995-01-01","2021-12-01"));xmaj=seq(xlim.val[1],xlim.val[2],"10 years");xmin=seq(xlim.val[1],xlim.val[2],"1 years")
ylim.val=c(-0.2,0.2);by.y=0.1;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(value~Date.mon,AMO.dat.melt,xlim=xlim.val,ylim=ylim.val,type="n",axes=F,ann=F)
abline(h=ymaj,v=xmaj,lty=1,col=adjustcolor("grey",0.5))
#with(AMO.dat.melt,lines(Date.mon,ma,col="red"))
with(subset(AMO.dat.melt,is.na(value)==F),shaded.range(Date.mon,rep(0,length(Date.mon)),ifelse(value>0,value,0),"indianred1",lty=1))
with(subset(AMO.dat.melt,is.na(value)==F),shaded.range(Date.mon,ifelse(value<0,value,0),rep(0,length(Date.mon)),"dodgerblue1",lty=1))
abline(h=0)
axis_fun(1,xmaj,xmin,NA)
axis_fun(2,ymaj,ymin,format(ymaj,nsmall=1));box(lwd=1)
abline(v=date.fun(c("1995-05-01","2019-05-01")),lty=2)
text(date.fun(date.fun("1995-05-01")+lubridate::ddays(4383)),ylim.val[2],"Study Period",cex=0.75,font=3)
mtext(side=2,line=2.5,"AMO Index")

ylim.val=c(-4,4);by.y=2;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(value~monCY.date,pdo,xlim=xlim.val,ylim=ylim.val,type="n",axes=F,ann=F)
abline(h=ymaj,v=xmaj,lty=1,col=adjustcolor("grey",0.5))
with(subset(pdo,is.na(value)==F),shaded.range(monCY.date,rep(0,length(monCY.date)),ifelse(value>0,value,0),"indianred1",lty=1))
with(subset(pdo,is.na(value)==F),shaded.range(monCY.date,ifelse(value<0,value,0),rep(0,length(monCY.date)),"dodgerblue1",lty=1))
abline(h=0)
axis_fun(1,xmaj,xmin,format(xmaj,"%Y"),line=-0.5)
axis_fun(2,ymaj,ymin,format(ymaj,nsmall=1));box(lwd=1)
abline(v=date.fun(c("1995-05-01","2019-05-01")),lty=2)
mtext(side=2,line=2.5,"PDO Index")
mtext(side=1,line=1.5,"Year")
dev.off()

# GAM ---------------------------------------------------------------------
# TN
dat.all.TN.GM2=merge(dat.all.TN.GM,sites.shp,"STATION")

WY.k=23
loc.k=200
m.TN<-bam(log(TN.GM)~
            s(WY,bs="cr",k=WY.k)+
            s(UTMX,UTMY,bs="ds",k=loc.k,m=c(1,0.5))+
            ti(UTMX,UTMY,WY,d=c(2,1),bs=c("ds","cr"),k=c(100,24)),
          data=dat.all.TN.GM2,
          nthreads = c(12),discrete=T)
summary(m.TN)
qq.gam(m.TN)
nvar=3;layout(matrix(1:nvar,1,nvar))
plot(m.TN,residuals=T,pch=21)
plot(m.TN,residuals=T,pch=21,select=1,shade=T)

nvar=4;layout(matrix(1:nvar,1,nvar))
gam.check(m.TN)
dev.off()

m.TN.sum=notidy_tidy_gam(m.TN)
m.TN.est=notidy_glance_gam(m.TN)
# write.csv(m.TN.sum,export.path("TN_gam_mod_sum.csv"),row.names=F)
# write.csv(m.TN.est,export.path("TN_gam_mod_est.csv"),row.names=F)

tmp=plot(m.TN,residuals=T)
# png(filename=paste0(plot.path,"GAM_TN_draw_base.png"),width=8,height=3,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(2,2,1,0.5),oma=c(2,2,0.25,0.25),xpd=F);
layout(matrix(1:3,1,3),widths=c(1,1,0.4))

crit=qnorm((1 - 0.95) / 2, lower.tail = FALSE)
ylim.val=c(-1.25,1);by.y=0.5;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(1996,2019);by.x=5;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/by.x)
plot(tmp[1][[1]]$fit~tmp[1][[1]]$x,ylim=ylim.val,xlim=xlim.val,ann=F,axes=F,type="n")
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(tmp[1][[1]],points(jitter(raw,0.5),p.resid,pch=19,col=adjustcolor("dodgerblue1",0.10)))
with(tmp[1][[1]],shaded.range(x,fit-(crit*se),fit+(crit*se),"grey",lty=1))
with(tmp[1][[1]],lines(x,fit,lwd=2))
abline(h=0)
axis_fun(1,xmaj,xmin,xmaj,line=-0.5);axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=3,adj=0,"Total Nitrogen")
mtext(side=3,adj=1,"s(WY)")
mtext(side=1,line=2,"WY")
mtext(side=2,line=2.5,"Effect")

bbox.lims=bbox(regions2)
tmp.ma=with(tmp[2][[1]],matrix(fit,nrow=length(y),ncol=length(x)))
dat1=list()
dat1$x=tmp[2][[1]]$x
dat1$y=tmp[2][[1]]$y
dat1$z=tmp.ma
r=raster(dat1)

brk=20
breaks.val=classInt::classIntervals(tmp[2][[1]]$fit[is.na(tmp[2][[1]]$fit)==F],style="equal",n=brk)
pal=hcl.colors(n=brk,alpha=0.75)
plot(ENP,ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],border=NA)
image(r,add=T,breaks=breaks.val$brks,col = pal)
plot(shore,add=T,lwd=0.1)
plot(ENP,add=T)
plot(sites.shp2,add=T,cex=0.5,pch=19,col=adjustcolor("red",0.25))
plot(rasterToContour(r),col=adjustcolor("white",0.5),add=T)
box(lwd=1)
mtext(side=3,adj=0,"ti(UTMX,UTMY)")

legend_image=as.raster(matrix(rev(pal),ncol=1))
par(xpd=NA,mar=c(2,1,1,0))
plot(c(0,1),c(0,1),type = 'n', axes = F,ann=F)
rasterImage(legend_image, 0, 0.25, 0.3,0.75)
leg.labs=with(breaks.val,c(format(round(min(brks),1),nsmall=1),"0.0",format(round(max(brks),1),nsmall=1)))
text(x=0.3, y = seq(0.25,0.75,length.out=3), labels = leg.labs,cex=1,pos=4)
text(0.15,0.76,"Effect",pos=3,xpd=NA)
dev.off()

# TP
dat.all.TP.GM2=merge(dat.all.TP.GM,sites.shp,"STATION")
head(dat.all.TP.GM2)
gc()
WY.k=23
loc.k=500
m.TP<-bam(log(TP.GM)~
            s(WY,bs="cr",k=WY.k)+
            s(UTMX,UTMY,bs="ds",k=loc.k,m=c(1,0.5))+
            ti(UTMX,UTMY,WY,d=c(2,1),bs=c("ds","cr"),m = list(c(1, 0.5), NA),k=c(85,24)),
          data=dat.all.TP.GM2,
          nthreads = 6,discrete=T)

summary(m.TP)
qq.gam(m.TP)
nvar=3;layout(matrix(1:nvar,1,nvar))
plot(m.TP,residuals=T,pch=21)

nvar=4;layout(matrix(1:nvar,1,nvar))
gam.check(m.TP)
dev.off()

m.TP.sum=notidy_tidy_gam(m.TP)
m.TP.est=notidy_glance_gam(m.TP)
# write.csv(m.TP.sum,export.path("TP_gam_mod_sum.csv"),row.names=F)
# write.csv(m.TP.est,export.path("TP_gam_mod_est.csv"),row.names=F)


tmp=plot(m.TP,residuals=T)
# png(filename=paste0(plot.path,"GAM_TP_draw_base.png"),width=8,height=3,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(2,2,1,0.5),oma=c(2,2,0.25,0.25),xpd=F);
layout(matrix(1:3,1,3),widths=c(1,1,0.4))

crit=qnorm((1 - 0.95) / 2, lower.tail = FALSE)
ylim.val=c(-1.25,1);by.y=0.5;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(1996,2019);by.x=5;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/by.x)
plot(tmp[1][[1]]$fit~tmp[1][[1]]$x,ylim=ylim.val,xlim=xlim.val,ann=F,axes=F,type="n")
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(tmp[1][[1]],points(jitter(raw,0.5),p.resid,pch=19,col=adjustcolor("dodgerblue1",0.10)))
with(tmp[1][[1]],shaded.range(x,fit-(crit*se),fit+(crit*se),"grey",lty=1))
with(tmp[1][[1]],lines(x,fit,lwd=2))
abline(h=0)
axis_fun(1,xmaj,xmin,xmaj,line=-0.5);axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=3,adj=0,"Total Phosphorus")
mtext(side=3,adj=1,"s(WY)")
mtext(side=1,line=2,"WY")
mtext(side=2,line=2.5,"Effect")

bbox.lims=bbox(regions2)
tmp.ma=with(tmp[2][[1]],matrix(fit,nrow=length(y),ncol=length(x)))
dat1=list()
dat1$x=tmp[2][[1]]$x
dat1$y=tmp[2][[1]]$y
dat1$z=tmp.ma
r=raster(dat1)

brk=20
breaks.val=classInt::classIntervals(tmp[2][[1]]$fit[is.na(tmp[2][[1]]$fit)==F],style="equal",n=brk)
pal=hcl.colors(n=brk,alpha=0.75)
plot(ENP,ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],border=NA)
image(r,add=T,breaks=breaks.val$brks,col = pal)
plot(shore,add=T,lwd=0.1)
plot(ENP,add=T)
plot(sites.shp2,add=T,cex=0.5,pch=19,col=adjustcolor("red",0.25))
plot(rasterToContour(r),col=adjustcolor("white",0.5),add=T)
box(lwd=1)
mtext(side=3,adj=0,"ti(UTMX,UTMY)")

legend_image=as.raster(matrix(rev(pal),ncol=1))
par(xpd=NA,mar=c(2,1,1,0))
plot(c(0,1),c(0,1),type = 'n', axes = F,ann=F)
rasterImage(legend_image, 0, 0.25, 0.3,0.75)
leg.labs=with(breaks.val,c(format(round(min(brks),1),nsmall=1),"0.0",format(round(max(brks),1),nsmall=1)))
text(x=0.3, y = seq(0.25,0.75,length.out=3), labels = leg.labs,cex=1,pos=4)
text(0.15,0.76,"Effect",pos=3,xpd=NA)
dev.off()
# Static Trend and GM maps ------------------------------------------------

cols.val=c("white","red")
# tiff(filename=paste0(plot.path,"TrendMaps_v1.tiff"),width=6.5,height=6,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
# png(filename=paste0(plot.path,"TrendMaps.png"),width=6.5,height=6,units="in",res=200,type="windows",bg="white")
par(family="serif",oma=c(0.25,0.25,0.25,0.25),mar=c(0.1,0.1,0.1,0.1),xpd=F)
layout(matrix(c(1:20),5,4,byrow=T),widths = c(1,0.4,1,0.4))
bbox.lims=bbox(region.mask)
{
# TN
plot(shore,col="cornsilk",border="grey",bg="lightblue",ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.05)
n=10
#int=classIntervals(values(tps.TN.trend)[is.na(values(tps.TN.trend))==F],n=n,style="equal")
int.bks=c(min(values(tps.TN.trend),na.rm=T),0+min(values(tps.TN.trend),na.rm=T)/2,0,0+max(values(tps.TN.trend),na.rm=T)/2,max(values(tps.TN.trend),na.rm=T)) #int$brks;#int.bks[2]=0
pal=hcl.colors(length(int.bks)-1, "viridis", rev = F,alpha=0.7)
#pal=hcl.colors(length(int$brks)-1, "GnBu", rev = TRUE,alpha=0.7)
image(tps.TN.trend,add=T,breaks=int.bks,col = pal)
plot(ENP,add=T,bg=NA,lwd=0.5)
plot(regions2,lty=2,add=T,border="grey80",lwd=0.5)
plot(sites.shp2,add=T,cex=0.40,pch=21,bg=adjustcolor("grey",0.5),col="grey",lwd=0.1)
plot(sites.shp.TN.trend,add=T,pch=21,cex=0.5,bg=adjustcolor(cols.val,0.5)[sites.shp.TN.trend$stat.sig],col=NA);
box(lwd=1)
mapmisc::scaleBar(utm17,"bottomright",bty="n",cex=1,seg.len=4)
#
plot(0:1,0:1,ann=F,axes=F,type="n")
# legend_image=as.raster(matrix(pal,ncol=1))
# text(x=0.25, y = c(0.83,0.47), labels = c(round(min(int.bks),3),round(max(int.bks),3)),cex=0.6,adj=0,pos=4)
# rasterImage(legend_image,0.15,0.45,0.25,0.85)
# text(x=0.15,y=0.95,"TN Thiel-Sen Slope\n(mg N L\u207B\u00B9 Yr\u207B\u00B9)",adj=0,cex=0.75)
int.bks.vals=format(round(int.bks,3),nsmall=3)
labs=c(paste0("< ",int.bks.vals[2]),paste(int.bks.vals[2:3],int.bks.vals[3:4],sep=" - "),paste(paste0(">",int.bks.vals[4])))
n.bks=length(int.bks)-1
bx.val= seq(0.45,0.85,(0.85-0.45)/n.bks)
rect(0.15,bx.val[1:n.bks],0.25,bx.val[2:(n.bks+1)],col=rev(pal),lty=0)
text(x=0.25, y = bx.val[2:(n.bks+1)]-c(mean(diff(bx.val[2:(n.bks+1)]))/2), labels = rev(labs),cex=0.6,adj=0,pos=4)
text(x=0.15,y=0.95,"TN Thiel-Sen Slope\n(mg N L\u207B\u00B9 Yr\u207B\u00B9)",adj=0,cex=0.75)
legend(0.5,0.4,legend=c("Insufficent Data","\u03C1 > 0.05","\u03C1 < 0.05"),
       pch=21,lty=c(NA),lwd=c(0.1),
       col=c("grey"),pt.bg=c("grey",cols.val),
       pt.cex=1,ncol=1,cex=0.6,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=1,title.adj = 0,title="Kendall's trend \u03C1-value")

plot(shore,col="cornsilk",border="grey",bg="lightblue",ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.05,xpd=F)
int=classIntervals(values(tps.TN.GM)[is.na(values(tps.TN.GM))==F],n=n,style="equal")
int.bks=int$brks
pal=hcl.colors(length(int$brks)-1, "Inferno", rev = F,alpha=0.7)
image(tps.TN.GM,add=T,breaks=int.bks,col = pal)
plot(ENP,add=T,bg=NA,lwd=0.5)
plot(regions2,lty=2,add=T,border="grey80",lwd=0.5)
plot(sites.shp.TN.GM,add=T,pch=21,cex=0.25,bg=adjustcolor("white",0.5),col=NA)
box(lwd=1)
plot(0:1,0:1,ann=F,axes=F,type="n")
legend_image=as.raster(matrix(pal,ncol=1))
text(x=0.25, y = c(0.83,0.47), labels = c("< 0.10","1.20"),cex=0.6,adj=0,pos=4)
rasterImage(legend_image,0.15,0.45,0.25,0.85)
text(x=0.15,y=0.95,"Average Annual GM TN\n(mg N L\u207B\u00B9)",adj=0,cex=0.75)

# DIN
plot(shore,col="cornsilk",border="grey",bg="lightblue",ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.05,xpd=F)
int.bks=c(min(values(tps.DIN.trend),na.rm=T),0+min(values(tps.DIN.trend),na.rm=T)/2,0,0+max(values(tps.DIN.trend),na.rm=T)/2,max(values(tps.DIN.trend),na.rm=T)) #int$brks;#int.bks[2]=0
pal=hcl.colors(length(int.bks)-1, "viridis", rev = F,alpha=0.7)
image(tps.DIN.trend,add=T,breaks=int.bks,col = pal)
plot(ENP,add=T,bg=NA,lwd=0.5)
plot(regions2,lty=2,add=T,border="grey80",lwd=0.5)
plot(sites.shp2,add=T,cex=0.40,pch=21,bg=adjustcolor("grey",0.5),col="grey",lwd=0.1)
plot(sites.shp.DIN.trend,add=T,pch=21,cex=0.5,bg=adjustcolor(cols.val,0.5)[sites.shp.DIN.trend$stat.sig],col=NA);
box(lwd=1)
#
plot(0:1,0:1,ann=F,axes=F,type="n")
int.bks.vals=format(round(int.bks,3),nsmall=3)
labs=c(paste0("< ",int.bks.vals[2]),paste(int.bks.vals[2:3],int.bks.vals[3:4],sep=" - "),paste(paste0(">",int.bks.vals[4])))
n.bks=length(int.bks)-1
bx.val= seq(0.45,0.85,(0.85-0.45)/n.bks)
rect(0.15,bx.val[1:n.bks],0.25,bx.val[2:(n.bks+1)],col=rev(pal),lty=0)
text(x=0.25, y = bx.val[2:(n.bks+1)]-c(mean(diff(bx.val[2:(n.bks+1)]))/2), labels = rev(labs),cex=0.6,adj=0,pos=4)
text(x=0.15,y=0.95,"DIN Thiel-Sen Slope\n(mg N L\u207B\u00B9 Yr\u207B\u00B9)",adj=0,cex=0.75)

plot(shore,col="cornsilk",border="grey",bg="lightblue",ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.05,xpd=F)
int=classIntervals(values(tps.DIN.GM)[is.na(values(tps.DIN.GM))==F],n=n,style="equal")
int.bks=int$brks
pal=hcl.colors(length(int$brks)-1, "Inferno", rev = F,alpha=0.7)
image(tps.DIN.GM,add=T,breaks=int.bks,col = pal)
plot(ENP,add=T,bg=NA,lwd=0.5)
plot(regions2,lty=2,add=T,border="grey80",lwd=0.5)
plot(sites.shp.DIN.GM,add=T,pch=21,cex=0.25,bg=adjustcolor("white",0.5),col=NA)
box(lwd=1)
plot(0:1,0:1,ann=F,axes=F,type="n")
legend_image=as.raster(matrix(pal,ncol=1))
text(x=0.25, y = c(0.83,0.47), labels = c("< 0.10","1.00"),cex=0.6,adj=0,pos=4)
rasterImage(legend_image,0.15,0.45,0.25,0.85)
text(x=0.15,y=0.95,"Average Annual GM DIN\n(mg N L\u207B\u00B9)",adj=0,cex=0.75)

# TP
plot(shore,col="cornsilk",border="grey",bg="lightblue",ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.05,xpd=F)
int.bks=c(min(values(tps.TP.trend),na.rm=T),0+min(values(tps.TP.trend),na.rm=T)/2,0,0+max(values(tps.TP.trend),na.rm=T)/2,max(values(tps.TP.trend),na.rm=T)) #int$brks;#int.bks[2]=0
pal=hcl.colors(length(int.bks)-1, "viridis", rev = F,alpha=0.7)
image(tps.TP.trend,add=T,breaks=int.bks,col = pal)
plot(ENP,add=T,bg=NA,lwd=0.5)
plot(regions2,lty=2,add=T,border="grey80",lwd=0.5)
plot(sites.shp2,add=T,cex=0.40,pch=21,bg=adjustcolor("grey",0.5),col="grey",lwd=0.1)
plot(sites.shp.TP.trend,add=T,pch=21,cex=0.5,bg=adjustcolor(cols.val,0.5)[sites.shp.TP.trend$stat.sig],col=NA);
box(lwd=1)
plot(0:1,0:1,ann=F,axes=F,type="n")
int.bks.vals=format(round(int.bks,3),nsmall=3)
labs=c(paste0("< ",int.bks.vals[2]),paste(int.bks.vals[2:3],int.bks.vals[3:4],sep=" - "),paste(paste0(">",int.bks.vals[4])))
n.bks=length(int.bks)-1
bx.val= seq(0.45,0.85,(0.85-0.45)/n.bks)
rect(0.15,bx.val[1:n.bks],0.25,bx.val[2:(n.bks+1)],col=rev(pal),lty=0)
text(x=0.25, y = bx.val[2:(n.bks+1)]-c(mean(diff(bx.val[2:(n.bks+1)]))/2), labels = rev(labs),cex=0.6,adj=0,pos=4)
text(x=0.15,y=0.95,"TP Thiel-Sen Slope\n(\u03BCg P L\u207B\u00B9 Yr\u207B\u00B9)",adj=0,cex=0.75)

plot(shore,col="cornsilk",border="grey",bg="lightblue",ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.05,xpd=F)
int=classIntervals(values(tps.TP.GM)[is.na(values(tps.TP.GM))==F],n=n,style="equal")
int.bks=int$brks
pal=hcl.colors(length(int$brks)-1, "Inferno", rev = F,alpha=0.7)
image(tps.TP.GM,add=T,breaks=int.bks,col = pal)
plot(ENP,add=T,bg=NA,lwd=0.5)
plot(regions2,lty=2,add=T,border="grey80",lwd=0.5)
plot(sites.shp.TP.GM,add=T,pch=21,cex=0.25,bg=adjustcolor("white",0.5),col=NA)
box(lwd=1)
plot(0:1,0:1,ann=F,axes=F,type="n")
legend_image=as.raster(matrix(pal,ncol=1))
text(x=0.25, y = c(0.83,0.47), labels = c("< 2","> 40"),cex=0.6,adj=0,pos=4)
rasterImage(legend_image,0.15,0.45,0.25,0.85)
text(x=0.15,y=0.95,"Average Annual GM TP\n(\u03BCg P L\u207B\u00B9)",adj=0,cex=0.75)

# SRP
plot(shore,col="cornsilk",border="grey",bg="lightblue",ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.05,xpd=F)
int.bks=c(min(values(tps.SRP.trend),na.rm=T),0+min(values(tps.SRP.trend),na.rm=T)/2,0,0+max(values(tps.SRP.trend),na.rm=T)/2,max(values(tps.SRP.trend),na.rm=T)) #int$brks;#int.bks[2]=0
pal=hcl.colors(length(int.bks)-1, "viridis", rev = F,alpha=0.7)
image(tps.SRP.trend,add=T,breaks=int.bks,col = pal)
plot(ENP,add=T,bg=NA,lwd=0.5)
plot(regions2,lty=2,add=T,border="grey80",lwd=0.5)
plot(sites.shp2,add=T,cex=0.40,pch=21,bg=adjustcolor("grey",0.5),col="grey",lwd=0.1)
plot(sites.shp.SRP.trend,add=T,pch=21,cex=0.5,bg=adjustcolor(cols.val,0.5)[sites.shp.SRP.trend$stat.sig],col=NA);
box(lwd=1)
plot(0:1,0:1,ann=F,axes=F,type="n")
int.bks.vals=format(round(int.bks,3),nsmall=3)
labs=c(paste0("< ",int.bks.vals[2]),paste(int.bks.vals[2:3],int.bks.vals[3:4],sep=" - "),paste(paste0(">",int.bks.vals[4])))
n.bks=length(int.bks)-1
bx.val= seq(0.45,0.85,(0.85-0.45)/n.bks)
rect(0.15,bx.val[1:n.bks],0.25,bx.val[2:(n.bks+1)],col=rev(pal),lty=0)
text(x=0.25, y = bx.val[2:(n.bks+1)]-c(mean(diff(bx.val[2:(n.bks+1)]))/2), labels = rev(labs),cex=0.6,adj=0,pos=4)
text(x=0.15,y=0.95,"SRP Thiel-Sen Slope\n(\u03BCg P L\u207B\u00B9 Yr\u207B\u00B9)",adj=0,cex=0.75)

plot(shore,col="cornsilk",border="grey",bg="lightblue",ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.05,xpd=F)
int=classIntervals(values(tps.SRP.GM)[is.na(values(tps.SRP.GM))==F],n=n,style="equal")
int.bks=int$brks
pal=hcl.colors(length(int$brks)-1, "Inferno", rev = F,alpha=0.7)
image(tps.SRP.GM,add=T,breaks=int.bks,col = pal)
plot(ENP,add=T,bg=NA,lwd=0.5)
plot(regions2,lty=2,add=T,border="grey80",lwd=0.5)
plot(sites.shp.SRP.GM,add=T,pch=21,cex=0.25,bg=adjustcolor("white",0.5),col=NA)
box(lwd=1)
plot(0:1,0:1,ann=F,axes=F,type="n")
legend_image=as.raster(matrix(pal,ncol=1))
text(x=0.25, y = c(0.83,0.47), labels = c("< 1","12"),cex=0.6,adj=0,pos=4)
rasterImage(legend_image,0.15,0.45,0.25,0.85)
text(x=0.15,y=0.95,"Average Annual GM SRP\n(\u03BCg P L\u207B\u00B9)",adj=0,cex=0.75)

# Chla
plot(shore,col="cornsilk",border="grey",bg="lightblue",ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.05,xpd=F)
int.bks=c(min(values(tps.Chla.trend),na.rm=T),0+min(values(tps.Chla.trend),na.rm=T)/2,0,0+max(values(tps.Chla.trend),na.rm=T)/2,max(values(tps.Chla.trend),na.rm=T)) #int$brks;#int.bks[2]=0
pal=hcl.colors(length(int.bks)-1, "viridis", rev = F,alpha=0.7)
image(tps.Chla.trend,add=T,breaks=int.bks,col = pal)
plot(ENP,add=T,bg=NA,lwd=0.5)
plot(regions2,lty=2,add=T,border="grey80",lwd=0.5)
plot(sites.shp2,add=T,cex=0.40,pch=21,bg=adjustcolor("grey",0.5),col="grey",lwd=0.1)
plot(sites.shp.Chla.trend,add=T,pch=21,cex=0.5,bg=adjustcolor(cols.val,0.5)[sites.shp.Chla.trend$stat.sig],col=NA);
box(lwd=1)
plot(0:1,0:1,ann=F,axes=F,type="n")
int.bks.vals=format(round(int.bks,3),nsmall=3)
labs=c(paste0("< ",int.bks.vals[2]),paste(int.bks.vals[2:3],int.bks.vals[3:4],sep=" - "),paste(paste0(">",int.bks.vals[4])))
n.bks=length(int.bks)-1
bx.val= seq(0.45,0.85,(0.85-0.45)/n.bks)
rect(0.15,bx.val[1:n.bks],0.25,bx.val[2:(n.bks+1)],col=rev(pal),lty=0)
text(x=0.25, y = bx.val[2:(n.bks+1)]-c(mean(diff(bx.val[2:(n.bks+1)]))/2), labels = rev(labs),cex=0.6,adj=0,pos=4)
text(x=0.15,y=0.95,"Chl-a Thiel-Sen Slope\n(\u03BCg L\u207B\u00B9 Yr\u207B\u00B9)",adj=0,cex=0.75)

plot(shore,col="cornsilk",border="grey",bg="lightblue",ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.05,xpd=F)
int=classIntervals(values(tps.Chla.GM)[is.na(values(tps.Chla.GM))==F],n=n,style="equal")
int.bks=int$brks
pal=hcl.colors(length(int$brks)-1, "Inferno", rev = F,alpha=0.7)
image(tps.Chla.GM,add=T,breaks=int.bks,col = pal)
plot(ENP,add=T,bg=NA,lwd=0.5)
plot(regions2,lty=2,add=T,border="grey80",lwd=0.5)
plot(sites.shp.Chla.GM,add=T,pch=21,cex=0.25,bg=adjustcolor("white",0.5),col=NA)
box(lwd=1)
plot(0:1,0:1,ann=F,axes=F,type="n")
legend_image=as.raster(matrix(pal,ncol=1))
text(x=0.25, y = c(0.83,0.47), labels = c("< 0.2","1.3"),cex=0.6,adj=0,pos=4)
rasterImage(legend_image,0.15,0.45,0.25,0.85)
text(x=0.15,y=0.95,"Average Annual GM Chl-a\n(\u03BCg L\u207B\u00B9)",adj=0,cex=0.75)
}
dev.off()

## ggplot version
library(ggplot2)
library(ggmap)
library(ggspatial)
library(ggsn)
library(cowplot)

theme_map <- function(...) {
  theme_minimal() +
    theme(
      text = element_text(family = "serif", color = "#22211d"),
      axis.line = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      # panel.grid.minor = element_line(color = "#ebebe5", size = 0.2),
      # panel.grid.major = element_line(color = "#ebebe5", size = 0.2),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      # plot.background = element_rect(fill = "lightblue", color = NA),
      panel.background = element_rect(fill = "lightblue", color = NA),
      # legend.background = element_rect(fill = "#f5f5f2", color = NA),
      plot.background = element_blank(),
      # panel.background = element_blank(),
      legend.background = element_blank(),
      panel.border = element_blank(),
      plot.title=element_text(size=12),
      plot.subtitle = element_text(color = "grey50",size=8),
      plot.caption = element_text(hjust = 0),
      ...
    )
}

##
shore2=sf::st_as_sf(shore)
ENP2=sf::st_as_sf(ENP)
regions2_gg=sf::st_as_sf(regions2)
sites.shp2_gg=sf::st_as_sf(sites.shp2)
sites.shp.TN.trend_gg=sf::st_as_sf(sites.shp.TN.trend)

bbox.lims=bbox(region.mask)

TN.trend.gg=ggplot()+
  geom_sf(data=shore2,fill="cornsilk",colour="grey",size=0.1)+
  layer_spatial(tps.TN.trend,alpha=0.75)+scale_fill_viridis(na.value=NA)+
  geom_sf(data=ENP2,fill=NA,colour="black")+
  geom_sf(data=regions2_gg,fill=NA,colour="grey80")+
  geom_sf(data=sites.shp.TN.trend_gg,aes(colour=stat.sig),size=2,shape=19,alpha=0.5)+
  scale_colour_manual(name = "Kendall's trend \u03C1-value",
                      labels=c("not-sig"="\u03C1 > 0.05","sig"="\u03C1 < 0.05"),
                      values = c("not-sig" = "grey","sig"="red"))+
  theme_map()+
  coord_sf(xlim=c(bbox.lims[1,1],bbox.lims[1,2]),ylim=c(bbox.lims[2,1],bbox.lims[2,2]))+
  labs(fill="TN Thiel-Sen Slope\n(mg N L\u207B\u00B9 Yr\u207B\u00B9)")

TN.conc.gg=ggplot()+
  geom_sf(data=shore2,fill="cornsilk",colour="grey",size=0.1)+
  layer_spatial(tps.TN.GM,alpha=0.75)+scale_fill_viridis(na.value=NA,option="A")+
  geom_sf(data=ENP2,fill=NA,colour="black")+
  geom_sf(data=regions2_gg,fill=NA,colour="grey80")+
  geom_sf(data=sf::st_as_sf(sites.shp.TN.GM),colour="white",size=1,alpha=0.5,shape=19)+
  theme_map()+
  coord_sf(xlim=c(bbox.lims[1,1],bbox.lims[1,2]),ylim=c(bbox.lims[2,1],bbox.lims[2,2]))+
  labs(fill="Average Annual GM TN\n(mg N L\u207B\u00B9)")

DIN.trend.gg=ggplot()+
  geom_sf(data=shore2,fill="cornsilk",colour="grey",size=0.1)+
  layer_spatial(tps.DIN.trend,alpha=0.75)+scale_fill_viridis(na.value=NA)+
  geom_sf(data=ENP2,fill=NA,colour="black")+
  geom_sf(data=regions2_gg,fill=NA,colour="grey80")+
  geom_sf(data=sf::st_as_sf(sites.shp.DIN.trend),aes(colour=stat.sig),size=2,shape=19,alpha=0.5)+
  scale_colour_manual(name = "Kendall's trend \u03C1-value",
                      values = c("not-sig" = "grey","sig"="red"),guide=F)+
  theme_map()+
  coord_sf(xlim=c(bbox.lims[1,1],bbox.lims[1,2]),ylim=c(bbox.lims[2,1],bbox.lims[2,2]))+
  labs(fill="DIN Thiel-Sen Slope\n(mg N L\u207B\u00B9 Yr\u207B\u00B9)")

DIN.conc.gg=ggplot()+
  geom_sf(data=shore2,fill="cornsilk",colour="grey",size=0.1)+
  layer_spatial(tps.DIN.GM,alpha=0.75)+scale_fill_viridis(na.value=NA,option="A")+
  geom_sf(data=ENP2,fill=NA,colour="black")+
  geom_sf(data=regions2_gg,fill=NA,colour="grey80")+
  geom_sf(data=sf::st_as_sf(sites.shp.DIN.GM),colour="white",size=1,alpha=0.5,shape=19)+
  theme_map()+
  coord_sf(xlim=c(bbox.lims[1,1],bbox.lims[1,2]),ylim=c(bbox.lims[2,1],bbox.lims[2,2]))+
  labs(fill="Average Annual GM DIN\n(mg N L\u207B\u00B9)")


comboplot=plot_grid(
  TN.trend.gg,TN.conc.gg,
  DIN.trend.gg,
  DIN.conc.gg,
  ncol=2,
  nrow=2)
comboplot
TN.trend.gg
ggsave(paste0(plot.path,"TrendMaps_vgg.tiff"),comboplot,device="tiff",height =7,width=6.5,units="in")
ggsave(paste0(plot.path,"TNtrendgg_vgg.tiff"),TN.trend.gg,device="tiff",height =4,width=6.5,units="in")


# Regional Summaries ------------------------------------------------------
dev.off()
library(dunn.test)
library(rcompanion)
###
## TN
levels.var=c("ENP","Coastal_Mangroves","FLBay","Shelf","Keys")
levels.var.labs=c("ENP","Mangrove Fringe","Florida Bay","West Florida Shelf","Keys")
dat.all.TN.GM2=merge(dat.all.TN.GM,all.sites.region,"STATION",all.x=T)
dat.all.TN.GM2$Region=factor(dat.all.TN.GM2$Region,levels=levels.var)

## DIN
dat.all.DIN.GM2=merge(dat.all.DIN.GM,all.sites.region,"STATION",all.x=T)
dat.all.DIN.GM2$Region=factor(dat.all.DIN.GM2$Region,levels=levels.var)

## TP
dat.all.TP.GM2=merge(dat.all.TP.GM,all.sites.region,"STATION",all.x=T)
dat.all.TP.GM2$Region=factor(dat.all.TP.GM2$Region,levels=c("ENP","Coastal_Mangroves","FLBay","Shelf","Keys"))

## SRP
dat.all.SRP.GM2=merge(dat.all.SRP.GM,all.sites.region,"STATION",all.x=T)
dat.all.SRP.GM2$Region=factor(dat.all.SRP.GM2$Region,levels=c("ENP","Coastal_Mangroves","FLBay","Shelf","Keys"))

## Chl-a
dat.all.Chla.GM2=merge(dat.all.Chla.GM,all.sites.region,"STATION",all.x=T)
dat.all.Chla.GM2$Region=factor(dat.all.Chla.GM2$Region,levels=c("ENP","Coastal_Mangroves","FLBay","Shelf","Keys"))
boxplot(Chla.GM~Region,dat.all.Chla.GM2,outline=F)


cols=c("white",adjustcolor(wesanderson::wes_palette("Zissou1",4,"continuous"),0.5))
levels.var.labs=c("ENP\n"," Mangrove\nFringe"," Florida\nBay","W. Florida\nShelf"," Keys\n")
# png(filename=paste0(plot.path,"RegionComp.png"),width=6.5,height=5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,3.5,0.5,0.75),oma=c(3,2,1,0.5));
layout(matrix(c(1:6),3,2,byrow=T))

ylim.val=c(0,2);by.y=0.5;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
x=boxplot(TN.GM~Region,dat.all.TN.GM2,outline=F,ylim=ylim.val,axes=F,ann=F,col=cols,boxwex=0.5)
TN.DT=with(dat.all.TN.GM2,dunn.test(TN.GM, Region))
TN.DT.ltr=cldList(P.adjusted ~ comparison,data=TN.DT,threshold = 0.05)
TN.DT.ltr$Letter=toupper(TN.DT.ltr$Letter)
TN.DT.ltr=TN.DT.ltr[order(match(TN.DT.ltr$Group,levels.var)),]
text(1:5,x$stats[5,],TN.DT.ltr$Letter,pos=3)
axis_fun(2,ymaj,ymin,format(ymaj))
axis_fun(1,1:5,1:5,NA);box(lwd=1)
mtext(side=2,line=2.5,"TN (mg N L\u207B\u00B9)")

ylim.val=c(0,2.1);by.y=0.5;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
x=boxplot(DIN.GM~Region,dat.all.DIN.GM2,outline=F,ylim=ylim.val,axes=F,ann=F,col=cols,boxwex=0.5)
DIN.DT=with(dat.all.DIN.GM2,dunn.test(DIN.GM, Region))
DIN.DT.ltr=cldList(P.adjusted ~ comparison,data=DIN.DT,threshold = 0.05)
DIN.DT.ltr$Letter=toupper(DIN.DT.ltr$Letter)
DIN.DT.ltr=DIN.DT.ltr[order(match(DIN.DT.ltr$Group,levels.var)),]
text(1:5,x$stats[5,],DIN.DT.ltr$Letter,pos=3)
axis_fun(2,ymaj,ymin,format(ymaj))
axis_fun(1,1:5,1:5,NA);box(lwd=1)
mtext(side=2,line=2.5,"DIN (mg N L\u207B\u00B9)")

ylim.val=c(0,50);by.y=10;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
x=boxplot(TP.GM~Region,dat.all.TP.GM2,outline=F,ylim=ylim.val,axes=F,ann=F,col=cols,boxwex=0.5)
TP.DT=with(dat.all.TP.GM2,dunn.test(TP.GM, Region))
TP.DT.ltr=cldList(P.adjusted ~ comparison,data=TP.DT,threshold = 0.05)
TP.DT.ltr$Letter=toupper(TP.DT.ltr$Letter)
TP.DT.ltr=TP.DT.ltr[order(match(TP.DT.ltr$Group,levels.var)),]
text(1:5,x$stats[5,],TP.DT.ltr$Letter,pos=3)
axis_fun(2,ymaj,ymin,format(ymaj))
axis_fun(1,1:5,1:5,NA);box(lwd=1)
mtext(side=2,line=2.5,"TP (\u03BCg P L\u207B\u00B9)")

ylim.val=c(0,10);by.y=2;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
x=boxplot(SRP.GM~Region,dat.all.SRP.GM2,outline=F,ylim=ylim.val,axes=F,ann=F,col=cols,boxwex=0.5)
SRP.DT=with(dat.all.SRP.GM2,dunn.test(SRP.GM, Region))
SRP.DT.ltr=cldList(P.adjusted ~ comparison,data=SRP.DT,threshold = 0.05)
SRP.DT.ltr$Letter=toupper(SRP.DT.ltr$Letter)
SRP.DT.ltr=SRP.DT.ltr[order(match(SRP.DT.ltr$Group,levels.var)),]
text(1:5,x$stats[5,],SRP.DT.ltr$Letter,pos=3)
axis_fun(2,ymaj,ymin,format(ymaj))
axis_fun(1,1:5,1:5,levels.var.labs,line=0.3,cex=0.9);box(lwd=1)
mtext(side=2,line=2.5,"SRP (\u03BCg P L\u207B\u00B9)")

ylim.val=c(0,3);by.y=1;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
x=boxplot(Chla.GM~Region,dat.all.Chla.GM2,outline=F,ylim=ylim.val,axes=F,ann=F,col=cols,boxwex=0.5)
Chla.DT=with(dat.all.Chla.GM2,dunn.test(Chla.GM, Region))
Chla.DT.ltr=cldList(P.adjusted ~ comparison,data=Chla.DT,threshold = 0.05)
Chla.DT.ltr$Letter=toupper(Chla.DT.ltr$Letter)
Chla.DT.ltr=Chla.DT.ltr[order(match(Chla.DT.ltr$Group,levels.var)),]
text(1:5,x$stats[5,],Chla.DT.ltr$Letter,pos=3)
axis_fun(2,ymaj,ymin,format(ymaj))
axis_fun(1,1:5,1:5,levels.var.labs,line=0.3,cex=0.9);box(lwd=1)
mtext(side=2,line=2.5,"Chl-a (\u03BCg L\u207B\u00B9)")
mtext(side=1,line=1,outer=T,"Region")
dev.off()

# Hydroperiod -------------------------------------------------------------

stage.sites=data.frame(SITE=c("NP201","NESRS1","NESRS2","NP203","P33","P36","NP-112","TSB","R127","TSH"),
                       Region=c(rep("UpperSRS",3),rep("MidSRS",2),"LowSRS",rep("UpperTS",2),"MidTS","LowTS"),
                       min.GndElev.NAVD88=c(5.32,4.27,4.19,3.01,3.81,1.65,1.98,1.29,0.00,-0.15),
                       NAVD88.to.NGVD29=c(-1.51,-1.53,-1.54,-1.51,-1.51,-1.51,-1.51,-1.51,-1.57,-1.56))
stage.sites$min.GndElev.NGVD29=with(stage.sites,min.GndElev.NAVD88-NAVD88.to.NGVD29)
stage.sites

stg.dbkeys=data.frame(SITE=c("NP201","NESRS1","NESRS2","NP203","P33","P36","NP-112","TSB","TSH"),
                      DBKEY=c("06719","01140","01218","G6154","06717","06718","H2427","H2442","07090"))


# download R127 from EDEN webpage
# https://sofia.usgs.gov/eden/station.php?stn_name=R127
R127=read.csv(paste0(data.path,"EDEN/1604536339_water_level.csv"),skip=4)
R127$SITE="R127"
R127$Date=date.fun(as.character(R127$Date))
R127$Data.Value=R127$R127.Daily.median.water.Level.feet.NAVD88.-(-1.57)
R127=R127[,c("SITE","Date","Data.Value")]

stg.dat=data.frame()
for(i in 1:nrow(stg.dbkeys)){
  tmp=DBHYDRO_daily(dates[1],dates[2],stg.dbkeys$DBKEY[i])
  tmp$DBKEY=as.character(stg.dbkeys$DBKEY[i])
  stg.dat=rbind(stg.dat,tmp)
  print(i)
}
stg.dat=merge(stg.dat,stg.dbkeys,"DBKEY")

stg.dat=stg.dat[,c("SITE","Date","Data.Value")]
stg.dat=rbind(stg.dat,R127)
stg.dat=merge(stg.dat,stage.sites[,c("SITE","Region","min.GndElev.NGVD29")],"SITE")
stg.dat$WL.GT.Gnd=with(stg.dat,ifelse(Data.Value>min.GndElev.NGVD29,1,0))
stg.dat$Depth.cm=with(stg.dat,ft.to.m(Data.Value-min.GndElev.NGVD29)*100)
stg.dat$WY=WY(stg.dat$Date)

stg.dat.WY.HP=ddply(stg.dat,c("SITE","Region","WY"),summarise,HP.days=sum(WL.GT.Gnd,na.rm=T),N.val=N.obs(SITE),meanDepth=mean(Depth.cm,na.rm=T))
stg.dat.WY.HP=subset(stg.dat.WY.HP,WY%in%seq(1996,2019,1))
stg.dat.WY.HP$HP.freq=with(stg.dat.WY.HP,HP.days/N.val)
subset(stg.dat.WY.HP,N.val!=365)

plot(HP.freq~WY,subset(stg.dat.WY.HP,SITE=="NP-112"))
plot(meanDepth~WY,subset(stg.dat.WY.HP,SITE=="NP-112"))
plot(HP.freq~WY,subset(stg.dat.WY.HP,SITE=="R127"))
plot(meanDepth~WY,subset(stg.dat.WY.HP,SITE=="R127"))
stg.dat.WY.HP.region=ddply(stg.dat.WY.HP,c("WY","Region"),summarise,meanHP.f=mean(HP.freq,na.rm=T),meanZ=mean(meanDepth,na.rm=T))

ddply(stg.dat.WY.HP.region,c("Region"),summarise,
      est=as.numeric(cor.test(WY,meanZ,method="kendall")$estimate),
      pval=cor.test(WY,meanZ,method="kendall")$p.value)



# Discharge Data ----------------------------------------------------------
SRS.dbkeys=data.frame(SITE=c(rep("S12A",3),rep("S12B",3),rep("S12C",3),rep("S12D",3),rep("S333",2),rep("S334",2),"S355A","S355B_P","S355B","S356","S335"),
                      DBKEY=c(c("FE771","P0796","03620"),c("FE772","P0950","03626"),c("FE773","P0951","03632"),c("FE774","P0952","03638"),c("15042","91487"),c("FB752","91488"),"MQ895","AM173","MQ896","64136","91489"),
                      Priority=paste0("P",c(rep(c(1,2,3),4),c(1,2),c(1,2),rep(1,5))),
                      WQSite=c(rep("S12A",3),rep("S12B",3),rep("S12C",3),rep("S12D",3),rep("S333",2),rep("S334",2),"S355A","S355B","S355B","S356-334","S356-334"))


srs.flow=data.frame()
for(i in 1:nrow(SRS.dbkeys)){
  tmp=DBHYDRO_daily(dates[1],dates[2],SRS.dbkeys$DBKEY[i])
  tmp$DBKEY=SRS.dbkeys$DBKEY[i]
  srs.flow=rbind(tmp,srs.flow)
  print(i)
}
srs.flow$WY=WY(srs.flow$Date)
srs.flow=merge(srs.flow,SRS.dbkeys,"DBKEY")
flow.xtab=data.frame(cast(srs.flow,Date+WY+SITE+WQSite~Priority,value="Data.Value",fun.aggregate=function(x) ifelse(sum(x,na.rm=T)==0,NA,sum(x,na.rm=T))))
flow.xtab$fflow.cfs=with(flow.xtab,ifelse(is.na(P1)==T&is.na(P2)==T,P3,ifelse(is.na(P2)==T&is.na(P3)==T,P1,ifelse(is.na(P1)==T&is.na(P3)==T,P2,P1))));#final flow value for analysis

srs.flow.da.xtab=cast(flow.xtab,Date+WY~SITE,value="fflow.cfs",fun.aggregate=function(x) mean(cfs.to.km3d(ifelse(x<0,NA,x)),na.rm=T))
srs.flow.da.xtab$S355B=rowSums(srs.flow.da.xtab[,c("S355B","S355B_P")],na.rm=T)
srs.flow.da.xtab=srs.flow.da.xtab[,c("Date","WY","S12A","S12B", "S12C", "S12D", "S333", "S334", "S355A", "S355B","S356","S335")]

plot(S12A~Date,srs.flow.da.xtab)
srs.flow.da.xtab[is.na(srs.flow.da.xtab)]<-0
srs.flow.da.xtab$minS356_S335=apply(srs.flow.da.xtab[,c('S356','S335')],1,min); # with(srs.flow.da.xtab,ifelse(S356>S335,S335,S356));
srs.flow.da.xtab$m15.ESRS.Q=rowSums(srs.flow.da.xtab[,c("S333","S355A","S355B","minS356_S335")],na.rm=T)
# srs.flow.da.xtab$S333R.m15=with(srs.flow.da.xtab,S333*(1-ifelse(m15.ESRS.Q>0,ifelse(S334/m15.ESRS.Q>1,1,S334/m15.ESRS.Q),0)))
# srs.flow.da.xtab$S355A.adj.m15=with(srs.flow.da.xtab,S355A*(1-ifelse(m15.ESRS.Q>0,ifelse(S334/m15.ESRS.Q>1,1,S334/m15.ESRS.Q),0)))
# srs.flow.da.xtab$S355B.adj.m15=with(srs.flow.da.xtab,S355B*(1-ifelse(m15.ESRS.Q>0,ifelse(S334/m15.ESRS.Q>1,1,S334/m15.ESRS.Q),0)))
# srs.flow.da.xtab$S356.adj.m15=with(srs.flow.da.xtab,minS356_S335*(1-ifelse(m15.ESRS.Q>0,ifelse(S334/m15.ESRS.Q>1,1,S334/m15.ESRS.Q),0)))
srs.flow.da.xtab$TFlow.m15=rowSums(srs.flow.da.xtab[,c("S12A","S12B","S12C","S12D","S333","S355A","S355B","minS356_S335")],na.rm=T)
srs.flow.da.xtab.WY=ddply(subset(srs.flow.da.xtab,WY%in%seq(1996,2019,1)),"WY",summarise,TFlow.km3=sum(TFlow.m15,na.rm=T))

plot(TFlow.km3~WY,srs.flow.da.xtab.WY,type="b")
with(srs.flow.da.xtab.WY,cor.test(WY,TFlow.km3,method="kendall"))

# https://www.sfwmd.gov/sites/default/files/documents/ag_item_4_TSCB_Method3.pdf
# https://www.sfwmd.gov/sites/default/files/documents/TS%26CB_TP_Compliance_1st_qtr_2020_Data_Report_method1%262%263_0.pdf

# TSCB
TS.dbkeys=data.frame(SITE=c("S332","S175","S332D","S332D","S332DX1","S328","G737","S18C"),
                     DBKEY=c("91486","15752","TA413","91485","91484","AN558","AN674","15760"),
                     Priority=c("P1","P1","P1","P2","P1","P1","P1","P1"))
TS.flow=data.frame()
for(i in 1:nrow(TS.dbkeys)){
  tmp=DBHYDRO_daily(dates[1],dates[2],TS.dbkeys$DBKEY[i])
  tmp$DBKEY=TS.dbkeys$DBKEY[i]
  TS.flow=rbind(tmp,TS.flow)
  print(i)
}
TS.flow$WY=WY(TS.flow$Date)
TS.flow=merge(TS.flow,TS.dbkeys,"DBKEY")
TS.flow$Date=date.fun(TS.flow$Date)

ts.flow.xtab=data.frame(cast(TS.flow,Date+WY+SITE~Priority,value="Data.Value",fun.aggregate=function(x) ifelse(sum(x,na.rm=T)==0,NA,sum(x,na.rm=T))))
ts.flow.xtab$P3=NA
ts.flow.xtab$fflow.cfs=with(ts.flow.xtab,ifelse(is.na(P1)==T&is.na(P2)==T,P3,ifelse(is.na(P2)==T&is.na(P3)==T,P1,ifelse(is.na(P1)==T&is.na(P3)==T,P2,P1))));#final flow value for analysis
ts.flow.xtab$fflow.cfs=with(ts.flow.xtab,ifelse(SITE%in%c("S332","S175")&Date>date.fun("1999-08-30"),NA,fflow.cfs))

ts.flow.da.xtab=cast(ts.flow.xtab,Date+WY~SITE,value="fflow.cfs",fun.aggregate=function(x) mean(cfs.to.km3d(ifelse(x<0,NA,x)),na.rm=T))
ts.flow.da.xtab[is.na(ts.flow.da.xtab)]<-0
ts.flow.da.xtab$DBasin.Flow=with(ts.flow.da.xtab,(S332D-S332DX1-S328))
ts.flow.da.xtab$DBasin.Flow=with(ts.flow.da.xtab,ifelse(DBasin.Flow<0,0,DBasin.Flow))
ts.flow.da.xtab$TFlow=with(ts.flow.da.xtab,DBasin.Flow+S328+G737+S18C+S332+S175)

ts.flow.da.xtab.WY=ddply(subset(ts.flow.da.xtab,WY%in%seq(1996,2019,1)),"WY",summarise,TFlow.km3=sum(TFlow,na.rm=T))

plot(TFlow.km3~WY,ts.flow.da.xtab.WY,type="b")
with(ts.flow.da.xtab.WY,cor.test(WY,TFlow.km3,method="kendall"))


##
SRS.hydro=merge(subset(stg.dat.WY.HP.region,Region%in%c("UpperSRS","MidSRS","LowSRS")),
                srs.flow.da.xtab.WY,"WY")
range(SRS.hydro$meanHP,na.rm=T)
ddply(SRS.hydro,"Region",summarise,mean.val=mean(meanHP.f))
range(SRS.hydro$TFlow.km3)

plot(meanHP.f~TFlow.km3,SRS.hydro,type="n")
with(subset(SRS.hydro,Region=="UpperSRS"),points(TFlow.km3,meanHP.f,pch=21,bg="indianred1"))
with(subset(SRS.hydro,Region=="MidSRS"),points(TFlow.km3,meanHP.f,pch=21,bg="dodgerblue1"))
with(subset(SRS.hydro,Region=="LowSRS"),points(TFlow.km3,meanHP.f,pch=21,bg="forestgreen"))

ddply(SRS.hydro,c("Region"),summarise,
      est=as.numeric(cor.test(TFlow.km3,meanHP.f,method="spearman")$estimate),
      pval=cor.test(TFlow.km3,meanHP.f,method="spearman")$p.value)

srs.HP.Q=lm(meanHP.f~TFlow.km3+WY+Region,SRS.hydro)
layout(matrix(1:4,2,2));plot(srs.HP.Q)
gvlma::gvlma(srs.HP.Q)

TS.hydro=merge(subset(stg.dat.WY.HP.region,Region%in%paste0(c("Upper","Mid","Low"),"TS")),
               ts.flow.da.xtab.WY,"WY")

range(TS.hydro$meanHP.f,na.rm=T)
ddply(TS.hydro,"Region",summarise,mean.val=mean(meanHP.f))
range(TS.hydro$TFlow.km3)

ddply(TS.hydro,c("Region"),summarise,
      est=as.numeric(cor.test(TFlow.km3,meanHP.f,method="spearman")$estimate),
      pval=cor.test(TFlow.km3,meanHP.f,method="spearman")$p.value)


### hydro figure
# png(filename=paste0(plot.path,"Hydro_Trend.png"),width=6.5,height=4,units="in",res=200,bg="white")
par(family="serif",oma=c(0.5,0.5,0.5,1),mar=c(0.1,0.1,0.1,0.1),xpd=F)
layout(matrix(c(1:6),2,3,byrow=T))

SRS.stg.site=c("NP-201","NESRS1","NESRS2","NP-203","NP-P33","NP-P36")
TS.stg.site=c("NP-112","NP-TSB","NP-TSH")

bbox.lims=bbox(subset(sloughs,NAME=="Shark River"))
plot(ENP.shore,col="white",border="grey",lwd=0.05,ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)])
plot(sloughs.clp,col="grey90",border=NA,add=T)
plot(canal,col="grey",add=T)
plot(ENP,bg=NA,lwd=0.5,add=T)
plot(subset(structures,NAME%in%c(paste0("S12",LETTERS[1:4]),"S333","S334","S356","S355A","S355B")),add=T,pch=21,bg="black")
tmp=subset(wmd.mon,STATION%in%SRS.stg.site&ACTIVITY_S=="Stage")
tmp=tmp[order(match(tmp$STATION,SRS.stg.site)),]
plot(tmp,pch=21,bg=c(rep("indianred1",3),rep("dodgerblue1",2),"forestgreen"),add=T,lwd=0.01)
text(tmp,"SITE",pos=1,cex=0.5,halo=T)
mapmisc::scaleBar(utm17,"bottomleft",bty="n",cex=1,seg.len=4)
legend("topleft",legend=c("Discharge","Upper Slough Stage","Mid Slough Stage","Lower Slough Stage"),
       pch=c(21),lty=c(NA),lwd=c(0.01),
       pt.bg=c("black","indianred1","dodgerblue1","forestgreen"),
       pt.cex=1,ncol=1,cex=0.7,bg="white",box.lty=0,y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=1,title.adj = 0,title=" Monitoring Locations")
box(lwd=1)

xlim.val=c(1996,2019);by.x=10;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/by.x)
ylim.val=c(0.5,1);by.y=0.1;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
par(family="serif",mar=c(3,4,0.5,0.5),xpd=F)
plot(meanHP.f~WY,stg.dat.WY.HP.region,type="n",ylim=ylim.val,xlim=xlim.val,ann=F,axes=F)
abline(h=ymaj,v=xmaj,lty=3,col="grey",lwd=0.5)
with(subset(stg.dat.WY.HP.region,Region=="UpperSRS"),pt_line(WY,meanHP.f,1,"indianred1",1,21,"indianred1"))
with(subset(stg.dat.WY.HP.region,Region=="MidSRS"),pt_line(WY,meanHP.f,1,"dodgerblue1",1,21,"dodgerblue1"))
with(subset(stg.dat.WY.HP.region,Region=="LowSRS"),pt_line(WY,meanHP.f,1,"forestgreen",1,21,"forestgreen"))
axis_fun(1,xmaj,xmin,xmaj,line=-0.5)
axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=2,line=2.5,"Hydroperiod (Prop of WY)",cex=0.8)

ylim.val=c(0,3);by.y=1;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(TFlow.km3~WY,srs.flow.da.xtab.WY,type="n",ylim=ylim.val,xlim=xlim.val,ann=F,axes=F)
abline(h=ymaj,v=xmaj,lty=3,col="grey",lwd=0.5)
with(srs.flow.da.xtab.WY,pt_line(WY,TFlow.km3,1,"black",1.25,21,"black"))
axis_fun(1,xmaj,xmin,xmaj,line=-0.5)
axis_fun(2,ymaj,ymin,format(ymaj,nsmall=1));box(lwd=1)
mtext(side=2,line=2.5,"Discharge (km\u00B3 WY\u207B\u00B9)",cex=0.8)
text(xlim.val[2]+2,ylim.val[2],"A",xpd=NA,cex=1.5)

par(mar=c(0.1,0.1,0.1,0.1))
bbox.lims=bbox(subset(sloughs,NAME=="Taylor Sloug"))
plot(ENP.shore,col="white",border="grey",lwd=0.05,ylim=c(bbox.lims[c(2)],2818533),xlim=bbox.lims[c(1,3)])
plot(sloughs.clp,col="grey90",border=NA,add=T)
plot(canal,col="grey",add=T)
plot(ENP,bg=NA,lwd=0.5,add=T)
plot(subset(structures,NAME%in%c("S332","S332D","S332DX1","S328","G737","S18C","S175")),add=T,pch=21,bg="black")
tmp=subset(wmd.mon,STATION%in%TS.stg.site&ACTIVITY_S=="Stage")
tmp=tmp[order(match(tmp$STATION,TS.stg.site)),]
plot(tmp,pch=21,bg=c(rep("indianred1",2),"forestgreen"),add=T,lwd=0.01)
text(tmp,"SITE",pos=2,cex=0.5,halo=T)
EDEN.R127=SpatialPointsDataFrame(data.frame(UTMX=c(539626.9),UTMY=c(2804113.2)),data=data.frame(SITE="R127"),proj4string=CRS(SRS_string="EPSG:26917"))
text(EDEN.R127,"SITE",pos=2,cex=0.5,halo=T)
plot(EDEN.R127,add=T,pch=21,bg="dodgerblue1",lwd=0.01)
mapmisc::scaleBar(utm17,"bottomleft",bty="n",cex=1,seg.len=4)
box(lwd=1)

xlim.val=c(1996,2019);by.x=10;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/by.x)
ylim.val=c(0.5,1);by.y=0.1;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
par(family="serif",mar=c(3,4,0.5,0.5),xpd=F)
plot(meanHP.f~WY,stg.dat.WY.HP.region,type="n",ylim=ylim.val,xlim=xlim.val,ann=F,axes=F)
abline(h=ymaj,v=xmaj,lty=3,col="grey",lwd=0.5)
with(subset(stg.dat.WY.HP.region,Region=="UpperTS"),pt_line(WY,meanHP.f,1,"indianred1",1,21,"indianred1"))
with(subset(stg.dat.WY.HP.region,Region=="MidTS"),pt_line(WY,meanHP.f,1,"dodgerblue1",1,21,"dodgerblue1"))
with(subset(stg.dat.WY.HP.region,Region=="LowTS"),pt_line(WY,meanHP.f,1,"forestgreen",1,21,"forestgreen"))
axis_fun(1,xmaj,xmin,xmaj,line=-0.5)
axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=2,line=2.5,"Hydroperiod (Prop of WY)",cex=0.8)
mtext(side=1,line=1.75,"Water Year",cex=0.8)

ylim.val=c(0,1);by.y=0.2;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(TFlow.km3~WY,ts.flow.da.xtab.WY,type="n",ylim=ylim.val,xlim=xlim.val,ann=F,axes=F)
abline(h=ymaj,v=xmaj,lty=3,col="grey",lwd=0.5)
with(ts.flow.da.xtab.WY,pt_line(WY,TFlow.km3,1,"black",1.25,21,"black"))
axis_fun(1,xmaj,xmin,xmaj,line=-0.5)
axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=2,line=2.5,"Discharge (km\u00B3 WY\u207B\u00B9)",cex=0.8)
mtext(side=1,line=1.75,"Water Year",cex=0.8)
text(xlim.val[2]+2,ylim.val[2],"B",xpd=NA,cex=1.5)
dev.off()
# END ---------------------------------------------------------------------
