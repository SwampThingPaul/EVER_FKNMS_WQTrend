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
library(reshape2)
library(openxlsx)
library(classInt)
library(zoo)

#GIS Libraries
library(sp)
library(rgdal)
library(PROJ)
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
library(DHARMa)
library(marginaleffects)

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

wgs84=CRS("+init=epsg:4326")
utm17=CRS("+init=epsg:26917")

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

tmap_mode("view")

## for GAM periods of change
#  download.file("https://gist.github.com/gavinsimpson/e73f011fdaaab4bb5a30/raw/82118ee30c9ef1254795d2ec6d356a664cc138ab/Deriv.R",
#               paste0(wd,"/src/Deriv.R"))
source(paste0(wd,"/src/Deriv.R"))

dev.data.val=function(model,smooth.id,n=400,var.names,nc=nc){
  sm.term=smooths(model)
  
  df.res <- df.residual(model)
  crit.t <- qt(0.025, df.res, lower.tail = FALSE)
  
  newd <- gratia:::derivative_data(model,
                                   id = smooth.id, n = n,
                                   offset = NULL, order = 1L,
                                   type = "central", eps = 1e-7)
  terms.fit=predict(model,newd,type="terms",se.fit = T,nthreads=nc,discrete=T)
  tmp.fit=terms.fit$fit
  colnames(tmp.fit)=paste("fit",var.names,sep=".")
  tmp.SE=terms.fit$se.fit
  colnames(tmp.SE)=paste("SE",var.names,sep=".")
  newd=cbind(newd,tmp.fit,tmp.SE)
  newd$upper.CI=newd[,paste0("fit.",var.names[smooth.id])] + (crit.t* newd[,paste0("SE.",var.names[smooth.id])])
  newd$lower.CI=newd[,paste0("fit.",var.names[smooth.id])] - (crit.t* newd[,paste0("SE.",var.names[smooth.id])])
  
  model.d <- derivatives(model,newd,term=sm.term[smooth.id],type = "central",interval="confidence",ncores=nc)
  dsig <- signifD(newd[,paste0("fit.",var.names[smooth.id])],
                  d=model.d$derivative,# ver >0.8.1  d=model.d$.derivative,
                  upper=model.d$upper, # upper=model.d$.upper_ci,
                  lower=model.d$lower # lower=model.d$.lower_ci
  )
  dsig.incr=unlist(dsig$incr)
  dsig.decr=unlist(dsig$decr)
  
  newd = merge(newd,model.d,by.x=names(model$model)[smooth.id+1],by.y="data",all.x=T)
  newd=cbind(newd,dsig.incr,dsig.decr)
  return(newd)
}


# load(paste0(export.path,"regiontrend_reanalysis.RData"))

# GIS ---------------------------------------------------------------------
ogrListLayers(paste(gen.GIS,"/SFER_GIS_Geodatabase.gdb",sep=""))
ogrListLayers(paste(gen.GIS,"/AHED_release/AHED_20171102.gdb",sep=""))

structures=spTransform(readOGR(paste(gen.GIS,"/AHED_release/AHED_20171102.gdb",sep=""),"STRUCTURE"),utm17)
wmd.mon=spTransform(readOGR(paste0(gen.GIS,"/SFWMD_Monitoring/20240423"),"DBHYDRO_SITE_STATION"),utm17)
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

shore2=spTransform(readOGR(paste0(gen.GIS,"/FWC"),"FWC_Shoreline"),utm17)
shore2=gSimplify(shore2,500)

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

# bbox(spTransform(regions2,wgs84)); #for metadata; 


# Data --------------------------------------------------------------------
WYs=seq(1996,2019,1)

#from Caccia and Boyer & Rookery Bay SERC report
mdls=data.frame(param=c("TN","TP","TOC","DOC","NOX","NO2","NH4","SRP","CHLA"),
                mdl=c(0.05,0.0012,0.04,0.06,0.0024,0.0003,0.0057,0.0022,0.1))
mdls$mdl.uM=with(mdls,ifelse(param%in%c("TN","NOX","NO2","NH4"),(mdl*1000)/N.mw,
                             ifelse(param%in%c("TP","SRP"),(mdl*1000)/P.mw,
                                    ifelse(param%in%c("TOC","DOC"),(mdl*1000)/C.mw,mdl))))

## FIU SERC ----------------------------------------------------------------
# http://serc.fiu.edu/wqmnetwork/
serc=read.xlsx(paste0(data.path,"FIU_SERC/WQFloridaKeys&Shelf (ppm) UPDATED 6-6-2020.xlsx"),sheet=3)
serc$DATE=date.fun(convertToDate(serc$DATE))
serc=subset(serc,DATE>date.fun("1995-01-01"))

serc.sites=ddply(serc,c("STATION"),summarise,LATDEC=mean(LATDEC,na.rm=T),LONDEC=mean(LONDEC,na.rm=T),n.val=N.obs(SITE),min.date=min(DATE,na.rm=T),max.date=max(DATE,na.rm=T))
serc.sites=subset(serc.sites,is.na(LONDEC)==F)
serc.sites.shp=spTransform(SpatialPointsDataFrame(serc.sites[,c("LONDEC","LATDEC")],data=serc.sites,proj4string=wgs84),utm17)
serc.sites.shp=cbind(serc.sites.shp,coordinates(serc.sites.shp))
plot(coordinates(serc.sites.shp))
# tm_shape(serc.sites.shp)+tm_dots()+tm_shape(ENP_FLB)+tm_dots(col="red",alpha=0.5)+tm_shape(lter)+tm_dots(col="yellow",alpha=0.5)
serc.sites.region=data.frame(sf::st_intersection(sf::st_as_sf(serc.sites.shp),sf::st_as_sf(regions)))[,c("STATION","ESTUARY","SEGMENT_NA","Region")]

serc.sites.region=subset(serc.sites.region,is.na(Region)==F)
serc.sites.region$source="SERC"
serc.sites.region=rbind(serc.sites.region,
                        data.frame(STATION=c(c(316,508,509,506,253,254,504),c(311,303,507,505),500,502,501,503),
                                   ESTUARY=c(rep("Florida Keys",13),"Florida Bay","Florida Bay"),
                                   SEGMENT_NA=c(rep("Lower Keys",7),rep("Back Bay",4),"Upper Keys","Middle Keys","East Central Florida Bay","Southern Florida Bay"),
                                   Region=c(rep("Keys",13),"FLBay","FLBay"),
                                   source="SERC"))
serc.sites.region=merge(serc.sites.region,serc.sites.shp@data[,c("STATION","LONDEC.1","LATDEC.1")],"STATION",all.x=T)
names(serc.sites.region)
colnames(serc.sites.region)=c("STATION", "ESTUARY", "SEGMENT_NA", "Region", "source", "UTMX","UTMY")
subset(serc.sites.region,is.na(Region)==T)

# Surface data only
names(serc)
serc2=serc[,c("DATE","STATION","TN-S","NOX-S","NH4-S","TP-S","SRP-S","CHLA-S","TOC-S","SAL-S")]
serc2=rename(serc2,c("TN-S"="TN","NOX-S"="NOx","NH4-S"="NH4","TP-S"="TP","SRP-S"="SRP","CHLA-S"="Chla","TOC-S"="TOC","SAL-S"="SALINITY"))
serc2[,c("TN","NOx","NH4","TP","SRP","Chla","TOC","SALINITY")]=sapply(serc2[,c("TN","NOx","NH4","TP","SRP","Chla","TOC","SALINITY")],as.numeric)
serc2$TN=with(serc2,ifelse(TN<0.05,0.05/2,TN))
serc2$NOx=with(serc2,ifelse(NOx<0.0024,0.0024/2,NOx))
serc2$NH4=with(serc2,ifelse(NH4<0.0057,0.0057/2,NH4))
serc2$DIN=with(serc2,ifelse(is.na(NH4)==T|is.na(NOx)==T,NA,NH4+NOx))
serc2$TP=with(serc2,ifelse(TP<0.0012,0.0012/2,TP))
serc2$SRP=with(serc2,ifelse(SRP<0.0022,0.0022/2,SRP))
serc2$TOC=with(serc2,ifelse(TOC< 0.04, 0.04/2,TOC))
serc2$DOC=NA # with(serc2,ifelse(DOC< 0.06, 0.06/2,DOC))
serc2=subset(serc2,is.na(TN)==F)
serc2$WY=WY(serc2$DATE)
serc2$season=FL.Hydroseason(serc2$DATE)

chk=ddply(serc2,c("STATION"),summarise,n.val=N.obs(TN),min.date=min(DATE,na.rm=T),max.date=max(DATE,na.rm=T))
chk[chk$max.date>date.fun("2011-07-30"),"STATION"]

serc2=ddply(serc2,c("STATION","DATE","WY","season"),summarise,
            TN=mean(TN,na.rm=T),DIN=mean(DIN,na.rm=T),
            TP=mean(TP,na.rm=T),SRP=mean(SRP,na.rm=T),
            Chla=mean(Chla,na.rm=T),TOC=mean(TOC,na.rm=T),
            DOC=NA,SALINITY=mean(SALINITY,na.rm=T))
serc2=subset(serc2,WY%in%WYs)
head(serc2)
serc2$source="SERC"

## LTER --------------------------------------------------------------------
fce.boyer=read.lter("knb-lter-fce/","1055/9/","acc2be4a77ab1e55b740efdef27648bd",na.string.val=c("-9999","-9999.00","-9999.000","-9999.0"))
fce.boyer=rename(fce.boyer,c("LO.9999EC"="LONDEC"))
fce.boyer.sites=ddply(fce.boyer,c("SITE","LATDEC","LONDEC"),summarise,n.val=N.obs(SITE))
fce.boyer.sites=subset(fce.boyer.sites,n.val>300)
fce.boyer.sites=subset(fce.boyer.sites,is.na(LATDEC)==F|is.na(LONDEC)==F)
fce.boyer.sites=spTransform(SpatialPointsDataFrame(fce.boyer.sites[,c("LONDEC","LATDEC")],data=fce.boyer.sites,proj4string=wgs84),utm17)
# only TS/PH9,10,11

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
fce.wq$TOC=(fce.wq$TOC/1000)*C.mw
range(fce.wq$DOC,na.rm=T)

fce.wq$TN=with(fce.wq,ifelse(TN<0.05,0.05/2,TN))
fce.wq$N.N=with(fce.wq,ifelse(N.N<0.0024,0.0024/2,N.N))
fce.wq$NH4=with(fce.wq,ifelse(NH4<0.0057,0.0057/2,NH4))
fce.wq$DIN=with(fce.wq,ifelse(is.na(NH4)==T|is.na(N.N)==T,NA,NH4+N.N))
fce.wq$TP=with(fce.wq,ifelse(TP<0.0012,0.0012/2,TP))
fce.wq$SRP=with(fce.wq,ifelse(SRP<0.0022,0.0022/2,SRP))

unique(fce.wq$SITENAME)
fce.wq$SITENAME=with(fce.wq,ifelse(substr(SITENAME,1,5)=="TS/Ph",paste0("TS/PH",substr(SITENAME,6,7)),as.character(SITENAME)))

fce.wq$WY=WY(fce.wq$DATE)
fce.wq=subset(fce.wq,WY%in%WYs)
fce.wq$season=FL.Hydroseason(fce.wq$DATE)
fce.wq$Chla=NA
fce.wq$STATION=fce.wq$SITENAME
fce.wq$source="LTER"
head(fce.wq)

fce.wq=fce.wq[,names(serc2)]
lter=cbind(lter,coordinates(lter))
plot(coordinates(lter))
LTER.sites.region=data.frame(sf::st_intersection(sf::st_as_sf(lter),sf::st_as_sf(regions)))[,c("SITE","ESTUARY","SEGMENT_NA","Region")]
LTER.sites.region=subset(LTER.sites.region,is.na(Region)==F)
LTER.sites.region$source="LTER"
LTER.sites.region=rbind(subset(LTER.sites.region,Region=="Coastal_Mangroves"),
                        data.frame(SITE=c(paste("SRS",c("1a","1b","1c","1d",2:3),sep="-"),paste("TS/Ph",c("1a","1b","2","2a","2b",4,3,5,"6a","6b","7a","7b",8:11),sep="-")),
                                   ESTUARY=c(rep("SRS",6),rep("TS",13),rep("FLBay",3)),
                                   SEGMENT_NA=c(rep("SRS_FW",5),rep("SRS_est",1),rep("TS_FW",6),rep("TS_est",7),rep("FLBay",3)),
                                   Region=c(rep("ENP",19),rep("FLBay",3)),
                                   source="LTER"))
LTER.sites.region=merge(LTER.sites.region,lter@data[,c("SITE","coords.x1","coords.x2")],"SITE",all.x=T)
names(LTER.sites.region)
colnames(LTER.sites.region)=c("STATION", "ESTUARY", "SEGMENT_NA", "Region", "source", "UTMX","UTMY")

spl=strsplit(as.character(LTER.sites.region$STATION),"-")
LTER.sites.region$STATION=paste0(toupper(sapply(spl,"[",1)),(sapply(spl,"[",2)))
LTER.sites.region

LTER.sites.region[LTER.sites.region$STATION=="TS/PH2",c("UTMX","UTMY")]=LTER.sites.region[LTER.sites.region$STATION=="TS/PH2a",c("UTMX","UTMY")]
plot(LTER.sites.region[,c("UTMX","UTMY")])

## SFWMD -------------------------------------------------------------------
dates=date.fun(c(paste(min(WYs)-1,"05-01",sep="-"),paste(max(WYs),"05-01",sep="-")))
# ENP_FLB$SITE
# ENP_FLB$STATION

ENP_FLB$STATION
wmd.dat=read.csv(paste0(export.path,"wmd_dat.csv"))
wmd.dat$Date.EST=date.fun(wmd.dat$Date.EST)


ddply(wmd.dat,c("param","Method"),summarise,N.val=N.obs(HalfMDL))

# use max?
ddply(wmd.dat,"param",summarise,min.val=min(MDL,na.rm=T),max.val=max(MDL,na.rm=T))
ddply(wmd.dat,"param",summarise,
      min.val=min(abs(Value[Value<0]),na.rm=T),
      max.val=max(abs(Value[Value<0]),na.rm=T))
mdls

wmd.dat.xtab=reshape2::dcast(wmd.dat,Station.ID+Date.EST~param,value.var="HalfMDL",mean)
wmd.dat.xtab$SALINITY=with(wmd.dat.xtab,ifelse(is.na(SALINITY)==T,SalinityCalc(SPC,TempC),SALINITY))
wmd.dat.xtab$TN=with(wmd.dat.xtab,TN_Combine(NOx,TKN,TN))
wmd.dat.xtab$DIN=with(wmd.dat.xtab,NH4+NOx)

# subset(wmd.dat.xtab,SRP==0)
wmd.dat.xtab$SRP[wmd.dat.xtab$SRP==0]<-NA

wmd.dat.xtab$WY=WY(wmd.dat.xtab$Date.EST)
wmd.dat.xtab=subset(wmd.dat.xtab,WY%in%WYs)
wmd.dat.xtab$season=FL.Hydroseason(wmd.dat.xtab$Date.EST)
wmd.dat.xtab$STATION=wmd.dat.xtab$Station.ID
wmd.dat.xtab$DATE=wmd.dat.xtab$Date.EST

# For analysis purposes, S332 monitoring was terminated and moved to S332D (i.e. inflow to Taylor slough shifted)
wmd.dat.xtab$STATION=with(wmd.dat.xtab,ifelse(STATION=="S332DX","S332D",as.character(STATION)))

# remove obvious outliers
wmd.dat.xtab$TN=with(wmd.dat.xtab,ifelse(TN>10,NA,TN))
wmd.dat.xtab$NH4=with(wmd.dat.xtab,ifelse(NH4>100,NA,NH4))
wmd.dat.xtab$DIN=with(wmd.dat.xtab,ifelse(DIN>10,NA,TN))
wmd.dat.xtab$Chla=with(wmd.dat.xtab,ifelse(Chla>100,NA,Chla))
wmd.dat.xtab$DOC=NA;#double check dbhydro
wmd.dat.xtab$source="DBHydro"

vars=c("STATION", "DATE", "WY", "season", "TN", "DIN", "TP", "SRP","Chla", "TOC", "DOC", "SALINITY", "source")

wmd.dat.melt=melt(wmd.dat.xtab[,vars],id.vars = c('STATION',"DATE","season","WY","source"))
wmd.dat.melt=subset(wmd.dat.melt,is.na(value)==F)

subset(wmd.dat,Station.ID=='FLAB35'&Date.EST==date.fun("2009-10-13"))
subset(wmd.dat.melt,Station.ID=='FLAB35'&Date.EST==date.fun("2009-10-13")&variable=="Chla")

## data from SFWMD to fill early 2000s gap
wmd.dat2=read.csv(paste0(data.path,"SFWMD/CWQMN_thru_Apr2021_For_Analysis.csv"))
names(wmd.dat2)

wmd.names=c("DATE", "month", "Yr", "WY.MON", "WY", "Season", 
            "Station.Num", "Station.ID", "Station.ID.2", "Region", "ZSI", 
            "Temp", "spc_uScm", "Sal", "pH", "DO.mgl", 
            "DO_persat", "Turb", "Chla", "TN.mgL", "NOx.uM", 
            "NO2.uM", "NH4.uM", "TN.uM", "DIN.uM", "TON.uM", 
            "TP.uM", "SRP.uM", "TOC.uM")
colnames(wmd.dat2)=wmd.names
head(wmd.dat2)
summary(wmd.dat2)

wmd.names=c("DATE", "month", "Yr", "WY.MON", "WY", "Season", 
            "Station.Num", "Station.ID", "Station.ID.2", "Region", "ZSI", 
            "Temp", "spc_uScm", "Sal", "pH", "DO.mgl", 
            "DO_persat", "Turb", "Chla", "TN.mgL", "NOx.uM", 
            "NO2.uM", "NH4.uM", "TN.uM", "DIN.uM", "TON.uM", 
            "TP.uM", "SRP.uM", "TOC.uM")
colnames(wmd.dat2)=wmd.names
head(wmd.dat2)
summary(wmd.dat2)

subset(wmd.dat2,Turb<=0)
subset(wmd.dat2,Chla<=0)
subset(wmd.dat2,NOx.uM<=0)
subset(wmd.dat2,NH4.uM<=0)
subset(wmd.dat2,TOC.uM<=0)
subset(wmd.dat2,SRP.uM<=0)
subset(wmd.dat2,TP.uM<=0)

wmd.dat2$Turb[wmd.dat2$Turb<=0]=NA
wmd.dat2$Chla[wmd.dat2$Chla<=0]=NA
wmd.dat2$NOx.uM[wmd.dat2$NOx.uM==0]=NA
wmd.dat2$NH4.uM[wmd.dat2$NH4.uM==0]=NA
wmd.dat2$TOC.uM[wmd.dat2$TOC.uM<=0]=NA
wmd.dat2$TP.uM[wmd.dat2$TP.uM==0]=NA
wmd.dat2$SRP.uM[wmd.dat2$SRP.uM==0]=NA

wmd.dat2$DATE=date.fun(wmd.dat2$DATE,form="%m/%d/%Y")
wmd.dat2$Date.EST=wmd.dat2$DATE#date.fun(wmd.dat2$DATE,form="%m/%d/%Y")
wmd.dat2$DO=wmd.dat2$DO.mgl
wmd.dat2$NH4=with(wmd.dat2,ifelse(NH4.uM*N.mw/1000<0.0008,0.0008/2,round(NH4.uM*N.mw/1000,4)))
wmd.dat2$NOx=with(wmd.dat2,ifelse(NOx.uM*N.mw/1000<0.0098,0.0098/2,round(NOx.uM*N.mw/1000,4)))
wmd.dat2$DIN=with(wmd.dat2,NOx+NH4)
wmd.dat2$TKN=NA
wmd.dat2$TN=with(wmd.dat2,ifelse(TN.uM*N.mw/1000<0.02,0.02/2,round(TN.uM*N.mw/1000,4)))
wmd.dat2$TOC=with(wmd.dat2,ifelse(TOC.uM*C.mw/1000<0.05,0.05/2,round(TOC.uM*C.mw/1000,4)))
wmd.dat2$TP=with(wmd.dat2,ifelse(TP.uM*P.mw/1000<0.002,0.002/2,round(TP.uM*P.mw/1000,4)))
wmd.dat2$SRP=with(wmd.dat2,ifelse(SRP.uM*P.mw/1000<0.002,0.002/2,round(SRP.uM*P.mw/1000,4)))
wmd.dat2$DOC=NA
wmd.dat2$WY=WY(wmd.dat2$Date.EST)
wmd.dat2$source="CWQMN"
wmd.dat2$SALINITY=with(wmd.dat2,SalinityCalc(spc_uScm,Temp))

wmd.dat2$season=FL.Hydroseason(wmd.dat2$Date.EST)
wmd.dat2$STATION=wmd.dat2$Station.ID

vars=c("STATION", "DATE", "WY", "season", "TN", "DIN", "TP", "SRP","Chla", "TOC", "DOC", "SALINITY", "source")

wmd.dat2=wmd.dat2[,vars]
head(wmd.dat2)
head(wmd.dat.xtab)

wmd.dat2.melt=melt(wmd.dat2,id.vars = c('STATION',"DATE","season","WY","source"))
subset(wmd.dat2.melt,value==0)
# unique(subset(wmd.dat2.melt,value==0)$variable)
# wmd.dat2.melt$value=with(wmd.dat2.melt,ifelse(value==0,NA,value))
wmd.dat2.melt=subset(wmd.dat2.melt,is.na(value)==F)

dat.comb=rbind(wmd.dat2.melt,wmd.dat.melt)
dat.comb=reshape2::dcast(dat.comb,STATION+DATE+season+WY+variable~source,value.var = "value",mean)
dat.comb$diff_val=with(dat.comb,CWQMN-DBHydro)

ddply(dat.comb,"variable",summarise,min.dff=min(diff_val,na.rm=T),max.diff=max(diff_val,na.rm=T))

range(dat.comb$diff_val,na.rm=T)
boxplot(dat.comb$diff_val)
subset(dat.comb,diff_val>20)
subset(dat.comb,diff_val<20)

plot(DBHydro~CWQMN,dat.comb);abline(0,1,col="red")

subset(wmd.dat2,STATION=="FLAB35"&DATE==date.fun("2009-10-13"))
subset(wmd.dat.xtab,STATION=='FLAB35'&DATE==date.fun("2009-10-13"))
subset(dat.comb,STATION=="FLAB35"&DATE==date.fun("2009-10-13"))
subset(dat.comb,STATION=="FLAB13"&variable=="Chla")

plot(CWQMN~DATE,subset(dat.comb,STATION=="FLAB13"&variable=="Chla"))
points(DBHydro~DATE,subset(dat.comb,STATION=="FLAB13"&variable=="Chla"),
       pch=21,bg="blue")
dat.comb$value.f=with(dat.comb,ifelse(is.na(DBHydro)==T,CWQMN,DBHydro))

plot(value.f~DATE,subset(dat.comb,STATION=="FLAB13"&variable=="Chla"))
points(DBHydro~DATE,subset(dat.comb,STATION=="FLAB13"&variable=="Chla"),
       pch=21,bg="blue")

wmd.dat.xtab.combo=reshape2::dcast(dat.comb,STATION+DATE+WY+season~variable,value.var = "value.f",mean,na.rm=T)
wmd.dat.xtab.combo$source="WMD"
rm(dat.comb)

ENP_FLB=cbind(ENP_FLB,coordinates(ENP_FLB))
# wmd.sites.region=spatialEco::point.in.poly(ENP_FLB,regions)@data[,c("STATION","ESTUARY","SEGMENT_NA","Region")]
wmd.sites.region=data.frame(sf::st_intersection(sf::st_as_sf(ENP_FLB),sf::st_as_sf(regions)))[,c("STATION","ESTUARY","SEGMENT_NA","Region")]
wmd.sites.region=subset(wmd.sites.region,is.na(Region)==F)
wmd.sites.region$source="WMD"
wmd.sites.region=rbind(wmd.sites.region,
                       data.frame(STATION=c(paste0("S12",LETTERS[1:4]),"S333","SRS1B","SRS1C","NP201","NE1","SRS2",paste0("P",33:36),"G-3273","RG1","CR2","S332D","S332","TSB","S177","S197","S18C","EP","P37"),
                                  ESTUARY=c(rep("SRS",17),rep("TS",8)),
                                  SEGMENT_NA=c(rep("SRS_FW",17),rep("TS_FW",8)),
                                  Region=rep("ENP",25),
                                  source="WMD"))
wmd.sites.region=merge(wmd.sites.region,ENP_FLB@data[,c("STATION","coords.x1","coords.x2")],"STATION",all.x=T)
names(wmd.sites.region)
colnames(wmd.sites.region)=c("STATION", "ESTUARY", "SEGMENT_NA", "Region", "source", "UTMX","UTMY")

nrow(wmd.dat.xtab.combo)
unique(wmd.dat.xtab.combo$STATION)
unique(wmd.sites.region$STATION)
# unique(subset(wmd.dat.xtab.combo,!(STATION%in%wmd.sites.region$STATION))$STATION)
wmd.dat.xtab.combo=subset(wmd.dat.xtab.combo,STATION%in%wmd.sites.region$STATION)

## NOAA --------------------------------------------------------------------
noaa.mdls=data.frame(variable=c("Chla", "NH4.uM","NO2.uM", "NO3.uM", "NOx.uM", "SRP.uM","Si.uM"),
                     MDL=c(0.10,0.028,0.010,0.007,0.015,0.008,0.014))

noaa.dat=read.xlsx(paste0(data.path,"NOAA/AOML_SFP_regional_WQ_surface_v14.xlsx"))
noaa.dat$Date=date.fun(convertToDate(noaa.dat$Date))
unique(noaa.dat$Station)

# Station Name Cleanup
noaa.dat$Station=trimws(noaa.dat$Station)
noaa.dat$Station=toupper(noaa.dat$Station)
noaa.dat[nchar(noaa.dat$Station)>10,"Station"]=c("9.69","9.80","NapleBlue")

noaa.dat$Station=gsub("-","_",noaa.dat$Station)
noaa.dat$Station=gsub(" ","_",noaa.dat$Station)

# ddply(subset(noaa.dat,Station%in%c("21.5/LK","21.5")),c("Station"),summarise,
#       min.lat=min(Latitude),max.lat=max(Latitude),mean.lat=mean(Latitude))
# View(subset(noaa.dat,Station%in%c("21.5/LK","21.5")))
noaa.dat[noaa.dat$Station=="21.5/LK","Station"]=c("21.5")
# ddply(subset(noaa.dat,Station%in%c("21/LK","21")),c("Station"),summarise,
#       min.lat=min(Latitude),max.lat=max(Latitude),mean.lat=mean(Latitude))
noaa.dat[noaa.dat$Station=="21/LK","Station"]=c("21")

noaa.dat$STATION=paste("NOAA",noaa.dat$Station,sep="_")

# noaa.dat$Depth=as.numeric(noaa.dat$Depth)

noaa.vars=c("Record.#", "Longitude", "Latitude", "Station", "Depth", "Date", 
            "GMT", "SST", "SSS", "Chla", "Phaeo", "NH4.uM", "SRP.uM", "NOx.uM", 
            "NO2.uM", "NO3.uM", "Si.uM",'STATION')
colnames(noaa.dat)=noaa.vars

noaa.dat[,c("Chla", "Phaeo", "NH4.uM")]=sapply(noaa.dat[,c("Chla", "Phaeo", "NH4.uM")],as.numeric)

vars1=c("STATION","Date")
params.vars=c("Chla", "NH4.uM", "SRP.uM", "NOx.uM","NO2.uM", "NO3.uM","SSS")
noaa.dat.melt=melt(noaa.dat[,c(vars1,params.vars)],id.vars = vars1)

noaa.dat.melt=merge(noaa.dat.melt,noaa.mdls,"variable",all.x=T)
noaa.dat.melt
subset(noaa.dat.melt,variable=="SSS")
noaa.dat.melt$HalfMDL.uM=with(noaa.dat.melt,ifelse(variable=="SSS",value,ifelse(value<MDL,MDL/2,value)))

noaa.dat.xtab=dcast(noaa.dat.melt,STATION+Date~variable,value.var = "HalfMDL.uM",mean)
noaa.dat.xtab$NOx.uM.calc=with(noaa.dat.xtab,NO2.uM+NO3.uM)
plot(NOx.uM~NOx.uM.calc,noaa.dat.xtab);abline(0,1)

noaa.dat.xtab$NH4=(noaa.dat.xtab$NH4.uM/1000)*N.mw
noaa.dat.xtab$SRP=(noaa.dat.xtab$SRP.uM/1000)*P.mw
noaa.dat.xtab$NOx=(noaa.dat.xtab$NOx.uM/1000)*N.mw
noaa.dat.xtab$DIN=with(noaa.dat.xtab, NH4+NOx)
noaa.dat.xtab$SALINITY=noaa.dat.xtab$SSS

noaa.dat.xtab$DATE=noaa.dat.xtab$Date
noaa.dat.xtab$season=FL.Hydroseason(noaa.dat.xtab$Date)
noaa.dat.xtab$WY=WY(noaa.dat.xtab$Date)

noaa.dat.xtab[,c("TKN","TN","TOC","TP")]=NA
vars=c("STATION","DATE","WY","season","TN","DIN","TP","SRP","Chla","TOC","SALINITY")
noaa.dat.xtab=noaa.dat.xtab[,vars]

noaa.dat.xtab$WY=WY(noaa.dat.xtab$DATE)
noaa.dat.xtab$source="NOAA"
noaa.dat.xtab$season=FL.Hydroseason(noaa.dat.xtab$DATE)
names(wmd.dat.xtab.combo)
names(noaa.dat.xtab)
noaa.dat.xtab=noaa.dat.xtab[,names(wmd.dat.xtab.combo)]

## Locations
noaa.sites=ddply(noaa.dat,c('STATION'),summarise,
                 LAT2=mean(Latitude),LONG2=mean(Longitude),N.val=N.obs(STATION))
unique(noaa.sites$STATION)
noaa.sites.shp=SpatialPointsDataFrame(noaa.sites[,c("LONG2","LAT2")],
                                      data=noaa.sites,proj4string = wgs84)
noaa.sites.shp=spTransform(noaa.sites.shp,utm17)
noaa.sites.shp=cbind(noaa.sites.shp,coordinates(noaa.sites.shp))

noaa.sites.region=data.frame(sf::st_intersection(sf::st_as_sf(noaa.sites.shp),sf::st_as_sf(regions)))[,c("STATION","ESTUARY","SEGMENT_NA","Region")]
noaa.sites.region=subset(noaa.sites.region,is.na(Region)==F)
noaa.sites.region$source="NOAA"
noaa.sites.region

noaa.sites.region=merge(noaa.sites.region,noaa.sites.shp@data[,c("STATION","LONG2.1","LAT2.1")],"STATION",all.x=T)
names(noaa.sites.region)
colnames(noaa.sites.region)=c("STATION", "ESTUARY", "SEGMENT_NA", "Region", "source", "UTMX","UTMY")
plot(noaa.sites.region[,c("UTMX","UTMY")])

noaa.sites.shp=subset(noaa.sites.shp,STATION%in%noaa.sites.region$STATION)
# write.csv(noaa.sites.region,paste0(export.path,"noaa_sites.csv"),row.names = F)

noaa.dat.xtab=subset(noaa.dat.xtab,STATION%in%noaa.sites.region$STATION)

## Combine datasets --------------------------------------------------------
names(serc2)
names(fce.wq)
names(wmd.dat.xtab.combo)
names(noaa.dat.xtab)

vars=c("STATION","DATE","WY","season","TN","DIN","TP","SRP","Chla","TOC","SALINITY","source")
dat.all=rbind(serc2[,vars],fce.wq[,vars],wmd.dat.xtab.combo[,vars],noaa.dat.xtab[,vars])
# write.csv(dat.all,paste0(export.path,"SERC_FCE_WMD_data.csv"),row.names = F)
rm(serc2,fce.wq,wmd.dat.xtab.combo,noaa.dat.xtab)

all.sites.region=rbind(serc.sites.region,LTER.sites.region,wmd.sites.region,noaa.sites.region)
rm(serc.sites.region,LTER.sites.region,wmd.sites.region,noaa.sites.region)

## for fun - sampling density
all.sites.region.shp=SpatialPointsDataFrame(all.sites.region[,c("UTMX","UTMY")],
                                            data=all.sites.region,
                                            proj4string = utm17)
library(SpatialKDE)
band_width=15000
samp.grid=create_grid_rectangular(sf::st_as_sf(all.sites.region.shp),cell_size=2000,side_offset = band_width)
# samp.grid=create_grid_hexagonal(sf::st_as_sf(all.sites.region.shp),cell_size=2000,side_offset = band_width)
kde.val=kde(sf::st_as_sf(all.sites.region.shp),band_width = band_width, kernel = "quartic", grid = samp.grid)
plot(kde.val,main=NA,border=NA)

samp.raster=create_raster(sf::st_as_sf(all.sites.region.shp),cell_size=2000,side_offset = band_width)
kde.val2=kde(sf::st_as_sf(all.sites.region.shp),band_width = band_width, kernel = "quartic", grid = samp.raster)
plot(kde.val2)
### 
# dat.all2=dat.all;# prior analysis removed SFWMD SRS2 site
dat.all2=merge(dat.all,all.sites.region,c("STATION",'source'))


### Reversal screening ------------------------------------------------------
dat.all2$TPReversal=with(dat.all2,ifelse(is.na(SRP)==T|is.na(TP)==T,0,ifelse(SRP>(TP*1.3),1,0)));# Reversals identified as 1 reversals consistent with TP rule evaluation
dat.all2$TNReversal=with(dat.all2,ifelse(is.na(DIN)==T|is.na(TN)==T,0,ifelse(DIN>(TN*1.3),1,0)));

par(family="serif",oma=c(1,1,1,1),mar=c(4,4,1,1))
layout(matrix(1:2,1,2,byrow=F))
plot(TN~DIN,dat.all2,ylab="TN (mg L\u207B\u00b9)",xlab="DIN (mg L\u207B\u00b9)",pch=21,bg=ifelse(TNReversal==1,"dodgerblue1",NA),col=adjustcolor("grey",0.8));abline(0,1,col="dodgerblue1")
plot(TP~SRP,dat.all2,ylab="TP (mg L\u207B\u00b9)",xlab="SRP (mg L\u207B\u00b9)",pch=21,bg=ifelse(TPReversal==1,"red",NA),col=adjustcolor("grey",0.8));abline(0,1,col="red")

plot(TN~DIN,dat.all2,ylab="TN (mg L\u207B\u00b9)",xlab="DIN (mg L\u207B\u00b9)",pch=21,bg=ifelse(DIN>TN,"dodgerblue1",NA),col=adjustcolor("grey",0.8));abline(0,1,col="dodgerblue1")
plot(TP~SRP,dat.all2,ylab="TP (mg L\u207B\u00b9)",xlab="SRP (mg L\u207B\u00b9)",pch=21,bg=ifelse(SRP>TP,"red",NA),col=adjustcolor("grey",0.8));abline(0,1,col="red")
dev.off()

## check and clean
sum(dat.all2$TNReversal,na.rm=T)
sum(dat.all2$TPReversal,na.rm=T)
dat.all2$TN=with(dat.all2,ifelse(TNReversal==1,NA,TN))
dat.all2$DIN=with(dat.all2,ifelse(TNReversal==1,NA,DIN))
dat.all2$TP=with(dat.all2,ifelse(TPReversal==1,NA,TP))
dat.all2$SRP=with(dat.all2,ifelse(TPReversal==1,NA,SRP))

dat.all2$TP.mM=with(dat.all2,TP/P.mw)
dat.all2$TN.mM=with(dat.all2,TN/N.mw)
dat.all2$SRP.mM=with(dat.all2,SRP/P.mw)
dat.all2$DIN.mM=with(dat.all2,DIN/N.mw)
dat.all2$TOC.mM=with(dat.all2,TOC/C.mw)

dat.all2$TN_TP=with(dat.all2,TN.mM/TP.mM)
dat.all2$TOC_TP=with(dat.all2,TOC.mM/TP.mM)
dat.all2$TOC_TN=with(dat.all2,TOC.mM/TN.mM)
dat.all2$DIN_SRP=with(dat.all2,DIN.mM/SRP.mM)
dat.all2$TOC_SRP=with(dat.all2,TOC.mM/SRP.mM)
dat.all2$TOC_DIN=with(dat.all2,TOC.mM/DIN.mM)


## from https://www.alxyon.com/en/library/plant-nutrition/c-n-p-ratio-or-redfield-ratio.html
NO3=seq(0.1,50,length.out=50)
PO4=seq(0.01,2.5,length.out=50)
red.explore=expand.grid(NO3=NO3,PO4=PO4)
red.explore$NO3.mM=red.explore$NO3/N.mw
red.explore$PO4.mM=red.explore$PO4/P.mw
red.explore$N_P=with(red.explore,NO3.mM/PO4.mM)

plot(NO3.mM~PO4.mM,red.explore,pch=21,bg=ifelse(red.explore$N_P<5,"blue",ifelse(red.explore$N_P>20,"green","grey")),col=NA)
abline(0,16,col="red")
abline(0,20,col="red",lty=2)
abline(0,5,col="red",lty=2)

head(dat.all2)

plot(Chla~TOC,dat.all2,log="y")

## Season Screen -----------------------------------------------------------
head(dat.all2)
idvars=c("STATION", "DATE", "WY", "season","ESTUARY", "SEGMENT_NA", "Region", "source","UTMX","UTMY")
dat.all2.melt=melt(dat.all2,id.vars=idvars)
dat.all2.melt=subset(dat.all2.melt,is.na(value)==F);# remove NAs carried over from xtab
unique(subset(dat.all2.melt,variable=="TOC"&Region=="ENP"&is.na(value)==F)$STATION)

samp.size=dcast(dat.all2.melt,STATION+WY+ESTUARY+Region+source+variable~season,value.var = "value",fun.aggregate = function(x)N.obs(x))
samp.size$TSamp=rowSums(samp.size[,c("A_Wet","B_Dry")],na.rm=T)
samp.size$sea.screen=with(samp.size, ifelse(A_Wet>0&B_Dry>0&TSamp>=4,1,0))

#double check Everglades TOC
subset(samp.size,variable=="TOC"&Region=="ENP")
subset(samp.size,variable=="Chla"&Region=="Shelf")
subset(samp.size,variable=="Chla"&Region=="Shelf"&source=="NOAA")

# number of consecutive years https://stackoverflow.com/a/23786245/5213091
consec.WY=ddply(samp.size,c("STATION","ESTUARY","Region","source","variable"),summarise,
                max.consec=max(rle(cumsum(c(1,diff(WY)>1)))$lengths),
                min.consec=min(rle(cumsum(c(1,diff(WY)>1)))$lengths))
consec.WY$ceonsec.diff=with(consec.WY,max.consec-min.consec)
consec.WY$consec.screen=with(consec.WY,ifelse(max.consec>=5,1,0))
consec.WY

vars_join1=c("STATION","WY","ESTUARY","Region","source","variable","sea.screen")
vars_by1=c("STATION","WY","ESTUARY","Region","source","variable")
vars_join2=c("STATION","ESTUARY","Region","source","variable","consec.screen")
vars_by2=c("STATION","ESTUARY","Region","source","variable")
dat.all2.melt=merge(dat.all2.melt,samp.size[,vars_join1],vars_by1)
dat.all2.melt=merge(dat.all2.melt,consec.WY[,vars_join2],vars_by2)

# change TP and SRP to ug/L
dat.all2.melt$value=with(dat.all2.melt,ifelse(variable%in%c("TP","SRP"),value*1000,value))

dat.all.GM=ddply(subset(dat.all2.melt,sea.screen==1&consec.screen==1),
                 c("STATION","WY","variable","source","ESTUARY","Region","UTMX","UTMY"),
                 summarise,
                 GM=exp(mean(log(value),na.rm=T)),N.val=N.obs(value))
head(dat.all.GM)
range(dat.all.GM$WY)
dat.all=subset(dat.all,WY%in%WYs)
# sanity spot checks
# subset(dat.all.GM,STATION==200&variable=="TN")
# subset(dat.all.GM,STATION==200&variable=="DIN")
# subset(dat.all.GM,STATION=="S12A"&variable=="DIN")
# subset(dat.all.GM,STATION=="SRS2"&variable=="TP")

tm_shape(serc.sites.shp)+tm_dots()+
  tm_shape(ENP_FLB)+tm_dots(col="red",alpha=0.5)+
  tm_shape(lter)+tm_dots(col="yellow",alpha=0.5)+
  tm_shape(noaa.sites.shp)+tm_dots(col="blue",alpha=0.5)


unique(dat.all.GM$variable)
dat.all.GM2=subset(dat.all.GM,
                   variable%in%c("TN", "DIN", "TP", "SRP", "Chla","TOC","TN_TP","DIN_SRP","SALINITY"))
dat.all.GM2.xtab=reshape2::dcast(subset(dat.all.GM2,WY%in%WYs),STATION+ESTUARY+Region+source+WY~variable,value.var = "GM",mean)
head(dat.all.GM2.xtab)
# write.csv(dat.all.GM2.xtab,paste0(export.path,"WYGeomean.csv"),row.names=F)

dat.all.GM=subset(dat.all.GM,WY%in%WYs)
# Data Analyses  ----------------------------------------------------------
library(dunn.test)
library(rcompanion)
library(flextable)

levels.var=c("ENP","Coastal_Mangroves","FLBay","Shelf","Keys")
levels.var.labs=c("ENP","Mangrove Fringe","Florida Bay","West Florida Shelf","Keys")
dat.all.GM$Region.f=factor(dat.all.GM$Region,levels=levels.var)
dat.all.GM=merge(dat.all.GM,data.frame(Region=levels.var,Region.txt=levels.var.labs),"Region",all.x=T)
dat.all.GM$Region.txt=factor(dat.all.GM$Region.txt,levels=levels.var.labs)
## Regional Comparison -----------------------------------------------------

region.sum=ddply(subset(dat.all.GM,variable%in%c("TN", "DIN", "TP", "SRP", "Chla","TOC")),
                 c("variable","Region.txt"),summarise,
                 N.val=N.obs(GM),
                 min.val=min(GM,na.rm=T),
                 med.val=median(GM,na.rm=T),
                 mean.val=mean(GM,na.rm=T),
                 max.val=max(GM,na.rm=T))

region.sum$Region.txt=with(region.sum,ifelse(Region.txt=="ENP","ENP (Marsh)",trimws(as.character(Region.txt))))

vars=c( "variable","Region.txt", "N.val", "min.val", "med.val", "mean.val","max.val")
region.sum[,vars]%>%
  flextable()%>%
  colformat_num(j=3,big.mark="")%>%
  colformat_double(j=4:7,i=1:10,digits=3,na_str = " ")%>%
  colformat_double(j=4:7,i=11:20,digits=2,na_str = " ")%>%
  colformat_double(j=4:7,i=21:30,digits=2,na_str = " ")%>%
  merge_v(j=1)%>%
  fix_border_issues()%>%
  valign(j=1,valign="top")%>%
  set_header_labels(
    "param.vals"="Parameter\n(Units)",
    "Region.txt"="Region",
    "N.val"="N",
    "min.val"="Minimum",
    "med.val"="Median",
    "mean.val"="Mean",
    "max.val"="Maximum"
  )%>%
  hline(i=c(5,10,15,20,25))%>%
  align(j=3:7,align="center",part="header")%>%
  padding(padding=1,part="all")%>%
  font(fontname="Times New Roman",part="all")%>%
  bold(part="header")%>%
  autofit()  #%>%print("docx")


region.sum=ddply(subset(dat.all.GM,variable%in%c("TN_TP", "TOC_TP", "TOC_TN")),
                 c("variable","Region.txt"),summarise,
                 N.val=N.obs(GM),
                 min.val=min(GM,na.rm=T),
                 med.val=median(GM,na.rm=T),
                 mean.val=mean(GM,na.rm=T),
                 max.val=max(GM,na.rm=T))

region.sum$Region.txt=with(region.sum,ifelse(Region.txt=="ENP","ENP (Marsh)",trimws(as.character(Region.txt))))

vars=c( "variable","Region.txt", "N.val", "min.val", "med.val", "mean.val","max.val")
region.sum[,vars]%>%
  flextable()%>%
  colformat_num(j=3,big.mark="")%>%
  colformat_double(j=4:7,digits=0,na_str = " ",big.mark="")%>%
  compose(j="Region.txt",i=~Region.txt=="West Florida Shelf",value=as_paragraph('W. Florida Shelf'))%>%
  compose(j="variable",i=~variable=="TN_TP",value=as_paragraph('TN:TP\n(molar ratio)'))%>%
  compose(j="variable",i=~variable=="TOC_TP",value=as_paragraph('TOC:TP\n(molar ratio)'))%>%
  compose(j="variable",i=~variable=="TOC_TN",value=as_paragraph('TOC:TN\n(molar ratio)'))%>%
  merge_v(j=1)%>%
  fix_border_issues()%>%
  valign(j=1,valign="top")%>%
  set_header_labels(
    "variable"="Parameter\n(Units)",
    "Region.txt"="Region",
    "N.val"="N",
    "min.val"="Minimum",
    "med.val"="Median",
    "mean.val"="Mean",
    "max.val"="Maximum"
  )%>%
  hline(i=c(5,10))%>%
  align(j=3:7,align="center",part="header")%>%
  padding(padding=1,part="all")%>%
  footnote(j="variable",i=~variable%in%c("TOC_TP","TOC_TN"),ref_symbols = " 1 ",part="body",
                                           value=as_paragraph("Limited availability of TOC measurements in ENP (Marsh) restricted sample size"))%>%
  font(fontname="Times New Roman",part="all")%>%
  bold(part="header")%>%
  autofit()#%>%print("docx")


## KW test
ddply(subset(dat.all.GM,variable%in%c("TN","DIN","TP","SRP","Chla","TOC","TN_TP","TOC_TP","TOC_TN","DIN_SRP","TOC_SRP","TOC_DIN")),"variable",summarise,
      chisq=as.numeric(kruskal.test(GM~Region.f)$statistic),
      df=as.numeric(kruskal.test(GM~Region.f)$parameter),
      pval=as.numeric(kruskal.test(GM~Region.f)$p.value))


cols=c("white",adjustcolor(wesanderson::wes_palette("Zissou1",4,"continuous"),0.5))
levels.var.labs=c("ENP\n(Marsh)","Mangrove\nFringe","Florida\nBay","W. Florida\nShelf","Keys\n")
# png(filename=paste0(plot.path,"revised/Fig2_RegionComp.png"),width=6.5,height=5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,3.5,0.5,0.75),oma=c(3,1,1,0.5));
layout(matrix(c(1:6),3,2,byrow=T))

ylim.val=c(0,2);by.y=0.5;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
x=boxplot(GM~Region.f,subset(dat.all.GM,variable=="TN"),outline=F,ylim=ylim.val,axes=F,ann=F,col=cols,boxwex=0.5)
TN.DT=with(subset(dat.all.GM,variable=="TN"),dunn.test(GM, Region.f))
TN.DT.ltr=cldList(P.adjusted ~ comparison,data=TN.DT,threshold = 0.05)
TN.DT.ltr$Letter=toupper(TN.DT.ltr$Letter)
TN.DT.ltr=TN.DT.ltr[order(match(TN.DT.ltr$Group,levels.var)),]
text(1:5,x$stats[5,],TN.DT.ltr$Letter,pos=3)
axis_fun(2,ymaj,ymin,format(ymaj))
axis_fun(1,1:5,1:5,NA);box(lwd=1)
mtext(side=2,line=2.5,"TN (mg N L\u207B\u00B9)")
mtext(side=3,adj=1,line=-1.25,"A ",font=2)

tmp=data.frame(com=TN.DT$comparisons,z.val=TN.DT$Z)
tmp$comp1=sapply(strsplit(tmp$com,"-"),"[",1)
tmp$comp2=sapply(strsplit(tmp$com,"-"),"[",2)
# tapply(subset(dat.all.GM,variable=="TN")$GM, INDEX = subset(dat.all.GM,variable=="TN")$Region.f, FUN = mean)
tapply(tmp$z.val, INDEX = tmp$comp1, FUN = median)

ylim.val=c(0,2.1);by.y=0.5;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
x=boxplot(GM~Region.f,subset(dat.all.GM,variable=="DIN"),outline=F,ylim=ylim.val,axes=F,ann=F,col=cols,boxwex=0.5)
DIN.DT=with(subset(dat.all.GM,variable=="DIN"),dunn.test(GM, Region.f))
DIN.DT.ltr=cldList(P.adjusted ~ comparison,data=DIN.DT,threshold = 0.05)
DIN.DT.ltr$Letter=toupper(DIN.DT.ltr$Letter)
DIN.DT.ltr=DIN.DT.ltr[order(match(DIN.DT.ltr$Group,levels.var)),]
text(1:5,x$stats[5,],DIN.DT.ltr$Letter,pos=3)
axis_fun(2,ymaj,ymin,format(ymaj))
axis_fun(1,1:5,1:5,NA);box(lwd=1)
mtext(side=2,line=2.5,"DIN (mg N L\u207B\u00B9)")
mtext(side=3,adj=1,line=-1.25,"B ",font=2)

ylim.val=c(0,60);by.y=20;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
x=boxplot(GM~Region.f,subset(dat.all.GM,variable=="TP"),outline=F,ylim=ylim.val,axes=F,ann=F,col=cols,boxwex=0.5)
TP.DT=with(subset(dat.all.GM,variable=="TP"),dunn.test(GM, Region.f))
TP.DT.ltr=cldList(P.adjusted ~ comparison,data=TP.DT,threshold = 0.05)
TP.DT.ltr$Letter=toupper(TP.DT.ltr$Letter)
TP.DT.ltr=TP.DT.ltr[order(match(TP.DT.ltr$Group,levels.var)),]
text(1:5,x$stats[5,],TP.DT.ltr$Letter,pos=3)
axis_fun(2,ymaj,ymin,format(ymaj))
axis_fun(1,1:5,1:5,NA);box(lwd=1)
mtext(side=2,line=2.5,"TP (\u03BCg P L\u207B\u00B9)")
mtext(side=3,adj=1,line=-1.25,"C ",font=2)

ylim.val=c(0,10);by.y=2;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
x=boxplot(GM~Region.f,subset(dat.all.GM,variable=="SRP"),outline=F,ylim=ylim.val,axes=F,ann=F,col=cols,boxwex=0.5)
SRP.DT=with(subset(dat.all.GM,variable=="SRP"),dunn.test(GM, Region.f))
SRP.DT.ltr=cldList(P.adjusted ~ comparison,data=SRP.DT,threshold = 0.05)
SRP.DT.ltr$Letter=toupper(SRP.DT.ltr$Letter)
SRP.DT.ltr=SRP.DT.ltr[order(match(SRP.DT.ltr$Group,levels.var)),]
text(1:5,x$stats[5,],SRP.DT.ltr$Letter,pos=3)
axis_fun(2,ymaj,ymin,format(ymaj))
# axis_fun(1,1:5,1:5,levels.var.labs,line=0.3,cex=0.9);box(lwd=1)
axis_fun(1,1:5,1:5,NA,line=0.3,cex=0.9);box(lwd=1)
mtext(side=2,line=2.5,"SRP (\u03BCg P L\u207B\u00B9)")
mtext(side=3,adj=1,line=-1.25,"D ",font=2)

ylim.val=c(0,7);by.y=2;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
x=boxplot(GM~Region.f,subset(dat.all.GM,variable=="Chla"),outline=F,ylim=ylim.val,axes=F,ann=F,col=cols,boxwex=0.5)
Chla.DT=with(subset(dat.all.GM,variable=="SRP"),dunn.test(GM, Region.f))
Chla.DT.ltr=cldList(P.adjusted ~ comparison,data=Chla.DT,threshold = 0.05)
Chla.DT.ltr$Letter=toupper(Chla.DT.ltr$Letter)
Chla.DT.ltr=Chla.DT.ltr[order(match(Chla.DT.ltr$Group,levels.var)),]
text(1:5,x$stats[5,],Chla.DT.ltr$Letter,pos=3)
axis_fun(2,ymaj,ymin,format(ymaj))
axis_fun(1,1:5,1:5,levels.var.labs,line=0.3,cex=0.8);box(lwd=1)
mtext(side=2,line=2.5,"Chl-a (\u03BCg L\u207B\u00B9)")
mtext(side=1,line=1,outer=T,"Region")
mtext(side=3,adj=1,line=-1.25,"E ",font=2)

ylim.val=c(0,30);by.y=10;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
x=boxplot(GM~Region.f,subset(dat.all.GM,variable=="TOC"),outline=F,ylim=ylim.val,axes=F,ann=F,col=cols,boxwex=0.5)
TOC.DT=with(subset(dat.all.GM,variable=="SRP"),dunn.test(GM, Region.f))
TOC.DT.ltr=cldList(P.adjusted ~ comparison,data=TOC.DT,threshold = 0.05)
TOC.DT.ltr$Letter=toupper(TOC.DT.ltr$Letter)
TOC.DT.ltr=TOC.DT.ltr[order(match(TOC.DT.ltr$Group,levels.var)),]
text(1:5,x$stats[5,],TOC.DT.ltr$Letter,pos=3)
axis_fun(2,ymaj,ymin,format(ymaj))
axis_fun(1,1:5,1:5,levels.var.labs,line=0.3,cex=0.8);box(lwd=1)
mtext(side=2,line=2.5,"TOC (mg C L\u207B\u00B9)")
mtext(side=1,line=1,outer=T,"Region")
mtext(side=3,adj=1,line=-1.25,"F ",font=2)
dev.off()


# png(filename=paste0(plot.path,"revised/Figx_RegionComp_stoich.png"),width=6.5,height=5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,6,0.5,0.75),oma=c(3,1,1,0.5));
layout(matrix(c(1:6),3,2,byrow=T))

# ylim.val=c(0,900);by.y=200;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
ylim.val=c(10,1500);by.y=200;ymaj=log.scale.fun(ylim.val,"major");ymin=log.scale.fun(ylim.val,"minor")
x=boxplot(GM~Region.f,subset(dat.all.GM,variable=="TN_TP"),outline=F,
          ylim=ylim.val,axes=F,ann=F,col=cols,boxwex=0.5,log="y");
abline(h=16,lty=2,col="indianred1")
NP.DT=with(subset(dat.all.GM,variable=="TN_TP"),dunn.test(GM, Region.f))
NP.DT.ltr=cldList(P.adjusted ~ comparison,data=NP.DT,threshold = 0.05)
NP.DT.ltr$Letter=toupper(NP.DT.ltr$Letter)
NP.DT.ltr=NP.DT.ltr[order(match(NP.DT.ltr$Group,levels.var)),]
text(1:5,x$stats[5,],NP.DT.ltr$Letter,pos=3)
axis_fun(2,ymaj,ymin,format(ymaj))
axis_fun(1,1:5,1:5,NA);box(lwd=1)
mtext(side=2,line=3,"TN:TP\n(molar ratio)")
mtext(side=3,adj=1,line=-1.25,"A ",font=2)

# ylim.val=c(0,3000);by.y=1000;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
ylim.val=c(1,8000);ymaj=log.scale.fun(ylim.val,"major");ymin=log.scale.fun(ylim.val,"minor")
x=boxplot(GM~Region.f,subset(dat.all.GM,variable=="DIN_SRP"),outline=F,
          ylim=ylim.val,axes=F,ann=F,col=cols,boxwex=0.5,log="y");
abline(h=16,lty=2,col="indianred1")
NP.DT=with(subset(dat.all.GM,variable=="DIN_SRP"),dunn.test(GM, Region.f))
NP.DT.ltr=cldList(P.adjusted ~ comparison,data=NP.DT,threshold = 0.05)
NP.DT.ltr$Letter=toupper(NP.DT.ltr$Letter)
NP.DT.ltr=NP.DT.ltr[order(match(NP.DT.ltr$Group,levels.var)),]
text(1:5,x$stats[5,],NP.DT.ltr$Letter,pos=3)
axis_fun(2,ymaj,ymin,format(ymaj,scientific = F))
axis_fun(1,1:5,1:5,NA);box(lwd=1)
mtext(side=2,line=3,"DIN:SRP\n(molar ratio)")
mtext(side=3,adj=1,line=-1.25,"B ",font=2)

ylim.val=c(100,20000);ymaj=log.scale.fun(ylim.val,"major");ymin=log.scale.fun(ylim.val,"minor")
x=boxplot(GM~Region.f,subset(dat.all.GM,variable=="TOC_TP"),outline=F,
          ylim=ylim.val,axes=F,ann=F,col=cols,boxwex=0.5,log="y");
abline(h=106,lty=2,col="indianred1")
NP.DT=with(subset(dat.all.GM,variable=="TOC_TP"),dunn.test(GM, Region.f))
NP.DT.ltr=cldList(P.adjusted ~ comparison,data=NP.DT,threshold = 0.05)
NP.DT.ltr$Letter=toupper(NP.DT.ltr$Letter)
NP.DT.ltr=NP.DT.ltr[order(match(NP.DT.ltr$Group,levels.var)),]
text(1:5,x$stats[5,],NP.DT.ltr$Letter,pos=3)
axis_fun(2,ymaj,ymin,format(ymaj/1000))
axis_fun(1,1:5,1:5,NA);box(lwd=1)
mtext(side=2,line=3,"OC:TP\n(x10\u00B3 molar ratio)")
mtext(side=3,adj=1,line=-1.25,"C ",font=2)

ylim.val=c(100,100000);ymaj=log.scale.fun(ylim.val,"major");ymin=log.scale.fun(ylim.val,"minor")
x=boxplot(GM~Region.f,subset(dat.all.GM,variable=="TOC_SRP"),outline=F,
          ylim=ylim.val,axes=F,ann=F,col=cols,boxwex=0.5,log="y");
abline(h=106,lty=2,col="indianred1")
NP.DT=with(subset(dat.all.GM,variable=="TOC_SRP"),dunn.test(GM, Region.f))
NP.DT.ltr=cldList(P.adjusted ~ comparison,data=NP.DT,threshold = 0.05)
NP.DT.ltr$Letter=toupper(NP.DT.ltr$Letter)
NP.DT.ltr=NP.DT.ltr[order(match(NP.DT.ltr$Group,levels.var)),]
text(1:5,x$stats[5,],NP.DT.ltr$Letter,pos=3)
axis_fun(2,ymaj,ymin,format(ymaj/1000))
axis_fun(1,1:5,1:5,NA);box(lwd=1)
mtext(side=2,line=3,"OC:SRP\n(x10\u00B3 molar ratio)")
mtext(side=3,adj=1,line=-1.25,"D ",font=2)

ylim.val=c(1,100);ymaj=log.scale.fun(ylim.val,"major");ymin=log.scale.fun(ylim.val,"minor")
x=boxplot(GM~Region.f,subset(dat.all.GM,variable=="TOC_TN"),outline=F,
          ylim=ylim.val,axes=F,ann=F,col=cols,boxwex=0.5,log="y");
abline(h=6.6,lty=2,col="indianred1")
NP.DT=with(subset(dat.all.GM,variable=="TOC_TN"),dunn.test(GM, Region.f))
NP.DT.ltr=cldList(P.adjusted ~ comparison,data=NP.DT,threshold = 0.05)
NP.DT.ltr$Letter=toupper(NP.DT.ltr$Letter)
NP.DT.ltr=NP.DT.ltr[order(match(NP.DT.ltr$Group,levels.var)),]
text(1:5,x$stats[5,],NP.DT.ltr$Letter,pos=3)
axis_fun(2,ymaj,ymin,format(ymaj))
axis_fun(1,1:5,1:5,levels.var.labs,line=0.3,cex=0.8);box(lwd=1)
mtext(side=2,line=3,"OC:TN\n(molar ratio)")
mtext(side=3,adj=1,line=-1.25,"E ",font=2)

ylim.val=c(1,2000);ymaj=log.scale.fun(ylim.val,"major");ymin=log.scale.fun(ylim.val,"minor")
x=boxplot(GM~Region.f,subset(dat.all.GM,variable=="TOC_DIN"),outline=F,
          ylim=ylim.val,axes=F,ann=F,col=cols,boxwex=0.5,log="y");
abline(h=6.6,lty=2,col="indianred1")
NP.DT=with(subset(dat.all.GM,variable=="TOC_DIN"),dunn.test(GM, Region.f))
NP.DT.ltr=cldList(P.adjusted ~ comparison,data=NP.DT,threshold = 0.05)
NP.DT.ltr$Letter=toupper(NP.DT.ltr$Letter)
NP.DT.ltr=NP.DT.ltr[order(match(NP.DT.ltr$Group,levels.var)),]
text(1:5,x$stats[5,],NP.DT.ltr$Letter,pos=3)
axis_fun(2,ymaj,ymin,format(ymaj))
axis_fun(1,1:5,1:5,levels.var.labs,line=0.3,cex=0.8);box(lwd=1)
mtext(side=2,line=3,"OC:DIN\n(molar ratio)")
mtext(side=3,adj=1,line=-1.25,"F ",font=2)


mtext(side=1,line=1,outer=T,"Region")


dev.off()

#C:N:P 106:16:1
# png(filename=paste0(plot.path,"revised/Figx_RegionComp_stoich2.png"),width=3.5,height=5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,3,0.5,0.5),oma=c(3,3,0.5,0.5));
layout(matrix(c(1:3),3,1,byrow=T))

# ylim.val=c(0,900);by.y=200;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
ylim.val=c(10,1500);by.y=200;ymaj=log.scale.fun(ylim.val,"major");ymin=log.scale.fun(ylim.val,"minor")
x=boxplot(GM~Region.f,subset(dat.all.GM,variable=="TN_TP"),outline=F,
          ylim=ylim.val,axes=F,ann=F,col=cols,boxwex=0.5,log="y");
abline(h=16/1,lty=2,col="indianred1")
NP.DT=with(subset(dat.all.GM,variable=="TN_TP"),dunn.test(GM, Region.f))
NP.DT.ltr=cldList(P.adjusted ~ comparison,data=NP.DT,threshold = 0.05)
NP.DT.ltr$Letter=toupper(NP.DT.ltr$Letter)
NP.DT.ltr=NP.DT.ltr[order(match(NP.DT.ltr$Group,levels.var)),]
text(1:5,x$stats[5,],NP.DT.ltr$Letter,pos=3)
axis_fun(2,ymaj,ymin,format(ymaj))
axis_fun(1,1:5,1:5,NA);box(lwd=1)
mtext(side=2,line=3,"TN:TP\n(molar ratio)")
mtext(side=3,adj=1,line=-1.25,"A ",font=2)

ylim.val=c(100,20000);ymaj=log.scale.fun(ylim.val,"major");ymin=log.scale.fun(ylim.val,"minor")
x=boxplot(GM~Region.f,subset(dat.all.GM,variable=="TOC_TP"),outline=F,
          ylim=ylim.val,axes=F,ann=F,col=cols,boxwex=0.5,log="y");
abline(h=106/1,lty=2,col="indianred1")
NP.DT=with(subset(dat.all.GM,variable=="TOC_TP"),dunn.test(GM, Region.f))
NP.DT.ltr=cldList(P.adjusted ~ comparison,data=NP.DT,threshold = 0.05)
NP.DT.ltr$Letter=toupper(NP.DT.ltr$Letter)
NP.DT.ltr=NP.DT.ltr[order(match(NP.DT.ltr$Group,levels.var)),]
text(1:5,x$stats[5,],NP.DT.ltr$Letter,pos=3)
axis_fun(2,ymaj,ymin,format(ymaj/1000))
axis_fun(1,1:5,1:5,NA);box(lwd=1)
mtext(side=2,line=3,"TOC:TP\n(x10\u00B3 molar ratio)")
mtext(side=3,adj=1,line=-1.25,"B ",font=2)

ylim.val=c(1,100);ymaj=log.scale.fun(ylim.val,"major");ymin=log.scale.fun(ylim.val,"minor")
x=boxplot(GM~Region.f,subset(dat.all.GM,variable=="TOC_TN"),outline=F,
          ylim=ylim.val,axes=F,ann=F,col=cols,boxwex=0.5,log="y");
abline(h=106/16,lty=2,col="indianred1")
NP.DT=with(subset(dat.all.GM,variable=="TOC_TN"),dunn.test(GM, Region.f))
NP.DT.ltr=cldList(P.adjusted ~ comparison,data=NP.DT,threshold = 0.05)
NP.DT.ltr$Letter=toupper(NP.DT.ltr$Letter)
NP.DT.ltr=NP.DT.ltr[order(match(NP.DT.ltr$Group,levels.var)),]
text(1:5,x$stats[5,],NP.DT.ltr$Letter,pos=3)
axis_fun(2,ymaj,ymin,format(ymaj))
axis_fun(1,1:5,1:5,levels.var.labs,line=0.3,cex=0.8);box(lwd=1)
mtext(side=2,line=3,"TOC:TN\n(molar ratio)")
mtext(side=3,adj=1,line=-1.25,"C ",font=2)
mtext(side=1,line=3,"Region")
dev.off()

# png(filename=paste0(plot.path,"revised/Figx_RegionComp_NP.png"),width=5,height=3,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,3.5,0.25,0.75),oma=c(3,1,0.5,0.5));
ylim.val=c(0,900);by.y=200;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)

x=boxplot(GM~Region.f,subset(dat.all.GM,variable=="TN_TP"),outline=F,
          ylim=ylim.val,axes=F,ann=F,col=cols,boxwex=0.5);
abline(h=16,lty=2,col="indianred1")
NP.DT=with(subset(dat.all.GM,variable=="TN_TP"),dunn.test(GM, Region.f))
NP.DT.ltr=cldList(P.adjusted ~ comparison,data=NP.DT,threshold = 0.05)
NP.DT.ltr$Letter=toupper(NP.DT.ltr$Letter)
NP.DT.ltr=NP.DT.ltr[order(match(NP.DT.ltr$Group,levels.var)),]
text(1:5,x$stats[5,],NP.DT.ltr$Letter,pos=3)
axis_fun(2,ymaj,ymin,format(ymaj))
axis_fun(1,1:5,1:5,levels.var.labs,line=0.3,cex=0.8);box(lwd=1)
mtext(side=2,line=2.5,"TN:TP (molar ratio)")
mtext(side=1,line=1,outer=T,"Region")
dev.off()



# Chla-TOC  ---------------------------------------------------------------
# analysis based on Isles et al (2021) Trade-offs Between Light and Nutrient 
## Availability Across Gradients of Dissolved Organic Carbon Lead to Spatially
## and Temporally Variable Responses of Lake Phytoplankton Biomass to Browning.
## Ecosystems 24:1837–1852. doi: 10.1007/s10021-021-00619-7

### Future analysis

library(ggplot2)

ggplot(dat.all.GM2.xtab,aes(y=Chla,x=TOC))+
  geom_point(fill="indianred1",alpha=0.5,shape=21)+
  scale_y_continuous(trans="log")+
  facet_wrap(~Region)+
  theme_bw()+theme(text=element_text(family="serif"))

fitG =
  function(x,y,mu,sig,scale){
    
    f = function(p){
      d = p[3]*dnorm(x,mean=p[1],sd=p[2])
      sum((d-y)^2,na.rm=T)
    }
    
    optim(c(mu,sig,scale),f)
  }
p=fitG(x=dat.all.GM2.xtab$TOC,y=dat.all.GM2.xtab$Chla,
       mu=mean(dat.all.GM2.xtab$TOC,na.rm=T),sig=sd(dat.all.GM2.xtab$TOC,na.rm=T),scale=0.001)

plot(Chla~TOC,dat.all.GM2.xtab,log="y")
lines(dat.all.GM2.xtab$TOC,p$par[3]*dnorm(dat.all.GM2.xtab$TOC,p$par[1],p$par[2]),col="red")

###
dat.scrn.xtab=dcast(subset(dat.all2.melt,sea.screen==1&consec.screen==1),STATION+ESTUARY+Region+source+DATE+WY~variable,value.var = "value",mean)
dat.scrn.xtab$month=as.numeric(format(dat.scrn.xtab$DATE,"%m"))
dat.scrn.xtab=merge(dat.scrn.xtab,data.frame(month=1:12,season=c(rep("EarlyDry",2),rep("LateDry",3),rep("Wet",5),rep("EarlyDry",2))),"month",all.x=T)


plot(Chla~TOC,subset(dat.scrn.xtab,season=="EarlyDry"),log="y")
points(Chla~TOC,subset(dat.scrn.xtab,season=="LateDry"),pch=21,bg="green")
points(Chla~TOC,subset(dat.scrn.xtab,season=="Wet"),pch=21,bg="blue")

ggplot(dat.scrn.xtab,aes(y=Chla,x=TOC,fill=season))+
  geom_point(alpha=0.5,shape=21)+
  scale_y_continuous(trans="log")+
  facet_wrap(~Region)+
  theme_bw()+theme(text=element_text(family="serif"))

ggplot(dat.scrn.xtab,aes(y=Chla,x=TOC,fill=Region))+
  geom_point(alpha=0.25,shape=21)+
  scale_y_continuous(trans="log")+
  facet_wrap(~season)+
  theme_bw()+theme(text=element_text(family="serif"))

library(quantreg)
fit1.nlrq <- nlrq(Chla ~ a+b*TOC, data=dat.scrn.xtab, start =list(a=1,b=5), tau=0.5, trace=T)
summary(fit1.nlrq)
fit2.nlrq <- nlrq(Chla ~ a+(b*TOC)+(c*TOC**2), data=dat.scrn.xtab, start =list(a=1,b=5,c=2.5), tau=0.5, trace=T)
summary(fit2.nlrq)
fit3.nlrq <- nlrq(Chla ~ a*exp(-((TOC-b)/c)**2), data=dat.scrn.xtab, start =list(a=2,b=30,c=20), tau=0.5, trace=T)
summary(fit3.nlrq)



x.val=seq(0.01,100,length.out=100)
pred.fit1.nlrq=predict(fit1.nlrq,data.frame(TOC=x.val))
pred.fit2.nlrq=predict(fit2.nlrq,data.frame(TOC=x.val))
pred.fit3.nlrq=predict(fit3.nlrq,data.frame(TOC=x.val))

plot(Chla~TOC,dat.scrn.xtab,log="y",xlim=c(0,50),col=adjustcolor("grey",0.5))
lines(x.val,pred.fit1.nlrq,col="red",lwd=2)
lines(x.val,pred.fit2.nlrq,col="green",lwd=2)
lines(x.val,pred.fit3.nlrq,col="dodgerblue",lwd=2)






# Trend Analysis ----------------------------------------------------------
Ncheck=ddply(dat.all.GM,c("STATION","variable"),summarise,N.val=N.obs(WY))
subset(Ncheck,N.val<3)

ann.trend=ddply(subset(dat.all.GM,WY%in%seq(1996,2019,1)),c("STATION","variable"),summarise,
                est=as.numeric(cor.test(WY,GM,method="kendall")$estimate),
                pval=cor.test(WY,GM,method="kendall")$p.value,
                sen.slope=as.numeric(zyp::zyp.sen(GM~WY)$coefficients[2]),
                N.WY=N.obs(WY))
subset(ann.trend,pval<0.05)
ann.trend$stat.sig=with(ann.trend,ifelse(pval<0.05,"sig","not-sig"))
ann.trend$stat.sig=with(ann.trend,ifelse(is.na(stat.sig)==T,"Empty",as.character(stat.sig)))
ann.trend$stat.sig=as.factor(ann.trend$stat.sig)



AGM.all=ddply(subset(dat.all.GM,WY%in%seq(1996,2019,1)),c("STATION","variable"),summarise,
              mean.GM=mean(GM,na.rm=T),
              N.val=N.obs(GM),
              SE.val=SE(GM),
              var.val=var(GM,na.rm=T),# sample variance s^2; use this one (variance within site)
              pop.var=var.val*((N.val-1)/N.val), # population variance sigma^2 
              sd.val=sd(GM,na.rm=T),
              CV.val=cv.per(GM)*100,
              min.GM=min(GM,na.rm=T))

vars=c("TN","DIN","TP","SRP","Chla","TOC")
col.vars=c("STATION", "variable", "est", "pval", "sen.slope", "N.WY")
# write.csv(subset(ann.trend,variable%in%vars)[,col.vars],paste0(export.path,"TrendRslt_TableS2.csv"),row.names=F)

## Spatial -----------------------------------------------------------------
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

noaa.sites.shp2=cbind(data.frame(STATION=noaa.sites.shp$STATION),coordinates(noaa.sites.shp))
colnames(noaa.sites.shp2)=c("STATION","UTMX","UTMY")

sites.shp=rbind(serc.shp,lter.shp,wmd.shp,noaa.sites.shp2)
rm(serc.shp,lter.shp,wmd.shp,noaa.sites.shp2)
sites.shp.TableS1=merge(sites.shp,all.sites.region[,c("STATION","Region","source")],'STATION',all.x=T)

# write.csv(sites.shp.TableS1,paste0(export.path,"Sites_TableS1.csv"),row.names = F)

length(unique(sites.shp.TableS1$STATION))

### Quick Data inventory ----------------------------------------------------
dat.inv=dat.all.GM# merge(dat.all.GM,sites.shp.TableS1[,c("STATION","Region")],"STATION",all.x=T)
unique(dat.inv$Region)
unique(subset(dat.inv,is.na(Region)==T)$STATION)

dat.inv2=ddply(subset(dat.inv,is.na(Region)==F),
               c("Region","WY","variable"),summarise,N.val=N.obs(GM))
fill=data.frame(expand.grid(Region=unique(dat.inv2$Region),
                            WY=seq(1996,2019,1),
                            variable=unique(dat.inv2$variable)))
dat.inv2=merge(dat.inv2,fill,c("Region","WY","variable"),all.y=T)
dat.inv2$N.val[is.na(dat.inv2$N.val)==T]=0

dat.inv2=merge(dat.inv2,data.frame(Region=c("Keys", "Shelf", "FLBay", "ENP", "Coastal_Mangroves"),
                                   Region.plot=c(5,4,3,1,2)),"Region")
subset(dat.inv2,Region=="Shelf"&variable=="TP")

bks=c(0,1,10,25,50,100,200)
dat.inv2$count.cat=as.factor(findInterval(dat.inv2$N.val,
                                          bks))#,rightmost.closed = T,left.open = T))
cols.vir=rev(viridis::inferno(6))
cols.vir=rev(viridis::plasma(6))
cols.vir[1]="white"
dat.inv2$cols=cols.vir[dat.inv2$count.cat]

ylim.val=c(1,5);by.y=1;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/by.y)
ylim.val2=ylim.val
ylim.val=c(ylim.val[1]-0.5,ylim.val[2]+0.5)
xlim.val=c(1992,2019);by.x=5;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/by.x)
xlim.val=c(xlim.val[1]-0.5,xlim.val[2]+0.5)

par(family="serif",mar=c(2,1,1,0.75),oma=c(1,5,0.5,1));
layout(matrix(c(1:2),1,2,byrow=T),widths=c(1,0.5))
  
plot(Region.plot~WY,subset(dat.inv2,variable=='TP'),xlim=xlim.val,ylim=rev(ylim.val),type="n",xaxs="i",yaxs="i",axes=F,ann=F)
for(i in seq(ylim.val2[1],ylim.val2[2],1)){
  tmp=subset(dat.inv2,variable=='DIN'&Region.plot==i)  
  rect(seq(xlim.val[1],xlim.val[2],1),i+0.5,2020.5,i-0.5,col=tmp$cols,border=adjustcolor("grey",0.3),lwd=0.1)
}
axis_fun(1,xmaj,xmin,xmaj,line=-0.75,cex=0.80)
axis_fun(2,rev(ymaj),ymin,rev(c("ENP", "Coastal Mangroves", "FLBay", 
                                "Shelf", "Keys")),cex=0.75)
box(lwd=1)
mtext(side=2,line=4,"Region")
mtext(side=1,line=1.75,"WY",cex=0.95)

plot(0:1,0:1,axes=F,ann=F,type="n")
legend("center",legend=c("No Samples","1 - 10", "10 - 25","25 - 50","50 - 100","100 - 200"),
       lty=0,col="black",pch=22,
       pt.bg=cols.vir,lwd=0.1,
       pt.cex=1.5,ncol=1,cex=0.9,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0,yjust=1,
       title.adj = 0,title="Number of Sites\nDIN")


#####
# write.csv(dat.all.GM,paste0(export.path,"20230103_GM_allparam.csv"),row.names = F)
#####  
sites.shp2=SpatialPointsDataFrame(sites.shp[,c("UTMX","UTMY")],data=sites.shp,proj4string = utm17)

sites.shp.TN.trend=merge(sites.shp2,subset(ann.trend,variable=="TN"),"STATION",all.y=T)
sites.shp.TP.trend=merge(sites.shp2,subset(ann.trend,variable=="TP"),"STATION",all.y=T)
sites.shp.DIN.trend=merge(sites.shp2,subset(ann.trend,variable=="DIN"),"STATION",all.y=T)
sites.shp.SRP.trend=merge(sites.shp2,subset(ann.trend,variable=="SRP"),"STATION",all.y=T)
sites.shp.TOC.trend=merge(sites.shp2,subset(ann.trend,variable=="TOC"),"STATION",all.y=T)
sites.shp.Chla.trend=merge(sites.shp2,subset(ann.trend,variable=="Chla"),"STATION",all.y=T)
sites.shp.NP.trend=merge(sites.shp2,subset(ann.trend,variable=="TN_TP"),"STATION",all.y=T)

sites.shp.TN.GM=merge(sites.shp2,subset(AGM.all,variable=="TN"),"STATION",all.y=T)
sites.shp.TP.GM=merge(sites.shp2,subset(AGM.all,variable=="TP"),"STATION",all.y=T)
sites.shp.DIN.GM=merge(sites.shp2,subset(AGM.all,variable=="DIN"),"STATION",all.y=T)
sites.shp.SRP.GM=merge(sites.shp2,subset(AGM.all,variable=="SRP"),"STATION",all.y=T)
sites.shp.TOC.GM=merge(sites.shp2,subset(AGM.all,variable=="TOC"),"STATION",all.y=T)
sites.shp.Chla.GM=merge(sites.shp2,subset(AGM.all,variable=="Chla"),"STATION",all.y=T)
sites.shp.NP.GM=merge(sites.shp2,subset(AGM.all,variable=="TN_TP"),"STATION",all.y=T)

## Spline ------------------------------------------------------------------
#thin plate spline https://rspatial.org/raster/analysis/4-interpolation.html

region.buf.r=raster(region.mask)
res(region.buf.r)=1000

### Spatial Trend -----------------------------------------------------------
# TN
m=Tps(coordinates(sites.shp.TN.trend),sites.shp.TN.trend$sen.slope)
tps=interpolate(region.buf.r,m)
tps.TN.trend=mask(tps,region.mask)
# plot(tps.TN.trend)

# DIN
m=Tps(coordinates(sites.shp.DIN.trend),sites.shp.DIN.trend$sen.slope)
tps=interpolate(region.buf.r,m)
tps.DIN.trend=mask(tps,region.mask)
# plot(tps.DIN.trend)

# TP
m=Tps(coordinates(sites.shp.TP.trend),sites.shp.TP.trend$sen.slope)
tps=interpolate(region.buf.r,m)
tps.TP.trend=mask(tps,region.mask)
# plot(tps.TP.trend)

# SRP
m=Tps(coordinates(sites.shp.SRP.trend),sites.shp.SRP.trend$sen.slope)
tps=interpolate(region.buf.r,m)
tps.SRP.trend=mask(tps,region.mask)
# plot(tps.SRP.trend)

# Chla
m=Tps(coordinates(sites.shp.Chla.trend),sites.shp.Chla.trend$sen.slope)
tps=interpolate(region.buf.r,m)
tps.Chla.trend=mask(tps,region.mask)
# plot(tps.Chla.trend)

# plot(tps.Chla.trend)
# plot(sites.shp.Chla.trend,pch=ifelse(sites.shp.Chla.trend$sen.slope<0,25,21),
#      bg=ifelse(sites.shp.Chla.trend$sen.slope<0,"red","black"),
#      add=T)
# plot(GM~WY,subset(dat.all.GM,variable=="Chla"&STATION=="FLAB13"))

# TOC
m=Tps(coordinates(sites.shp.TOC.trend),sites.shp.TOC.trend$sen.slope)
tps=interpolate(region.buf.r,m)
tps.TOC.trend=mask(tps,region.mask)
# plot(tps.TOC.trend)

# NP
m=Tps(coordinates(sites.shp.NP.trend),sites.shp.NP.trend$sen.slope)
tps=interpolate(region.buf.r,m)
tps.NP.trend=mask(tps,region.mask)
# plot(tps.NP.trend)

### Spatial AGM -------------------------------------------------------------
# TN
sites.shp.TN.GM=subset(sites.shp.TN.GM,is.na(mean.GM)==F)
m.GM=Tps(coordinates(sites.shp.TN.GM),sites.shp.TN.GM$mean.GM)
tps.GM=interpolate(region.buf.r,m.GM)
tps.TN.GM=mask(tps.GM,region.mask)
# plot(tps.TN.GM)

m.GM.var=Tps(coordinates(sites.shp.TN.GM),sites.shp.TN.GM$var.val)
tps.GM.var=interpolate(region.buf.r,m.GM.var)
tps.TN.GM.var=mask(tps.GM.var,region.mask)

m.GM.CV=Tps(coordinates(sites.shp.TN.GM),sites.shp.TN.GM$CV.val)
tps.GM.CV=interpolate(region.buf.r,m.GM.CV)
tps.TN.GM.CV=mask(tps.GM.CV,region.mask)

# DIN
sites.shp.DIN.GM=subset(sites.shp.DIN.GM,is.na(mean.GM)==F)
m.GM=Tps(coordinates(sites.shp.DIN.GM),sites.shp.DIN.GM$mean.GM)
tps.GM=interpolate(region.buf.r,m.GM)
tps.DIN.GM=mask(tps.GM,region.mask)
# plot(tps.DIN.GM)

m.GM.var=Tps(coordinates(sites.shp.DIN.GM),sites.shp.DIN.GM$var.val)
tps.GM.var=interpolate(region.buf.r,m.GM.var)
tps.DIN.GM.var=mask(tps.GM.var,region.mask)

m.GM.CV=Tps(coordinates(sites.shp.DIN.GM),sites.shp.DIN.GM$CV.val)
tps.GM.CV=interpolate(region.buf.r,m.GM.CV)
tps.DIN.GM.CV=mask(tps.GM.CV,region.mask)

# TP
sites.shp.TP.GM=subset(sites.shp.TP.GM,is.na(mean.GM)==F)
m.GM=Tps(coordinates(sites.shp.TP.GM),sites.shp.TP.GM$mean.GM)
tps.GM=interpolate(region.buf.r,m.GM)
tps.TP.GM=mask(tps.GM,region.mask)
# plot(tps.TP.GM)

m.GM.var=Tps(coordinates(sites.shp.TP.GM),sites.shp.TP.GM$var.val)
tps.GM.var=interpolate(region.buf.r,m.GM.var)
tps.TP.GM.var=mask(tps.GM.var,region.mask)
# plot(tps.TP.GM.var)

m.GM.CV=Tps(coordinates(sites.shp.TP.GM),sites.shp.TP.GM$CV.val)
tps.GM.CV=interpolate(region.buf.r,m.GM.CV)
tps.TP.GM.CV=mask(tps.GM.CV,region.mask)
# plot(tps.TP.GM.CV)

# SRP
sites.shp.SRP.GM=subset(sites.shp.SRP.GM,is.na(mean.GM)==F)
m.GM=Tps(coordinates(sites.shp.SRP.GM),sites.shp.SRP.GM$mean.GM)
tps.GM=interpolate(region.buf.r,m.GM)
tps.SRP.GM=mask(tps.GM,region.mask)
# plot(tps.SRP.GM)

m.GM.var=Tps(coordinates(sites.shp.SRP.GM),sites.shp.SRP.GM$var.val)
tps.GM.var=interpolate(region.buf.r,m.GM.var)
tps.SRP.GM.var=mask(tps.GM.var,region.mask)
# plot(tps.SRP.GM.var)

m.GM.CV=Tps(coordinates(sites.shp.SRP.GM),sites.shp.SRP.GM$CV.val)
tps.GM.CV=interpolate(region.buf.r,m.GM.CV)
tps.SRP.GM.CV=mask(tps.GM.CV,region.mask)
# plot(tps.SRP.GM.CV)

# Chla
sites.shp.Chla.GM=subset(sites.shp.Chla.GM,is.na(mean.GM)==F)
m.GM=Tps(coordinates(sites.shp.Chla.GM),sites.shp.Chla.GM$mean.GM)
tps.GM=interpolate(region.buf.r,m.GM)
tps.Chla.GM=mask(tps.GM,region.mask)
# plot(tps.Chla.GM)

m.GM.var=Tps(coordinates(sites.shp.Chla.GM),sites.shp.Chla.GM$var.val)
tps.GM.var=interpolate(region.buf.r,m.GM.var)
tps.Chla.GM.var=mask(tps.GM.var,region.mask)
# plot(tps.Chla.GM.var)

m.GM.CV=Tps(coordinates(sites.shp.Chla.GM),sites.shp.Chla.GM$CV.val)
tps.GM.CV=interpolate(region.buf.r,m.GM.CV)
tps.Chla.GM.CV=mask(tps.GM.CV,region.mask)
# plot(tps.Chla.GM.CV)

# TOC
sites.shp.TOC.GM=subset(sites.shp.TOC.GM,is.na(mean.GM)==F)
m.GM=Tps(coordinates(sites.shp.TOC.GM),sites.shp.TOC.GM$mean.GM)
tps.GM=interpolate(region.buf.r,m.GM)
tps.TOC.GM=mask(tps.GM,region.mask)
# plot(tps.TOC.GM)

m.GM.var=Tps(coordinates(sites.shp.TOC.GM),sites.shp.TOC.GM$var.val)
tps.GM.var=interpolate(region.buf.r,m.GM.var)
tps.TOC.GM.var=mask(tps.GM.var,region.mask)
# plot(tps.TOC.GM.var)

m.GM.CV=Tps(coordinates(sites.shp.TOC.GM),sites.shp.TOC.GM$CV.val)
tps.GM.CV=interpolate(region.buf.r,m.GM.CV)
tps.TOC.GM.CV=mask(tps.GM.CV,region.mask)
# plot(tps.TOC.GM.CV)

# NP
sites.shp.NP.GM=subset(sites.shp.NP.GM,is.na(mean.GM)==F)
m.GM=Tps(coordinates(sites.shp.NP.GM),sites.shp.NP.GM$mean.GM)
tps.GM=interpolate(region.buf.r,m.GM)
tps.NP.GM=mask(tps.GM,region.mask)
# plot(tps.NP.GM)

m.GM.var=Tps(coordinates(sites.shp.NP.GM),sites.shp.NP.GM$var.val)
tps.GM.var=interpolate(region.buf.r,m.GM.var)
tps.NP.GM.var=mask(tps.GM.var,region.mask)
# plot(tps.NP.GM.var)

m.GM.CV=Tps(coordinates(sites.shp.NP.GM),sites.shp.NP.GM$CV.val)
tps.GM.CV=interpolate(region.buf.r,m.GM.CV)
tps.NP.GM.CV=mask(tps.GM.CV,region.mask)
# plot(tps.NP.GM.CV)



### Trend maps --------------------------------------------------------------
cols.val=c("grey","white","red")
# png(filename=paste0(plot.path,"revised/TrendMaps_TP_TN_TNTP.png"),width=6.5,height=4.5,units="in",res=200,type="windows",bg="white")
par(family="serif",oma=c(0.25,0.25,0.25,0.25),mar=c(0.1,0.1,0.1,0.1),xpd=F)
layout(matrix(c(1:12),3,4,byrow=T),widths = c(1,0.5,1,0.5))
bbox.lims=bbox(region.mask)
leg.pos.x=c(0.25,0.35)
leg.pos.y=c(0.4,0.8)
{
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
  mapmisc::scaleBar(utm17,"bottomright",bty="n",cex=1,seg.len=4)
  mtext(side=3,adj=0,line=-1.25," A",font=2)
  
  plot(0:1,0:1,ann=F,axes=F,type="n")
  leg.fun(format(round(int.bks,3),nsmall=3),pal,
          "TP Thiel-Sen Slope\n(\u03BCg L\u207B\u00B9 Yr\u207B\u00B9)",
          top.val=leg.pos.y[2],bot.val=leg.pos.y[1],x.min = leg.pos.x[1],x.max = leg.pos.x[2],
          leg.type = "categorical",title.cex = 0.75)
  legend(0.5,leg.pos.y[1],legend=c("Insufficent Data","\u03C1 > 0.05","\u03C1 < 0.05"),
         pch=21,lty=c(NA),lwd=c(0.1),
         col=c("grey"),pt.bg=c("grey",cols.val[2:3]),
         pt.cex=1,ncol=1,cex=0.6,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=1,title.adj = 0,title="Kendall's trend \u03C1-value")
  
  plot(shore,col="cornsilk",border="grey",bg="lightblue",ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.05,xpd=F)
  n=10
  int=classIntervals(values(tps.TP.GM)[is.na(values(tps.TP.GM))==F],n=n,style="equal")
  int.bks=int$brks
  pal=hcl.colors(length(int$brks)-1, "Inferno", rev = F,alpha=0.7)
  image(tps.TP.GM,add=T,breaks=int.bks,col = pal)
  plot(ENP,add=T,bg=NA,lwd=0.5)
  plot(regions2,lty=2,add=T,border="grey80",lwd=0.5)
  plot(sites.shp.TP.GM,add=T,pch=21,cex=0.25,bg=adjustcolor("white",0.5),col=NA)
  box(lwd=1)
  mtext(side=3,adj=0,line=-1.25," B",font=2)
  
  plot(0:1,0:1,ann=F,axes=F,type="n")
  range(int.bks)
  txt.val=c("< 3","40")
  leg.fun(format(round(int.bks,1),nsmall=3),pal,leg.txt=txt.val,
          "Average\nAnnual GM TP\n(\u03BCg L\u207B\u00B9)",
          top.val=leg.pos.y[2],bot.val=leg.pos.y[1],x.min = leg.pos.x[1],x.max = leg.pos.x[2],
          leg.type = "continuous",title.cex = 0.75)

  
  }
{
  # TN
  plot(shore,col="cornsilk",border="grey",bg="lightblue",ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.05,xpd=F)
  int.bks=c(min(values(tps.TN.trend),na.rm=T),0+min(values(tps.TN.trend),na.rm=T)/2,0,0+max(values(tps.TN.trend),na.rm=T)/2,max(values(tps.TN.trend),na.rm=T)) #int$brks;#int.bks[2]=0
  pal=hcl.colors(length(int.bks)-1, "viridis", rev = F,alpha=0.7)
  image(tps.TN.trend,add=T,breaks=int.bks,col = pal)
  plot(ENP,add=T,bg=NA,lwd=0.5)
  plot(regions2,lty=2,add=T,border="grey80",lwd=0.5)
  plot(sites.shp2,add=T,cex=0.40,pch=21,bg=adjustcolor("grey",0.5),col="grey",lwd=0.1)
  plot(sites.shp.TN.trend,add=T,pch=21,cex=0.5,bg=adjustcolor(cols.val,0.5)[sites.shp.TN.trend$stat.sig],col=NA);
  box(lwd=1)
  # mapmisc::scaleBar(utm17,"bottomright",bty="n",cex=1,seg.len=4)
  mtext(side=3,adj=0,line=-1.25," C",font=2)
  
  plot(0:1,0:1,ann=F,axes=F,type="n")
  leg.fun(format(round(int.bks,3),nsmall=3),pal,
          "TN Thiel-Sen Slope\n(mg L\u207B\u00B9 Yr\u207B\u00B9)",
          top.val=leg.pos.y[2],bot.val=leg.pos.y[1],x.min = leg.pos.x[1],x.max = leg.pos.x[2],
          leg.type = "categorical",title.cex = 0.75)
  # legend(0.5,0.4,legend=c("Insufficent Data","\u03C1 > 0.05","\u03C1 < 0.05"),
  #        pch=21,lty=c(NA),lwd=c(0.1),
  #        col=c("grey"),pt.bg=c("grey",cols.val[2:3]),
  #        pt.cex=1,ncol=1,cex=0.6,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=1,title.adj = 0,title="Kendall's trend \u03C1-value")
  
  plot(shore,col="cornsilk",border="grey",bg="lightblue",ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.05,xpd=F)
  n=10
  int=classIntervals(values(tps.TN.GM)[is.na(values(tps.TN.GM))==F],n=n,style="equal")
  int.bks=int$brks
  pal=hcl.colors(length(int$brks)-1, "Inferno", rev = F,alpha=0.7)
  image(tps.TN.GM,add=T,breaks=int.bks,col = pal)
  plot(ENP,add=T,bg=NA,lwd=0.5)
  plot(regions2,lty=2,add=T,border="grey80",lwd=0.5)
  plot(sites.shp.TN.GM,add=T,pch=21,cex=0.25,bg=adjustcolor("white",0.5),col=NA)
  box(lwd=1)
  mtext(side=3,adj=0,line=-1.25," D",font=2)
  
  plot(0:1,0:1,ann=F,axes=F,type="n")
  range(int.bks)
  txt.val=c("< 0.10","1.18")
  leg.fun(format(round(int.bks,1),nsmall=3),pal,leg.txt=txt.val,
          "Average\nAnnual GM TN\n(mg L\u207B\u00B9)",
          top.val=leg.pos.y[2],bot.val=leg.pos.y[1],x.min = leg.pos.x[1],x.max = leg.pos.x[2],
          leg.type = "continuous",title.cex = 0.75)
  
  
}
{
  # TN:TP
  plot(shore,col="cornsilk",border="grey",bg="lightblue",ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.05,xpd=F)
  int.bks=c(min(values(tps.NP.trend),na.rm=T),0+min(values(tps.NP.trend),na.rm=T)/2,0,0+max(values(tps.NP.trend),na.rm=T)/2,max(values(tps.NP.trend),na.rm=T)) #int$brks;#int.bks[2]=0
  pal=hcl.colors(length(int.bks)-1, "viridis", rev = F,alpha=0.7)
  image(tps.NP.trend,add=T,breaks=int.bks,col = pal)
  plot(ENP,add=T,bg=NA,lwd=0.5)
  plot(regions2,lty=2,add=T,border="grey80",lwd=0.5)
  plot(sites.shp2,add=T,cex=0.40,pch=21,bg=adjustcolor("grey",0.5),col="grey",lwd=0.1)
  plot(sites.shp.NP.trend,add=T,pch=21,cex=0.5,bg=adjustcolor(cols.val,0.5)[sites.shp.TN.trend$stat.sig],col=NA);
  box(lwd=1)
  # mapmisc::scaleBar(utm17,"bottomright",bty="n",cex=1,seg.len=4)
  mtext(side=3,adj=0,line=-1.25," E",font=2)
  
  plot(0:1,0:1,ann=F,axes=F,type="n")
  leg.fun(format(round(int.bks,3),nsmall=3),pal,
          "N:P Thiel-Sen Slope\n(Yr\u207B\u00B9)",
          top.val=leg.pos.y[2],bot.val=leg.pos.y[1],x.min = leg.pos.x[1],x.max = leg.pos.x[2],
          leg.type = "categorical",title.cex = 0.75)
  # legend(0.5,0.4,legend=c("Insufficent Data","\u03C1 > 0.05","\u03C1 < 0.05"),
  #        pch=21,lty=c(NA),lwd=c(0.1),
  #        col=c("grey"),pt.bg=c("grey",cols.val[2:3]),
  #        pt.cex=1,ncol=1,cex=0.6,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=1,title.adj = 0,title="Kendall's trend \u03C1-value")
  
  plot(shore,col="cornsilk",border="grey",bg="lightblue",ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.05,xpd=F)
  n=10
  int=classIntervals(values(tps.NP.GM)[is.na(values(tps.NP.GM))==F],n=n,style="equal")
  int.bks=int$brks
  pal=hcl.colors(length(int$brks)-1, "Inferno", rev = F,alpha=0.7)
  image(tps.NP.GM,add=T,breaks=int.bks,col = pal)
  plot(ENP,add=T,bg=NA,lwd=0.5)
  plot(regions2,lty=2,add=T,border="grey80",lwd=0.5)
  plot(sites.shp.NP.GM,add=T,pch=21,cex=0.25,bg=adjustcolor("white",0.5),col=NA)
  box(lwd=1)
  mtext(side=3,adj=0,line=-1.25," F",font=2)
  
  plot(0:1,0:1,ann=F,axes=F,type="n")
  range(int.bks)
  txt.val=c("3","600")
  leg.fun(format(round(int.bks,1),nsmall=3),pal,leg.txt=txt.val,
          "Average\nAnnual GM N:P\n(unitless)",
          top.val=leg.pos.y[2],bot.val=leg.pos.y[1],x.min = leg.pos.x[1],x.max = leg.pos.x[2],
          leg.type = "continuous",title.cex = 0.75)
  
  
}
dev.off()

# png(filename=paste0(plot.path,"revised/TrendMaps_SRP_DIN.png"),width=6.5,height=3.25,units="in",res=200,type="windows",bg="white")
# can add DIN:SRP
par(family="serif",oma=c(0.25,0.25,0.25,0.25),mar=c(0.1,0.1,0.1,0.1),xpd=F)
layout(matrix(c(1:8),2,4,byrow=T),widths = c(1,0.5,1,0.5))
bbox.lims=bbox(region.mask)
leg.pos.x=c(0.25,0.35)
leg.pos.y=c(0.4,0.8)
{
  # SRP
  plot(shore,col="cornsilk",border="grey",bg="lightblue",ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.05,xpd=F)
  int.bks=c(min(values(tps.SRP.trend),na.rm=T),0+min(values(tps.SRP.trend),na.rm=T)/2,0,0+max(values(tps.SRP.trend),na.rm=T)/2,max(values(tps.SRP.trend),na.rm=T)) #int$brks;#int.bks[2]=0
  pal=hcl.colors(length(int.bks)-1, "viridis", rev = F,alpha=0.7)
  image(tps.SRP.trend,add=T,breaks=int.bks,col = pal)
  plot(ENP,add=T,bg=NA,lwd=0.5)
  plot(regions2,lty=2,add=T,border="grey80",lwd=0.5)
  plot(sites.shp2,add=T,cex=0.40,pch=21,bg=adjustcolor("grey",0.5),col="grey",lwd=0.1)
  plot(sites.shp.SRP.trend,add=T,pch=21,cex=0.5,bg=adjustcolor(cols.val,0.5)[sites.shp.TP.trend$stat.sig],col=NA);
  box(lwd=1)
  mapmisc::scaleBar(utm17,"bottomright",bty="n",cex=1,seg.len=4)
  mtext(side=3,adj=0,line=-1.25," A",font=2)
  
  plot(0:1,0:1,ann=F,axes=F,type="n")
  leg.fun(format(round(int.bks,3),nsmall=3),pal,
          "SRP Thiel-Sen Slope\n(\u03BCg L\u207B\u00B9 Yr\u207B\u00B9)",
          top.val=leg.pos.y[2],bot.val=leg.pos.y[1],x.min = leg.pos.x[1],x.max = leg.pos.x[2],
          leg.type = "categorical",title.cex = 0.75)
  legend(0.5,leg.pos.y[1],legend=c("Insufficent Data","\u03C1 > 0.05","\u03C1 < 0.05"),
         pch=21,lty=c(NA),lwd=c(0.1),
         col=c("grey"),pt.bg=c("grey",cols.val[2:3]),
         pt.cex=1,ncol=1,cex=0.6,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=1,title.adj = 0,title="Kendall's trend \u03C1-value")
  
  plot(shore,col="cornsilk",border="grey",bg="lightblue",ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.05,xpd=F)
  n=10
  int=classIntervals(values(tps.SRP.GM)[is.na(values(tps.SRP.GM))==F],n=n,style="equal")
  int.bks=int$brks
  pal=hcl.colors(length(int$brks)-1, "Inferno", rev = F,alpha=0.7)
  image(tps.SRP.GM,add=T,breaks=int.bks,col = pal)
  plot(ENP,add=T,bg=NA,lwd=0.5)
  plot(regions2,lty=2,add=T,border="grey80",lwd=0.5)
  plot(sites.shp.SRP.GM,add=T,pch=21,cex=0.25,bg=adjustcolor("white",0.5),col=NA)
  box(lwd=1)
  mtext(side=3,adj=0,line=-1.25," B",font=2)
  
  plot(0:1,0:1,ann=F,axes=F,type="n")
  range(int.bks)
  txt.val=c("< 1","10")
  leg.fun(format(round(int.bks,1),nsmall=3),pal,leg.txt=txt.val,
          "Average\nAnnual GM SRP\n(\u03BCg L\u207B\u00B9)",
          top.val=leg.pos.y[2],bot.val=leg.pos.y[1],x.min = leg.pos.x[1],x.max = leg.pos.x[2],
          leg.type = "continuous",title.cex = 0.75)
  
  
}
{
  # TN
  plot(shore,col="cornsilk",border="grey",bg="lightblue",ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.05,xpd=F)
  int.bks=c(min(values(tps.DIN.trend),na.rm=T),0+min(values(tps.DIN.trend),na.rm=T)/2,0,0+max(values(tps.DIN.trend),na.rm=T)/2,max(values(tps.DIN.trend),na.rm=T)) #int$brks;#int.bks[2]=0
  pal=hcl.colors(length(int.bks)-1, "viridis", rev = F,alpha=0.7)
  image(tps.DIN.trend,add=T,breaks=int.bks,col = pal)
  plot(ENP,add=T,bg=NA,lwd=0.5)
  plot(regions2,lty=2,add=T,border="grey80",lwd=0.5)
  plot(sites.shp2,add=T,cex=0.40,pch=21,bg=adjustcolor("grey",0.5),col="grey",lwd=0.1)
  plot(sites.shp.DIN.trend,add=T,pch=21,cex=0.5,bg=adjustcolor(cols.val,0.5)[sites.shp.DIN.trend$stat.sig],col=NA);
  box(lwd=1)
  # mapmisc::scaleBar(utm17,"bottomright",bty="n",cex=1,seg.len=4)
  mtext(side=3,adj=0,line=-1.25," C",font=2)
  
  plot(0:1,0:1,ann=F,axes=F,type="n")
  leg.fun(format(round(int.bks,3),nsmall=3),pal,
          "DIN Thiel-Sen Slope\n(mg L\u207B\u00B9 Yr\u207B\u00B9)",
          top.val=leg.pos.y[2],bot.val=leg.pos.y[1],x.min = leg.pos.x[1],x.max = leg.pos.x[2],
          leg.type = "categorical",title.cex = 0.75)
  # legend(0.5,0.4,legend=c("Insufficent Data","\u03C1 > 0.05","\u03C1 < 0.05"),
  #        pch=21,lty=c(NA),lwd=c(0.1),
  #        col=c("grey"),pt.bg=c("grey",cols.val[2:3]),
  #        pt.cex=1,ncol=1,cex=0.6,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=1,title.adj = 0,title="Kendall's trend \u03C1-value")
  
  plot(shore,col="cornsilk",border="grey",bg="lightblue",ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.05,xpd=F)
  n=10
  int=classIntervals(values(tps.DIN.GM)[is.na(values(tps.DIN.GM))==F],n=n,style="equal")
  int.bks=int$brks
  pal=hcl.colors(length(int$brks)-1, "Inferno", rev = F,alpha=0.7)
  image(tps.DIN.GM,add=T,breaks=int.bks,col = pal)
  plot(ENP,add=T,bg=NA,lwd=0.5)
  plot(regions2,lty=2,add=T,border="grey80",lwd=0.5)
  plot(sites.shp.DIN.GM,add=T,pch=21,cex=0.25,bg=adjustcolor("white",0.5),col=NA)
  box(lwd=1)
  mtext(side=3,adj=0,line=-1.25," D",font=2)
  
  plot(0:1,0:1,ann=F,axes=F,type="n")
  range(int.bks)
  txt.val=c("< 0.10","1.00")
  leg.fun(format(round(int.bks,1),nsmall=3),pal,leg.txt=txt.val,
          "Average\nAnnual GM DIN\n(mg L\u207B\u00B9)",
          top.val=leg.pos.y[2],bot.val=leg.pos.y[1],x.min = leg.pos.x[1],x.max = leg.pos.x[2],
          leg.type = "continuous",title.cex = 0.75)
  
  
}
dev.off()

# png(filename=paste0(plot.path,"revised/TrendMaps_Chla_TOC.png"),width=6.5,height=3.25,units="in",res=200,type="windows",bg="white")
par(family="serif",oma=c(0.25,0.25,0.25,0.25),mar=c(0.1,0.1,0.1,0.1),xpd=F)
layout(matrix(c(1:8),2,4,byrow=T),widths = c(1,0.5,1,0.5))
bbox.lims=bbox(region.mask)
leg.pos.x=c(0.25,0.35)
leg.pos.y=c(0.4,0.8)
{
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
  mapmisc::scaleBar(utm17,"bottomright",bty="n",cex=1,seg.len=4)
  mtext(side=3,adj=0,line=-1.25," A",font=2)
  
  plot(0:1,0:1,ann=F,axes=F,type="n")
  leg.fun(format(round(int.bks,3),nsmall=3),pal,
          "Chl-a Thiel-Sen Slope\n(\u03BCg L\u207B\u00B9 Yr\u207B\u00B9)",
          top.val=leg.pos.y[2],bot.val=leg.pos.y[1],x.min = leg.pos.x[1],x.max = leg.pos.x[2],
          leg.type = "categorical",title.cex = 0.75)
  legend(0.5,leg.pos.y[1],legend=c("Insufficent Data","\u03C1 > 0.05","\u03C1 < 0.05"),
         pch=21,lty=c(NA),lwd=c(0.1),
         col=c("grey"),pt.bg=c("grey",cols.val[2:3]),
         pt.cex=1,ncol=1,cex=0.6,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=1,title.adj = 0,title="Kendall's trend \u03C1-value")
  
  plot(shore,col="cornsilk",border="grey",bg="lightblue",ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.05,xpd=F)
  n=10
  int=classIntervals(values(tps.Chla.GM)[is.na(values(tps.Chla.GM))==F],n=n,style="equal")
  int.bks=int$brks
  pal=hcl.colors(length(int$brks)-1, "Inferno", rev = F,alpha=0.7)
  image(tps.Chla.GM,add=T,breaks=int.bks,col = pal)
  plot(ENP,add=T,bg=NA,lwd=0.5)
  plot(regions2,lty=2,add=T,border="grey80",lwd=0.5)
  plot(sites.shp.Chla.GM,add=T,pch=21,cex=0.25,bg=adjustcolor("white",0.5),col=NA)
  box(lwd=1)
  mtext(side=3,adj=0,line=-1.25," B",font=2)
  
  plot(0:1,0:1,ann=F,axes=F,type="n")
  range(int.bks)
  txt.val=c("< 0.05","6.5")
  leg.fun(format(round(int.bks,1),nsmall=3),pal,leg.txt=txt.val,
          "Average\nAnnual GM Chl-a\n(\u03BCg L\u207B\u00B9)",
          top.val=leg.pos.y[2],bot.val=leg.pos.y[1],x.min = leg.pos.x[1],x.max = leg.pos.x[2],
          leg.type = "continuous",title.cex = 0.75)
  
  
}
{
  # TOC
  plot(shore,col="cornsilk",border="grey",bg="lightblue",ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.05,xpd=F)
  int.bks=c(min(values(tps.TOC.trend),na.rm=T),0+min(values(tps.TOC.trend),na.rm=T)/2,0,0+max(values(tps.TOC.trend),na.rm=T)/2,max(values(tps.TOC.trend),na.rm=T)) #int$brks;#int.bks[2]=0
  pal=hcl.colors(length(int.bks)-1, "viridis", rev = F,alpha=0.7)
  image(tps.TOC.trend,add=T,breaks=int.bks,col = pal)
  plot(ENP,add=T,bg=NA,lwd=0.5)
  plot(regions2,lty=2,add=T,border="grey80",lwd=0.5)
  plot(sites.shp2,add=T,cex=0.40,pch=21,bg=adjustcolor("grey",0.5),col="grey",lwd=0.1)
  plot(sites.shp.TOC.trend,add=T,pch=21,cex=0.5,bg=adjustcolor(cols.val,0.5)[sites.shp.TOC.trend$stat.sig],col=NA);
  box(lwd=1)
  # mapmisc::scaleBar(utm17,"bottomright",bty="n",cex=1,seg.len=4)
  mtext(side=3,adj=0,line=-1.25," C",font=2)
  
  plot(0:1,0:1,ann=F,axes=F,type="n")
  leg.fun(format(round(int.bks,3),nsmall=3),pal,
          "TOC Thiel-Sen Slope\n(mg L\u207B\u00B9 Yr\u207B\u00B9)",
          top.val=leg.pos.y[2],bot.val=leg.pos.y[1],x.min = leg.pos.x[1],x.max = leg.pos.x[2],
          leg.type = "categorical",title.cex = 0.75)
  # legend(0.5,0.4,legend=c("Insufficent Data","\u03C1 > 0.05","\u03C1 < 0.05"),
  #        pch=21,lty=c(NA),lwd=c(0.1),
  #        col=c("grey"),pt.bg=c("grey",cols.val[2:3]),
  #        pt.cex=1,ncol=1,cex=0.6,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=1,title.adj = 0,title="Kendall's trend \u03C1-value")
  
  plot(shore,col="cornsilk",border="grey",bg="lightblue",ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.05,xpd=F)
  n=10
  int=classIntervals(values(tps.TOC.GM)[is.na(values(tps.TOC.GM))==F],n=n,style="equal")
  int.bks=int$brks
  pal=hcl.colors(length(int$brks)-1, "Inferno", rev = F,alpha=0.7)
  image(tps.TOC.GM,add=T,breaks=int.bks,col = pal)
  plot(ENP,add=T,bg=NA,lwd=0.5)
  plot(regions2,lty=2,add=T,border="grey80",lwd=0.5)
  plot(sites.shp.TOC.GM,add=T,pch=21,cex=0.25,bg=adjustcolor("white",0.5),col=NA)
  box(lwd=1)
  mtext(side=3,adj=0,line=-1.25," D",font=2)
  
  plot(0:1,0:1,ann=F,axes=F,type="n")
  range(int.bks)
  txt.val=c("< 0.5","20")
  leg.fun(format(round(int.bks,1),nsmall=3),pal,leg.txt=txt.val,
          "Average\nAnnual GM TOC\n(mg L\u207B\u00B9)",
          top.val=leg.pos.y[2],bot.val=leg.pos.y[1],x.min = leg.pos.x[1],x.max = leg.pos.x[2],
          leg.type = "continuous",title.cex = 0.75)
  
  
}
dev.off()


# workspace image ---------------------------------------------------------
## Save workspace image to this point
# save.image(file=paste0(export.path,"regiontrend_reanalysis.RData"))
# load(paste0(export.path,"regiontrend_reanalysis.RData"))

# TP and TN - SALINITY GAM ------------------------------------------------
library(sf)
library(EVERSpatDat)
nad83.pro=st_crs("EPSG:4269")
utm17=st_crs("EPSG:26917")
wgs84=st_crs("EPSG:4326")

dec.month=function(date){
  yr=as.numeric(format(date,"%Y"))
  leap=(yr%%4==0)&((yr%%100!=0)|(yr%%400 == 0))
  
  month_x <- month.abb[as.numeric(format(date,"%m"))]
  N_DAYS_IN_MONTHS <- c(
    Jan = 31L, Feb = 28L, Mar = 31L,
    Apr = 30L, May = 31L, Jun = 30L,
    Jul = 31L, Aug = 31L, Sep = 30L,
    Oct = 31L, Nov = 30L, Dec = 31L
  )
  
  n_days <- N_DAYS_IN_MONTHS[month_x]
  n_days[month_x == "Feb" & leap==T] <- 29L
  
  decmon=as.numeric(format(date,"%m")) + as.numeric(format(date,"%d"))/n_days
  return(decmon)
}

####
region.mask.sf=st_as_sf(region.mask)
regions2.sf=st_as_sf(regions2)
ENP.sf=st_as_sf(ENP)


data("sfwmd_bound")

region.shore=st_intersection(region.mask.sf,sfwmd_bound)

###
tmp=ddply(subset(dat.all2.melt,sea.screen==1&consec.screen==1&variable=="SALINITY"),
          c("STATION"),summarise,N.sal=N.obs(value))

sal.sites=merge(tmp,sites.shp.TableS1,"STATION",all.y=T)
sal.sites.shp=st_as_sf(sal.sites,coords = c("UTMX","UTMY"),crs=utm17)

plot(st_geometry(sal.sites.shp),
     pch=ifelse(is.na(sal.sites.shp$N.sal)==T,4,21),
     bg="dodgerblue1")

da.dat=dcast(subset(dat.all2.melt,sea.screen==1&consec.screen==1),
             STATION+DATE+UTMX+UTMY+Region~variable,value.var="value",mean,na.rm=T)|>
  mutate(month=as.numeric(format(DATE,"%m")),
         CY=as.numeric(format(DATE,"%Y")),
         DOY=as.numeric(format(DATE,"%j")),
         DOWY=hydro.day(DATE,"Fed"),
         WY=WY(DATE,'Fed'),
         decMonth=dec.month(DATE),
         decWY=decimal.WY(DATE,"Fed"),
         decCY=lubridate::decimal_date(DATE)
  )
da.dat[is.na(da.dat$SALINITY)==F&da.dat$SALINITY>=100,]=NA; # error in data (salinity too high, outside the distribution of the data)

da.dat.WY=ddply(da.dat,c("WY","STATION","UTMX","UTMY","Region"), summarise,
                mean.TP = mean(TP,na.rm=T),
                mean.TN = mean(TN,na.rm=T),
                mean.sal = mean(SALINITY,na.rm=T))
da.dat.WY=subset(da.dat.WY,is.na(WY)==F)
## TP exploring
# library(fitdistrplus)
# library(MASS)
# 
# descdist(subset(da.dat,is.na(TP)==F)$TP,discrete = F)
# tmp=fitdist(subset(da.dat,is.na(TP)==F)$TP,"lnorm")
# plot(tmp)
# 
# descdist(subset(da.dat.WY,is.na(mean.TP)==F)$mean.TP,discrete = F)
# tmp=fitdist(subset(da.dat.WY,is.na(mean.TP)==F)$mean.TP,"lnorm")
# plot(tmp)


## GAM
library(parallel)
detectCores()
nc <- 5
cl <- makeCluster(nc)

clusterEvalQ(cl,c(library(mgcv)))

da.dat$Region.f=as.factor(da.dat$Region)
da.dat.WY$Region.f=as.factor(da.dat.WY$Region)


## TP GAM ------------------------------------------------------------------
# load(file=paste0(export.path,"TPSal_GAM.Rdata"))
# load(file=paste0(export.path,"TPSal_GAM_summary.Rdata"))
ctrl = gam.control(trace = T)
# TP.sal.WY.m.re<-bam(log(mean.TP)~
#                       s(Region.f,bs="re",k=5)+
#                       s(mean.sal,k=20)+
#                       s(WY,bs="cr",k=25)+
#                       s(UTMX,UTMY,bs="ds",m=c(1,0.5))+
#                       ti(UTMX,UTMY,WY,d=c(2,1),bs=c("ds","cr"),m = list(c(1, 0.5), NA)),
#                     data=da.dat.WY, nthreads = c(10,1), discrete=T,control = ctrl,chunk.size = 1000,samfrac=0.1
# )
# TP.sal.WY.m.re.sum=summary(TP.sal.WY.m.re);TP.sal.WY.m.re.sum
# layout(matrix(1:4,2,2,byrow=T));gam.check(TP.sal.WY.m.re,pch=21,col="lightblue",bg="grey")
# abline(0,1,lty=2)
# dev.off();plot(TP.sal.WY.m.re,pages=1)
# acf(residuals(TP.sal.WY.m.re))
# pacf(residuals(TP.sal.WY.m.re))
# 
# length(unique(da.dat.WY$WY))
# nrow(unique(da.dat.WY[,c("UTMX","UTMY")]))

TP.sal.WY.mi<-bam((mean.TP)~
                    s(mean.sal)+
                    s(WY,bs="cr")+
                    s(UTMX,UTMY,bs="ds",m=c(1,0.5))+
                    ti(UTMX,UTMY,WY,d=c(2,1),bs=c("ds","cr"),m = list(c(1, 0.5), NA)),
                  data=da.dat.WY, nthreads = c(10,1), discrete=T,
                  family=Gamma(link="log"),
                  control = ctrl,chunk.size = 1000,samfrac=0.1)
TP.sal.WY.mi.sum=summary(TP.sal.WY.mi);TP.sal.WY.mi.sum
layout(matrix(1:4,2,2,byrow=T));gam.check(TP.sal.WY.mi,pch=21,col="lightblue",bg="grey");abline(0,1,lty=2)

## got idea from here to speed up model fitting https://m-clark.github.io/posts/2019-10-20-big-mixed-models/#function-arguments
## and here https://stackoverflow.com/a/75352185/5213091
TP.sal.WY.m<-bam(log(mean.TP)~
                   s(mean.sal,k=10)+
                   s(WY,bs="cr",k=26)+
                   s(UTMX,UTMY,bs="ds",m=c(1,0.5),k=250)+
                   ti(UTMX,UTMY,WY,d=c(2,1),bs=c("ds","cr"),m = list(c(1, 0.5), NA),k=c(200,26)),
                 data=da.dat.WY, nthreads = c(10,1), discrete=T,
                 control = ctrl,chunk.size = 1000,samfrac=0.1)

## took ~30 minutes to fit
# save(TP.sal.WY.m,file=paste0(export.path,"TPSal_GAM.Rdata"))
TP.sal.WY.m.sum=summary(TP.sal.WY.m);TP.sal.WY.m.sum
# save(TP.sal.WY.m.sum,file=paste0(export.path,"TPSal_GAM_summary.Rdata"))

layout(matrix(1:4,2,2,byrow=T));gam.check(TP.sal.WY.m,pch=21,col="lightblue",bg="grey")
abline(0,1,lty=2)

dev.off();plot(TP.sal.WY.m,pages=1)
acf(residuals(TP.sal.WY.m))
pacf(residuals(TP.sal.WY.m))
# testResiduals(simulateResiduals(TP.sal.WY.m))

# residual check
pred.org=(predict(TP.sal.WY.m,type="terms"))
TP.sal.partial.resids<-pred.org+residuals(TP.sal.WY.m)
ncol(TP.sal.partial.resids)
layout(matrix(1:4,2:2))
for(i in 1:ncol(TP.sal.partial.resids)){hist(TP.sal.partial.resids[,i],main=smooths(TP.sal.WY.m)[i])}  


### TP Derivative -----------------------------------------------------------
smooths(TP.sal.WY.m)
dev.vars=c("Sal","WY","UTM","UTMWY")
TP.sal.WY.m.sal.d=dev.data.val(TP.sal.WY.m,1,n=400,var.names=dev.vars)
TP.sal.WY.m.WY.d=dev.data.val(TP.sal.WY.m,2,n=400,var.names=dev.vars)

# png(filename=paste0(plot.path,"TP_sal_gam.png"),width=6.5,height=4,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(2,2.75,1,0.5),oma=c(1,1,0.25,0.25),xpd=F);
layout(matrix(c(1:4),2,2,byrow=T))# ,widths=c(1,1,1.1,0.4))

ylim.val=c(-0.5,0.5);by.y=0.25;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(0,45);by.x=10;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
plot(fit.Sal~mean.sal,TP.sal.WY.m.sal.d,ylim=ylim.val,xlim=xlim.val,ann=F,axes=F,type="n")
abline(h=ymaj,v=xmaj,lty=3,col="grey")
abline(h=0)
# points(jitter(TP.sal.WY.m$model$mean.sal,0.5),TP.sal.partial.resids[,1],pch=19,col=adjustcolor("dodgerblue1",0.5))
with(TP.sal.WY.m.sal.d,shaded.range(mean.sal,lower.CI,upper.CI,"grey",lty=1))
lines(fit.Sal~mean.sal,TP.sal.WY.m.sal.d,lwd=2)
lines(dsig.incr ~ mean.sal, data = TP.sal.WY.m.sal.d, col = "red", lwd = 3,lty=1)
lines(dsig.decr ~ mean.sal, data = TP.sal.WY.m.sal.d, col = "blue", lwd = 3,lty=1)
axis_fun(1,xmaj,xmin,xmaj,line=-0.5);axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
# mtext(side=3,adj=0,paste0("Model 1 (Chla~",paste(smooths(chla.m1),collapse=" + ")))
# mtext(side=3,adj=0,"Model 1 (with STG)")
mtext(side=3,adj=1,"s(Salinity)")
mtext(side=1,line=1.75,"Salinity (ppt)",cex=0.9)
mtext(side=2,line=2.5,"TP Effect")
legend("topleft",legend=c("Sig. Increase","Sig. Decrease"),
       pch=c(NA),lty=c(1),lwd=c(1.5),
       col=c("red","blue"),pt.bg=NA,
       pt.cex=1.5,ncol=1,cex=0.75,bty="n",y.intersp=0.75,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=1)

range(TP.sal.WY.m.WY.d$fit.WY)
range(TP.sal.WY.m.WY.d[c("upper.CI","lower.CI")])
ylim.val=c(-1,0.5);by.y=0.5;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(1990,2023);by.x=10;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/by.x)
plot(fit.WY~WY,TP.sal.WY.m.WY.d,ylim=ylim.val,xlim=xlim.val,ann=F,axes=F,type="n")
abline(h=ymaj,v=xmaj,lty=3,col="grey")
abline(h=0)
# points(jitter(TP.sal.WY.m$model$WY,0.5),TP.sal.partial.resids[,2],pch=19,col=adjustcolor("dodgerblue1",0.5))
with(TP.sal.WY.m.WY.d,shaded.range(WY,lower.CI,upper.CI,"grey",lty=1))
lines(fit.WY~WY,TP.sal.WY.m.WY.d,lwd=2)
lines(dsig.incr ~ WY, data = TP.sal.WY.m.WY.d, col = "red", lwd = 3,lty=1)
lines(dsig.decr ~ WY, data = TP.sal.WY.m.WY.d, col = "blue", lwd = 3,lty=1)
axis_fun(1,xmaj,xmin,xmaj,line=-0.5);axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=3,adj=1,"s(WY)")
mtext(side=1,line=1.75,"Water Year",cex=0.9)

reg.ext=extent(region.mask.sf)
pdat.sp.UTM.TP.sal.m=expand.grid(
  mean.sal=TP.sal.WY.m$var.summary$mean.sal[2],
  WY=TP.sal.WY.m$var.summary$WY[2],
  UTMX=seq(reg.ext[1],reg.ext[2],by=2000),
  UTMY=seq(reg.ext[3],reg.ext[4],by=2000)
)
UTM.pred=predict(TP.sal.WY.m,pdat.sp.UTM.TP.sal.m,type="terms",se.fit = F,nthreads=6,discrete=T)
colnames(UTM.pred)=paste0("fit.",dev.vars)
pdat.sp.UTM.TP.sal.m=cbind(pdat.sp.UTM.TP.sal.m,UTM.pred)

tmp.ma=with(pdat.sp.UTM.TP.sal.m,matrix(fit.UTM,nrow=length(unique(UTMX)),ncol=length(unique(UTMY))))
dat1=list()
dat1$x=unique(pdat.sp.UTM.TP.sal.m$UTMX)
dat1$y=unique(pdat.sp.UTM.TP.sal.m$UTMY)
dat1$z=round(tmp.ma,2)
r=raster(dat1)
TPSal.GAM.UTM <- mask(r, region.mask.sf)

par(mar=c(0.5,0.5,2.5,0.5))
bbox.lims=st_bbox(st_buffer(region.mask.sf,-1000))
b1=seq(-1,1.5,0.5)
b2=seq(-1.1,1.5,length.out=50)
pal2=viridis(length(b2)-1,alpha=0.75)# colorRampPalette(c("blue","grey90","red"))(length(b2)-1)
plot(st_geometry(region.mask.sf),ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],border=NA)
image(TPSal.GAM.UTM,add=T,breaks=b2,col = pal2)
plot(rasterToContour(TPSal.GAM.UTM,levels=b1,nlevels=length(b1)),col=adjustcolor("black",0.5),lwd=1,add=T)
plot(st_geometry(region.mask.sf),add=T)
plot(st_geometry(subset(sal.sites.shp,N.sal>1)),add=T,pch=19,bg=NA,col=adjustcolor("grey2",0.5),cex=0.1)
plot(st_geometry(ENP.sf),add=T,bg=NA,lwd=0.5)
plot(st_geometry(regions2.sf),lty=2,add=T,border="white",lwd=0.5)
mapmisc::scaleBar(crs=region.mask.sf,"bottomleft",bty="n",cex=1,seg.len=4,outer=F)
box(lwd=1)
mtext(side=3,adj=0,"s(UTMX,UTMY)")

plot(0:1,0:1,ann=F,axes=F,type="n")
leg.fun(round(b2,1),pal2,leg.title="s(UTMX,UTMY)\nEffect",leg.type = "continuous",
        x.max=0.6,x.min=0.4)
legend(0.5,0.15,legend=c("Spatial Effect Contour","Monitoring Locations"),
       pch=c(NA,19),lty=c(1,NA),lwd=c(1,0.01),
       col=c(adjustcolor("black",0.5),"grey2"),pt.bg=c(NA,NA),
       pt.cex=1.5,ncol=1,cex=0.8,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=1)

dev.off()


## marginaleffects
# tmp=plot_slopes(TP.sal.WY.m, variables = "mean.sal", condition = "mean.sal")
# tmp$plot$dat
# lines(estimate~mean.sal,tmp$plot$dat,col="red",lwd=3)

## marginaleffects plot_slopes/slopes is using the 1st derivative...gratia is faster 


## response variable plotting
## https://ecogambler.netlify.app/blog/interpreting-gams/
## https://www.maths.ed.ac.uk/~swood34/mgcv/check-select.pdf
crit.t <- qt(0.025, df.residual(TP.sal.WY.m), lower.tail = FALSE)

sal.resp.dat=predict(TP.sal.WY.m,TP.sal.WY.m.sal.d,type="response",se=T)|>
  mutate(
    resp.UCI=exp(fit+(crit.t*se.fit)),
    resp.LCI=exp(fit-(crit.t*se.fit)),
    fit=exp(fit)
  )|>
  as.data.frame()
TP.sal.WY.m.sal.d=cbind(TP.sal.WY.m.sal.d,sal.resp.dat)

# png(filename=paste0(plot.path,"TP_sal_gam_salResponsePlot.png"),width=5,height=3,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(2,2.75,1,0.5),oma=c(1,1,0.25,0.25),xpd=F);

ylim.val=c(0,30);by.y=10;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(0,45);by.x=10;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
plot(fit~mean.sal,TP.sal.WY.m.sal.d,ylim=ylim.val,xlim=xlim.val,ann=F,axes=F,type="n")
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(TP.sal.WY.m.sal.d,shaded.range(mean.sal,resp.LCI,resp.UCI,"grey",lty=0))
lines(fit~mean.sal,TP.sal.WY.m.sal.d)
lines(resp.UCI~mean.sal,TP.sal.WY.m.sal.d,lty=2)
lines(resp.LCI~mean.sal,TP.sal.WY.m.sal.d,lty=2)

axis_fun(1,xmaj,xmin,xmaj,line=-0.5);axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=2,line=2,"Annual Mean TP (\u03BCg L\u207B\u00B9)")
mtext(side=1,line=1.75,"Annual Mean Salinity",cex=0.9)
legend("bottomleft",legend=c("Expected Response","95% CI"),
       pch=c(NA),lty=c(1,2),lwd=c(1),
       col=c("black"),pt.bg=NA,
       pt.cex=1.5,ncol=1,cex=0.75,bty="n",y.intersp=0.75,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=1)
dev.off()


# png(filename=paste0(plot.path,"TP_sal_gam_sal_SlopePlot.png"),width=5,height=3,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(2,2.75,1,0.5),oma=c(1,1,0.25,0.25),xpd=F);

ylim.val=c(-0.06,0.15);by.y=0.03;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(0,45);by.x=10;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
plot(derivative~mean.sal,TP.sal.WY.m.sal.d,ylim=ylim.val,xlim=xlim.val,ann=F,axes=F,type="n")
abline(h=ymaj,v=xmaj,lty=3,col="grey")
abline(h=0)
with(TP.sal.WY.m.sal.d,shaded.range(mean.sal,lower,upper,"grey",lty=0))
lines(derivative~mean.sal,TP.sal.WY.m.sal.d)
lines(upper~mean.sal,TP.sal.WY.m.sal.d,lty=2)
lines(lower~mean.sal,TP.sal.WY.m.sal.d,lty=2)

axis_fun(1,xmaj,xmin,xmaj,line=-0.5);axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=2,line=2.5,"Slope")
mtext(side=1,line=1.75,"Annual Mean Salinity",cex=0.9)
legend("topright",legend=c("1st derivative of linear predictor","95% CI"),
       pch=c(NA),lty=c(1,2),lwd=c(1),
       col=c("black"),pt.bg=NA,
       pt.cex=1.5,ncol=1,cex=0.75,bty="n",y.intersp=0.75,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=1)
dev.off()

## TN GAM ------------------------------------------------------------------
# load(file=paste0(export.path,"TNSal_GAM.Rdata"))
# load(file=paste0(export.path,"TNSal_GAM_summary.Rdata"))
ctrl = gam.control(trace = T)

TN.sal.WY.mi<-bam((mean.TN)~
                   s(mean.sal)+
                   s(WY,bs="cr")+
                   s(UTMX,UTMY,bs="ds",m=c(1,0.5))+
                   ti(UTMX,UTMY,WY,d=c(2,1),bs=c("ds","cr"),m = list(c(1, 0.5), NA)),
                 data=da.dat.WY, nthreads = c(10,1), discrete=T,
                 family=Gamma(link="log"),
                 control = ctrl,chunk.size = 1000,samfrac=0.1)
TN.sal.WY.mi.sum=summary(TN.sal.WY.mi);TN.sal.WY.mi.sum
layout(matrix(1:4,2,2,byrow=T));gam.check(TN.sal.WY.mi,pch=21,col="lightblue",bg="grey");abline(0,1,lty=2)

dev.off();plot(TN.sal.WY.mi,pages=1)
acf(residuals(TN.sal.WY.mi))
pacf(residuals(TN.sal.WY.mi))

TN.sal.WY.m<-bam(mean.TN~
                   s(mean.sal,k=20)+
                   s(WY,bs="cr",k=26)+
                   s(UTMX,UTMY,bs="ds",m=c(1,0.5),k=250)+
                   ti(UTMX,UTMY,WY,d=c(2,1),bs=c("ds","cr"),m = list(c(1, 0.5), NA),k=c(200,26)),
                 data=da.dat.WY, nthreads = c(10,1), discrete=T,
                 family=Gamma(link="log"),
                 control = ctrl,chunk.size = 1000,samfrac=0.1)

## took ~63 minutes to fit
# save(TN.sal.WY.m,file=paste0(export.path,"TNSal_GAM.Rdata"))
TN.sal.WY.m.sum=summary(TN.sal.WY.m);TN.sal.WY.m.sum
# save(TN.sal.WY.m.sum,file=paste0(export.path,"TNSal_GAM_summary.Rdata"))
layout(matrix(1:4,2,2,byrow=T));gam.check(TN.sal.WY.m,pch=21,col="lightblue",bg="grey");abline(0,1,lty=2)

dev.off();plot(TN.sal.WY.m,pages=1)
acf(residuals(TN.sal.WY.m))
pacf(residuals(TN.sal.WY.m))
# testResiduals(simulateResiduals(TN.sal.WY.m))

# residual check
pred.org=(predict(TN.sal.WY.m,type="terms"))
TN.sal.partial.resids<-pred.org+residuals(TN.sal.WY.m)
ncol(TN.sal.partial.resids)
layout(matrix(1:4,2:2))
for(i in 1:ncol(TN.sal.partial.resids)){hist(TN.sal.partial.resids[,i],main=smooths(TN.sal.WY.m)[i])}  

### TN Derivative -----------------------------------------------------------
smooths(TN.sal.WY.m)
dev.vars=c("Sal","WY","UTM","UTMWY")
TN.sal.WY.m.sal.d=dev.data.val(TN.sal.WY.m,1,n=400,var.names=dev.vars)
TN.sal.WY.m.WY.d=dev.data.val(TN.sal.WY.m,2,n=400,var.names=dev.vars)

# png(filename=paste0(plot.path,"TN_sal_gam.png"),width=6.5,height=4,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(2,2.75,1,0.5),oma=c(1,1,0.25,0.25),xpd=F);
layout(matrix(c(1:4),2,2,byrow=T))# ,widths=c(1,1,1.1,0.4))

ylim.val=c(-0.25,0.5);by.y=0.25;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(0,45);by.x=10;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
plot(fit.Sal~mean.sal,TN.sal.WY.m.sal.d,ylim=ylim.val,xlim=xlim.val,ann=F,axes=F,type="n")
abline(h=ymaj,v=xmaj,lty=3,col="grey")
abline(h=0)
# points(jitter(TN.sal.WY.m$model$mean.sal,0.5),TN.sal.partial.resids[,1],pch=19,col=adjustcolor("dodgerblue1",0.5))
with(TN.sal.WY.m.sal.d,shaded.range(mean.sal,lower.CI,upper.CI,"grey",lty=1))
lines(fit.Sal~mean.sal,TN.sal.WY.m.sal.d,lwd=2)
lines(dsig.incr ~ mean.sal, data = TN.sal.WY.m.sal.d, col = "red", lwd = 3,lty=1)
lines(dsig.decr ~ mean.sal, data = TN.sal.WY.m.sal.d, col = "blue", lwd = 3,lty=1)
axis_fun(1,xmaj,xmin,xmaj,line=-0.5);axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
# mtext(side=3,adj=0,paste0("Model 1 (Chla~",paste(smooths(chla.m1),collapse=" + ")))
# mtext(side=3,adj=0,"Model 1 (with STG)")
mtext(side=3,adj=1,"s(Salinity)")
mtext(side=1,line=1.75,"Salinity (ppt)",cex=0.9)
mtext(side=2,line=2.5,"TN Effect")
legend("topleft",legend=c("Sig. Increase","Sig. Decrease"),
       pch=c(NA),lty=c(1),lwd=c(1.5),
       col=c("red","blue"),pt.bg=NA,
       pt.cex=1.5,ncol=1,cex=0.75,bty="n",y.intersp=0.75,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=1)

range(TN.sal.WY.m.WY.d$fit.WY)
range(TN.sal.WY.m.WY.d[c("upper.CI","lower.CI")])
ylim.val=c(-1,1);by.y=0.5;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(1990,2023);by.x=10;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/by.x)
plot(fit.WY~WY,TN.sal.WY.m.WY.d,ylim=ylim.val,xlim=xlim.val,ann=F,axes=F,type="n")
abline(h=ymaj,v=xmaj,lty=3,col="grey")
abline(h=0)
# points(jitter(TN.sal.WY.m$model$WY,0.5),TN.sal.partial.resids[,2],pch=19,col=adjustcolor("dodgerblue1",0.5))
with(TN.sal.WY.m.WY.d,shaded.range(WY,lower.CI,upper.CI,"grey",lty=1))
lines(fit.WY~WY,TN.sal.WY.m.WY.d,lwd=2)
lines(dsig.incr ~ WY, data = TN.sal.WY.m.WY.d, col = "red", lwd = 3,lty=1)
lines(dsig.decr ~ WY, data = TN.sal.WY.m.WY.d, col = "blue", lwd = 3,lty=1)
axis_fun(1,xmaj,xmin,xmaj,line=-0.5);axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=3,adj=1,"s(WY)")
mtext(side=1,line=1.75,"Water Year",cex=0.9)

reg.ext=extent(region.mask.sf)
pdat.sp.UTM.TN.sal.m=expand.grid(
  mean.sal=TN.sal.WY.m$var.summary$mean.sal[2],
  WY=TN.sal.WY.m$var.summary$WY[2],
  UTMX=seq(reg.ext[1],reg.ext[2],by=2000),
  UTMY=seq(reg.ext[3],reg.ext[4],by=2000)
)
UTM.pred=predict(TN.sal.WY.m,pdat.sp.UTM.TN.sal.m,type="terms",se.fit = F,nthreads=6,discrete=T)
colnames(UTM.pred)=paste0("fit.",dev.vars)
pdat.sp.UTM.TN.sal.m=cbind(pdat.sp.UTM.TN.sal.m,UTM.pred)

tmp.ma=with(pdat.sp.UTM.TN.sal.m,matrix(fit.UTM,nrow=length(unique(UTMX)),ncol=length(unique(UTMY))))
dat1=list()
dat1$x=unique(pdat.sp.UTM.TN.sal.m$UTMX)
dat1$y=unique(pdat.sp.UTM.TN.sal.m$UTMY)
dat1$z=round(tmp.ma,2)
r=raster(dat1)
TNSal.GAM.UTM <- mask(r, region.mask.sf)

par(mar=c(0.5,0.5,2.5,0.5))
bbox.lims=st_bbox(st_buffer(region.mask.sf,-1000))
b1=seq(-1,1.6,0.5)
b2=seq(-1.1,1.6,length.out=50)
pal2=viridis(length(b2)-1,alpha=0.75)# colorRampPalette(c("blue","grey90","red"))(length(b2)-1)
plot(st_geometry(region.mask.sf),ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],border=NA)
image(TNSal.GAM.UTM,add=T,breaks=b2,col = pal2)
plot(rasterToContour(TNSal.GAM.UTM,levels=b1,nlevels=length(b1)),col=adjustcolor("black",0.5),lwd=1,add=T)
plot(st_geometry(region.mask.sf),add=T)
plot(st_geometry(subset(sal.sites.shp,N.sal>1)),add=T,pch=19,bg=NA,col=adjustcolor("grey2",0.5),cex=0.1)
plot(st_geometry(ENP.sf),add=T,bg=NA,lwd=0.5)
plot(st_geometry(regions2.sf),lty=2,add=T,border="white",lwd=0.5)
mapmisc::scaleBar(crs=region.mask.sf,"bottomleft",bty="n",cex=1,seg.len=4,outer=F)
box(lwd=1)
mtext(side=3,adj=0,"s(UTMX,UTMY)")

plot(0:1,0:1,ann=F,axes=F,type="n")
leg.fun(round(b2,1),pal2,leg.title="s(UTMX,UTMY)\nEffect",leg.type = "continuous",
        x.max=0.6,x.min=0.4)
legend(0.5,0.15,legend=c("Spatial Effect Contour","Monitoring Locations"),
       pch=c(NA,19),lty=c(1,NA),lwd=c(1,0.01),
       col=c(adjustcolor("black",0.5),"grey2"),pt.bg=c(NA,NA),
       pt.cex=1.5,ncol=1,cex=0.8,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=1)

dev.off()

## response variable plotting
## https://ecogambler.netlify.app/blog/interpreting-gams/
## https://www.maths.ed.ac.uk/~swood34/mgcv/check-select.pdf
crit.t <- qt(0.025, df.residual(TN.sal.WY.m), lower.tail = FALSE)

sal.resp.dat=predict(TN.sal.WY.m,TN.sal.WY.m.sal.d,type="response",se=T)|>
  mutate(
    resp.UCI=exp(fit+(crit.t*se.fit)),
    resp.LCI=exp(fit-(crit.t*se.fit)),
    fit=fit
  )|>
  as.data.frame()
TN.sal.WY.m.sal.d=cbind(TN.sal.WY.m.sal.d,sal.resp.dat)

# png(filename=paste0(plot.path,"TN_sal_gam_salResponsePlot.png"),width=5,height=3,units="in",res=200,type="cairo",bg="white")
par(family="serif",mar=c(2,2.75,1,0.5),oma=c(1,1,0.25,0.25),xpd=F);

ylim.val=c(1,1.3);by.y=0.1;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(0,45);by.x=10;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
plot(fit~mean.sal,TN.sal.WY.m.sal.d,ylim=ylim.val,xlim=xlim.val,ann=F,axes=F,type="n")
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(TN.sal.WY.m.sal.d,shaded.range(mean.sal,resp.LCI,resp.UCI,"grey",lty=0))
lines(fit~mean.sal,TN.sal.WY.m.sal.d)
lines(resp.UCI~mean.sal,TN.sal.WY.m.sal.d,lty=2)
lines(resp.LCI~mean.sal,TN.sal.WY.m.sal.d,lty=2)

axis_fun(1,xmaj,xmin,xmaj,line=-0.5);axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=2,line=2,"Annual Mean TN (\u03BCg L\u207B\u00B9)")
mtext(side=1,line=1.75,"Annual Mean Salinity",cex=0.9)
legend("bottomleft",legend=c("Expected Response","95% CI"),
       pch=c(NA),lty=c(1,2),lwd=c(1),
       col=c("black"),pt.bg=NA,
       pt.cex=1.5,ncol=1,cex=0.75,bty="n",y.intersp=0.75,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=1)
dev.off()


# png(filename=paste0(plot.path,"TN_sal_gam_sal_SlopePlot.png"),width=5,height=3,units="in",res=200,type="cairo",bg="white")
par(family="serif",mar=c(2,2.75,1,0.5),oma=c(1,1,0.25,0.25),xpd=F);

ylim.val=c(-0.06,0.15);by.y=0.03;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(0,45);by.x=10;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
plot(derivative~mean.sal,TN.sal.WY.m.sal.d,ylim=ylim.val,xlim=xlim.val,ann=F,axes=F,type="n")
abline(h=ymaj,v=xmaj,lty=3,col="grey")
abline(h=0)
with(TN.sal.WY.m.sal.d,shaded.range(mean.sal,lower,upper,"grey",lty=0))
lines(derivative~mean.sal,TN.sal.WY.m.sal.d)
lines(upper~mean.sal,TN.sal.WY.m.sal.d,lty=2)
lines(lower~mean.sal,TN.sal.WY.m.sal.d,lty=2)

axis_fun(1,xmaj,xmin,xmaj,line=-0.5);axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=2,line=2.5,"Slope")
mtext(side=1,line=1.75,"Annual Mean Salinity",cex=0.9)
legend("topright",legend=c("1st derivative of linear predictor","95% CI"),
       pch=c(NA),lty=c(1,2),lwd=c(1),
       col=c("black"),pt.bg=NA,
       pt.cex=1.5,ncol=1,cex=0.75,bty="n",y.intersp=0.75,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=1)
dev.off()



## Combined GAM plots ------------------------------------------------------
reg.ext=extent(region.mask.sf)
pdat.sp.UTM.TP.sal.m=expand.grid(
  mean.sal=TP.sal.WY.m$var.summary$mean.sal[2],
  WY=TP.sal.WY.m$var.summary$WY[2],
  UTMX=seq(reg.ext[1],reg.ext[2],by=2000),
  UTMY=seq(reg.ext[3],reg.ext[4],by=2000)
)
UTM.pred=predict(TP.sal.WY.m,pdat.sp.UTM.TP.sal.m,type="terms",se.fit = F,nthreads=6,discrete=T)
colnames(UTM.pred)=paste0("fit.",dev.vars)
pdat.sp.UTM.TP.sal.m=cbind(pdat.sp.UTM.TP.sal.m,UTM.pred)

tmp.ma=with(pdat.sp.UTM.TP.sal.m,matrix(fit.UTM,nrow=length(unique(UTMX)),ncol=length(unique(UTMY))))
dat1=list()
dat1$x=unique(pdat.sp.UTM.TP.sal.m$UTMX)
dat1$y=unique(pdat.sp.UTM.TP.sal.m$UTMY)
dat1$z=round(tmp.ma,2)
r=raster(dat1)
TPSal.GAM.UTM <- mask(r, region.mask.sf)

reg.ext=extent(region.mask.sf)
pdat.sp.UTM.TN.sal.m=expand.grid(
  mean.sal=TN.sal.WY.m$var.summary$mean.sal[2],
  WY=TN.sal.WY.m$var.summary$WY[2],
  UTMX=seq(reg.ext[1],reg.ext[2],by=2000),
  UTMY=seq(reg.ext[3],reg.ext[4],by=2000)
)
UTM.pred=predict(TN.sal.WY.m,pdat.sp.UTM.TN.sal.m,type="terms",se.fit = F,nthreads=6,discrete=T)
colnames(UTM.pred)=paste0("fit.",dev.vars)
pdat.sp.UTM.TN.sal.m=cbind(pdat.sp.UTM.TN.sal.m,UTM.pred)

tmp.ma=with(pdat.sp.UTM.TN.sal.m,matrix(fit.UTM,nrow=length(unique(UTMX)),ncol=length(unique(UTMY))))
dat1=list()
dat1$x=unique(pdat.sp.UTM.TN.sal.m$UTMX)
dat1$y=unique(pdat.sp.UTM.TN.sal.m$UTMY)
dat1$z=round(tmp.ma,2)
r=raster(dat1)
TNSal.GAM.UTM <- mask(r, region.mask.sf)

b1=seq(-1.1,1.6,0.5)
b2=seq(-1.1,1.6,length.out=50)
pal2=viridis(length(b2)-1,alpha=0.75)# colorRampPalette(c("blue","grey90","red"))(length(b2)-1)


# png(filename=paste0(plot.path,"TPTN_sal_gam.png"),width=7.5,height=4,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(2,2.75,1,0.5),oma=c(1,1,0.25,0.25),xpd=F);
layout(cbind(matrix(c(1:6),2,3,byrow=T),7),widths=c(1,1,1.5,1))

ylim.val=c(-0.5,0.5);by.y=0.25;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(0,45);by.x=10;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
plot(fit.Sal~mean.sal,TP.sal.WY.m.sal.d,ylim=ylim.val,xlim=xlim.val,ann=F,axes=F,type="n")
abline(h=ymaj,v=xmaj,lty=3,col="grey")
abline(h=0)
with(TP.sal.WY.m.sal.d,shaded.range(mean.sal,lower.CI,upper.CI,"grey",lty=1))
lines(fit.Sal~mean.sal,TP.sal.WY.m.sal.d,lwd=2)
lines(dsig.incr ~ mean.sal, data = TP.sal.WY.m.sal.d, col = "red", lwd = 3,lty=1)
lines(dsig.decr ~ mean.sal, data = TP.sal.WY.m.sal.d, col = "blue", lwd = 3,lty=1)
axis_fun(1,xmaj,xmin,xmaj,line=-0.5);axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=3,adj=1,"s(Salinity)")
mtext(side=1,line=1.75,"Salinity (ppt)",cex=0.9)
mtext(side=2,line=2.5,"TP Effect")
mtext(side=3,line=-1.25,adj=1,"A ",font=2)

range(TP.sal.WY.m.WY.d$fit.WY)
range(TP.sal.WY.m.WY.d[c("upper.CI","lower.CI")])
ylim.val=c(-1,0.5);by.y=0.5;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(1990,2023);by.x=10;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/by.x)
plot(fit.WY~WY,TP.sal.WY.m.WY.d,ylim=ylim.val,xlim=xlim.val,ann=F,axes=F,type="n")
abline(h=ymaj,v=xmaj,lty=3,col="grey")
abline(h=0)
with(TP.sal.WY.m.WY.d,shaded.range(WY,lower.CI,upper.CI,"grey",lty=1))
lines(fit.WY~WY,TP.sal.WY.m.WY.d,lwd=2)
lines(dsig.incr ~ WY, data = TP.sal.WY.m.WY.d, col = "red", lwd = 3,lty=1)
lines(dsig.decr ~ WY, data = TP.sal.WY.m.WY.d, col = "blue", lwd = 3,lty=1)
axis_fun(1,xmaj,xmin,xmaj,line=-0.5);axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=3,adj=1,"s(WY)")
mtext(side=1,line=1.75,"Water Year",cex=0.9)
mtext(side=3,line=-1.25,adj=1,"B ",font=2)

par(mar=c(0,0.5,1,0.5))
bbox.lims=st_bbox(st_buffer(region.mask.sf,-1000))
plot(st_geometry(region.mask.sf),ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],border=NA)
image(TPSal.GAM.UTM,add=T,breaks=b2,col = pal2)
plot(rasterToContour(TPSal.GAM.UTM,levels=b1,nlevels=length(b1)),col=adjustcolor("black",0.5),lwd=1,add=T)
plot(st_geometry(region.mask.sf),add=T)
plot(st_geometry(subset(sal.sites.shp,N.sal>1)),add=T,pch=19,bg=NA,col=adjustcolor("grey2",0.5),cex=0.1)
plot(st_geometry(ENP.sf),add=T,bg=NA,lwd=0.5)
plot(st_geometry(region.shore),lty=1,add=T,border="white",lwd=0.5)
# mapmisc::scaleBar(crs=region.mask.sf,"bottomleft",bty="n",cex=1,seg.len=4,outer=F)
box(lwd=1)
mtext(side=3,adj=1,"s(UTMX,UTMY)")
mtext(side=3,line=-1.25,adj=1,"C ",font=2)
## TN 
par(mar=c(2,2.75,1,0.5))
ylim.val=c(-0.25,0.5);by.y=0.25;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(0,45);by.x=10;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
plot(fit.Sal~mean.sal,TN.sal.WY.m.sal.d,ylim=ylim.val,xlim=xlim.val,ann=F,axes=F,type="n")
abline(h=ymaj,v=xmaj,lty=3,col="grey")
abline(h=0)
with(TN.sal.WY.m.sal.d,shaded.range(mean.sal,lower.CI,upper.CI,"grey",lty=1))
lines(fit.Sal~mean.sal,TN.sal.WY.m.sal.d,lwd=2)
lines(dsig.incr ~ mean.sal, data = TN.sal.WY.m.sal.d, col = "red", lwd = 3,lty=1)
lines(dsig.decr ~ mean.sal, data = TN.sal.WY.m.sal.d, col = "blue", lwd = 3,lty=1)
axis_fun(1,xmaj,xmin,xmaj,line=-0.5);axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=1,line=1.75,"Salinity (ppt)",cex=0.9)
mtext(side=2,line=2.5,"TN Effect")
mtext(side=3,line=-1.25,adj=1,"D ",font=2)

range(TN.sal.WY.m.WY.d$fit.WY)
range(TN.sal.WY.m.WY.d[c("upper.CI","lower.CI")])
ylim.val=c(-1,1);by.y=0.5;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(1990,2023);by.x=10;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/by.x)
plot(fit.WY~WY,TN.sal.WY.m.WY.d,ylim=ylim.val,xlim=xlim.val,ann=F,axes=F,type="n")
abline(h=ymaj,v=xmaj,lty=3,col="grey")
abline(h=0)
with(TN.sal.WY.m.WY.d,shaded.range(WY,lower.CI,upper.CI,"grey",lty=1))
lines(fit.WY~WY,TN.sal.WY.m.WY.d,lwd=2)
lines(dsig.incr ~ WY, data = TN.sal.WY.m.WY.d, col = "red", lwd = 3,lty=1)
lines(dsig.decr ~ WY, data = TN.sal.WY.m.WY.d, col = "blue", lwd = 3,lty=1)
axis_fun(1,xmaj,xmin,xmaj,line=-0.5);axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
# mtext(side=3,adj=1,"s(WY)")
mtext(side=1,line=1.75,"Water Year",cex=0.9)
mtext(side=3,line=-1.25,adj=1,"E ",font=2)

par(mar=c(0,0.5,1,0.5))
plot(st_geometry(region.mask.sf),ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],border=NA)
image(TNSal.GAM.UTM,add=T,breaks=b2,col = pal2)
plot(rasterToContour(TNSal.GAM.UTM,levels=b1,nlevels=length(b1)),col=adjustcolor("black",0.5),lwd=1,add=T)
plot(st_geometry(region.mask.sf),add=T)
plot(st_geometry(subset(sal.sites.shp,N.sal>1)),add=T,pch=19,bg=NA,col=adjustcolor("grey2",0.5),cex=0.1)
plot(st_geometry(ENP.sf),add=T,bg=NA,lwd=0.5)
plot(st_geometry(region.shore),lty=1,add=T,border="white",lwd=0.5)
mapmisc::scaleBar(crs=region.mask.sf,"bottomright",bty="n",cex=1,seg.len=4,outer=F)
box(lwd=1)
mtext(side=3,line=-1.25,adj=1,"F ",font=2)

par(mar=c(0,0.5,1,0))
plot(0:1,0:1,ann=F,axes=F,type="n")
leg.fun(round(b2,1),pal2,leg.title="s(UTMX,UTMY)\nEffect",leg.type = "continuous",
        x.max=0.6,x.min=0.4)
legend(0.5,0.15,legend=c("Spatial Effect Contour","Monitoring Locations",
                         "Sig. Increase","Sig. Decrease"),
       pch=c(NA,19,NA,NA),lty=c(1,NA,1,1),lwd=c(1,0.01,2,2),
       col=c(adjustcolor("black",0.5),"grey2","red","blue"),pt.bg=c(NA,NA,NA,NA),
       pt.cex=1.5,ncol=1,cex=0.8,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=1)
dev.off()




## look into the modelsummary r-package
## Just testing it out
# library(modelsummary)
# # 
# modelsummary(TP.sal.WY.mi,output = "flextable")
# mods = list(
#     TP = TP.sal.WY.mi,
#     TN = TN.sal.WY.mi
# )
#   
# modelsummary(mods,output = "flextable")

# parametric
TN.p.tab=TN.sal.WY.m.sum$p.table
TN.s.tab=TN.sal.WY.m.sum$s.table

TN.p.tab=TN.p.tab|>
  as.data.frame()|>
  mutate(term=row.names(TN.p.tab))
# TN.p.tab$term="Intercept"
row.names(TN.p.tab)=1

TN.s.tab=TN.s.tab|>
  as.data.frame()|>
  mutate(term=row.names(TN.s.tab))
row.names(TN.s.tab)=1:nrow(TN.s.tab)

TN.sum.tab=data.frame(term=c(row.names(TN.sal.WY.m.sum$p.table),row.names(TN.sal.WY.m.sum$s.table)))|>
  merge(TN.p.tab,"term",all.x=T)|>
  merge(TN.s.tab,"term",all.x=T)|>
  mutate(R2=TN.sal.WY.m.sum$r.sq,
         dev.exp=TN.sal.WY.m.sum$dev.expl,
         sm.crit=as.numeric(TN.sal.WY.m.sum$sp.criterion),
         sm.method=names(TN.sal.WY.m.sum$sp.criterion),
         scale=TN.sal.WY.m.sum$scale,
         n.val=TN.sal.WY.m.sum$n,
         Pred.var="TN"
         )

TP.p.tab=TP.sal.WY.m.sum$p.table
TP.s.tab=TP.sal.WY.m.sum$s.table

TP.p.tab=TP.p.tab|>
  as.data.frame()|>
  mutate(term=row.names(TP.p.tab))
# TP.p.tab$term="Intercept"
row.names(TP.p.tab)=1

TP.s.tab=TP.s.tab|>
  as.data.frame()|>
  mutate(term=row.names(TP.s.tab))
row.names(TP.s.tab)=1:nrow(TP.s.tab)

TP.sum.tab=data.frame(term=c(row.names(TP.sal.WY.m.sum$p.table),row.names(TP.sal.WY.m.sum$s.table)))|>
  merge(TP.p.tab,"term",all.x=T)|>
  merge(TP.s.tab,"term",all.x=T)|>
  mutate(R2=TP.sal.WY.m.sum$r.sq,
         dev.exp=TP.sal.WY.m.sum$dev.expl,
         sm.crit=as.numeric(TP.sal.WY.m.sum$sp.criterion),
         sm.method=names(TP.sal.WY.m.sum$sp.criterion),
         scale=TP.sal.WY.m.sum$scale,
         n.val=TP.sal.WY.m.sum$n,
         Pred.var="TP"
  )

GAM.sum.tab=rbind(TP.sum.tab,TN.sum.tab)           
GAM.sum.tab=GAM.sum.tab[,c(16,1:15)]
GAM.sum.tab=GAM.sum.tab|>
  mutate(df.comb=trimws(ifelse(is.na(edf)==T,NA,paste0(format(round(edf,2)),"\n(",format(round(Ref.df,2),trim=T),")"))),
         r2.comb=paste0(round(R2,2),"\n(",round(dev.exp,2),")"),
         sm.comb=paste0(round(sm.crit,2),"\n(",sm.method,")"),
         scale.comb=paste0(round(scale,2),"\n(",n.val,")"))
GAM.sum.tab=rename(GAM.sum.tab,c("Pr(>|t|)"="p_tval","p-value"="p_value"))

library(flextable)
vars=c("Pred.var", "term", "Estimate", "Std. Error", "t value", "p_tval",
       "df.comb","F","p_value","r2.comb","sm.comb","scale.comb")
GAM.sum.tab[,vars]|>
  flextable()|>
  colformat_double(j=c(3:6,8:9),digits=2,na_str=" --- ")|>
  colformat_char(j=7,na_str=" --- ")|>
  merge_v(j=c(1,10:12))|>
  valign(j=c(1,10:12),valign="top")|>
  fix_border_issues()|>
  compose(j="Pred.var",i=~Pred.var=="TP",value=as_paragraph("(1) Annual\nMean\nTotal\nPhosphorus"))|>
  compose(j="Pred.var",i=~Pred.var=="TN",value=as_paragraph("(2) Annual\nMean\nTotal\nNitrogen"))|>
  compose(j="term",i=~term=="(Intercept)",value=as_paragraph("Intercept"))|>
  compose(j="term",i=~term=="s(mean.sal)",value=as_paragraph("s(salinity)"))|>
  italic(j="term",i=~term!="(Intercept)")|>
  italic(j="p_tval",i=~p_tval<0.05)|>
  italic(j="p_value",i=~p_value<0.05)|>
  # compose(j="p_tval",i=~p_tval<0.05,value=as_paragraph("< 0.05"))|>
  compose(j="p_tval",i=~p_tval<0.01,value=as_paragraph("< 0.01"))|>
  compose(j="p_value",i=~p_value<0.01,value=as_paragraph("< 0.01"))|>
  set_header_labels(
    "Pred.var"="Predicted\nVariable",
    "term" = "Term",
    "Std. Error" = "Standard\nError",
    "t value"="t-value",
    "p_tval" = "\u03C1-value",
    "df.comb" = "edf\n(Ref df)",
    "F" = "F-value",
    "p_value" = "\u03C1-value",
    "r2.comb" = "Adj R\u00B2\n(Dev. Exp.)",
    "sm.comb" = "Smooth\nSelection\nCriterion",
    "scale.comb" = "Scale\n(n)"
  )|>
  align(j=3:12,align="center",part="all")|>
  hline(i=c(5))|>
  padding(padding=0.75,part="all")|>
  font(fontname="Times New Roman",part="all")|>
  fontsize(size=10,part="all")|>
  bold(part="header")|>
  width(width=c(0.7,1.2,rep(0.6,10)))|>print("docx")




colformat_num(j=3,big.mark="")%>%
  colformat_double(j=4:7,i=1:10,digits=3,na_str = " ")%>%
  colformat_double(j=4:7,i=11:20,digits=2,na_str = " ")%>%
  colformat_double(j=4:7,i=21:30,digits=2,na_str = " ")%>%
  merge_v(j=1)%>%
  fix_border_issues()%>%
  valign(j=1,valign="top")%>%
  set_header_labels(
    "param.vals"="Parameter\n(Units)",
    "Region.txt"="Region",
    "N.val"="N",
    "min.val"="Minimum",
    "med.val"="Median",
    "mean.val"="Mean",
    "max.val"="Maximum"
  )%>%
  hline(i=c(5,10,15,20,25))%>%
  align(j=3:7,align="center",part="header")%>%
  padding(padding=1,part="all")%>%
  font(fontname="Times New Roman",part="all")%>%
  bold(part="header")%>%
  autofit()  #%>%print("docx")

# Extra Analyses ----------------------------------------------------------



## Change analysis ---------------------------------------------------------
library(trend)


tmp.dat1=subset(dat.all.GM,WY%in%seq(1996,2019,1))
sites.vals=unique(tmp.dat1$STATION)
TP.ptrslt=data.frame()
for(i in 1:length(sites.vals)){
  tmp.dat=subset(tmp.dat1,variable=="TP"&STATION==sites.vals[i])
  pt.rslt=pettitt.test(tmp.dat$GM)
  yr.val=tmp.dat$WY[as.numeric(pt.rslt$estimate)]
  rslt=data.frame(STATION=sites.vals[i],p.val=pt.rslt$p.value,stat=pt.rslt$statistic,Yr.chnge=yr.val)
  TP.ptrslt=rbind(TP.ptrslt,rslt)
  }
subset(TP.ptrslt,p.val<0.05)

# Quantmod data analysis (TN & TP) ----------------------------------------




# PCA analysis ------------------------------------------------------------
# dat.all2 #daily data
dat.all.GM

library(REdaS)
library(vegan)

idvars.val=c("STATION","Region.f","WY")
params=c("TOC",'Chla',"TN","TP","TN_TP")

dat.all.GM.xtab=dcast(subset(dat.all.GM,variable%in%params),STATION+Region.f+WY~variable,value.var = "GM",mean)

pca.dat=na.omit(dat.all.GM.xtab[,c(idvars.val,"TN","TP","TOC","TN_TP","Chla")])
pca.dat1=pca.dat[,c("TN","TP","TOC")]


KMOS(pca.dat1)
bart_spher(pca.dat1)

pca1=rda(pca.dat1,scale=T)

eig <- pca1$CA$eig
variance <- eig*100/sum(eig)
cumvar <- cumsum(variance)
eig.pca <- data.frame(eig = eig, variance = variance,cumvariance = cumvar)
eig.pca

# library(eigenprcomp)
# boot_pr_comp(as.matrix(eig.pca))

plot(pca1)
scrs=scores(pca1,display=c("sites","species"),choices=c(1,2,3));

pca.dat$PCA1=scrs$sites[,1]
reg.vals=levels.var

x.vals.sp=seq(1996,2019,0.1)
cols2=c("grey",cols[2:5])
plot(PCA1~WY,pca.dat)
for(i in 2:length(reg.vals)){
k1=loess(PCA1~WY,subset(pca.dat,Region.f==reg.vals[i]))
k1.pred=predict(k1,data.frame(WY=x.vals.sp))
lines(k1.pred~x.vals.sp,col=cols2[i],lwd=2)
}

plot(Chla~PCA1,pca.dat)
abline(0:1)
for(i in 2:length(reg.vals)){
  k1=lm(Chla~PCA1,subset(pca.dat,Region.f==reg.vals[i]))
  xvals=with(subset(pca.dat,Region.f==reg.vals[i]),seq(min(PCA1),max(PCA1),length.out=20))
  pred.k1=predict(k1,data.frame(PCA1=xvals))
  lines(pred.k1~xvals,col=cols2[i],lwd=4)
}



