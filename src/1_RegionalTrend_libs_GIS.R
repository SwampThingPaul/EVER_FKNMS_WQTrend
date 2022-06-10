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

# utm17=wkt(CRS(SRS_string="EPSG:26917"))
# wgs84=wkt(CRS(SRS_string="EPSG:4326"))
# https://gis.stackexchange.com/questions/387072/r-spcrs-returns-na

# utm17=sf::st_crs(26917)[[2]]
# utm17=sp::CRS(utm17)
# 
# wgs84=sf::st_crs(4326)[[2]]
# wgs84=sp::CRS(wgs84)

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


# for metadata
bbox(spTransform(regions2,wgs84))
