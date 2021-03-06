
## Run 1_RegionalTrend_libs_GIS.R first

# Data --------------------------------------------------------------------
WYs=seq(1996,2019,1)

#from Caccia and Boyer & Rookery Bay SERC report
mdls=data.frame(param=c("TN","TP","TOC","DOC","NOX","NO2","NH4","SRP","CHLA"),
                mdl=c(0.05,0.0012,0.04,0.06,0.0024,0.0003,0.0057,0.0022,0.1))
mdls$mdl.uM=with(mdls,ifelse(param%in%c("TN","NOX","NO2","NH4"),(mdl*1000)/N.mw,
                             ifelse(param%in%c("TP","SRP"),(mdl*1000)/P.mw,
                                    ifelse(param%in%c("TOC","DOC"),(mdl*1000)/C.mw,mdl))))


# FIU SERC ----------------------------------------------------------------
# http://serc.fiu.edu/wqmnetwork/
serc=read.xlsx(paste0(data.path,"FIU_SERC/WQFloridaKeys&Shelf (ppm) UPDATED 6-6-2020.xlsx"),sheet=3)
serc$DATE=date.fun(convertToDate(serc$DATE))
serc=subset(serc,DATE>date.fun("1995-01-01"))
# serc=rename(serc,c("TN-S"="TN"))

serc.sites=ddply(serc,c("STATION"),summarise,LATDEC=mean(LATDEC,na.rm=T),LONDEC=mean(LONDEC,na.rm=T),n.val=N.obs(SITE),min.date=min(DATE,na.rm=T),max.date=max(DATE,na.rm=T))
serc.sites=subset(serc.sites,is.na(LONDEC)==F)
serc.sites.shp=spTransform(SpatialPointsDataFrame(serc.sites[,c("LONDEC","LATDEC")],data=serc.sites,proj4string=wgs84),utm17)

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
serc2=serc[,c("DATE","STATION","TN-S","NOX-S","NH4-S","TP-S","SRP-S","CHLA-S","TOC-S")]
serc2=rename(serc2,c("TN-S"="TN","DIN-S"="DIN","NOX-S"="NOx","NH4-S"="NH4","TP-S"="TP","SRP-S"="SRP","CHLA-S"="Chla","TOC-S"="TOC"))
serc2[,c("TN","NOx","NH4","TP","SRP","Chla","TOC")]=sapply(serc2[,c("TN","NOx","NH4","TP","SRP","Chla","TOC")],as.numeric)
serc2$TN=with(serc2,ifelse(TN<0.05,0.05/2,TN))
serc2$NOx=with(serc2,ifelse(NOx<0.0024,0.0024/2,NOx))
serc2$NH4=with(serc2,ifelse(NH4<0.0057,0.0057/2,NH4))
serc2$DIN=with(serc2,ifelse(is.na(NH4)==T|is.na(NOx)==T,NA,NH4+NOx))
serc2$TP=with(serc2,ifelse(TP<0.0012,0.0012/2,TP))
serc2$SRP=with(serc2,ifelse(SRP<0.0022,0.0022/2,SRP))
serc2$TOC=with(serc2,ifelse(TOC< 0.04, 0.04/2,TOC))
serc2$DOC=with(serc2,ifelse(DOC< 0.06, 0.06/2,DOC))
serc2=subset(serc2,is.na(TN)==F)
serc2$WY=WY(serc2$DATE)
serc2$season=FL.Hydroseason(serc2$DATE)

chk=ddply(serc2,c("STATION"),summarise,n.val=N.obs(TN),min.date=min(DATE,na.rm=T),max.date=max(DATE,na.rm=T))
chk[chk$max.date>date.fun("2011-07-30"),"STATION"]

head(serc2)
serc2=ddply(serc2,c("STATION","DATE","WY","season"),summarise,
            TN=mean(TN,na.rm=T),DIN=mean(DIN,na.rm=T),
            TP=mean(TP,na.rm=T),SRP=mean(SRP,na.rm=T),
            Chla=mean(Chla,na.rm=T),TOC=mean(TOC,na.rm=T))
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
fce.boyer.sites=spTransform(SpatialPointsDataFrame(fce.boyer.sites[,c("LONDEC","LATDEC")],data=fce.boyer.sites,proj4string=wgs84),utm17)
# only TS/PH9,10,11
# 
# tm_shape(serc.sites.shp)+tm_dots()+
#   tm_shape(ENP_FLB)+tm_dots(col="red",alpha=0.5)+
#   tm_shape(lter)+tm_dots(col="yellow",alpha=0.5)+
#   tm_shape(fce.boyer.sites)+tm_dots(col="green",alpha=0.5)

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
fce.wq$TOC=(fce.wq$TOC/1000)*C.mw
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
LTER.sites.region=rbind(subset(LTER.sites.region,Region=="Coastal_Mangroves"),
                        data.frame(SITE=c(paste("SRS",c("1a","1b","1c","1d",2:3),sep="-"),paste("TS/Ph",c("1a","1b","2","2a","2b",4,3,5,"6a","6b","7a","7b",8:11),sep="-")),
                                   ESTUARY=c(rep("SRS",6),rep("TS",13),rep("FLBay",3)),
                                   SEGMENT_NA=c(rep("SRS_FW",5),rep("SRS_est",1),rep("TS_FW",6),rep("TS_est",7),rep("FLBay",3)),
                                   Region=c(rep("ENP",19),rep("FLBay",3)),
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

# params=data.frame(Test.Number=c(21,20,18,80,25,23,61,179,112,178,100),
#                   param=c("TKN","NH4","NOx","TN","TP","SRP","Chla","Chla","Chla","Chla","TOC"))
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

# use max?
ddply(wmd.dat,"param",summarise,min.val=min(MDL,na.rm=T),max.val=max(MDL,na.rm=T))
ddply(wmd.dat,"param",summarise,
      min.val=min(abs(Value[Value<0]),na.rm=T),
      max.val=max(abs(Value[Value<0]),na.rm=T))
mdls

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
wmd.dat.xtab$NH4=with(wmd.dat.xtab,ifelse(NH4>100,NA,NH4))
wmd.dat.xtab$DIN=with(wmd.dat.xtab,ifelse(DIN>10,NA,TN))
wmd.dat.xtab$Chla=with(wmd.dat.xtab,ifelse(Chla>100,NA,Chla))
wmd.dat.xtab$source="DBHydro"

subset(wmd.dat.xtab,Station.ID=="FLAB35"&Date.EST==date.fun("2009-10-13"))

plot(TN~Date.EST,wmd.dat.xtab,type="n")
with(subset(wmd.dat.xtab,Station.ID=="S332"),points(Date.EST,TN,pch=21,bg="red"))
with(subset(wmd.dat.xtab,Station.ID=="S332DX"),points(Date.EST,TN,pch=21,bg="green"))
with(subset(wmd.dat.xtab,Station.ID=="S332D"),points(Date.EST,TN,pch=21,bg="blue"))

wmd.dat.melt=melt(wmd.dat.xtab,id.vars = c("Station.ID",'STATION',"Date.EST","DATE","season","WY","source"))
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

subset(wmd.dat2,Turb<=0)
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
wmd.dat2$FLWY=WY(wmd.dat2$Date.EST)
wmd.dat2$source="CWQMN"
wmd.dat2$Sal=with(wmd.dat2,SalinityCalc(spc_uScm,Temp))

wmd.dat2$season=FL.Hydroseason(wmd.dat2$Date.EST)
wmd.dat2$STATION=wmd.dat2$Station.ID

wmd.dat2=wmd.dat2[,names(wmd.dat.xtab)]
head(wmd.dat2)
head(wmd.dat.xtab)

wmd.dat2.melt=melt(wmd.dat2,id.vars = c("Station.ID",'STATION',"Date.EST","DATE","season","WY","source"))
subset(wmd.dat2.melt,value==0)
# unique(subset(wmd.dat2.melt,value==0)$variable)
# wmd.dat2.melt$value=with(wmd.dat2.melt,ifelse(value==0,NA,value))
wmd.dat2.melt=subset(wmd.dat2.melt,is.na(value)==F)

dat.comb=rbind(wmd.dat2.melt,wmd.dat.melt)
dat.comb=reshape2::dcast(dat.comb,Station.ID+STATION+Date.EST+DATE+season+WY+variable~source,value.var = "value",mean)
dat.comb$diff_val=with(dat.comb,CWQMN-DBHydro)

ddply(dat.comb,"variable",summarise,min.dff=min(diff_val,na.rm=T),max.diff=max(diff_val,na.rm=T))

range(dat.comb$diff_val,na.rm=T)
boxplot(dat.comb$diff_val)
subset(dat.comb,diff_val>20)

subset(wmd.dat2,Station.ID=="FLAB35"&Date.EST==date.fun("2009-10-13"))
subset(wmd.dat,Station.ID=='FLAB35'&Date.EST==date.fun("2009-10-13"))
subset(dat.comb,Station.ID=="FLAB35"&Date.EST==date.fun("2009-10-13"))
subset(dat.comb,Station.ID=="FLAB13"&variable=="Chla")


plot(CWQMN~Date.EST,subset(dat.comb,Station.ID=="FLAB13"&variable=="Chla"))
points(DBHydro~Date.EST,subset(dat.comb,Station.ID=="FLAB13"&variable=="Chla"),
       pch=21,bg="blue")

head(dat.comb)
dat.comb$value.f=with(dat.comb,ifelse(is.na(DBHydro)==T,CWQMN,DBHydro))

plot(value.f~Date.EST,subset(dat.comb,Station.ID=="FLAB13"&variable=="Chla"))
points(DBHydro~Date.EST,subset(dat.comb,Station.ID=="FLAB13"&variable=="Chla"),
       pch=21,bg="blue")

wmd.dat.xtab.combo=reshape2::dcast(dat.comb,STATION+DATE+WY+season~variable,value.var = "value.f",mean,na.rm=T)
rm(dat.comb)
wmd.sites.region=spatialEco::point.in.poly(ENP_FLB,regions)@data[,c("STATION","ESTUARY","SEGMENT_NA","Region")]
wmd.sites.region=subset(wmd.sites.region,is.na(Region)==F)
wmd.sites.region$source="WMD"
wmd.sites.region=rbind(wmd.sites.region,
                       data.frame(STATION=c(paste0("S12",LETTERS[1:4]),"S333","SRS1B","SRS1C","NP201","NE1","SRS2",paste0("P",33:36),"G-3273","RG1","CR2","S332D","S332","TSB","S177","S197","S18C","EP","P37"),
                                  ESTUARY=c(rep("SRS",17),rep("TS",8)),
                                  SEGMENT_NA=c(rep("SRS_FW",17),rep("TS_FW",8)),
                                  Region=rep("ENP",25),
                                  source="WMD"))

nrow(wmd.dat.xtab.combo)
unique(wmd.dat.xtab.combo$STATION)
unique(wmd.sites.region$STATION)
wmd.dat.xtab.combo=subset(wmd.dat.xtab.combo,STATION%in%wmd.sites.region$STATION)
# NOAA --------------------------------------------------------------------
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

ddply(subset(noaa.dat,Station%in%c("21.5/LK","21.5")),c("Station"),summarise,
      min.lat=min(Latitude),max.lat=max(Latitude),mean.lat=mean(Latitude))
# View(subset(noaa.dat,Station%in%c("21.5/LK","21.5")))
noaa.dat[noaa.dat$Station=="21.5/LK","Station"]=c("21.5")

ddply(subset(noaa.dat,Station%in%c("21/LK","21")),c("Station"),summarise,
      min.lat=min(Latitude),max.lat=max(Latitude),mean.lat=mean(Latitude))
noaa.dat[noaa.dat$Station=="21/LK","Station"]=c("21")

# test=unique(noaa.dat$Station)
# sta.l=nchar(test)
# site.end=substr(test,ifelse(sta.l==1,1,sta.l),sta.l)
# test[site.end%in%LETTERS[1:3]==T]

noaa.dat$STATION=paste("NOAA",noaa.dat$Station,sep="_")

# noaa.dat$Depth=as.numeric(noaa.dat$Depth)

noaa.vars=c("Record.#", "Longitude", "Latitude", "Station", "Depth", "Date", 
            "GMT", "SST", "SSS", "Chla", "Phaeo", "NH4.uM", "SRP.uM", "NOx.uM", 
            "NO2.uM", "NO3.uM", "Si.uM",'STATION')
colnames(noaa.dat)=noaa.vars

noaa.dat[,c("Chla", "Phaeo", "NH4.uM")]=sapply(noaa.dat[,c("Chla", "Phaeo", "NH4.uM")],as.numeric)

vars1=c("STATION","Date")
params.vars=c("Chla", "NH4.uM", "SRP.uM", "NOx.uM","NO2.uM", "NO3.uM")
noaa.dat.melt=melt(noaa.dat[,c(vars1,params.vars)],id.vars = vars1)

noaa.dat.melt=merge(noaa.dat.melt,noaa.mdls,"variable",all.x=T)
noaa.dat.melt
noaa.dat.melt$HalfMDL.uM=with(noaa.dat.melt,ifelse(value<MDL,MDL/2,value))

noaa.dat.xtab=dcast(noaa.dat.melt,STATION+Date~variable,value.var = "HalfMDL.uM",mean)
noaa.dat.xtab$NOx.uM.calc=with(noaa.dat.xtab,NO2.uM+NO3.uM)
plot(NOx.uM~NOx.uM.calc,noaa.dat.xtab);abline(0,1)

noaa.dat.xtab$NH4=(noaa.dat.xtab$NH4.uM/1000)*N.mw
noaa.dat.xtab$SRP=(noaa.dat.xtab$SRP.uM/1000)*P.mw
noaa.dat.xtab$NOx=(noaa.dat.xtab$NOx.uM/1000)*N.mw
noaa.dat.xtab$DIN=with(noaa.dat.xtab, NH4+NOx)

noaa.dat.xtab$DATE=noaa.dat.xtab$Date
noaa.dat.xtab$season=FL.Hydroseason(noaa.dat.xtab$Date)
noaa.dat.xtab$WY=WY(noaa.dat.xtab$Date)

noaa.dat.xtab[,c("TKN","TN","TOC","TP")]=NA
noaa.dat.xtab[,vars]

names(wmd.dat.xtab.combo)
names(noaa.dat.xtab)

## Locations
noaa.sites=ddply(noaa.dat,c('STATION'),summarise,
                 LAT2=mean(Latitude),LONG2=mean(Longitude),N.val=N.obs(STATION))
unique(noaa.sites$STATION)
noaa.sites.shp=SpatialPointsDataFrame(noaa.sites[,c("LONG2","LAT2")],
                                      data=noaa.sites,proj4string = wgs84)
noaa.sites.shp=spTransform(noaa.sites.shp,utm17)

noaa.sites.region=spatialEco::point.in.poly(noaa.sites.shp,regions)@data[,c("STATION","ESTUARY","SEGMENT_NA","Region")]
noaa.sites.region=subset(noaa.sites.region,is.na(Region)==F)
noaa.sites.region$source="NOAA"
noaa.sites.region

noaa.sites.shp=subset(noaa.sites.shp,STATION%in%noaa.sites.region$STATION)
# write.csv(noaa.sites.region,paste0(export.path,"noaa_sites.csv"),row.names = F)

noaa.dat.xtab=subset(noaa.dat.xtab,STATION%in%noaa.sites.region$STATION)
# -------------------------------------------------------------------------
# Combine datasets
names(serc2)
names(fce.wq)
# names(wmd.dat.xtab)
names(wmd.dat.xtab.combo)
names(noaa.dat.xtab)

vars=c("STATION","DATE","WY","season","TN","DIN","TP","SRP","Chla","TOC")
dat.all=rbind(serc2[,vars],fce.wq[,vars],wmd.dat.xtab.combo[,vars],noaa.dat.xtab[,vars])
# write.csv(dat.all,paste0(export.path,"SERC_FCE_WMD_data.csv"),row.names = F)
rm(serc2,fce.wq,wmd.dat.xtab.combo,noaa.dat.xtab)

all.sites.region=rbind(serc.sites.region,LTER.sites.region,wmd.sites.region,noaa.sites.region)
rm(serc.sites.region,LTER.sites.region,wmd.sites.region,noaa.sites.region)
subset(all.sites.region,STATION=="S12C")
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
# Merge site BGChem data with region delineation
# Attribute all SRS2 data to FCE site
all.sites.region=subset(all.sites.region,!(STATION=="SRS2"&source=="WMD"))
dat.all2=merge(dat.all,subset(all.sites.region,!(STATION=="SRS2"&source=="WMD")),"STATION",all.x=T)
nrow(dat.all)
nrow(dat.all2)
head(dat.all2)

dat.all2$TP.mM=with(dat.all2,TP/P.mw)
dat.all2$TN.mM=with(dat.all2,TN/N.mw)
dat.all2$SRP.mM=with(dat.all2,SRP/P.mw)
dat.all2$DIN.mM=with(dat.all2,DIN/N.mw)

dat.all2$TN_TP=with(dat.all2,TN.mM/TP.mM)
dat.all2$DIN_SRP=with(dat.all2,DIN.mM/SRP.mM)

with(dat.all2,DIN.mM/TP.mM)
# write.csv(dat.all2,paste0(export.path,"20211117_alldata_daily.csv"),row.names=F)

### Season Screen
idvars=c("STATION", "DATE", "WY", "season","ESTUARY", "SEGMENT_NA", "Region", "source")
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

#double check Everglades TOC
subset(consec.WY,variable=="TOC"&Region=="ENP")


vars_join1=c("STATION","WY","ESTUARY","Region","source","variable","sea.screen")
vars_by1=c("STATION","WY","ESTUARY","Region","source","variable")
vars_join2=c("STATION","ESTUARY","Region","source","variable","consec.screen")
vars_by2=c("STATION","ESTUARY","Region","source","variable")
dat.all2.melt=merge(dat.all2.melt,samp.size[,vars_join1],vars_by1)
dat.all2.melt=merge(dat.all2.melt,consec.WY[,vars_join2],vars_by2)

# change TP and SRP to ug/L
dat.all2.melt$value=with(dat.all2.melt,ifelse(variable%in%c("TP","SRP"),value*1000,value))

dat.all.GM=ddply(subset(dat.all2.melt,sea.screen==1&consec.screen==1),
                 c("STATION","WY","variable"),
                 summarise,
                 GM=exp(mean(log(value),na.rm=T)),N.val=N.obs(value))

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
                   variable%in%c("TN", "DIN", "TP", "SRP", "Chla","TOC"))
dat.all.GM2.xtab=reshape2::dcast(subset(dat.all.GM2,WY%in%WYs),STATION+WY~variable,value.var = "GM",mean)
# write.csv(dat.all.GM2.xtab,paste0(export.path,"WYGeomean.csv"),row.names=F)

dat.all.GM2.xtab.info=merge(dat.all.GM2.xtab,all.sites.region[,c("STATION","Region","source")],all.x=T)
nrow(dat.all.GM2.xtab)
nrow(dat.all.GM2.xtab.info)
summary(dat.all.GM2.xtab.info)
# write.csv(dat.all.GM2.xtab.info,paste0(export.path,"20220521_WYGeomean.csv"),row.names=F)

unique(dat.all.GM2.xtab.info$source)
