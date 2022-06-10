

mdls$mdl.uM=with(mdls,ifelse(param%in%c("TN","NOX","NO2","NH4"),(mdl*1000)/N.mw,
                             ifelse(param%in%c("TP","SRP"),(mdl*1000)/P.mw,
                                    ifelse(param%in%c("TOC","DOC"),(mdl*1000)/C.mw,mdl))))

noaa.mdls=data.frame(variable=c("Chla", "NH4.uM","NO2.uM", "NO3.uM", "NOx.uM", "SRP.uM","Si.uM"),
           MDL=c(0.10,0.028,0.010,0.007,0.015,0.008,0.014))

noaa.dat=read.xlsx(paste0(data.path,"NOAA/AOML_SFP_regional_WQ_surface_v14.xlsx"))
noaa.dat$Date=date.fun(convertToDate(noaa.dat$Date))
# noaa.dat$Depth=as.numeric(noaa.dat$Depth)

names(noaa.dat)
noaa.vars=c("Record.#", "Longitude", "Latitude", "Station", "Depth", "Date", 
  "GMT", "SST", "SSS", "Chla", "Phaeo", "NH4.uM", "SRP.uM", "NOx.uM", 
  "NO2.uM", "NO3.uM", "Si.uM")
colnames(noaa.dat)=noaa.vars
noaa.dat$WY=WY(noaa.dat$Date)

noaa.dat[,c("Chla", "Phaeo", "NH4.uM")]=sapply(noaa.dat[,c("Chla", "Phaeo", "NH4.uM")],as.numeric)

vars1=c("Station","Date")
params.vars=c("Chla", "NH4.uM", "SRP.uM", "NOx.uM","NO2.uM", "NO3.uM")
noaa.dat.melt=melt(noaa.dat[,c(vars1,params.vars)],id.vars = vars1)

noaa.dat.melt=merge(noaa.dat.melt,noaa.mdls,"variable",all.x=T)
noaa.dat.melt
noaa.dat.melt$HalfMDL.uM=with(noaa.dat.melt,ifelse(value<MDL,MDL/2,value))

noaa.dat.xtab=dcast(noaa.dat.melt,Station+Date~variable,value.var = "HalfMDL.uM",mean)
noaa.dat.xtab$NOx.uM.calc=with(noaa.dat.xtab,NO2.uM+NO3.uM)
plot(NOx.uM~NOx.uM.calc,noaa.dat.xtab);abline(0,1)

noaa.dat.xtab$NH4=(noaa.dat.xtab$NH4.uM/1000)*N.mw
noaa.dat.xtab$SRP=(noaa.dat.xtab$SRP.uM/1000)*P.mw
noaa.dat.xtab$NOx=(noaa.dat.xtab$NOx.uM/1000)*N.mw
noaa.dat.xtab$DIN=with(noaa.dat.xtab, NH4+NOx)




## Locations
noaa.dat$LONG2=round(noaa.dat$Longitude,3)
noaa.dat$LAT2=round(noaa.dat$Latitude,3)
noaa.sites=ddply(noaa.dat,c('LONG2','LAT2','Station'),summarise,
                 N.val=N.obs(Date))
noaa.sites=ddply(noaa.dat,c('Station'),summarise,
                 LAT2=mean(Latitude),LONG2=mean(Longitude),N.val=N.obs(Station))
subset(noaa.sites,N.val>2)
subset(noaa.sites,Station==1)

noaa.sites.shp=SpatialPointsDataFrame(noaa.sites[,c("LONG2","LAT2")],
                                      data=noaa.sites,proj4string = wgs84)
noaa.sites.shp=spTransform(noaa.sites.shp,utm17)

plot(sites.shp2,pch=21,bg="grey60")
plot(noaa.sites.shp,add=T,pch=21,bg="red")
plot(regions2,add=T,border="grey40")

noaa.sites.shp.region=spatialEco::point.in.poly(noaa.sites.shp,regions)@data
subset(noaa.sites.shp.region,is.na(Region)==F)
