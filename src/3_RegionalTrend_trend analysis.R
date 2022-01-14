
## Run 2_RegionalTrend_DataPrep.R prior to this file

# Trend -------------------------------------------------------------------
Ncheck=ddply(dat.all.GM,c("STATION","variable"),summarise,N.val=N.obs(WY))
subset(Ncheck,N.val<3)

ann.trend=ddply(dat.all.GM,c("STATION","variable"),summarise,
                   est=as.numeric(cor.test(WY,GM,method="kendall")$estimate),
                   pval=cor.test(WY,GM,method="kendall")$p.value,
                   sen.slope=as.numeric(zyp::zyp.sen(GM~WY)$coefficients[2]),
                   N.WY=N.obs(WY))
subset(ann.trend,pval<0.05)
ann.trend$stat.sig=with(ann.trend,ifelse(pval<0.05,"sig","not-sig"))
ann.trend$stat.sig=with(ann.trend,ifelse(is.na(stat.sig)==T,"Empty",as.character(stat.sig)))
ann.trend$stat.sig=as.factor(ann.trend$stat.sig)

AGM.all=ddply(dat.all.GM,c("STATION","variable"),summarise,
                 mean.GM=mean(GM,na.rm=T),
                 N.val=N.obs(GM),
                 SE.val=SE(GM),
                 var.val=var(GM,na.rm=T),# sample variance s^2; use this one (variance within site)
                 pop.var=var.val*((N.val-1)/N.val), # population variance sigma^2 
                 sd.val=sd(GM,na.rm=T),
                 CV.val=cv.per(GM)*100)

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


# Spline ------------------------------------------------------------------
#thin plate spline https://rspatial.org/raster/analysis/4-interpolation.html

region.buf.r=raster(region.mask)
res(region.buf.r)=1000

# Spatial Trend -----------------------------------------------------------
# TN
m=Tps(coordinates(sites.shp.TN.trend),sites.shp.TN.trend$sen.slope)
tps=interpolate(region.buf.r,m)
tps.TN.trend=mask(tps,region.mask)
plot(tps.TN.trend)

# DIN
m=Tps(coordinates(sites.shp.DIN.trend),sites.shp.DIN.trend$sen.slope)
tps=interpolate(region.buf.r,m)
tps.DIN.trend=mask(tps,region.mask)

# TP
m=Tps(coordinates(sites.shp.TP.trend),sites.shp.TP.trend$sen.slope)
tps=interpolate(region.buf.r,m)
tps.TP.trend=mask(tps,region.mask)

# SRP
m=Tps(coordinates(sites.shp.SRP.trend),sites.shp.SRP.trend$sen.slope)
tps=interpolate(region.buf.r,m)
tps.SRP.trend=mask(tps,region.mask)

# Chla
m=Tps(coordinates(sites.shp.Chla.trend),sites.shp.Chla.trend$sen.slope)
tps=interpolate(region.buf.r,m)
tps.Chla.trend=mask(tps,region.mask)

# TOC
m=Tps(coordinates(sites.shp.TOC.trend),sites.shp.TOC.trend$sen.slope)
tps=interpolate(region.buf.r,m)
tps.TOC.trend=mask(tps,region.mask)

# NP
m=Tps(coordinates(sites.shp.NP.trend),sites.shp.NP.trend$sen.slope)
tps=interpolate(region.buf.r,m)
tps.NP.trend=mask(tps,region.mask)


# Spatial AGM -------------------------------------------------------------

# TN
sites.shp.TN.GM=subset(sites.shp.TN.GM,is.na(mean.GM)==F)
m.GM=Tps(coordinates(sites.shp.TN.GM),sites.shp.TN.GM$mean.GM)
tps.GM=interpolate(region.buf.r,m.GM)
tps.TN.GM=mask(tps.GM,region.mask)
plot(tps.TN.GM)

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

m.GM.var=Tps(coordinates(sites.shp.TP.GM),sites.shp.TP.GM$var.val)
tps.GM.var=interpolate(region.buf.r,m.GM.var)
tps.TP.GM.var=mask(tps.GM.var,region.mask)
plot(tps.TP.GM.var)

m.GM.CV=Tps(coordinates(sites.shp.TP.GM),sites.shp.TP.GM$CV.val)
tps.GM.CV=interpolate(region.buf.r,m.GM.CV)
tps.TP.GM.CV=mask(tps.GM.CV,region.mask)
plot(tps.TP.GM.CV)

# SRP
sites.shp.SRP.GM=subset(sites.shp.SRP.GM,is.na(mean.GM)==F)

m.GM=Tps(coordinates(sites.shp.SRP.GM),sites.shp.SRP.GM$mean.GM)
tps.GM=interpolate(region.buf.r,m.GM)
tps.SRP.GM=mask(tps.GM,region.mask)

m.GM.var=Tps(coordinates(sites.shp.SRP.GM),sites.shp.SRP.GM$var.val)
tps.GM.var=interpolate(region.buf.r,m.GM.var)
tps.SRP.GM.var=mask(tps.GM.var,region.mask)
plot(tps.SRP.GM.var)

m.GM.CV=Tps(coordinates(sites.shp.SRP.GM),sites.shp.SRP.GM$CV.val)
tps.GM.CV=interpolate(region.buf.r,m.GM.CV)
tps.SRP.GM.CV=mask(tps.GM.CV,region.mask)
plot(tps.SRP.GM.CV)

# Chla
sites.shp.Chla.GM=subset(sites.shp.Chla.GM,is.na(mean.GM)==F)

m.GM=Tps(coordinates(sites.shp.Chla.GM),sites.shp.Chla.GM$mean.GM)
tps.GM=interpolate(region.buf.r,m.GM)
tps.Chla.GM=mask(tps.GM,region.mask)

m.GM.var=Tps(coordinates(sites.shp.Chla.GM),sites.shp.Chla.GM$var.val)
tps.GM.var=interpolate(region.buf.r,m.GM.var)
tps.Chla.GM.var=mask(tps.GM.var,region.mask)
plot(tps.Chla.GM.var)

m.GM.CV=Tps(coordinates(sites.shp.Chla.GM),sites.shp.Chla.GM$CV.val)
tps.GM.CV=interpolate(region.buf.r,m.GM.CV)
tps.Chla.GM.CV=mask(tps.GM.CV,region.mask)
plot(tps.Chla.GM.CV)

# TOC
sites.shp.TOC.GM=subset(sites.shp.TOC.GM,is.na(mean.GM)==F)

m.GM=Tps(coordinates(sites.shp.TOC.GM),sites.shp.TOC.GM$mean.GM)
tps.GM=interpolate(region.buf.r,m.GM)
tps.TOC.GM=mask(tps.GM,region.mask)
plot(tps.TOC.GM)

m.GM.var=Tps(coordinates(sites.shp.TOC.GM),sites.shp.TOC.GM$var.val)
tps.GM.var=interpolate(region.buf.r,m.GM.var)
tps.TOC.GM.var=mask(tps.GM.var,region.mask)
plot(tps.TOC.GM.var)

m.GM.CV=Tps(coordinates(sites.shp.TOC.GM),sites.shp.TOC.GM$CV.val)
tps.GM.CV=interpolate(region.buf.r,m.GM.CV)
tps.TOC.GM.CV=mask(tps.GM.CV,region.mask)
plot(tps.TOC.GM.CV)

# NP
sites.shp.NP.GM=subset(sites.shp.NP.GM,is.na(mean.GM)==F)

m.GM=Tps(coordinates(sites.shp.NP.GM),sites.shp.NP.GM$mean.GM)
tps.GM=interpolate(region.buf.r,m.GM)
tps.NP.GM=mask(tps.GM,region.mask)
plot(tps.NP.GM)

m.GM.var=Tps(coordinates(sites.shp.NP.GM),sites.shp.NP.GM$var.val)
tps.GM.var=interpolate(region.buf.r,m.GM.var)
tps.NP.GM.var=mask(tps.GM.var,region.mask)
plot(tps.NP.GM.var)

m.GM.CV=Tps(coordinates(sites.shp.NP.GM),sites.shp.NP.GM$CV.val)
tps.GM.CV=interpolate(region.buf.r,m.GM.CV)
tps.NP.GM.CV=mask(tps.GM.CV,region.mask)
plot(tps.NP.GM.CV)
