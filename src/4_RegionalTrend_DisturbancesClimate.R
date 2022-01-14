
dev.off()
# hurricane analysis ------------------------------------------------------
library(HURDAT)

EVER_FKNMS=bind(gSimplify(subset(regions2,Region=="Coastal_Mangroves"),500),
                gSimplify(subset(regions2,Region=="FLBay"),500),
                gSimplify(subset(regions2,Region=="Shelf"),500),
                gSimplify(subset(regions2,Region=="Keys"),500),ENP
)
plot(EVER_FKNMS)

hur.dat=get_hurdat(basin = c("AL"))
range(hur.dat$DateTime)
hur.dat$date=date.fun(hur.dat$DateTime,tz="UTC")
hur.dat$Year=as.numeric(format(hur.dat$date,"%Y"))
hur.dat$WY=WY(hur.dat$DateTime)
hur.year=ddply(hur.dat,c("Key","Name"),summarise,
               date.start=min(date,na.rm=T),
               Year.start=min(Year,na.rm=T),
               WY.start=min(WY,na.rm=T),
               max.wind.ms=max(Wind*0.44704,na.rm=T),
               min.pres=min(Pressure,na.rm=T))
hur.dat.sp.pt=SpatialPointsDataFrame(coords=hur.dat[,c("Lon","Lat")],data=hur.dat,proj4string = CRS("+init=epsg:4269"))

#Convert point to line data
hur.dat2=hur.dat
hur.id=ddply(hur.dat,c("Key","Name"),summarise,N.val=length(which(Key!="NA")))
hur.dat2$dLat=with(hur.dat2,ave(Lat,Key,FUN=function(x)c(0,diff(x))))
hur.dat2$dLon=with(hur.dat2,ave(Lon,Key,FUN=function(x)c(0,diff(x))))
hur.dat2$dist.m=with(hur.dat2,sqrt((dLat^2)+(dLon^2)));#calculates distance between points

hur.id=ddply(hur.dat2,c("Key","Name"),summarise,N.val=N.obs(Key),max.dist=max(dist.m,na.rm=T))
hur.id=subset(hur.id,is.na(Key)==F)
hur.dat2=subset(hur.dat2,Key%in%subset(hur.id,N.val>1&max.dist>0)$Key);#need greater than one point and some distance between points to create a line
coordinates(hur.dat2)=c("Lon","Lat")
hur.dat2=hur.dat2[order(hur.dat2$Key,hur.dat2$DateTime),]
hur.dat2$WY=WY(hur.dat2$DateTime)
path=sp::split(hur.dat2,hur.dat2$Key)

##Convert spatial points to spatial lines for each hurricane
sp_lines=SpatialLinesDataFrame(SpatialLines(list(Lines(list(Line(path[[1]])),unique(path[[1]]@data$Key))),CRS("+init=epsg:4269")),data.frame(row.names=hur.id$Key,Key=hur.id$Key,Name=hur.id$Name))
pb=txtProgressBar(1,max=length(path),style=3)
for(i in 2:length(path)){
  tmp=SpatialLinesDataFrame(SpatialLines(list(Lines(list(Line(path[[i]])),unique(path[[i]]@data$Key))),CRS("+init=epsg:4269")),data.frame(row.names=hur.id$Key,Key=hur.id$Key,Name=hur.id$Name))
  #tmp=SpatialLines(list(Lines(list(Line(path[[i]])),unique(path[[i]]@data$Key))),CRS("+init=epsg:4269"))
  sp_lines=rbind(sp_lines,tmp)
  setTxtProgressBar(pb,i)
}
chk=data.frame(gIsValid(sp_lines,byid=T));#checks
chk$Key=rownames(chk)
colnames(chk)=c("geo.valid","Key")
subset(chk,geo.valid=="FALSE")

hur.tck=spTransform(sp_lines,utm17)
hur.tck=merge(hur.tck,hur.year,by.x=c("Key","Name"),by.y=c("Key","Name"))


EVER_FKNMS.buffer1=gBuffer(EVER_FKNMS,width=75*1000);
plot(EVER_FKNMS.buffer1,add=T)
EVER_FKNMS.buffer1=SpatialPolygonsDataFrame(EVER_FKNMS.buffer1,data.frame(row.names = "buffer",width.km=75,area="LOK"))
EVER_FKNMS.buffer1=spTransform(EVER_FKNMS.buffer1,utm17)

EVER_FKNMS.hurr=over(EVER_FKNMS.buffer1,hur.tck,returnList = T,byid=T)[[1]];#select only hurricanes that cross the 200 km buffer.
EVER_FKNMS.hurr=EVER_FKNMS.hurr[order(EVER_FKNMS.hurr$WY.start,EVER_FKNMS.hurr$Key),]
subset(EVER_FKNMS.hurr,Year.start%in%seq(1995,2019,1))

hurr.yr.N=ddply(subset(EVER_FKNMS.hurr,Year.start%in%seq(1995,2019,1)),"Year.start",summarise,
                N.val=N.obs(Year.start))

EVER_FKNMS.hurr=subset(EVER_FKNMS.hurr,Year.start%in%seq(1995,2019,1))
EVER_FKNMS.hurr.trk=subset(EVER_FKNMS.hurr,Key%in%EVER_FKNMS.hurr$Key)

plot(rep(1,nrow(EVER_FKNMS.hurr))~EVER_FKNMS.hurr$date.start)

plot(EVER_FKNMS.buffer1)
plot(EVER_FKNMS,add=T)
plot(subset(hur.tck,Key%in%EVER_FKNMS.hurr$Key),add=T)

EVER_FKNMS.hurr.shp=subset(hur.tck,Key%in%EVER_FKNMS.hurr$Key)
EVER_FKNMS.hurr.shp$Year.start=as.factor(EVER_FKNMS.hurr.shp$Year.start)

# Frequency analysis
grd=as.data.frame(spsample(EVER_FKNMS.buffer1,type="regular",cellsize=c(10000,10000)))
names(grd)       = c("UTMX", "UTMY")
coordinates(grd) = c("UTMX", "UTMY")
plot(grd)
plot(EVER_FKNMS.buffer1,add=T)
gridded(grd)=T;fullgrid(grd)=T
proj4string(grd)=wkt(utm17)

grd2=raster(grd)
values(grd2)=1:length(grd2)
grd2.2=rasterToPolygons(grd2)
# writeOGR(grd2.2,paste0(export.path,"GIS"),'SFLGrid',driver="ESRI Shapefile")

grd.hur=rasterize(EVER_FKNMS.hurr.shp,grd2,fun="count")
plot(grd.hur)

# png(filename=paste0(plot.path,"HurricaneFreq.png"),width=6.5,height=4,units="in",res=200,type="windows",bg="white")
par(family="serif",oma=c(0.25,0.25,0.25,0.25),mar=c(0.1,0.1,0.1,0.1),xpd=F)
layout(matrix(1:2,1,2,byrow=F),widths=c(1,0.3))
bbox.lims=bbox(EVER_FKNMS.buffer1)

plot(shore2,col="cornsilk",border="grey",bg="lightblue",ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.01)
plot(wca,add=T,col="grey90",border=NA)
plot(bcnp,add=T,col="grey90",border=NA)
plot(canal,add=T,col="lightblue",lwd=1)
b=1:4
pal=viridis::plasma(length(b)-1,alpha = 0.5)# hcl.colors(length(b)-1, "Spectral")
image(grd.hur,add=T,breaks=b,col=pal)
plot(grd,add=T,col=adjustcolor("grey",0.5),lwd=0.5)
plot(EVER_FKNMS,add=T,lwd=0.5)
box(lwd=1)
mapmisc::scaleBar(utm17,"bottomright",bty="n",cex=1,seg.len=4)

plot(0:1,0:1,ann=F,axes=F,type="n")
l.b=length(b)
labs=c(b[1],paste(b[2:(l.b-2)],b[3:(l.b-1)],sep=" - "),b[l.b])
n.bks=length(b)-1
top.val=0.8
bot.val=0.2
mid.v.val=bot.val+(top.val-bot.val)/2
x.max=0.6
x.min=0.3
mid.val=x.min+(x.max-x.min)/2
txt.offset.val=-0.01
# legend_image=as.raster(matrix(rev(pal2),ncol=1))
# rasterImage(legend_image,x.min,bot.val,x.max,top.val)
# text(x=x.max, y = c(bot.val,mid.v.val,top.val), labels = format(c(min(b2),0,max(b2))),cex=0.75,adj=0,pos=4,offset=0.5)
bx.val= seq(bot.val,top.val,(top.val-bot.val)/n.bks)
rect(x.min,bx.val[1:n.bks],x.max,bx.val[2:(n.bks+1)],col=rev(pal),lty=0)
text(y=bx.val[2:(n.bks+1)]-c(mean(diff(bx.val[2:(n.bks+1)]))/2), x = x.max, labels = rev(labs),cex=0.75,xpd=NA,pos=4,adj=0)
text(x=mid.val,y=top.val,"Number of storm tracks\n(1995 - 2019)",adj=0,cex=0.75,pos=3,xpd=NA)
dev.off()


cols=wesanderson::wes_palette("Zissou1",4,"continuous")
cols2=gray.colors(nrow(hurr.yr.N))
cols2[c(5,11)]=c("indianred1","indianred4")
# png(filename=paste0(plot.path,"HurricaneHist.png"),width=6.5,height=4,units="in",res=200,type="windows",bg="white")
par(family="serif",oma=c(0.25,0.25,0.25,0.25),mar=c(0.1,0.1,0.1,0.1),xpd=F)
# layout(matrix(c(1:2),1,2,byrow=T),widths = c(1,0.3))
# par(family="serif",mar=c(0.1,0.1,0.1,0.1),oma=c(0.5,0.5,0.5,0.5));
layout(matrix(c(1,1,2,2),2,2,byrow=F),widths=c(1,0.3),heights=c(0.75,0.75))
bbox.lims=bbox(EVER_FKNMS.buffer1)

plot(shore2,col="cornsilk",border="grey",bg="lightblue",ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.01)
plot(wca,add=T,col="grey90",border=NA)
plot(bcnp,add=T,col="grey90",border=NA)
plot(canal,add=T,col="lightblue",lwd=1)
plot(gSimplify(subset(regions2,Region=="Coastal_Mangroves"),500),add=T,col=adjustcolor(cols[1],0.5),lwd=0.5)
plot(gSimplify(subset(regions2,Region=="FLBay"),500),add=T,col=adjustcolor(cols[2],0.5),lwd=0.5)
plot(gSimplify(subset(regions2,Region=="Shelf"),500),add=T,col=adjustcolor(cols[3],0.5),lwd=0.5)
plot(gSimplify(subset(regions2,Region=="Keys"),500),add=T,col=adjustcolor(cols[4],0.5),lwd=0.5)
plot(ENP,add=T,bg=NA,lwd=1)
plot(EVER_FKNMS.hurr.shp,add=T,col=cols2[EVER_FKNMS.hurr.shp$Year.start],lwd=2)
plot(EVER_FKNMS.buffer1,add=T)
mapmisc::scaleBar(utm17,"bottomright",bty="n",cex=1,seg.len=4)
box(lwd=1)


plot(0:1,0:1,ann=F,axes=F,type="n")
legend("center",legend=c(paste(hurr.yr.N$Year.start,"(",hurr.yr.N$N.val,")")),
       pch=NA,lty=c(1),lwd=c(2),
       col=c(cols2),
       pt.cex=1.25,ncol=1,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=1,
       title.adj = 0.5,title="Year\n(Number of Hurricanes\nw/in 75 km)")

dev.off()


# Fire History ------------------------------------------------------------
all_fire=spTransform(readOGR("C:/Julian_LaCie/_GitHub/vistafire/data","all_fires_cost"),utm17)
sum(gIsValid(all_fire,byid=T)==F)
all_fire=gBuffer(all_fire, byid=TRUE, width=0)
sum(gIsValid(all_fire,byid=T)==F)
plot(all_fire)
# 
# # Geometry issues
# all_fire=readOGR(paste0(GIS.path,"/Fire"),"ENP_FIRES_1948-2019_utm")
# # cleangeo::clgeo_IsValid(all_fire)
# # report <- cleangeo::clgeo_CollectionReport(all_fire)
# # cleangeo::clgeo_SummaryReport(report)
# # test=cleangeo::clgeo_Clean(all_fire)
# 
# test=gIsValid(all_fire,byid=T,reason=T)
# test[test!="Valid Geometry"]
# test=gIsValid(all_fire,byid=T)
# test2=spTransform(all_fire,utm17)
# # test[test==F]
# # test[test==T]
# # gIsValid(all_fire[test==T,])
# # gIsValid(all_fire[test==F,])
# # 
# # fix.fire=all_fire[test==F,]
# # test.cln=cleangeo::clgeo_Clean(fix.fire)
# # 
# # gIsValid(test.cln)
# # fire.fix=bind(all_fire[test==T,],test.cln)
# # plot(fire.fix)
# # fire.fix=spTransform(fire.fix,utm17)
# 
# all_fire.smp=gSimplify(all_fire, tol = 0.0001)
# gIsValid(all_fire.smp)
# test=rmapshaper::ms_clip(all_fire,bbox=bbox(spTransform(gBuffer(ENP,width=1000),wgs84)))
# plot(test)
# 
# all_fire=spTransform(all_fire,utm17)
# plot(ENP)
# plot(all_fire,add=T)
# writeOGR(all_fire,paste0(GIS.path,"/Fire"),"ENP_FIRES_1948-2019_fix",driver="ESRI Shapefile")

range(all_fire$FireCalend)
all_fire1=subset(all_fire,FireCalend%in%seq(1995,2019,1))
head(all_fire1)
range(all_fire1$FireCalend)

all_fire_sub=raster::intersect(all_fire1,ENP)# gIntersection(all_fire1,ENP,byid = TRUE)
all_fire_sub$area.sqkm=area(all_fire_sub)*1e-6
plot(all_fire_sub)
all_fire_sub.sum=ddply(all_fire_sub@data,"FireCalend",summarise,TArea=sum(area.sqkm,na.rm=T),
                       per.area=(TArea/(gArea(raster::intersect(shore2,ENP))*1e-6))*100)

plot(all_fire_sub.sum$FireCalend,rep(1,length(all_fire_sub.sum$FireCalend)),pch=21,bg="grey",cex=all_fire_sub.sum$per.area/3)
plot(TArea~FireCalend,all_fire_sub.sum)
plot(per.area~FireCalend,all_fire_sub.sum)


grd.ENP=as.data.frame(spsample(ENP,type="regular",cellsize=c(100,100)))
names(grd.ENP)       = c("UTMX", "UTMY")
coordinates(grd.ENP) = c("UTMX", "UTMY")
plot(grd.ENP)
plot(ENP,add=T)
gridded(grd.ENP)=T;fullgrid(grd.ENP)=T
proj4string(grd.ENP)=wkt(utm17)

grd2.ENP=raster(grd.ENP)

grd.fire=rasterize(all_fire1,grd2.ENP,field=all_fire1@data$FireID,fun="count")
plot(grd.fire)
grd.fire <- mask(grd.fire, ENP)

# png(filename=paste0(plot.path,"FireFreq.png"),width=6.5,height=4,units="in",res=200,type="windows",bg="white")
par(family="serif",oma=c(0.25,0.25,0.25,0.25),mar=c(0.1,0.1,0.1,0.1),xpd=F)
layout(matrix(1:2,1,2,byrow=F),widths=c(1,0.3))
bbox.lims=bbox(ENP)

plot(shore2,col="cornsilk",border="grey",bg="lightblue",ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.01)
plot(wca,add=T,col="grey90",border=NA)
plot(bcnp,add=T,col="grey90",border=NA)
plot(canal,add=T,col="lightblue",lwd=1)
b=c(1,3,5,10)
pal=viridis::plasma(length(b)-1,alpha = 0.5)# hcl.colors(length(b)-1, "Spectral")
image(grd.fire,add=T,breaks=b,col=pal)
# plot(grd,add=T,col=adjustcolor("grey",0.5),lwd=0.5)
plot(EVER_FKNMS,add=T,lwd=0.5)
box(lwd=1)
mapmisc::scaleBar(utm17,"bottomright",bty="n",cex=1,seg.len=4)

plot(0:1,0:1,ann=F,axes=F,type="n")
l.b=length(b)
labs=c(paste(b[1],b[2],sep="-"),
       paste(b[2:(l.b-1)]+1,b[3:(l.b)],sep="-"))
n.bks=length(b)-1
top.val=0.8
bot.val=0.2
mid.v.val=bot.val+(top.val-bot.val)/2
x.max=0.5
x.min=0.4
mid.val=x.min+(x.max-x.min)/2
txt.offset.val=-0.01
# legend_image=as.raster(matrix(rev(pal2),ncol=1))
# rasterImage(legend_image,x.min,bot.val,x.max,top.val)
# text(x=x.max, y = c(bot.val,mid.v.val,top.val), labels = format(c(min(b2),0,max(b2))),cex=0.75,adj=0,pos=4,offset=0.5)
bx.val= seq(bot.val,top.val,(top.val-bot.val)/n.bks)
rect(x.min,bx.val[1:n.bks],x.max,bx.val[2:(n.bks+1)],col=rev(pal),lty=0)
text(y=bx.val[2:(n.bks+1)]-c(mean(diff(bx.val[2:(n.bks+1)]))/2), x = x.max, labels = rev(labs),cex=0.75,xpd=NA,pos=4,adj=0)
# text(y=bx.val, x = x.max, labels = rev(labs),cex=0.75,xpd=NA,pos=4,adj=0)
text(x=mid.val,y=top.val,"Number of fires\n(1995 - 2017)",adj=0,cex=0.75,pos=3,xpd=NA)
dev.off()

# Climate -----------------------------------------------------------------
## AMO
# vars=c('year',month.abb)
# row.count=length(seq(1856,2021,1))
# noaa.amo.path="https://psl.noaa.gov/data/correlation/amon.us.long.data"
# 
# # AMO.dat=read.table("https://psl.noaa.gov/data/correlation/amon.us.long.data",header=F,skip=1,col.names=vars,nrows=row.count,na.string="-99.990")
# AMO.dat=read.table("https://psl.noaa.gov/data/correlation/amon.sm.long.data",header=F,skip=1,col.names=vars,nrows=row.count,na.string="-99.990")
# AMO.dat.melt=melt(AMO.dat,id.vars="year")
# AMO.dat.melt=merge(AMO.dat.melt,data.frame(variable=month.abb,month=1:12))
# AMO.dat.melt$Date.mon=with(AMO.dat.melt,date.fun(paste(year,month,"01",sep="-")))
# AMO.dat.melt=AMO.dat.melt[order(AMO.dat.melt$Date.mon),c("Date.mon","value")]
# AMO.dat.melt$warm=with(AMO.dat.melt,ifelse(value>0,value,0))
# AMO.dat.melt$dry=with(AMO.dat.melt,ifelse(value<0,value,0))
# AMO.dat.melt$ma=with(AMO.dat.melt,c(rep(NA,120),zoo::rollapply(value,width=121,FUN=function(x)mean(x,na.rm=T))))
# head(AMO.dat.melt)
# tail(AMO.dat.melt)
# 
# layout(matrix(c(1:2),2,1,byrow=T))
# ylim.val=c(-0.4,0.4);by.y=0.5;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
# xlim.val=date.fun(c("1870-01-01","2016-12-01"));xmaj=seq(xlim.val[1],xlim.val[2],"20 years");xmin=seq(xlim.val[1],xlim.val[2],"1 years")
# 
# plot(value~Date.mon,AMO.dat.melt,xlim=xlim.val,ylim=ylim.val,type="n")
# #with(AMO.dat.melt,lines(Date.mon,ma,col="red"))
# with(subset(AMO.dat.melt,is.na(value)==F),shaded.range(Date.mon,rep(0,length(Date.mon)),ifelse(value>0,value,0),"indianred1",lty=1))
# with(subset(AMO.dat.melt,is.na(value)==F),shaded.range(Date.mon,ifelse(value<0,value,0),rep(0,length(Date.mon)),"dodgerblue1",lty=1))
# abline(h=0)
# 
# # PDO
# # https://www.ncdc.noaa.gov/teleconnections/pdo/
# # pdo=read.csv("https://www.ncdc.noaa.gov/teleconnections/pdo/data.csv",skip=1)
# nrow.val=length(seq(1854,2021,1))
# head.val=c("yr",month.abb)
# pdo=read.table(paste0(data.path,"NOAA/PDO/data.txt"),skip=2,header=F,col.names=head.val,nrows = nrow.val-1)
# head(pdo);tail(pdo)
# pdo$yr=as.numeric(pdo$yr)
# pdo=melt(pdo,id.var="yr")
# pdo$month.num=with(pdo,as.numeric(match(variable,month.abb)))
# pdo$monCY.date=with(pdo,date.fun(paste(yr,month.num,1,sep="-")))
# pdo$WY=WY(pdo$monCY.date)
# pdo$dec.WY=decimal.WY(pdo$monCY.date)
# pdo=pdo[order(pdo$dec.WY),]
# 
# pdo.WY.dat=ddply(pdo,"WY",summarise,mean.PDO=mean(as.numeric(value),na.rm=T),sd.PDO=sd(as.numeric(value),na.rm=T),N.val=N.obs(value))
# pdo.WY.dat$UCI=with(pdo.WY.dat,mean.PDO+qnorm(0.975)*sd.PDO/sqrt(N.val))
# pdo.WY.dat$LCI=with(pdo.WY.dat,mean.PDO-qnorm(0.975)*sd.PDO/sqrt(N.val))
# 
# plot(mean.PDO~WY,pdo.WY.dat)
# with(pdo.WY.dat,lines(WY,UCI))
# with(pdo.WY.dat,lines(WY,LCI))
# 
# # png(filename=paste0(plot.path,"ClimateIndex.png"),width=6.5,height=4,units="in",res=200,type="windows",bg="white")
# par(family="serif",mar=c(1,3,0.5,1.5),oma=c(2,1,0.25,0.25),xpd=F);
# layout(matrix(c(1:2),2,1,byrow=T))
# 
# xlim.val=date.fun(c("1870-01-01","2021-12-01"));xmaj=seq(xlim.val[1],xlim.val[2],"20 years");xmin=seq(xlim.val[1],xlim.val[2],"1 years")
# ylim.val=c(-0.4,0.4);by.y=0.2;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
# plot(value~Date.mon,AMO.dat.melt,xlim=xlim.val,ylim=ylim.val,type="n",axes=F,ann=F)
# abline(h=ymaj,v=xmaj,lty=1,col=adjustcolor("grey",0.5))
# #with(AMO.dat.melt,lines(Date.mon,ma,col="red"))
# with(subset(AMO.dat.melt,is.na(value)==F),shaded.range(Date.mon,rep(0,length(Date.mon)),ifelse(value>0,value,0),"indianred1",lty=1))
# with(subset(AMO.dat.melt,is.na(value)==F),shaded.range(Date.mon,ifelse(value<0,value,0),rep(0,length(Date.mon)),"dodgerblue1",lty=1))
# abline(h=0)
# axis_fun(1,xmaj,xmin,NA)
# axis_fun(2,ymaj,ymin,format(ymaj,nsmall=1));box(lwd=1)
# abline(v=date.fun(c("1995-05-01","2019-05-01")),lty=2)
# text(date.fun(date.fun("1995-05-01")+lubridate::ddays(4383)),ylim.val[2],"Study Period",cex=0.75,font=3)
# mtext(side=2,line=2.5,"AMO Index")
# 
# ylim.val=c(-4,4);by.y=2;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
# plot(value~monCY.date,pdo,xlim=xlim.val,ylim=ylim.val,type="n",axes=F,ann=F)
# abline(h=ymaj,v=xmaj,lty=1,col=adjustcolor("grey",0.5))
# with(subset(pdo,is.na(value)==F),shaded.range(monCY.date,rep(0,length(monCY.date)),ifelse(value>0,value,0),"indianred1",lty=1))
# with(subset(pdo,is.na(value)==F),shaded.range(monCY.date,ifelse(value<0,value,0),rep(0,length(monCY.date)),"dodgerblue1",lty=1))
# abline(h=0)
# axis_fun(1,xmaj,xmin,format(xmaj,"%Y"),line=-0.5)
# axis_fun(2,ymaj,ymin,format(ymaj,nsmall=1));box(lwd=1)
# abline(v=date.fun(c("1995-05-01","2019-05-01")),lty=2)
# mtext(side=2,line=2.5,"PDO Index")
# mtext(side=1,line=1.5,"Year")
# dev.off()
# 
# # png(filename=paste0(plot.path,"ClimateIndex_v2.png"),width=6.5,height=4,units="in",res=200,type="windows",bg="white")
# par(family="serif",mar=c(1,3,0.5,1.5),oma=c(2,1,0.25,0.25),xpd=F);
# layout(matrix(c(1:2),2,1,byrow=T))
# 
# xlim.val=date.fun(c("1995-01-01","2021-12-01"));xmaj=seq(xlim.val[1],xlim.val[2],"10 years");xmin=seq(xlim.val[1],xlim.val[2],"1 years")
# ylim.val=c(-0.2,0.2);by.y=0.1;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
# plot(value~Date.mon,AMO.dat.melt,xlim=xlim.val,ylim=ylim.val,type="n",axes=F,ann=F)
# abline(h=ymaj,v=xmaj,lty=1,col=adjustcolor("grey",0.5))
# #with(AMO.dat.melt,lines(Date.mon,ma,col="red"))
# with(subset(AMO.dat.melt,is.na(value)==F),shaded.range(Date.mon,rep(0,length(Date.mon)),ifelse(value>0,value,0),"indianred1",lty=1))
# with(subset(AMO.dat.melt,is.na(value)==F),shaded.range(Date.mon,ifelse(value<0,value,0),rep(0,length(Date.mon)),"dodgerblue1",lty=1))
# abline(h=0)
# axis_fun(1,xmaj,xmin,NA)
# axis_fun(2,ymaj,ymin,format(ymaj,nsmall=1));box(lwd=1)
# abline(v=date.fun(c("1995-05-01","2019-05-01")),lty=2)
# text(date.fun(date.fun("1995-05-01")+lubridate::ddays(4383)),ylim.val[2],"Study Period",cex=0.75,font=3)
# mtext(side=2,line=2.5,"AMO Index")
# 
# ylim.val=c(-4,4);by.y=2;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
# plot(value~monCY.date,pdo,xlim=xlim.val,ylim=ylim.val,type="n",axes=F,ann=F)
# abline(h=ymaj,v=xmaj,lty=1,col=adjustcolor("grey",0.5))
# with(subset(pdo,is.na(value)==F),shaded.range(monCY.date,rep(0,length(monCY.date)),ifelse(value>0,value,0),"indianred1",lty=1))
# with(subset(pdo,is.na(value)==F),shaded.range(monCY.date,ifelse(value<0,value,0),rep(0,length(monCY.date)),"dodgerblue1",lty=1))
# abline(h=0)
# axis_fun(1,xmaj,xmin,format(xmaj,"%Y"),line=-0.5)
# axis_fun(2,ymaj,ymin,format(ymaj,nsmall=1));box(lwd=1)
# abline(v=date.fun(c("1995-05-01","2019-05-01")),lty=2)
# mtext(side=2,line=2.5,"PDO Index")
# mtext(side=1,line=1.5,"Year")
# dev.off()


# https://origin.cpc.ncep.noaa.gov/products/analysis_monitoring/ensostuff/ONI_v5.php
ONI.dat=read.table("https://origin.cpc.ncep.noaa.gov/products/analysis_monitoring/ensostuff/detrend.nino34.ascii.txt",header=T)
ONI.dat$monCY.date=with(ONI.dat,date.fun(paste(YR,MON,01,sep="-")))
ONI.dat$mean.3mon=c(NA,round(rollmean(ONI.dat$ANOM,k=3,align="center"),1),NA)
ONI.dat$WY=WY(ONI.dat$monCY.date)
plot(mean.3mon~monCY.date,ONI.dat)

dates

xlim.val=dates;xmaj=seq(xlim.val[1],xlim.val[2],"5 years");xmin=seq(xlim.val[1],xlim.val[2],"1 years")

plot(mean.3mon~monCY.date,ONI.dat,xlim=xlim.val)
abline(h=0)

# H. Wilma October-2005
# Mustang Corner Fire May-2008
# Cold Snaps January 2010, December 2011
# Droughts October 2010-April 2011, May 2015-October 2015
# Flood November 2015-March 2016 
# H. Irma September 2017

events=data.frame(event=c("H. Wilma","Fire","Cold","Cold"),
                  Date.monCY=c("2005-10-01","2008-05-01","2010-01-01","2011-12-01"),
                  col=c("lightblue","indianred1","black","black"))

ylim.val=c(0,0.65)
plot(0:1~dates,type="n",ann=F,axes=F,xlim=xlim.val,ylim=ylim.val)
arrows(x0=date.fun("2005-10-01"),y0=0,y1=0.2,code=1,lwd=2,length=0.1)
text(date.fun("2005-10-01"),0.2,"H. Wilma",pos=3,cex=0.5,srt=45)
arrows(x0=date.fun("2017-09-10"),y0=0,y1=0.15,code=1,lwd=2,length=0.1)
text(date.fun("2017-09-10"),0.15,"H. Irma",pos=3,cex=0.5,srt=45)
# arrows(x0=date.fun("2017-07-2017"),y0=0,y1=0.20,code=1,lwd=2,length=0.1)
# text(date.fun("2017-07-2017"),0.20,"TS Emily",pos=3,cex=0.5,srt=90)
# arrows(x0=date.fun("2017-10-28"),y0=0,y1=0.25,code=1,lwd=2,length=0.1)
# text(date.fun("2017-10-28"),0.25,"TS Philippee",pos=3,cex=0.5,srt=90)
arrows(x0=date.fun("2008-05-01"),y0=0,y1=0.2,code=1,lwd=2,length=0.1,col="indianred1")
text(date.fun("2008-05-01"),0.2,"Fire",pos=3,cex=0.5,offset=0.1,col="indianred1")
arrows(x0=date.fun("2010-01-01"),y0=0,y1=0.2,code=1,lwd=2,length=0.1,col="dodgerblue1")
text(date.fun("2010-01-01"),0.2,"Cold",pos=3,cex=0.5,offset=0.1,col="dodgerblue1")
arrows(x0=date.fun("2011-12-01"),y0=0,y1=0.2,code=1,lwd=2,length=0.1,col="dodgerblue1")
text(date.fun("2011-12-01"),0.2,"Cold",pos=3,cex=0.5,offset=0.1,col="dodgerblue1")
text(xlim.val[1],0.1,"Pulse Events",xpd=NA)
xx=date.fun(c("2000-01-01","2004-12-01"));yy=c(0.35,0.35,0.45,0.45)
polygon(c(xx,rev(xx)),yy,col="khaki")
xx=date.fun(c("2006-01-01","2007-12-01"))
polygon(c(xx,rev(xx)),yy,col="khaki")
xx=date.fun(c("2010-10-01","2011-04-01"))
polygon(c(xx,rev(xx)),yy,col="khaki")
xx=date.fun(c("2015-05-01","2015-10-01"));
polygon(c(xx,rev(xx)),yy,col="khaki")
text(xlim.val[1],0.4,"Drought",xpd=NA)
# https://apps.sfwmd.gov/sfwmd/SFER/2017_sfer_final/v1/chapters/v1_ch3a.pdf
xx=date.fun(c("2015-10-01","2016-01-01"));yy=c(0.55,0.55,0.65,0.65)
polygon(c(xx,rev(xx)),yy,col="dodgerblue1")
# https://apps.sfwmd.gov/sfwmd/SFER/2019_sfer_final/v1/chapters/v1_ch3a.pdf
xx=date.fun(c("2017-06-05","2018-01-01"));yy=c(0.55,0.55,0.65,0.65)
polygon(c(xx,rev(xx)),yy,col="dodgerblue1")
text(xlim.val[1],0.6,"Flood",xpd=NA)
axis_fun(1,xmaj,xmin,format(xmaj,"%m-%Y"))

