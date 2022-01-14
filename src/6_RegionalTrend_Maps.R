
## Run 5_RegionalTrend_GAMs.R first

## 
# SamplingMap -------------------------------------------------------------
serc.shp2=SpatialPointsDataFrame(serc.shp[,c("UTMX","UTMY")],data=serc.shp,proj4string =utm17)
wmd.shp2=SpatialPointsDataFrame(wmd.shp[,c("UTMX","UTMY")],data=wmd.shp,proj4string =utm17)
lter.shp2=SpatialPointsDataFrame(lter.shp[,c("UTMX","UTMY")],data=lter.shp,proj4string =utm17)
# utm17=sf::st_crs(26917)[[2]]
# utm17=sp::CRS(utm17)

shore2=spTransform(readOGR(paste0(gen.GIS,"/FWC"),"FWC_Shoreline"),wkt(utm17))
shore2=gSimplify(shore2,250)
# roads.all=spTransform(readOGR(paste0(gen.GIS,"/FDOT"),"FDOT_Roads"),utm17)
lakes=spTransform(readOGR(paste0(gen.GIS,"/NHD"),"NHD100_Waterbody"),utm17)
wetland=subset(lakes,FTYPE%in%c("466"))

library(USAboundaries)

states.shp=us_boundaries(resolution ="low")
states.shp=as(states.shp,"Spatial")
states.shp=spTransform(states.shp,utm17)
attributes(states.shp)
SW.US=c("Florida","Georgia","Alabama","South Carolina")
FL.shp=us_boundaries(resolution ="low",states=SW.US)
FL.shp=as(FL.shp,"Spatial")
FL.shp=spTransform(FL.shp,utm17)
attributes(FL.shp)

cols=wesanderson::wes_palette("Zissou1",4,"continuous")
# png(filename=paste0(plot.path,"SamplingMap.png"),width=6.5,height=4,units="in",res=200,type="windows",bg="white")
# png(filename=paste0(plot.path,"SamplingMap.png"),width=6.5,height=4,units="in",res=200,type="windows",bg="white")
par(family="serif",oma=c(0.25,0.25,0.25,0.25),mar=c(0.1,0.1,0.1,0.1),xpd=F)
# layout(matrix(c(1:2),1,2,byrow=T),widths = c(1,0.3))
# par(family="serif",mar=c(0.1,0.1,0.1,0.1),oma=c(0.5,0.5,0.5,0.5));
layout(matrix(c(1,1,2,3),2,2,byrow=F),widths=c(1,0.3),heights=c(0.75,0.75))
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

AOI=raster::extent(gBuffer(region.mask,width=1000))
AOI.poly=as(AOI,"SpatialPolygons")
proj4string(AOI.poly)=utm17

bbox.lims=bbox(shore2)# bbox(gBuffer(region.mask,width=2000))
plot(states.shp,ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],col="grey75",border="grey",lwd=0.2,xpd=F)
plot(FL.shp,add=T,col=NA,border="white",xpd=F)
plot(lakes,add=T,border=NA,col=adjustcolor("grey95",0.5))
plot(wetland,add=T,border=NA,col=adjustcolor("grey40",0.5))
plot(canal,add=T,col="black",lwd=0.5)
# plot(roads.all,add=T,col="grey50",lwd=0.5,lty=1)
plot(AOI.poly,add=T,border=adjustcolor("red",0.5),lwd=2,lty=2)
box(lwd=1)
plot(0:1,0:1,ann=F,axes=F,type="n")
legend(0.5,0.9,legend=c(paste0("SERC (",length(serc.shp2$STATION),")"),
                        paste0("SFWMD (",length(wmd.shp2$STATION),")"),
                        paste0("FCE LTER (",length(lter.shp2$STATION),")")),
       pch=c(21,22,24),lty=c(NA),lwd=c(0.1),
       col=c("black"),pt.bg=c("white","grey","black"),
       pt.cex=1.25,ncol=1,cex=0.8,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=1,
       title.adj = 0,title="Monitoring\n(Number of Sites)")
legend(0.5,0.5,legend=c("Freshwater Everglades","Mangrove Fringe","Florida Bay","W. Florida Shelf","Keys"),
       pch=c(22),lty=c(NA),lwd=c(0.1),
       col=c("black"),pt.bg=c("white",adjustcolor(cols,0.5)),
       pt.cex=1.25,ncol=1,cex=0.8,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=1,
       title.adj = 0,title="Regions")
dev.off()


# Static Trend and GM maps ------------------------------------------------
sites.shp2=SpatialPointsDataFrame(sites.shp2[,c("UTMX","UTMY")],data=sites.shp2,proj4string=utm17)

cols.val=c("grey","white","red")
# png(filename=paste0(plot.path,"TrendMapsV2.png"),width=6.5,height=7,units="in",res=200,type="windows",bg="white")
par(family="serif",oma=c(0.25,0.25,0.25,0.25),mar=c(0.1,0.1,0.1,0.1),xpd=F)
layout(matrix(c(1:24),6,4,byrow=T),widths = c(1,0.4,1,0.5))
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
         col=c("grey"),pt.bg=c("grey",cols.val[2:3]),
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
  text(x=0.25, y = c(0.83,0.47), labels = c("< 2","40"),cex=0.6,adj=0,pos=4)
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
  plot(0:1,0:1,ann=F,axes=F,type="n")
  int.bks.vals=format(round(int.bks,3),nsmall=3)
  labs=c(paste0("< ",int.bks.vals[2]),paste(int.bks.vals[2:3],int.bks.vals[3:4],sep=" - "),paste(paste0(">",int.bks.vals[4])))
  n.bks=length(int.bks)-1
  bx.val= seq(0.45,0.85,(0.85-0.45)/n.bks)
  rect(0.15,bx.val[1:n.bks],0.25,bx.val[2:(n.bks+1)],col=rev(pal),lty=0)
  text(x=0.25, y = bx.val[2:(n.bks+1)]-c(mean(diff(bx.val[2:(n.bks+1)]))/2), labels = rev(labs),cex=0.6,adj=0,pos=4)
  text(x=0.15,y=0.95,"TOC Thiel-Sen Slope\n(mg C L\u207B\u00B9 Yr\u207B\u00B9)",adj=0,cex=0.75)
  
  plot(shore,col="cornsilk",border="grey",bg="lightblue",ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.05,xpd=F)
  int=classIntervals(values(tps.TOC.GM)[is.na(values(tps.TOC.GM))==F],n=n,style="equal")
  int.bks=int$brks
  pal=hcl.colors(length(int$brks)-1, "Inferno", rev = F,alpha=0.7)
  image(tps.TOC.GM,add=T,breaks=int.bks,col = pal)
  plot(ENP,add=T,bg=NA,lwd=0.5)
  plot(regions2,lty=2,add=T,border="grey80",lwd=0.5)
  plot(sites.shp.TOC.GM,add=T,pch=21,cex=0.25,bg=adjustcolor("white",0.5),col=NA)
  box(lwd=1)
  plot(0:1,0:1,ann=F,axes=F,type="n")
  legend_image=as.raster(matrix(pal,ncol=1))
  text(x=0.25, y = c(0.83,0.47), labels = c("< 0.5",">20"),cex=0.6,adj=0,pos=4)
  rasterImage(legend_image,0.15,0.45,0.25,0.85)
  text(x=0.15,y=0.95,"Average Annual GM TOC\n(mg C L\u207B\u00B9)",adj=0,cex=0.75)
}
dev.off()

# tm_shape(tps.TOC.trend)+tm_raster(palette="viridis",
#                                   breaks=c(-1.8, -0.9, 0, 0.20,0.5))+
#   tm_shape(sites.shp.TOC.trend)+tm_dots(col="stat.sig")
# subset(sites.shp.TOC.trend@data,sen.slope>0.12)

# png(filename=paste0(plot.path,"varianceMaps.png"),width=6.5,height=7,units="in",res=200,type="windows",bg="white")
par(family="serif",oma=c(0.25,0.25,0.25,0.25),mar=c(0.1,0.1,0.1,0.1),xpd=F)
layout(matrix(c(1:24),6,4,byrow=T),widths = c(1,0.5,1,0.5))
bbox.lims=bbox(region.mask)
n=10
CV.map.bks=seq(0,100,20);#c(0,5,10,20,40,80,100)
{
  # TN
  plot(shore,col="cornsilk",border="grey",bg="lightblue",ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.05,xpd=F)
  int=classIntervals(values(tps.TN.GM.var)[is.na(values(tps.TN.GM.var))==F],n=n,style="equal")
  int.bks=int$brks
  pal=hcl.colors(length(int$brks)-1, "Inferno", rev = F,alpha=0.7)
  image(tps.TN.GM.var,add=T,breaks=int.bks,col = pal)
  plot(ENP,add=T,bg=NA,lwd=0.5)
  plot(regions2,lty=2,add=T,border="grey80",lwd=0.5)
  box(lwd=1)
  plot(0:1,0:1,ann=F,axes=F,type="n")
  legend_image=as.raster(matrix(pal,ncol=1))
  text(x=0.25, y = c(0.83,0.47), labels = c(format(min(int.bks),digits=2),format(max(int.bks),digits=2)),cex=0.6,adj=0,pos=4)
  rasterImage(legend_image,0.15,0.45,0.25,0.85)
  text(x=0.15,y=0.95,"Variance Annual GM TN\n(mg N L\u207B\u00B9)\u00B2",adj=0,cex=0.75)
  
  plot(shore,col="cornsilk",border="grey",bg="lightblue",ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.05,xpd=F)
  int=classIntervals(values(tps.TN.GM.CV)[is.na(values(tps.TN.GM.CV))==F],n=n,style="equal")
  int.bks=CV.map.bks# c(0,25,50,75,100)# int$brks
  pal=hcl.colors(length(int.bks)-1, "Plasma", rev = F,alpha=0.7)
  image(tps.TN.GM.CV,add=T,breaks=int.bks,col = pal)
  plot(ENP,add=T,bg=NA,lwd=0.5)
  plot(regions2,lty=2,add=T,border="grey80",lwd=0.5)
  box(lwd=1)
  plot(0:1,0:1,ann=F,axes=F,type="n")
  int.bks.vals=format(round(int.bks,1),nsmall=0)
  n.bks=length(int.bks)
  labs=c(paste(int.bks.vals[1:n.bks-1],int.bks.vals[2:n.bks],sep=" - "))
  n.bks=n.bks-1
  bx.val= seq(0.45,0.85,(0.85-0.45)/n.bks)
  rect(0.15,bx.val[1:n.bks],0.25,bx.val[2:(n.bks+1)],col=rev(pal),lty=0)
  text(x=0.25, y = bx.val[2:(n.bks+1)]-c(mean(diff(bx.val[2:(n.bks+1)]))/2), labels = rev(labs),cex=0.6,adj=0,pos=4)
  text(x=0.15,y=0.95,"CV Annual GM TN\n(Percent)",adj=0,cex=0.75)
  
  # DIN
  plot(shore,col="cornsilk",border="grey",bg="lightblue",ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.05,xpd=F)
  int=classIntervals(values(tps.DIN.GM.var)[is.na(values(tps.DIN.GM.var))==F],n=n,style="equal")
  int.bks=int$brks
  pal=hcl.colors(length(int$brks)-1, "Inferno", rev = F,alpha=0.7)
  image(tps.DIN.GM.var,add=T,breaks=int.bks,col = pal)
  plot(ENP,add=T,bg=NA,lwd=0.5)
  plot(regions2,lty=2,add=T,border="grey80",lwd=0.5)
  box(lwd=1)
  plot(0:1,0:1,ann=F,axes=F,type="n")
  legend_image=as.raster(matrix(pal,ncol=1))
  text(x=0.25, y = c(0.83,0.47), labels = c(format(min(int.bks),digits=2),format(max(int.bks),digits=2)),cex=0.6,adj=0,pos=4)
  rasterImage(legend_image,0.15,0.45,0.25,0.85)
  text(x=0.15,y=0.95,"Variance Annual GM DIN\n(mg N L\u207B\u00B9)\u00B2",adj=0,cex=0.75)
  
  plot(shore,col="cornsilk",border="grey",bg="lightblue",ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.05,xpd=F)
  int=classIntervals(values(tps.DIN.GM.CV)[is.na(values(tps.DIN.GM.CV))==F],n=n,style="equal")
  int.bks=CV.map.bks# c(0,25,50,75,100)# int$brks
  pal=hcl.colors(length(int.bks)-1, "Plasma", rev = F,alpha=0.7)
  image(tps.DIN.GM.CV,add=T,breaks=int.bks,col = pal)
  plot(ENP,add=T,bg=NA,lwd=0.5)
  plot(regions2,lty=2,add=T,border="grey80",lwd=0.5)
  box(lwd=1)
  plot(0:1,0:1,ann=F,axes=F,type="n")
  int.bks.vals=format(round(int.bks,1),nsmall=0)
  n.bks=length(int.bks)
  labs=c(paste(int.bks.vals[1:n.bks-1],int.bks.vals[2:n.bks],sep=" - "))
  n.bks=n.bks-1
  bx.val= seq(0.45,0.85,(0.85-0.45)/n.bks)
  rect(0.15,bx.val[1:n.bks],0.25,bx.val[2:(n.bks+1)],col=rev(pal),lty=0)
  text(x=0.25, y = bx.val[2:(n.bks+1)]-c(mean(diff(bx.val[2:(n.bks+1)]))/2), labels = rev(labs),cex=0.6,adj=0,pos=4)
  text(x=0.15,y=0.95,"CV Annual GM DIN\n(Percent)",adj=0,cex=0.75)
  
  # TP
  plot(shore,col="cornsilk",border="grey",bg="lightblue",ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.05,xpd=F)
  int=classIntervals(values(tps.TP.GM.var)[is.na(values(tps.TP.GM.var))==F],n=n,style="equal")
  int.bks=int$brks
  pal=hcl.colors(length(int$brks)-1, "Inferno", rev = F,alpha=0.7)
  image(tps.TP.GM.var,add=T,breaks=int.bks,col = pal)
  plot(ENP,add=T,bg=NA,lwd=0.5)
  plot(regions2,lty=2,add=T,border="grey80",lwd=0.5)
  box(lwd=1)
  plot(0:1,0:1,ann=F,axes=F,type="n")
  legend_image=as.raster(matrix(pal,ncol=1))
  text(x=0.25, y = c(0.83,0.47), labels = c(format(min(int.bks),digits=2),format(max(int.bks),digits=2)),cex=0.6,adj=0,pos=4)
  rasterImage(legend_image,0.15,0.45,0.25,0.85)
  text(x=0.15,y=0.95,"Variance Annual GM TP\n(\u03BCg N L\u207B\u00B9)\u00B2",adj=0,cex=0.75)
  
  plot(shore,col="cornsilk",border="grey",bg="lightblue",ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.05,xpd=F)
  int=classIntervals(values(tps.TP.GM.CV)[is.na(values(tps.TP.GM.CV))==F],n=n,style="equal")
  int.bks=CV.map.bks# c(0,25,50,75,100)# int$brks
  pal=hcl.colors(length(int.bks)-1, "Plasma", rev = F,alpha=0.7)
  image(tps.TP.GM.CV,add=T,breaks=int.bks,col = pal)
  plot(ENP,add=T,bg=NA,lwd=0.5)
  plot(regions2,lty=2,add=T,border="grey80",lwd=0.5)
  box(lwd=1)
  plot(0:1,0:1,ann=F,axes=F,type="n")
  int.bks.vals=format(round(int.bks,1),nsmall=0)
  n.bks=length(int.bks)
  labs=c(paste(int.bks.vals[1:n.bks-1],int.bks.vals[2:n.bks],sep=" - "))
  n.bks=n.bks-1
  bx.val= seq(0.45,0.85,(0.85-0.45)/n.bks)
  rect(0.15,bx.val[1:n.bks],0.25,bx.val[2:(n.bks+1)],col=rev(pal),lty=0)
  text(x=0.25, y = bx.val[2:(n.bks+1)]-c(mean(diff(bx.val[2:(n.bks+1)]))/2), labels = rev(labs),cex=0.6,adj=0,pos=4)
  text(x=0.15,y=0.95,"CV Annual GM TP\n(Percent)",adj=0,cex=0.75)
  
  # SRP
  plot(shore,col="cornsilk",border="grey",bg="lightblue",ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.05,xpd=F)
  int=classIntervals(values(tps.SRP.GM.var)[is.na(values(tps.SRP.GM.var))==F],n=n,style="equal")
  int.bks=int$brks
  pal=hcl.colors(length(int$brks)-1, "Inferno", rev = F,alpha=0.7)
  image(tps.SRP.GM.var,add=T,breaks=int.bks,col = pal)
  plot(ENP,add=T,bg=NA,lwd=0.5)
  plot(regions2,lty=2,add=T,border="grey80",lwd=0.5)
  box(lwd=1)
  plot(0:1,0:1,ann=F,axes=F,type="n")
  legend_image=as.raster(matrix(pal,ncol=1))
  text(x=0.25, y = c(0.83,0.47), labels = c(format(min(int.bks),digits=2),format(max(int.bks),digits=2)),cex=0.6,adj=0,pos=4)
  rasterImage(legend_image,0.15,0.45,0.25,0.85)
  text(x=0.15,y=0.95,"Variance Annual GM SRP\n(\u03BCg N L\u207B\u00B9)\u00B2",adj=0,cex=0.75)
  
  plot(shore,col="cornsilk",border="grey",bg="lightblue",ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.05,xpd=F)
  int=classIntervals(values(tps.SRP.GM.CV)[is.na(values(tps.SRP.GM.CV))==F],n=n,style="equal")
  int.bks=CV.map.bks# c(0,25,50,75,100)# int$brks
  pal=hcl.colors(length(int.bks)-1, "Plasma", rev = F,alpha=0.7)
  image(tps.SRP.GM.CV,add=T,breaks=int.bks,col = pal)
  plot(ENP,add=T,bg=NA,lwd=0.5)
  plot(regions2,lty=2,add=T,border="grey80",lwd=0.5)
  box(lwd=1)
  plot(0:1,0:1,ann=F,axes=F,type="n")
  int.bks.vals=format(round(int.bks,1),nsmall=0)
  n.bks=length(int.bks)
  labs=c(paste(int.bks.vals[1:n.bks-1],int.bks.vals[2:n.bks],sep=" - "))
  n.bks=n.bks-1
  bx.val= seq(0.45,0.85,(0.85-0.45)/n.bks)
  rect(0.15,bx.val[1:n.bks],0.25,bx.val[2:(n.bks+1)],col=rev(pal),lty=0)
  text(x=0.25, y = bx.val[2:(n.bks+1)]-c(mean(diff(bx.val[2:(n.bks+1)]))/2), labels = rev(labs),cex=0.6,adj=0,pos=4)
  text(x=0.15,y=0.95,"CV Annual GM SRP\n(Percent)",adj=0,cex=0.75)
  
  # Chla
  plot(shore,col="cornsilk",border="grey",bg="lightblue",ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.05,xpd=F)
  int=classIntervals(values(tps.Chla.GM.var)[is.na(values(tps.Chla.GM.var))==F],n=n,style="equal")
  int.bks=int$brks
  pal=hcl.colors(length(int$brks)-1, "Inferno", rev = F,alpha=0.7)
  image(tps.Chla.GM.var,add=T,breaks=int.bks,col = pal)
  plot(ENP,add=T,bg=NA,lwd=0.5)
  plot(regions2,lty=2,add=T,border="grey80",lwd=0.5)
  box(lwd=1)
  plot(0:1,0:1,ann=F,axes=F,type="n")
  legend_image=as.raster(matrix(pal,ncol=1))
  text(x=0.25, y = c(0.83,0.47), labels = c(format(min(int.bks),digits=2),format(max(int.bks),digits=2)),cex=0.6,adj=0,pos=4)
  rasterImage(legend_image,0.15,0.45,0.25,0.85)
  text(x=0.15,y=0.95,"Variance Annual GM Chla\n(\u03BCg L\u207B\u00B9)\u00B2",adj=0,cex=0.75)
  
  plot(shore,col="cornsilk",border="grey",bg="lightblue",ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.05,xpd=F)
  int=classIntervals(values(tps.Chla.GM.CV)[is.na(values(tps.Chla.GM.CV))==F],n=n,style="equal")
  int.bks=CV.map.bks# c(0,25,50,75,100)# int$brks
  pal=hcl.colors(length(int.bks)-1, "Plasma", rev = F,alpha=0.7)
  image(tps.Chla.GM.CV,add=T,breaks=int.bks,col = pal)
  plot(ENP,add=T,bg=NA,lwd=0.5)
  plot(regions2,lty=2,add=T,border="grey80",lwd=0.5)
  box(lwd=1)
  plot(0:1,0:1,ann=F,axes=F,type="n")
  int.bks.vals=format(round(int.bks,1),nsmall=0)
  n.bks=length(int.bks)
  labs=c(paste(int.bks.vals[1:n.bks-1],int.bks.vals[2:n.bks],sep=" - "))
  n.bks=n.bks-1
  bx.val= seq(0.45,0.85,(0.85-0.45)/n.bks)
  rect(0.15,bx.val[1:n.bks],0.25,bx.val[2:(n.bks+1)],col=rev(pal),lty=0)
  text(x=0.25, y = bx.val[2:(n.bks+1)]-c(mean(diff(bx.val[2:(n.bks+1)]))/2), labels = rev(labs),cex=0.6,adj=0,pos=4)
  text(x=0.15,y=0.95,"CV Annual GM Chla\n(Percent)",adj=0,cex=0.75)
  
  # TOC
  plot(shore,col="cornsilk",border="grey",bg="lightblue",ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.05,xpd=F)
  int=classIntervals(values(tps.TOC.GM.var)[is.na(values(tps.TOC.GM.var))==F],n=n,style="equal")
  int.bks=int$brks
  pal=hcl.colors(length(int$brks)-1, "Inferno", rev = F,alpha=0.7)
  image(tps.TOC.GM.var,add=T,breaks=int.bks,col = pal)
  plot(ENP,add=T,bg=NA,lwd=0.5)
  plot(regions2,lty=2,add=T,border="grey80",lwd=0.5)
  box(lwd=1)
  plot(0:1,0:1,ann=F,axes=F,type="n")
  legend_image=as.raster(matrix(pal,ncol=1))
  text(x=0.25, y = c(0.83,0.47), labels = c(format(min(int.bks),digits=2),format(max(int.bks),digits=2)),cex=0.6,adj=0,pos=4)
  rasterImage(legend_image,0.15,0.45,0.25,0.85)
  text(x=0.15,y=0.95,"Variance Annual GM TOC\n(mg L\u207B\u00B9)\u00B2",adj=0,cex=0.75)
  
  plot(shore,col="cornsilk",border="grey",bg="lightblue",ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.05,xpd=F)
  int=classIntervals(values(tps.TOC.GM.CV)[is.na(values(tps.TOC.GM.CV))==F],n=n,style="equal")
  int.bks=CV.map.bks# c(0,25,50,75,100)# int$brks
  pal=hcl.colors(length(int.bks)-1, "Plasma", rev = F,alpha=0.7)
  image(tps.TOC.GM.CV,add=T,breaks=int.bks,col = pal)
  plot(ENP,add=T,bg=NA,lwd=0.5)
  plot(regions2,lty=2,add=T,border="grey80",lwd=0.5)
  box(lwd=1)
  plot(0:1,0:1,ann=F,axes=F,type="n")
  int.bks.vals=format(round(int.bks,1),nsmall=0)
  n.bks=length(int.bks)
  labs=c(paste(int.bks.vals[1:n.bks-1],int.bks.vals[2:n.bks],sep=" - "))
  n.bks=n.bks-1
  bx.val= seq(0.45,0.85,(0.85-0.45)/n.bks)
  rect(0.15,bx.val[1:n.bks],0.25,bx.val[2:(n.bks+1)],col=rev(pal),lty=0)
  text(x=0.25, y = bx.val[2:(n.bks+1)]-c(mean(diff(bx.val[2:(n.bks+1)]))/2), labels = rev(labs),cex=0.6,adj=0,pos=4)
  text(x=0.15,y=0.95,"CV Annual GM TOC\n(Percent)",adj=0,cex=0.75)
  
}
dev.off()

# png(filename=paste0(plot.path,"TrendMaps_TN.png"),width=7.5,height=2.25,units="in",res=200,type="windows",bg="white")
par(family="serif",oma=c(0.25,0.25,0.25,0.25),mar=c(0.1,0.1,0.1,0.1),xpd=F)
layout(matrix(c(1:4),1,4,byrow=T),widths = c(1,0.4,1,0.5))
bbox.lims=bbox(region.mask)
{
  # TN
  plot(shore,col="cornsilk",border="grey",bg="lightblue",ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.05)
  n=10
  int.bks=c(min(values(tps.TN.trend),na.rm=T),0+min(values(tps.TN.trend),na.rm=T)/2,0,0+max(values(tps.TN.trend),na.rm=T)/2,max(values(tps.TN.trend),na.rm=T)) #int$brks;#int.bks[2]=0
  pal=hcl.colors(length(int.bks)-1, "viridis", rev = F,alpha=0.7)
  image(tps.TN.trend,add=T,breaks=int.bks,col = pal)
  plot(ENP,add=T,bg=NA,lwd=0.5)
  plot(regions2,lty=2,add=T,border="grey80",lwd=0.5)
  plot(sites.shp2,add=T,cex=0.40,pch=21,bg=adjustcolor("grey",0.5),col="grey",lwd=0.1)
  plot(sites.shp.TN.trend,add=T,pch=21,cex=0.5,bg=adjustcolor(cols.val,0.5)[sites.shp.TN.trend$stat.sig],col=NA);
  box(lwd=1)
  mapmisc::scaleBar(utm17,"bottomright",bty="n",cex=1,seg.len=4)
  #
  plot(0:1,0:1,ann=F,axes=F,type="n")
  int.bks.vals=format(round(int.bks,3),nsmall=3)
  labs=c(paste0("< ",int.bks.vals[2]),paste(int.bks.vals[2:3],int.bks.vals[3:4],sep=" - "),paste(paste0(">",int.bks.vals[4])))
  n.bks=length(int.bks)-1
  bx.val= seq(0.45,0.85,(0.85-0.45)/n.bks)
  rect(0.15,bx.val[1:n.bks],0.25,bx.val[2:(n.bks+1)],col=rev(pal),lty=0)
  text(x=0.25, y = bx.val[2:(n.bks+1)]-c(mean(diff(bx.val[2:(n.bks+1)]))/2), labels = rev(labs),cex=0.6,adj=0,pos=4)
  text(x=0.15,y=0.95,"TN Thiel-Sen Slope\n(mg N L\u207B\u00B9 Yr\u207B\u00B9)",adj=0,cex=0.75)
  legend(0.5,0.4,legend=c("Insufficent Data","\u03C1 > 0.05","\u03C1 < 0.05"),
         pch=21,lty=c(NA),lwd=c(0.1),
         col=c("grey"),pt.bg=c("grey",cols.val[2:3]),
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
}
dev.off()

# png(filename=paste0(plot.path,"TrendMaps_TP.png"),width=7.5,height=2.25,units="in",res=200,type="windows",bg="white")
par(family="serif",oma=c(0.25,0.25,0.25,0.25),mar=c(0.1,0.1,0.1,0.1),xpd=F)
layout(matrix(c(1:4),1,4,byrow=T),widths = c(1,0.4,1,0.5))
bbox.lims=bbox(region.mask)
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
  plot(0:1,0:1,ann=F,axes=F,type="n")
  int.bks.vals=format(round(int.bks,3),nsmall=3)
  labs=c(paste0("< ",int.bks.vals[2]),paste(int.bks.vals[2:3],int.bks.vals[3:4],sep=" - "),paste(paste0(">",int.bks.vals[4])))
  n.bks=length(int.bks)-1
  bx.val= seq(0.45,0.85,(0.85-0.45)/n.bks)
  rect(0.15,bx.val[1:n.bks],0.25,bx.val[2:(n.bks+1)],col=rev(pal),lty=0)
  text(x=0.25, y = bx.val[2:(n.bks+1)]-c(mean(diff(bx.val[2:(n.bks+1)]))/2), labels = rev(labs),cex=0.6,adj=0,pos=4)
  text(x=0.15,y=0.95,"TP Thiel-Sen Slope\n(\u03BCg P L\u207B\u00B9 Yr\u207B\u00B9)",adj=0,cex=0.75)
  legend(0.5,0.4,legend=c("Insufficent Data","\u03C1 > 0.05","\u03C1 < 0.05"),
         pch=21,lty=c(NA),lwd=c(0.1),
         col=c("grey"),pt.bg=c("grey",cols.val[2:3]),
         pt.cex=1,ncol=1,cex=0.6,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=1,title.adj = 0,title="Kendall's trend \u03C1-value")
  
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
  
}
dev.off()

# png(filename=paste0(plot.path,"TrendMaps_NP.png"),width=7.5,height=2.25,units="in",res=200,type="windows",bg="white")
par(family="serif",oma=c(0.25,0.25,0.25,0.25),mar=c(0.1,0.1,0.1,0.1),xpd=F)
layout(matrix(c(1:4),1,4,byrow=T),widths = c(1,0.4,1,0.5))
bbox.lims=bbox(region.mask)
{
  plot(shore,col="cornsilk",border="grey",bg="lightblue",ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.05,xpd=F)
  int.bks=c(min(values(tps.NP.trend),na.rm=T),0+min(values(tps.NP.trend),na.rm=T)/2,0,0+max(values(tps.NP.trend),na.rm=T)/2,max(values(tps.NP.trend),na.rm=T)) #int$brks;#int.bks[2]=0
  pal=hcl.colors(length(int.bks)-1, "viridis", rev = F,alpha=0.7)
  image(tps.NP.trend,add=T,breaks=int.bks,col = pal)
  plot(ENP,add=T,bg=NA,lwd=0.5)
  plot(regions2,lty=2,add=T,border="grey80",lwd=0.5)
  plot(sites.shp2,add=T,cex=0.40,pch=21,bg=adjustcolor("grey",0.5),col="grey",lwd=0.1)
  plot(sites.shp.NP.trend,add=T,pch=21,cex=0.5,bg=adjustcolor(cols.val,0.5)[sites.shp.NP.trend$stat.sig],col=NA);
  box(lwd=1)
  mapmisc::scaleBar(utm17,"bottomright",bty="n",cex=1,seg.len=4)
  plot(0:1,0:1,ann=F,axes=F,type="n")
  int.bks.vals=format(round(int.bks,3),nsmall=3)
  labs=c(paste0("< ",int.bks.vals[2]),paste(int.bks.vals[2:3],int.bks.vals[3:4],sep=" - "),paste(paste0(">",int.bks.vals[4])))
  n.bks=length(int.bks)-1
  bx.val= seq(0.45,0.85,(0.85-0.45)/n.bks)
  rect(0.15,bx.val[1:n.bks],0.25,bx.val[2:(n.bks+1)],col=rev(pal),lty=0)
  text(x=0.25, y = bx.val[2:(n.bks+1)]-c(mean(diff(bx.val[2:(n.bks+1)]))/2), labels = rev(labs),cex=0.6,adj=0,pos=4)
  text(x=0.15,y=0.95,"N:P Thiel-Sen Slope\n(Yr\u207B\u00B9)",adj=0,cex=0.75)
  legend(0.5,0.4,legend=c("Insufficent Data","\u03C1 > 0.05","\u03C1 < 0.05"),
         pch=21,lty=c(NA),lwd=c(0.1),
         col=c("grey"),pt.bg=c("grey",cols.val[2:3]),
         pt.cex=1,ncol=1,cex=0.6,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=1,title.adj = 0,title="Kendall's trend \u03C1-value")
  
  plot(shore,col="cornsilk",border="grey",bg="lightblue",ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.05,xpd=F)
  int=classIntervals(values(tps.NP.GM)[is.na(values(tps.NP.GM))==F],n=n,style="equal")
  int.bks=int$brks
  pal=hcl.colors(length(int$brks)-1, "Inferno", rev = F,alpha=0.7)
  image(tps.NP.GM,add=T,breaks=int.bks,col = pal)
  plot(ENP,add=T,bg=NA,lwd=0.5)
  plot(regions2,lty=2,add=T,border="grey80",lwd=0.5)
  plot(sites.shp.NP.GM,add=T,pch=21,cex=0.25,bg=adjustcolor("white",0.5),col=NA)
  box(lwd=1)
  plot(0:1,0:1,ann=F,axes=F,type="n")
  legend_image=as.raster(matrix(pal,ncol=1))
  text(x=0.25, y = c(0.83,0.47), labels = c("< 5","> 550"),cex=0.6,adj=0,pos=4)
  rasterImage(legend_image,0.15,0.45,0.25,0.85)
  text(x=0.15,y=0.95,"Average Annual TN:TP\n(unitless)",adj=0,cex=0.75)
  
}
dev.off()

# png(filename=paste0(plot.path,"minNP_map.png"),width=6,height=3,units="in",res=200,type="windows",bg="white")
par(family="serif",oma=c(0.25,0.25,0.25,0.25),mar=c(0.1,0.1,0.1,0.1),xpd=F)
layout(matrix(c(1:2),1,2,byrow=T),widths = c(1,0.5))

plot(shore,col="cornsilk",border="grey",bg="lightblue",ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.05,xpd=F)
plot(ENP,add=T,bg=NA,lwd=0.5)
plot(regions2,lty=2,add=T,border="grey80",lwd=0.5)
bks=c(0,16,21,1000)
sites.shp.NP.GM$redfield=with(sites.shp.NP.GM@data, as.factor(findInterval(min.val,bks,left.open = FALSE,rightmost.closed = TRUE)))
cols.lim=c("red","yellow","green")
# cols.lim=with(sites.shp.NP.GM@data,ifelse(min.val<=16,"red","green"))
plot(sites.shp.NP.GM,add=T,pch=21,cex=0.8,
     bg=cols.lim[sites.shp.NP.GM@data$redfield],lwd=0.1)
box(lwd=1)
mapmisc::scaleBar(utm17,"bottomright",bty="n",cex=1,seg.len=4)
plot(0:1,0:1,ann=F,axes=F,type="n")
legend("center",legend=c("\u2264 16","16 - 21", "\u003E 21"),
       pch=21,lty=c(NA),lwd=c(0.1),
       col=c("grey"),pt.bg=c("red","yellow","green"),
       pt.cex=1.5,ncol=1,cex=1,bty="n",y.intersp=1.25,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=1,title.adj = 0.5,title="Minimum POR\nTN:TP")
dev.off()


# png(filename=paste0(plot.path,"varianceMap_TNTP.png"),width=6,height=3,units="in",res=200,type="windows",bg="white")
par(family="serif",oma=c(0.25,0.25,0.25,0.25),mar=c(0.1,0.1,0.1,0.1),xpd=F)
layout(matrix(c(1:2),1,2,byrow=T),widths = c(1,0.5))
bbox.lims=bbox(region.mask)
n=10
CV.map.bks=seq(0,100,20);#c(0,5,10,20,40,80,100)
{
  # TN
  plot(shore,col="cornsilk",border="grey",bg="lightblue",ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.05,xpd=F)
  int.bks=CV.map.bks
  pal=hcl.colors(length(int.bks)-1, "Inferno", rev = F,alpha=0.7)
  image(tps.NP.GM.CV,add=T,breaks=int.bks,col = pal)
  plot(ENP,add=T,bg=NA,lwd=0.5)
  plot(regions2,lty=2,add=T,border="grey80",lwd=0.5)
  box(lwd=1)
  mapmisc::scaleBar(utm17,"bottomright",bty="n",cex=1,seg.len=4)
  plot(0:1,0:1,ann=F,axes=F,type="n")
  int.bks.vals=format(round(int.bks,1),nsmall=0)
  n.bks=length(int.bks)
  labs=c(paste(int.bks.vals[1:n.bks-1],int.bks.vals[2:n.bks],sep=" - "))
  n.bks=n.bks-1
  bx.val= seq(0.45,0.85,(0.85-0.45)/n.bks)
  rect(0.15,bx.val[1:n.bks],0.25,bx.val[2:(n.bks+1)],col=rev(pal),lty=0)
  text(x=0.25, y = bx.val[2:(n.bks+1)]-c(mean(diff(bx.val[2:(n.bks+1)]))/2), labels = rev(labs),cex=0.6,adj=0,pos=4)
  text(x=0.15,y=0.95,"CV Annual GM TN:TP\n(Percent)",adj=0,cex=0.75)
}
dev.off()

# png(filename=paste0(plot.path,"TrendMaps_Chla.png"),width=7.5,height=2.25,units="in",res=200,type="windows",bg="white")
par(family="serif",oma=c(0.25,0.25,0.25,0.25),mar=c(0.1,0.1,0.1,0.1),xpd=F)
layout(matrix(c(1:4),1,4,byrow=T),widths = c(1,0.4,1,0.5))
bbox.lims=bbox(region.mask)
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
  plot(0:1,0:1,ann=F,axes=F,type="n")
  int.bks.vals=format(round(int.bks,3),nsmall=3)
  labs=c(paste0("< ",int.bks.vals[2]),paste(int.bks.vals[2:3],int.bks.vals[3:4],sep=" - "),paste(paste0(">",int.bks.vals[4])))
  n.bks=length(int.bks)-1
  bx.val= seq(0.45,0.85,(0.85-0.45)/n.bks)
  rect(0.15,bx.val[1:n.bks],0.25,bx.val[2:(n.bks+1)],col=rev(pal),lty=0)
  text(x=0.25, y = bx.val[2:(n.bks+1)]-c(mean(diff(bx.val[2:(n.bks+1)]))/2), labels = rev(labs),cex=0.6,adj=0,pos=4)
  text(x=0.15,y=0.95,"Chl-a Thiel-Sen Slope\n(\u03BCg L\u207B\u00B9 Yr\u207B\u00B9)",adj=0,cex=0.75)
  legend(0.5,0.4,legend=c("Insufficent Data","\u03C1 > 0.05","\u03C1 < 0.05"),
         pch=21,lty=c(NA),lwd=c(0.1),
         col=c("grey"),pt.bg=c("grey",cols.val[2:3]),
         pt.cex=1,ncol=1,cex=0.6,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=1,title.adj = 0,title="Kendall's trend \u03C1-value")
  
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

# png(filename=paste0(plot.path,"TrendMaps_TOC.png"),width=7.5,height=2.25,units="in",res=200,type="windows",bg="white")
par(family="serif",oma=c(0.25,0.25,0.25,0.25),mar=c(0.1,0.1,0.1,0.1),xpd=F)
layout(matrix(c(1:4),1,4,byrow=T),widths = c(1,0.4,1,0.5))
bbox.lims=bbox(region.mask)
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
  mapmisc::scaleBar(utm17,"bottomright",bty="n",cex=1,seg.len=4)
  
  plot(0:1,0:1,ann=F,axes=F,type="n")
  int.bks.vals=format(round(int.bks,3),nsmall=3)
  labs=c(paste0("< ",int.bks.vals[2]),paste(int.bks.vals[2:3],int.bks.vals[3:4],sep=" - "),paste(paste0(">",int.bks.vals[4])))
  n.bks=length(int.bks)-1
  bx.val= seq(0.45,0.85,(0.85-0.45)/n.bks)
  rect(0.15,bx.val[1:n.bks],0.25,bx.val[2:(n.bks+1)],col=rev(pal),lty=0)
  text(x=0.25, y = bx.val[2:(n.bks+1)]-c(mean(diff(bx.val[2:(n.bks+1)]))/2), labels = rev(labs),cex=0.6,adj=0,pos=4)
  text(x=0.15,y=0.95,"TOC Thiel-Sen Slope\n(mg C L\u207B\u00B9 Yr\u207B\u00B9)",adj=0,cex=0.75)
  legend(0.5,0.4,legend=c("Insufficent Data","\u03C1 > 0.05","\u03C1 < 0.05"),
         pch=21,lty=c(NA),lwd=c(0.1),
         col=c("grey"),pt.bg=c("grey",cols.val[2:3]),
         pt.cex=1,ncol=1,cex=0.6,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=1,title.adj = 0,title="Kendall's trend \u03C1-value")
  
  plot(shore,col="cornsilk",border="grey",bg="lightblue",ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.05,xpd=F)
  int=classIntervals(values(tps.TOC.GM)[is.na(values(tps.TOC.GM))==F],n=n,style="equal")
  int.bks=int$brks
  pal=hcl.colors(length(int$brks)-1, "Inferno", rev = F,alpha=0.7)
  image(tps.TOC.GM,add=T,breaks=int.bks,col = pal)
  plot(ENP,add=T,bg=NA,lwd=0.5)
  plot(regions2,lty=2,add=T,border="grey80",lwd=0.5)
  plot(sites.shp.TOC.GM,add=T,pch=21,cex=0.25,bg=adjustcolor("white",0.5),col=NA)
  box(lwd=1)
  plot(0:1,0:1,ann=F,axes=F,type="n")
  legend_image=as.raster(matrix(pal,ncol=1))
  text(x=0.25, y = c(0.83,0.47), labels = c("< 0.2",">20"),cex=0.6,adj=0,pos=4)
  rasterImage(legend_image,0.15,0.45,0.25,0.85)
  text(x=0.15,y=0.95,"Average Annual GM TOC\n(mg C L\u207B\u00B9)",adj=0,cex=0.75)
}
dev.off()