## Run 3_RegionalTrend_trend analysis.R first

# Pulse metrics -----------------------------------------------------------
## analysis idea from quantmod_sev42_ppt.R
dat.all2
dat.all3=subset(dat.all2,WY%in%seq(1995,2019,1))

dat.all2.melt

# WY.samples.count=ddply(dat.all2.melt,c("STATION","WY","variable"),summarise,N.samples=N.obs(value))
# WY.samples.count=merge(WY.samples.count,
#                        ddply(WY.samples.count,c("STATION","variable"),summarise,N.WY=N.obs(WY)),
#                        c("STATION","variable"),all.x=T)

WY.samples.count=ddply(subset(dat.all2.melt,sea.screen==1&consec.screen==1),c("STATION","WY","variable"),summarise,N.samples=N.obs(value))
WY.samples.count=ddply(WY.samples.count,c("STATION","variable"),summarise,N.WYs=N.obs(WY))

dat.all2.melt2=merge(dat.all2.melt,WY.samples.count,c("STATION","variable"))
dat.all2.melt2=subset(dat.all2.melt2,sea.screen==1&consec.screen==1&N.WYs>=4)

library(quantmod)


# TP ----------------------------------------------------------------------


dat.all2.melt2.TP=subset(dat.all2.melt2,variable=="TP")
station.vals=unique(dat.all2.melt2.TP$STATION)
pulse.metric.landscape=data.frame()
for(i in 1:length(station.vals)){
tmp=subset(dat.all2.melt2.TP,variable=="TP"&STATION== station.vals[i])
tmp=subset(tmp,is.na(value)==F)

threshold.peak=5
TPpeaks<-findPeaks(tmp$value, thresh=threshold.peak)

# plot(value~DATE,tmp)
# points(value~DATE,tmp[TPpeaks,],pch=21,bg="green",cex=2)

TPpeaks.data=tmp[TPpeaks-1,c("DATE","value")]

nYrs=as.numeric(diff(range(TPpeaks.data$DATE))/365)
samp.freq=length(tmp$value)/nYrs
#average peak magnitude
peak_mean<-mean(TPpeaks.data$value)

#peak standard deviation
peak_sd<-sd(TPpeaks.data$value)

#peak CV
peak_CV<-peak_sd/peak_mean
# peak_CV

# Standardized regression models for temporal change
# peak number vs. time: Has the number of peaks increased/decreased over time?
# get slope of (number of peaks per year for each year) vs. year (and p-value) to look for temporal change in number of peaks

tmp$peak=NA
tmp[TPpeaks-1,"peak"]=1

tmp.WY=ddply(tmp,c("WY"),summarise,Npeak=sum(peak,na.rm=T))
tmp.WY$scale.Npeak=scale(tmp.WY$Npeak)
tmp.WY$scale.WY=scale(tmp.WY$WY)

if(sum(tmp.WY$Npeak)<2){
  peak.number.slope=NA
  peak.number.p=NA
}else{
peak.number.lm<-lm(scale.Npeak~scale.WY,data=tmp.WY)
summary(peak.number.lm)

# kendall.rslt=with(tmp.WY,cor.test(Npeak,WY,method="kendall"))
lmsum.number<-summary(peak.number.lm)
peak.number.slope<-as.numeric(peak.number.lm$coefficients[2])
peak.number.slope
peak.number.p<-as.numeric(lmsum.number$coefficients[2,4])
peak.number.p
}
rslt.all=data.frame(
  STATION=station.vals[i],
  samp.freq.perYr=samp.freq,
  peak_mean=peak_mean,
  peak_sd=peak_sd,
  peak_CV=peak_CV,
  peak.number.slope=peak.number.slope,
  peak.number.p=peak.number.p
)
pulse.metric.landscape=rbind(pulse.metric.landscape,rslt.all)
print(i)
}

pulse.metric.landscape

pulse.metric.landscape$stat.sig=with(pulse.metric.landscape,ifelse(peak.number.p<0.05,"sig","not-sig"))
pulse.metric.landscape$stat.sig=with(pulse.metric.landscape,ifelse(is.na(stat.sig)==T,"Empty",as.character(stat.sig)))
pulse.metric.landscape$stat.sig=as.factor(pulse.metric.landscape$stat.sig)


pulse.metric.landscape.shp=merge(sites.shp2,pulse.metric.landscape,"STATION",all.y=T)

# Spline ------------------------------------------------------------------
#thin plate spline https://rspatial.org/raster/analysis/4-interpolation.html

region.buf.r=raster(region.mask)
res(region.buf.r)=1000

# Spatial Trend -----------------------------------------------------------
# Peak Mean
m=Tps(coordinates(pulse.metric.landscape.shp),pulse.metric.landscape.shp$peak_mean)
tps=interpolate(region.buf.r,m)
tps.peakmean.trend=mask(tps,region.mask)
plot(tps.peakmean.trend)

# Peak sd
m=Tps(coordinates(pulse.metric.landscape.shp),pulse.metric.landscape.shp$peak_sd)
tps=interpolate(region.buf.r,m)
tps.peaksd.trend=mask(tps,region.mask)
plot(tps.peaksd.trend)

# Peak CV
m=Tps(coordinates(pulse.metric.landscape.shp),pulse.metric.landscape.shp$peak_CV)
tps=interpolate(region.buf.r,m)
tps.peakCV.trend=mask(tps,region.mask)
plot(tps.peakCV.trend)

# Peak count slope
m=Tps(coordinates(pulse.metric.landscape.shp),pulse.metric.landscape.shp$peak.number.slope)
tps=interpolate(region.buf.r,m)
tps.slope.trend=mask(tps,region.mask)
plot(tps.slope.trend)


leg.fun=function(b,pal,leg.title,
                 top.val=0.8,bot.val=0.2,mid.v.val=NULL,
                 x.max=0.3,x.min=0.1,mid.val=NULL,
                 txt.offset.val=-0.01,txt.y=NULL,leg.txt=NULL,
                 txt.cex=0.75,txt.adj=0,txt.pos=4,txt.offset=0.5,
                 title.cex=0.8,title.pos=3,title.adj=0,
                 title.x=NULL,title.y=NULL,
                 leg.type=c("continuous","categorical"), ...){
  l.b=length(b)
  labs=c(paste0("< ",b[2]),paste(b[2:(l.b-2)],b[3:(l.b-1)],sep=" - "),paste(paste0(">",b[(l.b-1)])))
  n.bks=length(b)-1
  mid.v.val=if(is.null(mid.v.val)==T){bot.val+(top.val-bot.val)/2}else{mid.v.val}
  
  mid.val=if(is.null(mid.val)==T){x.min+(x.max-x.min)/2}else{mid.val}
  if(leg.type=="continuous"){
    legend_image=as.raster(matrix(rev(pal),ncol=1))
    rasterImage(legend_image,x.min,bot.val,x.max,top.val)
    txt.y=if(is.null(txt.y)==T){c(bot.val,top.val)}else(txt.y)
    leg.txt=if(is.null(leg.txt)==T){format(c(min(b),max(b)))}else(leg.txt)
    text(x=x.max, y = txt.y, labels =leg.txt,cex=txt.cex,adj=txt.adj,pos=txt.pos,offset=txt.offset, ...)
  }
  if(leg.type=="categorical"){
    bx.val= seq(bot.val,top.val,(top.val-bot.val)/n.bks)
    rect(x.min,bx.val[1:n.bks],x.max,bx.val[2:(n.bks+1)],col=rev(pal),lty=0)
    text(y=bx.val[2:(n.bks+1)]-c(mean(diff(bx.val[2:(n.bks+1)]))/2), x = x.max, 
         labels = rev(labs),cex=txt.cex,xpd=NA,pos=txt.pos,adj=txt.adj)
  }
  
  title.x=if(is.null(title.x)==T){mid.val}else{title.x}
  title.y=if(is.null(title.y)==T){top.val}else{title.y}
  text(x=title.x,y=title.y,leg.title,adj=title.adj,cex=title.cex,pos=title.pos,xpd=NA)
}




cols.val=c("grey","white","red")
# png(filename=paste0(plot.path,"PulseMetrics.png"),width=6,height=6.5,units="in",res=200,type="windows",bg="white")
par(family="serif",oma=c(0.25,0.25,0.25,0.25),mar=c(0.1,0.1,0.1,0.1),xpd=F)
layout(matrix(c(1:8),4,2,byrow=T),widths = c(1,0.4))
bbox.lims=bbox(region.mask)

plot(shore,col="cornsilk",border="grey",bg="lightblue",ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.05)
b=seq(10,50,5)
pal=viridis::plasma(length(b)-1,alpha = 0.8)# hcl.colors(length(b)-1, "Spectral")
image(tps.peakmean.trend,add=T,breaks=b,col = pal)
plot(ENP,add=T,bg=NA,lwd=0.5)
plot(regions2,lty=2,add=T,border="grey80",lwd=0.5)
# plot(pulse.metric.landscape.shp,add=T,cex=0.40,pch=21,bg=adjustcolor("grey",0.5),col="grey",lwd=0.1)
# plot(pulse.metric.landscape.shp,add=T,pch=21,cex=0.75,bg=adjustcolor(cols.val,0.5)[pulse.metric.landscape.shp$stat.sig],col=NA);
mapmisc::scaleBar(utm17,"bottomright",bty="n",cex=1,seg.len=4)
box(lwd=1)

plot(0:1,0:1,ann=F,axes=F,type="n")
leg.fun(b,
        pal,
        "TP Peak Mean\n(\u03BCg L\u207B\u00B9)",
        leg.txt=c(min(b),max(b)),
        leg.type="continuous",
        x.max=0.6,x.min=0.4)


plot(shore,col="cornsilk",border="grey",bg="lightblue",ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.05,xpd=F)
b=seq(0,25,2.5)
pal=viridis::plasma(length(b)-1,alpha = 0.8)# hcl.colors(length(b)-1, "Spectral")
image(tps.peaksd.trend,add=T,breaks=b,col = pal)
plot(ENP,add=T,bg=NA,lwd=0.5)
plot(regions2,lty=2,add=T,border="grey80",lwd=0.5)
# plot(pulse.metric.landscape.shp,add=T,cex=0.40,pch=21,bg=adjustcolor("grey",0.5),col="grey",lwd=0.1)
# plot(pulse.metric.landscape.shp,add=T,pch=21,cex=0.75,bg=adjustcolor(cols.val,0.5)[pulse.metric.landscape.shp$stat.sig],col=NA);
mapmisc::scaleBar(utm17,"bottomright",bty="n",cex=1,seg.len=4)
box(lwd=1)

plot(0:1,0:1,ann=F,axes=F,type="n")
leg.fun(b,
        pal,
        "TP Peak SD\n(\u03BCg L\u207B\u00B9)",
        leg.txt=c(min(b),max(b)),
        leg.type="continuous",
        x.max=0.6,x.min=0.4)

plot(shore,col="cornsilk",border="grey",bg="lightblue",ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.05,xpd=F)
b=seq(0,0.6,0.1)
pal=viridis::plasma(length(b)-1,alpha = 0.8)# hcl.colors(length(b)-1, "Spectral")
image(tps.peakCV.trend,add=T,breaks=b,col = pal)
plot(ENP,add=T,bg=NA,lwd=0.5)
plot(regions2,lty=2,add=T,border="grey80",lwd=0.5)
# plot(pulse.metric.landscape.shp,add=T,cex=0.40,pch=21,bg=adjustcolor("grey",0.5),col="grey",lwd=0.1)
# plot(pulse.metric.landscape.shp,add=T,pch=21,cex=0.75,bg=adjustcolor(cols.val,0.5)[pulse.metric.landscape.shp$stat.sig],col=NA);
mapmisc::scaleBar(utm17,"bottomright",bty="n",cex=1,seg.len=4)
box(lwd=1)

plot(0:1,0:1,ann=F,axes=F,type="n")
leg.fun(b*100,
        pal,
        "TP Peak CV\n(Percent)",
        leg.txt=c(min(b*100),max(b*100)),
        leg.type="continuous",
        x.max=0.6,x.min=0.4)

plot(shore,col="cornsilk",border="grey",bg="lightblue",ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.05,xpd=F)
b=round(seq(-0.5,0.25,0.05),1)
pal=c(colorRampPalette(c("blue","grey90"))(length(b[b<=0])-1),
       colorRampPalette(c("grey90",'red'))(length(b[b>0])))
# pal=viridis::plasma(length(b)-1,alpha = 0.8)# hcl.colors(length(b)-1, "Spectral")
image(tps.slope.trend,add=T,breaks=b,col = pal)
plot(ENP,add=T,bg=NA,lwd=0.5)
plot(regions2,lty=2,add=T,border="grey80",lwd=0.5)
# plot(pulse.metric.landscape.shp,add=T,cex=0.40,pch=21,bg=adjustcolor("grey",0.5),col="grey",lwd=0.1)
plot(pulse.metric.landscape.shp,add=T,pch=21,cex=0.75,bg=adjustcolor(cols.val,0.5)[pulse.metric.landscape.shp$stat.sig],col=NA);
mapmisc::scaleBar(utm17,"bottomright",bty="n",cex=1,seg.len=4)
box(lwd=1)

plot(0:1,0:1,ann=F,axes=F,type="n")
leg.fun(b,
        pal,
        "TP Peak\nTemporal Change (Yr\u207B\u00b9)",
        leg.txt=c(min(b),max(b)),
        leg.type="continuous",
        x.max=0.6,x.min=0.4)
dev.off()