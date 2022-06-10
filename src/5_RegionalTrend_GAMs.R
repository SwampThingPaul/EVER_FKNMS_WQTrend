
## Run 4_RegionalTrend_DisturbanceClimate.R first

# GAM ---------------------------------------------------------------------
# some ideas https://raw.githack.com/eco4cast/Statistical-Methods-Seminar-Series/main/bolker_mixedmodels/outputs/full_notes.html

gc()
# TN ----------------------------------------------------------------------
dat.all.TN.GM=merge(subset(dat.all.GM,variable=="TN"&WY%in%WYs),sites.shp,"STATION")
head(dat.all.TN.GM)
WY.k=23
loc.k=200
m.TN<-bam(log(GM)~
            s(WY,bs="cr",k=WY.k)+
            s(UTMX,UTMY,bs="ds",k=loc.k,m=c(1,0.5))+
            ti(UTMX,UTMY,WY,d=c(2,1),bs=c("ds","cr"),k=c(100,24)),
          data=dat.all.TN.GM,
          nthreads = c(12),discrete=T)
summary(m.TN)
qq.gam(m.TN)

pred.org=predict(m.TN,type="terms")
partial.resids.TN<-pred.org+residuals(m.TN)

hist(partial.resids.TN[,1])
shapiro.test(partial.resids.TN[,1])

hist(partial.resids.TN[,2])
shapiro.test(partial.resids.TN[,2])

hist(partial.resids.TN[,3])
shapiro.test(partial.resids.TN[,3])

nvar=3;layout(matrix(1:nvar,1,nvar))
plot(m.TN,residuals=T,pch=21)
plot(m.TN,residuals=T,pch=21,select=1,shade=T)
dev.off()
plot(m.TN,select=3)

nvar=4;layout(matrix(1:nvar,1,nvar))
gam.check(m.TN)
dev.off()

m.TN.sum=notidy_tidy_gam(m.TN)
m.TN.est=notidy_glance_gam(m.TN)
# save(m.TN,file=paste0(export.path,"TN_GAM.Rdata"))
# load(paste0(export.path,"TN_GAM.Rdata"))
# write.csv(m.TN.sum,paste0(export.path,"TN_gam_mod_sum.csv"),row.names=F)
# write.csv(m.TN.est,paste0(export.path,"TN_gam_mod_est.csv"),row.names=F)
m.TN.sum=read.csv(paste0(export.path,"TN_gam_mod_sum.csv"))
m.TN.est=read.csv(paste0(export.path,"TN_gam_mod_est.csv"))

library(flextable)
library(magrittr)
notidy_as_flextable_gam(x=NULL,data_t=m.TN.sum,data_g=m.TN.est,dig.num=2)%>%
  font(fontname="Times New Roman",part="all") %>%print("docx")
# as_flextable(m.TN)

reg.ext=extent(region.mask)
pdat.sp<-data.frame(expand.grid(
                UTMX=seq(reg.ext[1],reg.ext[2],by=4000),
                UTMY=seq(reg.ext[3],reg.ext[4],by=4000),
                WY=seq(1996,2019,0.25)
              ))
diff(unique(pdat.sp$UTMX))
diff(unique(pdat.sp$UTMY))

pred.mod=predict(m.TN,newdata=pdat.sp,type="terms",se=T,newdata.guaranteed=T)

pdat.TN=pdat.sp
pdat.TN$fit.WY=pred.mod$fit[,1]
pdat.TN$se.WY=pred.mod$se[,1]
pdat.TN$fit.UTM=pred.mod$fit[,2]
pdat.TN$se.UTM=pred.mod$se[,2]
pdat.TN$fit.UTMWY=pred.mod$fit[,3]
pdat.TN$se.UTMWY=pred.mod$se[,3]

df.res <- df.residual(m.TN)
crit.t <- qt(0.025, df.res, lower.tail = FALSE)
# crit.t <- qnorm((1 - 0.95) / 2, lower.tail = FALSE)
pdat.TN <- transform(pdat.TN,
                      WY.UCI = fit.WY + (crit.t * se.WY),
                      WY.LCI = fit.WY - (crit.t * se.WY),
                      UTM.UCI = fit.UTM + (crit.t * se.UTM),
                      UTM.LCI = fit.UTM - (crit.t * se.UTM),
                      UTMWY.UCI = fit.UTMWY + (crit.t * se.UTMWY),
                      UTMWY.LCI = fit.UTMWY - (crit.t * se.UTMWY))


tmp.ma=with(pdat.TN,matrix(fit.UTM,nrow=length(unique(UTMX)),ncol=length(unique(UTMY))))
dat1=list()
dat1$x=unique(pdat.TN$UTMX)
dat1$y=unique(pdat.TN$UTMY)
dat1$z=tmp.ma
r=raster(dat1)
TN.GAM.UTM <- mask(r, region.mask)# mask(r, gBuffer(region.mask,width=2000))
plot(TN.GAM.UTM)
plot(r) ## FIXED? 


par(family="serif",mar=c(2,3,1,1.5),oma=c(2,1.5,0.25,0.25));
layout(matrix(c(1:3),1,3),widths=c(1,1,0.4))

ylim.val=c(-1.25,1.25);by.y=0.5;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(1995,2019);by.x=10;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/by.x)
# pdat.TSS=pdat.TSS[order(pdat.TSS$WY),]
plot(fit.WY~WY,pdat.TN,ylim=ylim.val,xlim=xlim.val,ann=F,axes=F,type="n")
abline(h=ymaj,v=xmaj,lty=3,col="grey")
abline(h=0)
points(jitter(m.TN$model$WY,0.5),partial.resids.TN[,1],pch=19,col=adjustcolor("dodgerblue1",0.5))
with(pdat.TN,shaded.range(WY,WY.LCI,WY.UCI,"grey",lty=1))
lines(fit.WY~WY,pdat.TN,lwd=2)
axis_fun(1,xmaj,xmin,xmaj,line=-0.5);axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=3,adj=0,"Total Nitrogen")
mtext(side=3,adj=1,"s(WY)")
mtext(side=1,line=2,"WY")
mtext(side=2,line=2,"Effect")

# par(mar=c(0.1,0.1,0.1,0.1))
# par(mar=c(2,0.1,1,0.1))
bbox.lims=bbox(region.mask)
b2=seq(-2,2,0.1)
pal2=colorRampPalette(c("blue","grey90","red"))(length(b2)-1)
plot(ENP,ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],border=NA)
image(TN.GAM.UTM,add=T,breaks=b2,col = pal2)
contour(TN.GAM.UTM,breaks = b2, add = TRUE,col="white",lwd=0.4,labcex=0.5)
plot(shore,add=T,lwd=0.1)
plot(ENP,add=T)
plot(sites.shp2,pch=19,col=adjustcolor("black",0.5),cex=0.5,add=T)
mapmisc::scaleBar(utm17,"bottomright",bty="n",cex=1,seg.len=4,outer=F)
box(lwd=1)
mtext(side=3,adj=0,"ti(UTMX,UTMY)")

plot(0:1,0:1,ann=F,axes=F,type="n")
# b2=format(b2)
l.b=length(b2)
labs=c(paste0("< ",b2[2]),paste(b2[2:(l.b-2)],b2[3:(l.b-1)],sep=" - "),paste(paste0(">",b2[(l.b-1)])))
n.bks=length(b2)-1
top.val=0.8
bot.val=0.2
mid.v.val=bot.val+(top.val-bot.val)/2
x.max=0.3
x.min=0.1
mid.val=x.min+(x.max-x.min)/2
txt.offset.val=-0.01
legend_image=as.raster(matrix(rev(pal2),ncol=1))
rasterImage(legend_image,x.min,bot.val,x.max,top.val)
text(x=x.max, y = c(bot.val,mid.v.val,top.val), labels = format(c(min(b2),0,max(b2))),cex=0.75,adj=0,pos=4,offset=0.5)
# bx.val= seq(bot.val,top.val,(top.val-bot.val)/n.bks)
# rect(x.min,bx.val[1:n.bks],x.max,bx.val[2:(n.bks+1)],col=rev(pal2),lty=0)
# text(y=bx.val[2:(n.bks+1)]-c(mean(diff(bx.val[2:(n.bks+1)]))/2), x = x.max, labels = rev(labs),cex=0.75,xpd=NA,pos=4,adj=0)
text(x=mid.val,y=top.val,"s(UTMX,UTMY)\nEffect",adj=0,cex=0.8,pos=3,xpd=NA)
dev.off()

## supplemental material
# png(filename=paste0(plot.path,"GAM_TN_UTMWY.png"),width=12,height=7,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(0.1,0.1,0.1,0.1),oma=c(0.2,0.2,0.2,0.2))
layout(matrix(c(1:25),5,5,byrow=T))

b2=seq(-1,1,0.1)
pal2=colorRampPalette(c("blue","grey90","red"))(length(b2)-1)
WY.vals=seq(1996,2019,1)
for(i in 1:length(WY.vals)){
  tmp.ma=with(subset(pdat.TN,WY==WY.vals[i]),
              matrix(fit.UTMWY,nrow=length(unique(UTMX)),ncol=length(unique(UTMY))))
  
  dat1=list()
  dat1$x=unique(pdat.TN$UTMX)
  dat1$y=unique(pdat.TN$UTMY)
  dat1$z=tmp.ma
  r=raster(dat1)
  TN.GAM.UTMWY <- mask(r, region.mask)
  
  bbox.lims=bbox(region.mask)
  plot(ENP,ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],border=NA)
  image(TN.GAM.UTMWY,add=T,breaks=b2,col = pal2)
  contour(TN.GAM.UTMWY,breaks = b2, add = TRUE,col="white",lwd=0.4,labcex=0.5)
  plot(shore,add=T,lwd=0.1)
  plot(ENP,add=T)
  mtext(side=3,line=-1.5,adj=0,paste0(" WY",WY.vals[i]))
  box(lwd=1)
}
plot(0:1,0:1,ann=F,axes=F,type="n")
# b2=format(b2)
l.b=length(b2)
labs=c(paste0("< ",b2[2]),paste(b2[2:(l.b-2)],b2[3:(l.b-1)],sep=" - "),paste(paste0(">",b2[(l.b-1)])))
n.bks=length(b2)-1
top.val=0.7
bot.val=0.2
mid.v.val=bot.val+(top.val-bot.val)/2
x.max=0.50
x.min=0.25
mid.val=x.min+(x.max-x.min)/2
txt.offset.val=-0.01
legend_image=as.raster(matrix(rev(pal2),ncol=1))
rasterImage(legend_image,x.min,bot.val,x.max,top.val)
text(x=x.max, y = c(bot.val,mid.v.val,top.val), labels = format(c(min(b2),0,max(b2))),cex=0.75,adj=0,pos=4,offset=0.5)
# bx.val= seq(bot.val,top.val,(top.val-bot.val)/n.bks)
# rect(x.min,bx.val[1:n.bks],x.max,bx.val[2:(n.bks+1)],col=rev(pal2),lty=0)
# text(y=bx.val[2:(n.bks+1)]-c(mean(diff(bx.val[2:(n.bks+1)]))/2), x = x.max, labels = rev(labs),cex=0.75,xpd=NA,pos=4,adj=0)
text(x=mid.val,y=top.val,"ti(UTMX,UTMY,WY)\nEffect",adj=0,cex=0.8,pos=3,xpd=NA)
dev.off()
# TN period of change -----------------------------------------------------
## see https://fromthebottomoftheheap.net/2014/05/15/identifying-periods-of-change-with-gams/
tmpf <- tempfile()
download.file("https://gist.github.com/gavinsimpson/e73f011fdaaab4bb5a30/raw/82118ee30c9ef1254795d2ec6d356a664cc138ab/Deriv.R",
              tmpf)
source(tmpf)

sites.shp.tmp=subset(sites.shp,STATION%in%unique(subset(dat.all.GM,variable=="TN")$STATION))

# pdat.TN2=reshape::expand.grid.df(sites.shp.tmp,data.frame(WY=seq(1996,2019,1)))
pdat.TN2=reshape::expand.grid.df(sites.shp,data.frame(WY=seq(1996,2019,0.25)))
p2 <- predict(m.TN, newdata=pdat.TN2,type = "terms", se.fit = TRUE)
pdat.TN2$p2=p2$fit[,1]
pdat.TN2$se2 = p2$se.fit[,1]
pdat.TN2$loc.WY.p=p2$fit[,3]
pdat.TN2$loc.WY.se = p2$se.fit[,3]

df.res <- df.residual(m.TN)
crit.t <- qt(0.025, df.res, lower.tail = FALSE)
pdat.TN2 <- transform(pdat.TN2,
                     upper = p2 + (crit.t * se2),
                     lower = p2 - (crit.t * se2))
pdat.TN2=pdat.TN2[order(pdat.TN2$WY,pdat.TN2$STATION),]
gc()

m.TN.d <- derivatives(m.TN,newdata=pdat.TN2,
                      term='s(WY)',type = "central",interval="confidence",ncores=12)

m.TN.dsig <- signifD(pdat.TN2$p2,
                     d=m.TN.d$derivative,
                     m.TN.d$upper,m.TN.d$lower)
pdat.TN2$dsig.incr=unlist(m.TN.dsig$incr)
pdat.TN2$dsig.decr=unlist(m.TN.dsig$decr)

unique(subset(pdat.TN2,is.na(dsig.incr)==F)$WY)
unique(subset(pdat.TN2,is.na(dsig.decr)==F)$WY)



# TP ----------------------------------------------------------------------
dat.all.TP.GM=merge(subset(dat.all.GM,variable=="TP"&WY%in%WYs),sites.shp,"STATION")
head(dat.all.TP.GM)
gc()

WY.k=23
loc.k=500
m.TP<-bam(log(GM)~
            s(WY,bs="cr",k=WY.k)+
            s(UTMX,UTMY,bs="ds",k=loc.k,m=c(1,0.5))+
            ti(UTMX,UTMY,WY,d=c(2,1),bs=c("ds","cr"),k=c(90,24)),
          data=dat.all.TP.GM,
          nthreads = 12,discrete=T)

summary(m.TP)
qq.gam(m.TP)
nvar=3;layout(matrix(1:nvar,1,nvar))
plot(m.TP,residuals=T,pch=21)

nvar=4;layout(matrix(1:nvar,1,nvar))
gam.check(m.TP)
dev.off()

pred.org=predict(m.TP,type="terms")
partial.resids.TP<-pred.org+residuals(m.TP)

hist(partial.resids.TP[,1])
shapiro.test(partial.resids.TP[,1])

hist(partial.resids.TP[,2])
shapiro.test(partial.resids.TP[,2])

hist(partial.resids.TP[,3])
shapiro.test(partial.resids.TP[,3])



m.TP.sum=notidy_tidy_gam(m.TP)
m.TP.est=notidy_glance_gam(m.TP)
# save(m.TP,file=paste0(export.path,"TP_GAM.Rdata"))
# load(paste0(export.path,"TP_GAM.Rdata"))
# write.csv(m.TP.sum,paste0(export.path,"TP_gam_mod_sum.csv"),row.names=F)
# write.csv(m.TP.est,paste0(export.path,"TP_gam_mod_est.csv"),row.names=F)
m.TP.sum=read.csv(paste0(export.path,"TP_gam_mod_sum.csv"))
m.TP.est=read.csv(paste0(export.path,"TP_gam_mod_est.csv"))

notidy_as_flextable_gam(x=NULL,data_t=m.TP.sum,data_g=m.TP.est,dig.num=2)%>%
  font(fontname="Times New Roman",part="all")# %>%print("docx")

pred.mod=predict(m.TP,newdata=pdat.sp,type="terms",se=T,newdata.guaranteed=T)

pdat.TP=pdat.sp
pdat.TP$fit.WY=pred.mod$fit[,1]
pdat.TP$se.WY=pred.mod$se[,1]
pdat.TP$fit.UTM=pred.mod$fit[,2]
pdat.TP$se.UTM=pred.mod$se[,2]
pdat.TP$fit.UTMWY=pred.mod$fit[,3]
pdat.TP$se.UTMWY=pred.mod$se[,3]

df.res <- df.residual(m.TP)
crit.t <- qt(0.025, df.res, lower.tail = FALSE)
# crit.t <- qnorm((1 - 0.95) / 2, lower.tail = FALSE)
pdat.TP <- transform(pdat.TP,
                     WY.UCI = fit.WY + (crit.t * se.WY),
                     WY.LCI = fit.WY - (crit.t * se.WY),
                     UTM.UCI = fit.UTM + (crit.t * se.UTM),
                     UTM.LCI = fit.UTM - (crit.t * se.UTM),
                     UTMWY.UCI = fit.UTMWY + (crit.t * se.UTMWY),
                     UTMWY.LCI = fit.UTMWY - (crit.t * se.UTMWY))


tmp.ma=with(pdat.TP,matrix(fit.UTM,nrow=length(unique(UTMX)),ncol=length(unique(UTMY))))
dat1=list()
dat1$x=unique(pdat.TP$UTMX)
dat1$y=unique(pdat.TP$UTMY)
dat1$z=tmp.ma
r=raster(dat1)
TP.GAM.UTM <- mask(r, region.mask)# mask(r, gBuffer(region.mask,width=2000))
plot(TP.GAM.UTM)
plot(r) ## FIXED? 

par(family="serif",mar=c(2,3,1,1.5),oma=c(2,1.5,0.25,0.25));
layout(matrix(c(1:3),1,3),widths=c(1,1,0.4))

ylim.val=c(-1.25,1.25);by.y=0.5;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(1995,2019);by.x=10;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/by.x)
# pdat.TSS=pdat.TSS[order(pdat.TSS$WY),]
plot(fit.WY~WY,pdat.TP,ylim=ylim.val,xlim=xlim.val,ann=F,axes=F,type="n")
abline(h=ymaj,v=xmaj,lty=3,col="grey")
abline(h=0)
points(jitter(m.TP$model$WY,0.5),partial.resids.TP[,1],pch=19,col=adjustcolor("dodgerblue1",0.5))
with(pdat.TP,shaded.range(WY,WY.LCI,WY.UCI,"grey",lty=1))
lines(fit.WY~WY,pdat.TP,lwd=2)
axis_fun(1,xmaj,xmin,xmaj,line=-0.5);axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=3,adj=0,"Total Phosphorus")
mtext(side=3,adj=1,"s(WY)")
mtext(side=1,line=2,"WY")
mtext(side=2,line=2,"Effect")

# par(mar=c(0.1,0.1,0.1,0.1))
# par(mar=c(2,0.1,1,0.1))
bbox.lims=bbox(region.mask)
b2=seq(-2,2,0.1)
pal2=colorRampPalette(c("blue","grey90","red"))(length(b2)-1)
plot(ENP,ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],border=NA)
image(TP.GAM.UTM,add=T,breaks=b2,col = pal2)
contour(TP.GAM.UTM,breaks = b2, add = TRUE,col="white",lwd=0.4,labcex=0.5)
plot(shore,add=T,lwd=0.1)
plot(ENP,add=T)
plot(sites.shp2,pch=19,col=adjustcolor("black",0.5),cex=0.5,add=T)
mapmisc::scaleBar(utm17,"bottomright",bty="n",cex=1,seg.len=4,outer=F)
box(lwd=1)
mtext(side=3,adj=0,"ti(UTMX,UTMY)")

plot(0:1,0:1,ann=F,axes=F,type="n")
# b2=format(b2)
l.b=length(b2)
labs=c(paste0("< ",b2[2]),paste(b2[2:(l.b-2)],b2[3:(l.b-1)],sep=" - "),paste(paste0(">",b2[(l.b-1)])))
n.bks=length(b2)-1
top.val=0.8
bot.val=0.2
mid.v.val=bot.val+(top.val-bot.val)/2
x.max=0.3
x.min=0.1
mid.val=x.min+(x.max-x.min)/2
txt.offset.val=-0.01
legend_image=as.raster(matrix(rev(pal2),ncol=1))
rasterImage(legend_image,x.min,bot.val,x.max,top.val)
text(x=x.max, y = c(bot.val,mid.v.val,top.val), labels = format(c(min(b2),0,max(b2))),cex=0.75,adj=0,pos=4,offset=0.5)
# bx.val= seq(bot.val,top.val,(top.val-bot.val)/n.bks)
# rect(x.min,bx.val[1:n.bks],x.max,bx.val[2:(n.bks+1)],col=rev(pal2),lty=0)
# text(y=bx.val[2:(n.bks+1)]-c(mean(diff(bx.val[2:(n.bks+1)]))/2), x = x.max, labels = rev(labs),cex=0.75,xpd=NA,pos=4,adj=0)
text(x=mid.val,y=top.val,"s(UTMX,UTMY)\nEffect",adj=0,cex=0.8,pos=3,xpd=NA)
dev.off()

## supplemental material
# png(filename=paste0(plot.path,"GAM_TP_UTMWY.png"),width=12,height=7,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(0.1,0.1,0.1,0.1),oma=c(0.2,0.2,0.2,0.2))
layout(matrix(c(1:25),5,5,byrow=T))

b2=seq(-0.8,0.8,0.1)
pal2=colorRampPalette(c("blue","grey90","red"))(length(b2)-1)
WY.vals=seq(1996,2019,1)
for(i in 1:length(WY.vals)){
  tmp.ma=with(subset(pdat.TP,WY==WY.vals[i]),
              matrix(fit.UTMWY,nrow=length(unique(UTMX)),ncol=length(unique(UTMY))))
  
  dat1=list()
  dat1$x=unique(pdat.TP$UTMX)
  dat1$y=unique(pdat.TP$UTMY)
  dat1$z=tmp.ma
  r=raster(dat1)
  TP.GAM.UTMWY <- mask(r, region.mask)
  
  bbox.lims=bbox(region.mask)
  plot(ENP,ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],border=NA)
  image(TP.GAM.UTMWY,add=T,breaks=b2,col = pal2)
  contour(TP.GAM.UTMWY,breaks = b2, add = TRUE,col="white",lwd=0.4,labcex=0.5)
  plot(shore,add=T,lwd=0.1)
  plot(ENP,add=T)
  mtext(side=3,line=-1.5,adj=0,paste0(" WY",WY.vals[i]))
  box(lwd=1)
}
plot(0:1,0:1,ann=F,axes=F,type="n")
# b2=format(b2)
l.b=length(b2)
labs=c(paste0("< ",b2[2]),paste(b2[2:(l.b-2)],b2[3:(l.b-1)],sep=" - "),paste(paste0(">",b2[(l.b-1)])))
n.bks=length(b2)-1
top.val=0.7
bot.val=0.2
mid.v.val=bot.val+(top.val-bot.val)/2
x.max=0.50
x.min=0.25
mid.val=x.min+(x.max-x.min)/2
txt.offset.val=-0.01
legend_image=as.raster(matrix(rev(pal2),ncol=1))
rasterImage(legend_image,x.min,bot.val,x.max,top.val)
text(x=x.max, y = c(bot.val,mid.v.val,top.val), labels = format(c(min(b2),0,max(b2))),cex=0.75,adj=0,pos=4,offset=0.5)
# bx.val= seq(bot.val,top.val,(top.val-bot.val)/n.bks)
# rect(x.min,bx.val[1:n.bks],x.max,bx.val[2:(n.bks+1)],col=rev(pal2),lty=0)
# text(y=bx.val[2:(n.bks+1)]-c(mean(diff(bx.val[2:(n.bks+1)]))/2), x = x.max, labels = rev(labs),cex=0.75,xpd=NA,pos=4,adj=0)
text(x=mid.val,y=top.val,"ti(UTMX,UTMY,WY)\nEffect",adj=0,cex=0.8,pos=3,xpd=NA)
dev.off()


# TP Period of change -----------------------------------------------------

pdat.TP2=reshape::expand.grid.df(sites.shp,data.frame(WY=seq(1996,2019,0.25)))
p2 <- predict(m.TP, newdata=pdat.TP2,type = "terms", se.fit = TRUE)
pdat.TP2$p2=p2$fit[,1]
pdat.TP2$se2 = p2$se.fit[,1]
pdat.TP2$loc.WY.p=p2$fit[,3]
pdat.TP2$loc.WY.se = p2$se.fit[,3]

df.res <- df.residual(m.TP)
crit.t <- qt(0.025, df.res, lower.tail = FALSE)
pdat.TP2 <- transform(pdat.TP2,
                      upper = p2 + (crit.t * se2),
                      lower = p2 - (crit.t * se2))
pdat.TP2=pdat.TP2[order(pdat.TP2$WY,pdat.TP2$STATION),]

m.TP.d <- derivatives(m.TP,newdata=pdat.TP2,
                      term='s(WY)',type = "central",interval="confidence",ncores=12)

m.TP.dsig <- signifD(pdat.TP2$p2,
                     d=m.TP.d$derivative,
                     m.TP.d$upper,m.TP.d$lower)
pdat.TP2$dsig.incr=unlist(m.TP.dsig$incr)
pdat.TP2$dsig.decr=unlist(m.TP.dsig$decr)

unique(subset(pdat.TP2,is.na(dsig.incr)==F)$WY)
unique(subset(pdat.TP2,is.na(dsig.decr)==F)$WY)



# GAM Animation -----------------------------------------------------------
# reg.ext=extent(region.mask)
# pdat.sp<-data.frame(expand.grid(
#   UTMX=seq(reg.ext[1],reg.ext[2],by=4000),
#   UTMY=seq(reg.ext[3],reg.ext[4],by=4000),
#   WY=seq(1996,2019,1)
# ))

fit.TP <- predict(m.TP, pdat.sp)
pred.TP <- cbind(pdat.sp, Fitted = exp(fit.TP))

yrs=seq(1996,2019,1)
for(i in 1:length(yrs)){
  tmp=subset(pred.TP,WY==yrs[i])[,c("Fitted","UTMX","UTMY")]
  coordinates(tmp)<-~UTMX + UTMY
  gridded(tmp)<-TRUE
  # rasterDF<-raster::raster(tmp,layer=1,values=T)
  tmp=as(tmp,"RasterLayer")
  proj4string(tmp)<-utm17
  tmp.m=raster::mask(tmp,region.mask)
  assign(paste0("GAM.TP.",yrs[i]),tmp.m)
  print(i)
}

GAM.TP.stack2=stack(GAM.TP.1996,
                    GAM.TP.1997,
                    GAM.TP.1998,
                    GAM.TP.1999,
                    GAM.TP.2000,
                    GAM.TP.2001,
                    GAM.TP.2002,
                    GAM.TP.2003,
                    GAM.TP.2004,
                    GAM.TP.2005,
                    GAM.TP.2006,
                    GAM.TP.2007,
                    GAM.TP.2008,
                    GAM.TP.2009,
                    GAM.TP.2010,
                    GAM.TP.2011,
                    GAM.TP.2012,
                    GAM.TP.2013,
                    GAM.TP.2014,
                    GAM.TP.2015,
                    GAM.TP.2016,
                    GAM.TP.2017,
                    GAM.TP.2018,
                    GAM.TP.2019)
plot(GAM.TP.stack2)
bbox.lims=bbox(region.mask)
GAM.TP.map=tm_shape(shore2,bbox=bbox.lims)+tm_fill(col="cornsilk")+
  tm_shape(GAM.TP.stack2)+
  tm_raster(title="",palette="-viridis",
            breaks=c(0,5,10,15,20,30,40,50,61),
            labels=c("< 5","5 - 10","10 - 15","15 - 20","20 - 30","30 - 40","40 - 50","50 - 60"))+
  #tm_raster(title="",palette="-cividis",alpha=0.75,style="cont")+
  tm_shape(shore2)+
  tm_borders(col="grey30",lwd=0.1)+
  # tm_polygons(col="cornsilk",border.col="grey")+
  tm_facets(free.scales=FALSE,nrow=1,ncol=1)+
  tm_layout(panel.labels=paste0("WY",c(1996:2019)),
            fontfamily = "serif",bg.color="lightblue", 
            panel.label.size=0.8,legend.title.size=0.8,legend.text.size=0.6)+
  tm_legend(title="Annual GM TP\n(\u03BCg L\u207B\u00B9)",legend.outside=T)
GAM.TP.map
# tmap_animation(GAM.TP.map,filename="./Plots/TP_GAM.gif",delay=100,width=550,height=250,loop=T,dpi=150)

fit.TN <- predict(m.TN, pdat.sp)
pred.TN <- cbind(pdat.sp, Fitted = exp(fit.TN))

yrs=1996:2019
for(i in 1:length(yrs)){
  tmp=subset(pred.TN,WY==yrs[i])[,c("Fitted","UTMX","UTMY")]
  coordinates(tmp)<-~UTMX + UTMY
  gridded(tmp)<-TRUE
  # rasterDF<-raster::raster(tmp,layer=1,values=T)
  tmp=as(tmp,"RasterLayer")
  proj4string(tmp)<-utm17
  tmp.m=raster::mask(tmp,region.mask)
  assign(paste0("GAM.TN.",yrs[i]),tmp.m)
  print(i)
}

GAM.TN.stack2=stack(GAM.TN.1996,
                    GAM.TN.1997,
                    GAM.TN.1998,
                    GAM.TN.1999,
                    GAM.TN.2000,
                    GAM.TN.2001,
                    GAM.TN.2002,
                    GAM.TN.2003,
                    GAM.TN.2004,
                    GAM.TN.2005,
                    GAM.TN.2006,
                    GAM.TN.2007,
                    GAM.TN.2008,
                    GAM.TN.2009,
                    GAM.TN.2010,
                    GAM.TN.2011,
                    GAM.TN.2012,
                    GAM.TN.2013,
                    GAM.TN.2014,
                    GAM.TN.2015,
                    GAM.TN.2016,
                    GAM.TN.2017,
                    GAM.TN.2018,
                    GAM.TN.2019)
plot(GAM.TN.stack2)
range(GAM.TN.stack2)
bbox.lims=bbox(region.mask)
GAM.TN.map=tm_shape(shore2,bbox=bbox.lims)+tm_fill(col="cornsilk")+
  tm_shape(GAM.TN.stack2)+
  tm_raster(title="",palette="-viridis",
            breaks=c(0,0.1,0.2,0.4,0.6,1,1.5,2),
            labels=c("< 0.1","0.1 - 0.2","0.2 - 0.4","0.4 - 0.6","0.6 - 1.0","1.0 - 1.5","1.5 - 2.0"))+
  #tm_raster(title="",palette="-cividis",alpha=0.75,style="cont")+
  tm_shape(shore2)+
  tm_borders(col="grey30",lwd=0.1)+
  # tm_polygons(col="cornsilk",border.col="grey")+
  tm_facets(free.scales=FALSE,nrow=1,ncol=1)+
  tm_layout(panel.labels=paste0("WY",c(1996:2019)),
            fontfamily = "serif",bg.color="lightblue", 
            panel.label.size=0.8,legend.title.size=0.8,legend.text.size=0.6)+
  tm_legend(title="Annual GM TN\n(mg L\u207B\u00B9)",legend.outside=T)
GAM.TN.map
# tmap_animation(GAM.TN.map,filename="./Plots/TN_GAM.gif",delay=100,width=550,height=250,loop=T,dpi=150)




# -------------------------------------------------------------------------

# png(filename=paste0(plot.path,"GAM_effect_TN_TP.png"),width=6.5,height=4.25,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(2,1,1,0.5),oma=c(2,3,0.25,0.1),xpd=F);
layout(matrix(1:6,2,3,byrow = T),widths=c(1,1,0.4))

ylim.val=c(-1.25,1.25);by.y=0.5;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(1995,2019);by.x=10;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/by.x)
plot(fit.WY~WY,pdat.TN,ylim=ylim.val,xlim=xlim.val,ann=F,axes=F,type="n")
abline(h=ymaj,v=xmaj,lty=3,col="grey")
abline(h=0)
points(jitter(m.TN$model$WY,0.5),partial.resids.TN[,1],pch=19,col=adjustcolor("dodgerblue1",0.1))
with(pdat.TN,shaded.range(WY,WY.LCI,WY.UCI,"grey",lty=1))
lines(fit.WY~WY,pdat.TN,lwd=2)
axis_fun(1,xmaj,xmin,xmaj,line=-0.5);axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=3,adj=0,"Total Nitrogen")
mtext(side=3,adj=1,"s(WY)")
# mtext(side=1,line=2,"WY")
mtext(side=2,line=2.75,"Effect")

bbox.lims=bbox(region.mask)
b2=seq(-2,2,0.1)
pal2=colorRampPalette(c("blue","grey90","red"))(length(b2)-1)
plot(ENP,ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],border=NA)
image(TN.GAM.UTM,add=T,breaks=b2,col = pal2)
contour(TN.GAM.UTM,breaks = b2, add = TRUE,col="white",lwd=0.4,labcex=0.5)
plot(shore,add=T,lwd=0.1)
plot(ENP,add=T)
plot(sites.shp2,pch=19,col=adjustcolor("black",0.1),cex=0.5,add=T)
mapmisc::scaleBar(utm17,"bottomright",bty="n",cex=1,seg.len=4,outer=F)
box(lwd=1)
mtext(side=3,adj=0,"ti(UTMX,UTMY)")

plot(0:1,0:1,ann=F,axes=F,type="n")
# b2=format(b2)
l.b=length(b2)
labs=c(paste0("< ",b2[2]),paste(b2[2:(l.b-2)],b2[3:(l.b-1)],sep=" - "),paste(paste0(">",b2[(l.b-1)])))
n.bks=length(b2)-1
top.val=0.8
bot.val=0.2
mid.v.val=bot.val+(top.val-bot.val)/2
x.max=0.3
x.min=0.1
mid.val=x.min+(x.max-x.min)/2
txt.offset.val=-0.01
legend_image=as.raster(matrix(rev(pal2),ncol=1))
rasterImage(legend_image,x.min,bot.val,x.max,top.val)
text(x=x.max, y = c(bot.val,mid.v.val,top.val), labels = format(c(min(b2),0,max(b2))),cex=0.75,adj=0,pos=4,offset=0.5)
# bx.val= seq(bot.val,top.val,(top.val-bot.val)/n.bks)
# rect(x.min,bx.val[1:n.bks],x.max,bx.val[2:(n.bks+1)],col=rev(pal2),lty=0)
# text(y=bx.val[2:(n.bks+1)]-c(mean(diff(bx.val[2:(n.bks+1)]))/2), x = x.max, labels = rev(labs),cex=0.75,xpd=NA,pos=4,adj=0)
text(x=mid.val,y=top.val,"ti(UTMX,UTMY)\nEffect",adj=0,cex=0.8,pos=3,xpd=NA)


plot(fit.WY~WY,pdat.TP,ylim=ylim.val,xlim=xlim.val,ann=F,axes=F,type="n")
abline(h=ymaj,v=xmaj,lty=3,col="grey")
abline(h=0)
points(jitter(m.TP$model$WY,0.5),partial.resids.TP[,1],pch=19,col=adjustcolor("dodgerblue1",0.1))
with(pdat.TP,shaded.range(WY,WY.LCI,WY.UCI,"grey",lty=1))
lines(fit.WY~WY,pdat.TP,lwd=2)
axis_fun(1,xmaj,xmin,xmaj,line=-0.5);axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=3,adj=0,"Total Phosphorus")
mtext(side=3,adj=1,"s(WY)")
mtext(side=1,line=2,"WY")
mtext(side=2,line=2.75,"Effect")

b2=seq(-2,2,0.1)
pal2=colorRampPalette(c("blue","grey90","red"))(length(b2)-1)
plot(ENP,ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],border=NA)
image(TP.GAM.UTM,add=T,breaks=b2,col = pal2)
contour(TP.GAM.UTM,breaks = b2, add = TRUE,col="white",lwd=0.4,labcex=0.5)
plot(shore,add=T,lwd=0.1)
plot(ENP,add=T)
plot(sites.shp2,pch=19,col=adjustcolor("black",0.1),cex=0.5,add=T)
mapmisc::scaleBar(utm17,"bottomright",bty="n",cex=1,seg.len=4,outer=F)
box(lwd=1)
mtext(side=3,adj=0,"ti(UTMX,UTMY)")

plot(0:1,0:1,ann=F,axes=F,type="n")
# b2=format(b2)
l.b=length(b2)
labs=c(paste0("< ",b2[2]),paste(b2[2:(l.b-2)],b2[3:(l.b-1)],sep=" - "),paste(paste0(">",b2[(l.b-1)])))
n.bks=length(b2)-1
top.val=0.8
bot.val=0.2
mid.v.val=bot.val+(top.val-bot.val)/2
x.max=0.3
x.min=0.1
mid.val=x.min+(x.max-x.min)/2
txt.offset.val=-0.01
legend_image=as.raster(matrix(rev(pal2),ncol=1))
rasterImage(legend_image,x.min,bot.val,x.max,top.val)
text(x=x.max, y = c(bot.val,mid.v.val,top.val), labels = format(c(min(b2),0,max(b2))),cex=0.75,adj=0,pos=4,offset=0.5)
# bx.val= seq(bot.val,top.val,(top.val-bot.val)/n.bks)
# rect(x.min,bx.val[1:n.bks],x.max,bx.val[2:(n.bks+1)],col=rev(pal2),lty=0)
# text(y=bx.val[2:(n.bks+1)]-c(mean(diff(bx.val[2:(n.bks+1)]))/2), x = x.max, labels = rev(labs),cex=0.75,xpd=NA,pos=4,adj=0)
text(x=mid.val,y=top.val,"ti(UTMX,UTMY)\nEffect",adj=0,cex=0.8,pos=3,xpd=NA)
dev.off()

## TN and TP significant change plots
date_decimal_WY=function(decimal,tz.val="EST"){
  # from lubridate::date_decimal()
  WY=trunc(decimal)
  # start=as.Date(paste(WY-1,5,1,sep="-"),tz=tz.val)
  start= lubridate::make_datetime(WY-1,5,1,tz=tz.val)
  # end=as.Date(paste(WY,4,30,sep="-"),tz=tz.val)
  end= lubridate::make_datetime(WY,4,30,tz=tz.val)
  frac <- decimal - WY
  seconds=as.numeric(difftime(end, start, units = "secs"))
  end <- start+seconds*frac
  end <- date.fun(end)
  return(end)
}
date_decimal_WY(1996.25)

# png(filename=paste0(plot.path,"GAM_sigchange.png"),width=4.5,height=5.5,units="in",res=200,type="windows",bg="white")
layout(matrix(1:5,5,1),heights=c(0.8,1,1,1,0.5))
par(family="serif",mar=c(1,1,0.25,0.5),oma=c(1,4.5,0.75,1.5));

# H. Wilma October-2005
# Mustang Corner Fire May-2008
# Cold Snaps January 2010, December 2011
# Droughts October 2010-April 2011, May 2015-October 2015
# Flood November 2015-March 2016 
# H. Irma September 2017
ylim.val=c(0,1.1)
xlim.val=date.fun(c("1995-05-01","2018-05-01"));xmaj=seq(xlim.val[1],xlim.val[2],"5 years");xmin=seq(xlim.val[1],xlim.val[2],"1 years")
event.lab.cex=0.65
plot(0:1~dates,type="n",ann=F,axes=F,xlim=xlim.val,ylim=ylim.val)
arrows(x0=date.fun("2005-10-01"),y0=0,y1=0.15,code=1,lwd=2,length=0.05)
text(date.fun("2005-10-01"),0.1,"H. Wilma",pos=3,cex=event.lab.cex,srt=0)
arrows(x0=date.fun("2017-09-10"),y0=0,y1=0.15,code=1,lwd=2,length=0.05)
text(date.fun("2017-09-10"),0.1,"H. Irma",pos=3,cex=event.lab.cex,srt=0)
# arrows(x0=date.fun("2017-07-2017"),y0=0,y1=0.20,code=1,lwd=2,length=0.1)
# text(date.fun("2017-07-2017"),0.20,"TS Emily",pos=3,cex=0.5,srt=90)
# arrows(x0=date.fun("2017-10-28"),y0=0,y1=0.25,code=1,lwd=2,length=0.1)
# text(date.fun("2017-10-28"),0.25,"TS Philippee",pos=3,cex=0.5,srt=90)
arrows(x0=date.fun("2008-05-01"),y0=0,y1=0.2,code=1,lwd=2,length=0.05,col="indianred1")
text(date.fun("2008-05-01"),0.2,"Fire",pos=3,cex=event.lab.cex,offset=0.1,col="indianred1")
arrows(x0=date.fun("2010-01-01"),y0=0,y1=0.2,code=1,lwd=2,length=0.05,col="dodgerblue1")
text(date.fun("2010-01-01"),0.2,"Cold",pos=3,cex=event.lab.cex,offset=0.1,col="dodgerblue1")
arrows(x0=date.fun("2011-12-01"),y0=0,y1=0.2,code=1,lwd=2,length=0.05,col="dodgerblue1")
text(date.fun("2011-12-01"),0.2,"Cold",pos=3,cex=event.lab.cex,offset=0.1,col="dodgerblue1")
arrows(x0=date.fun("2015-12-01"),y0=0,y1=0.15,code=1,lwd=2,length=0.05,col="Forestgreen")
text(date.fun("2015-12-01"),0.15,"Seagrass\ndie-off",pos=3,cex=event.lab.cex,offset=0.1,col="Forestgreen")


with(subset(hur.year,Key%in%EVER_FKNMS.hurr$Key),
     points(date.start,rep(0.4,length(Key)),
            pch=21,bg=adjustcolor("grey",0.5),
            col=adjustcolor("grey",0.5),cex=1,lwd=0.1))
# with(hurr.yr.N,text(date.fun(paste(Year.start,"10-30",sep="-")),
#                     rep(0.4,length(N.val)),labels=N.val,cex=0.5))
with(all_fire_sub.sum,points(date.fun(paste(FireCalend,"10-30",sep="-")),
                             rep(0.6,length(FireCalend)),
                             pch=21,bg="grey",cex=per.area/3,lwd=0.1))

# text(xlim.val[1],0.1,"Pulse Events",xpd=NA)
xx=date.fun(c("2000-01-01","2004-12-01"));yy=c(0.75,0.75,0.85,0.85)
polygon(c(xx,rev(xx)),yy,col=adjustcolor("khaki",0.5),border="khaki")
xx=date.fun(c("2006-01-01","2007-12-01"))
polygon(c(xx,rev(xx)),yy,col=adjustcolor("khaki",0.5),border="khaki")
xx=date.fun(c("2010-10-01","2011-04-01"))
polygon(c(xx,rev(xx)),yy,col=adjustcolor("khaki",0.5),border="khaki")
xx=date.fun(c("2015-05-01","2015-10-01"));
polygon(c(xx,rev(xx)),yy,col=adjustcolor("khaki",0.5),border="khaki")
# text(xlim.val[1],0.4,"Drought",xpd=NA)
# https://apps.sfwmd.gov/sfwmd/SFER/2017_sfer_final/v1/chapters/v1_ch3a.pdf
xx=date.fun(c("2015-10-01","2016-01-01"));yy=c(0.95,0.95,1.05,1.05)
polygon(c(xx,rev(xx)),yy,col=adjustcolor("lightsteelblue",0.5),border="lightsteelblue")
# https://apps.sfwmd.gov/sfwmd/SFER/2019_sfer_final/v1/chapters/v1_ch3a.pdf
xx=date.fun(c("2017-06-05","2018-01-01"));yy=c(0.95,0.95,1.05,1.05)
polygon(c(xx,rev(xx)),yy,col=adjustcolor("lightsteelblue",0.5),border="lightsteelblue")
# text(xlim.val[1],0.6,"Flood",xpd=NA)
axis_fun(1,xmaj,xmin,NA)
axis_fun(2,
         c(0.1,0.4,0.6,0.8,1),
         c(0.1,0.4,0.6,0.8,1),
         c("Pulse\nEvents","Hurricanes","Fire","Drought","Flood"),cex=0.9)

box(lwd=1)
text(xlim.val[2]+lubridate::ddays(600),ylim.val[2],"A",xpd=NA,cex=1.5)

ylim.val=c(-2,3);by.y=2;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(mean.3mon~monCY.date,ONI.dat,type="n",ann=F,axes=F,xlim=xlim.val,ylim=ylim.val)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(subset(ONI.dat,is.na(mean.3mon)==F),shaded.range(monCY.date,rep(0,length(monCY.date)),ifelse(mean.3mon>0,mean.3mon,0),"indianred1",lty=1))
with(subset(ONI.dat,is.na(mean.3mon)==F),shaded.range(monCY.date,ifelse(mean.3mon<0,mean.3mon,0),rep(0,length(monCY.date)),"dodgerblue1",lty=1))
abline(h=c(0,-0.5,0.5),lty=c(1,2,2))
axis_fun(1,xmaj,xmin,NA)
axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
mtext(side=2,line=2,"ONI Index")
text(xlim.val[2]+lubridate::ddays(600),ylim.val[2],"B",xpd=NA,cex=1.5)
# mtext(side=3,adj=1,line=-1.25,"B",cex=1.5)

pdat.TN2$WY.monWY=date_decimal_WY(pdat.TN2$WY)
m.TN.tmp.WY=date.fun(paste((m.TN$model$WY)-1,"05-01",sep="-"))
tmp=m.TN.tmp.WY+lubridate::ddays(runif(length(m.TN.tmp.WY),-90,90))

# xlim.val=c(1996,2019);
xmaj.lab=WY(xmaj);xmin.lab=WY(xmin)
ylim.val=c(-1,1);by.y=1;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(p2 ~ WY.monWY, data = pdat.TN2,type="n",ylim=ylim.val,xlim=xlim.val,ann=F,axes=F)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
points(tmp,partial.resids.TN[,1],pch=21,bg=adjustcolor("grey",0.25),col=adjustcolor("grey",0.25))
lines(p2 ~ WY.monWY, data = pdat.TN2)
lines(upper ~ WY.monWY, data = pdat.TN2, lty = "dashed")
lines(lower ~ WY.monWY, data = pdat.TN2, lty = "dashed")
lines(dsig.incr ~ WY.monWY, data = pdat.TN2, col = "red", lwd = 3,lty=1)
lines(dsig.decr ~ WY.monWY, data = pdat.TN2, col = "blue", lwd = 3,lty=1)
abline(h=0)
axis_fun(1,xmaj,xmin,NA,line=-0.5)
axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
# mtext(side=1,line=2,"Water Year")
mtext(side=2,line=2.5,"Total Nitrogen\nWY Effect",cex=0.9)
# mtext(side=3,adj=0,"Total Nitrogen")
# legend("topleft",legend=c("Fitted (WY)","95% CI (WY)","Significant Change (Increasing)","Significant Change (Decreasing)","Predicted (ti(UTMX,UTMY,WY))"),
#        lty=c(1,2,1,1,0),col=c("black","black","red","blue",adjustcolor("grey",0.25)),
#        pch=c(NA,NA,NA,NA,21),pt.bg=c(NA,NA,NA,NA,adjustcolor("grey",0.25)),
#        pt.cex=1,ncol=3,cex=0.6,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=1)
# text(xlim.val[2]+1.5,ylim.val[2],"C",xpd=NA,cex=1.5)
text(xlim.val[2]+lubridate::ddays(600),ylim.val[2],"C",xpd=NA,cex=1.5)

pdat.TP2$WY.monWY=date_decimal_WY(pdat.TP2$WY)
m.TP.tmp.WY=date.fun(paste((m.TP$model$WY)-1,"05-01",sep="-"))
tmp=m.TP.tmp.WY+lubridate::ddays(runif(length(m.TP.tmp.WY),-90,90))

ylim.val=c(-1,1);by.y=1;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(p2 ~ WY.monWY, data = pdat.TP2,type="n",ylim=ylim.val,xlim=xlim.val,ann=F,axes=F)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
points(tmp,partial.resids.TP[,1],pch=21,bg=adjustcolor("grey",0.25),col=adjustcolor("grey",0.25))
lines(p2 ~ WY.monWY, data = pdat.TP2)
lines(upper ~ WY.monWY, data = pdat.TP2, lty = "dashed")
lines(lower ~ WY.monWY, data = pdat.TP2, lty = "dashed")
lines(dsig.incr ~ WY.monWY, data = pdat.TP2, col = "red", lwd = 3,lty=1)
lines(dsig.decr ~ WY.monWY, data = pdat.TP2, col = "blue", lwd = 3,lty=1)
abline(h=0)
axis_fun(1,xmaj,xmin,xmaj.lab,line=-0.5)
axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=1,line=2,"Water Year")
mtext(side=2,line=2.5,"Total Phosphorus\nWY Effect",cex=0.9)
# mtext(side=3,adj=0,"Total Phosphorus")
# text(xlim.val[2]+1.5,ylim.val[2],"D",xpd=NA,cex=1.5)
text(xlim.val[2]+lubridate::ddays(600),ylim.val[2],"D",xpd=NA,cex=1.5)

plot(0:1,0:1,ann=F,axes=F,type="n")
legend(0.5,0.5,legend=c("Fitted (WY)","95% CI (WY)","Significant Change (Increasing)","Significant Change (Decreasing)","Partial Residual"),
       lty=c(1,2,1,1,0),col=c("black","black","red","blue",adjustcolor("grey",0.25)),
       pch=c(NA,NA,NA,NA,21),pt.bg=c(NA,NA,NA,NA,adjustcolor("grey",0.25)),
       pt.cex=1,ncol=2,cex=0.8,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=1)

dev.off()