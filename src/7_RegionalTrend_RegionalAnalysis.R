
## Run 5_RegionalTrend_Maps.R first

# Regional Summaries ------------------------------------------------------
dev.off()
library(dunn.test)
library(rcompanion)
###
## TN
levels.var=c("ENP","Coastal_Mangroves","FLBay","Shelf","Keys")
levels.var.labs=c("ENP","Mangrove Fringe","Florida Bay","West Florida Shelf","Keys")
dat.all.TN.GM2=merge(subset(dat.all.GM,variable=="TN"),all.sites.region,"STATION",all.x=T)
dat.all.TN.GM2$Region=factor(dat.all.TN.GM2$Region,levels=levels.var)

## DIN
dat.all.DIN.GM2=merge(subset(dat.all.GM,variable=="DIN"),all.sites.region,"STATION",all.x=T)
dat.all.DIN.GM2$Region=factor(dat.all.DIN.GM2$Region,levels=levels.var)

## TP
dat.all.TP.GM2=merge(subset(dat.all.GM,variable=="TP"),all.sites.region,"STATION",all.x=T)
dat.all.TP.GM2$Region=factor(dat.all.TP.GM2$Region,levels=c("ENP","Coastal_Mangroves","FLBay","Shelf","Keys"))

## SRP
dat.all.SRP.GM2=merge(subset(dat.all.GM,variable=="SRP"),all.sites.region,"STATION",all.x=T)
dat.all.SRP.GM2$Region=factor(dat.all.SRP.GM2$Region,levels=c("ENP","Coastal_Mangroves","FLBay","Shelf","Keys"))

## Chl-a
dat.all.Chla.GM2=merge(subset(dat.all.GM,variable=="Chla"),all.sites.region,"STATION",all.x=T)
unique(subset(dat.all.Chla.GM2,is.na(Region)==T)$STATION)
dat.all.Chla.GM2$Region=factor(dat.all.Chla.GM2$Region,levels=c("ENP","Coastal_Mangroves","FLBay","Shelf","Keys"))
boxplot(GM~Region,dat.all.Chla.GM2,outline=F)

## TOC
dat.all.TOC.GM2=merge(subset(dat.all.GM,variable=="TOC"),all.sites.region,"STATION",all.x=T)
dat.all.TOC.GM2$Region=factor(dat.all.TOC.GM2$Region,levels=c("ENP","Coastal_Mangroves","FLBay","Shelf","Keys"))
boxplot(GM~Region,dat.all.TOC.GM2,outline=F)


## TNTP
dat.all.TNTP.GM2=merge(subset(dat.all.GM,variable=="TN_TP"),all.sites.region,"STATION",all.x=T)
dat.all.TNTP.GM2$Region=factor(dat.all.TNTP.GM2$Region,levels=c("ENP","Coastal_Mangroves","FLBay","Shelf","Keys"))
boxplot(GM~Region,dat.all.TNTP.GM2,outline=F);abline(h=20)

vars=c("STATION","WY","Region","param")
dat.all.GM$param=dat.all.GM$variable
region.para.melt=melt(
  merge(dat.all.GM,all.sites.region,"STATION",all.x=T)[,c(vars,"GM")],
  id.vars = vars)
head(region.para.melt)

region.para.melt=merge(region.para.melt,data.frame(Region=c("ENP","Coastal_Mangroves","FLBay","Shelf","Keys"),
                                                   Region.txt=c("ENP"," Mangrove Fringe"," Florida Bay","W. Florida Shelf"," Keys")),"Region",sort=F,all.x=T)
region.para.melt$Region.txt=factor(region.para.melt$Region.txt,levels=c("ENP"," Mangrove Fringe"," Florida Bay","W. Florida Shelf"," Keys"))

param.vals=c("TN (mg N L\u207B\u00B9)","DIN (mg N L\u207B\u00B9)","TP (\u03BCg P L\u207B\u00B9)","SRP (\u03BCg P L\u207B\u00B9)","Chl-a (\u03BCg L\u207B\u00B9)","TOC (mg C L\u207B\u00B9)")
region.para.melt=merge(region.para.melt,data.frame(param=c("TN","DIN","TP","SRP","Chla","TOC"),
                                                   param.vals=param.vals),"param")
region.para.melt$param.vals=factor(region.para.melt$param.vals,levels=param.vals)
region.sum=ddply(subset(region.para.melt,Region%in%c("ENP","Coastal_Mangroves","FLBay","Shelf","Keys")),
                 c("param.vals","Region.txt"),summarise,
                 N.val=N.obs(value),
                 min.val=min(value,na.rm=T),
                 med.val=median(value,na.rm=T),
                 mean.val=mean(value,na.rm=T),
                 max.val=max(value,na.rm=T))
library(flextable)
region.sum$Region.txt=with(region.sum,ifelse(Region.txt=="ENP","ENP (Marsh)",trimws(as.character(Region.txt))))
vars=c( "param.vals","Region.txt", "N.val", "min.val", "med.val", "mean.val","max.val")
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


cols=c("white",adjustcolor(wesanderson::wes_palette("Zissou1",4,"continuous"),0.5))
levels.var.labs=c("ENP\n(Marsh)","Mangrove\nFringe","Florida\nBay","W. Florida\nShelf","Keys\n")
# png(filename=paste0(plot.path,"RegionComp.png"),width=6.5,height=5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,3.5,0.5,0.75),oma=c(3,1,1,0.5));
layout(matrix(c(1:6),3,2,byrow=T))

ylim.val=c(0,2);by.y=0.5;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
x=boxplot(GM~Region,dat.all.TN.GM2,outline=F,ylim=ylim.val,axes=F,ann=F,col=cols,boxwex=0.5)
TN.DT=with(dat.all.TN.GM2,dunn.test(GM, Region))
TN.DT.ltr=cldList(P.adjusted ~ comparison,data=TN.DT,threshold = 0.05)
TN.DT.ltr$Letter=toupper(TN.DT.ltr$Letter)
TN.DT.ltr=TN.DT.ltr[order(match(TN.DT.ltr$Group,levels.var)),]
text(1:5,x$stats[5,],TN.DT.ltr$Letter,pos=3)
axis_fun(2,ymaj,ymin,format(ymaj))
axis_fun(1,1:5,1:5,NA);box(lwd=1)
mtext(side=2,line=2.5,"TN (mg N L\u207B\u00B9)")

ylim.val=c(0,2.1);by.y=0.5;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
x=boxplot(GM~Region,dat.all.DIN.GM2,outline=F,ylim=ylim.val,axes=F,ann=F,col=cols,boxwex=0.5)
DIN.DT=with(dat.all.DIN.GM2,dunn.test(GM, Region))
DIN.DT.ltr=cldList(P.adjusted ~ comparison,data=DIN.DT,threshold = 0.05)
DIN.DT.ltr$Letter=toupper(DIN.DT.ltr$Letter)
DIN.DT.ltr=DIN.DT.ltr[order(match(DIN.DT.ltr$Group,levels.var)),]
text(1:5,x$stats[5,],DIN.DT.ltr$Letter,pos=3)
axis_fun(2,ymaj,ymin,format(ymaj))
axis_fun(1,1:5,1:5,NA);box(lwd=1)
mtext(side=2,line=2.5,"DIN (mg N L\u207B\u00B9)")

ylim.val=c(0,60);by.y=20;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
x=boxplot(GM~Region,dat.all.TP.GM2,outline=F,ylim=ylim.val,axes=F,ann=F,col=cols,boxwex=0.5)
TP.DT=with(dat.all.TP.GM2,dunn.test(GM, Region))
TP.DT.ltr=cldList(P.adjusted ~ comparison,data=TP.DT,threshold = 0.05)
TP.DT.ltr$Letter=toupper(TP.DT.ltr$Letter)
TP.DT.ltr=TP.DT.ltr[order(match(TP.DT.ltr$Group,levels.var)),]
text(1:5,x$stats[5,],TP.DT.ltr$Letter,pos=3)
axis_fun(2,ymaj,ymin,format(ymaj))
axis_fun(1,1:5,1:5,NA);box(lwd=1)
mtext(side=2,line=2.5,"TP (\u03BCg P L\u207B\u00B9)")

ylim.val=c(0,10);by.y=2;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
x=boxplot(GM~Region,dat.all.SRP.GM2,outline=F,ylim=ylim.val,axes=F,ann=F,col=cols,boxwex=0.5)
SRP.DT=with(dat.all.SRP.GM2,dunn.test(GM, Region))
SRP.DT.ltr=cldList(P.adjusted ~ comparison,data=SRP.DT,threshold = 0.05)
SRP.DT.ltr$Letter=toupper(SRP.DT.ltr$Letter)
SRP.DT.ltr=SRP.DT.ltr[order(match(SRP.DT.ltr$Group,levels.var)),]
text(1:5,x$stats[5,],SRP.DT.ltr$Letter,pos=3)
axis_fun(2,ymaj,ymin,format(ymaj))
# axis_fun(1,1:5,1:5,levels.var.labs,line=0.3,cex=0.9);box(lwd=1)
axis_fun(1,1:5,1:5,NA,line=0.3,cex=0.9);box(lwd=1)
mtext(side=2,line=2.5,"SRP (\u03BCg P L\u207B\u00B9)")

ylim.val=c(0,7);by.y=2;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
x=boxplot(GM~Region,dat.all.Chla.GM2,outline=F,ylim=ylim.val,axes=F,ann=F,col=cols,boxwex=0.5)
Chla.DT=with(dat.all.Chla.GM2,dunn.test(GM, Region))
Chla.DT.ltr=cldList(P.adjusted ~ comparison,data=Chla.DT,threshold = 0.05)
Chla.DT.ltr$Letter=toupper(Chla.DT.ltr$Letter)
Chla.DT.ltr=Chla.DT.ltr[order(match(Chla.DT.ltr$Group,levels.var)),]
text(1:5,x$stats[5,],Chla.DT.ltr$Letter,pos=3)
axis_fun(2,ymaj,ymin,format(ymaj))
axis_fun(1,1:5,1:5,levels.var.labs,line=0.3,cex=0.8);box(lwd=1)
mtext(side=2,line=2.5,"Chl-a (\u03BCg L\u207B\u00B9)")
mtext(side=1,line=1,outer=T,"Region")

ylim.val=c(0,30);by.y=10;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
x=boxplot(GM~Region,dat.all.TOC.GM2,outline=F,ylim=ylim.val,axes=F,ann=F,col=cols,boxwex=0.5)
TOC.DT=with(dat.all.TOC.GM2,dunn.test(GM, Region))
TOC.DT.ltr=cldList(P.adjusted ~ comparison,data=TOC.DT,threshold = 0.05)
TOC.DT.ltr$Letter=toupper(TOC.DT.ltr$Letter)
TOC.DT.ltr=TOC.DT.ltr[order(match(TOC.DT.ltr$Group,levels.var)),]
text(1:5,x$stats[5,],TOC.DT.ltr$Letter,pos=3)
axis_fun(2,ymaj,ymin,format(ymaj))
axis_fun(1,1:5,1:5,levels.var.labs,line=0.3,cex=0.8);box(lwd=1)
mtext(side=2,line=2.5,"TOC (mg C L\u207B\u00B9)")
mtext(side=1,line=1,outer=T,"Region")

dev.off()

# png(filename=paste0(plot.path,"RegionComp_NP.png"),width=5,height=3,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,3.5,0.25,0.75),oma=c(3,1,0.5,0.5));
ylim.val=c(0,900);by.y=200;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)

x=boxplot(GM~Region,dat.all.TNTP.GM2,outline=F,
          ylim=ylim.val,axes=F,ann=F,col=cols,boxwex=0.5);
abline(h=16,lty=2,col="indianred1")
NP.DT=with(dat.all.TNTP.GM2,dunn.test(GM, Region))
NP.DT.ltr=cldList(P.adjusted ~ comparison,data=NP.DT,threshold = 0.05)
NP.DT.ltr$Letter=toupper(NP.DT.ltr$Letter)
NP.DT.ltr=NP.DT.ltr[order(match(NP.DT.ltr$Group,levels.var)),]
text(1:5,x$stats[5,],NP.DT.ltr$Letter,pos=3)
axis_fun(2,ymaj,ymin,format(ymaj))
axis_fun(1,1:5,1:5,levels.var.labs,line=0.3,cex=0.8);box(lwd=1)
mtext(side=2,line=2.5,"TN:TP (molar ratio)")
mtext(side=1,line=1,outer=T,"Region")
dev.off()