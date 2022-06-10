
## Run 7_RegionalTrend_RegionalAnalysis.R first

# WYs=seq(1996,2019,1)
# dates=date.fun(c(paste(min(WYs)-1,"05-01",sep="-"),paste(max(WYs),"05-01",sep="-")))

# Hydroperiod -------------------------------------------------------------

stage.sites=data.frame(SITE=c("NP201","NESRS1","NESRS2","NP203","P33","P36","NP-112","TSB","R127","TSH"),
                       Region=c(rep("UpperSRS",3),rep("MidSRS",2),"LowSRS",rep("UpperTS",2),"MidTS","LowTS"),
                       min.GndElev.NAVD88=c(5.32,4.27,4.19,3.01,3.81,1.65,1.98,1.29,0.00,-0.15),
                       NAVD88.to.NGVD29=c(-1.51,-1.53,-1.54,-1.51,-1.51,-1.51,-1.51,-1.51,-1.57,-1.56))
stage.sites$min.GndElev.NGVD29=with(stage.sites,min.GndElev.NAVD88-NAVD88.to.NGVD29)
stage.sites

stg.dbkeys=data.frame(SITE=c("NP201","NESRS1","NESRS2","NP203","P33","P36","NP-112","TSB","TSH"),
                      DBKEY=c("06719","01140","01218","G6154","06717","06718","H2427","H2442","07090"))


# download R127 from EDEN webpage
# https://sofia.usgs.gov/eden/station.php?stn_name=R127
R127=read.csv(paste0(data.path,"EDEN/1604536339_water_level.csv"),skip=4)
R127$SITE="R127"
R127$Date=date.fun(as.character(R127$Date))
R127$Data.Value=R127$R127.Daily.median.water.Level.feet.NAVD88.-(-1.57)
R127=R127[,c("SITE","Date","Data.Value")]

stg.dat=data.frame()
for(i in 1:nrow(stg.dbkeys)){
  tmp=DBHYDRO_daily(dates[1],dates[2],stg.dbkeys$DBKEY[i])
  tmp$DBKEY=as.character(stg.dbkeys$DBKEY[i])
  stg.dat=rbind(stg.dat,tmp)
  print(i)
}
stg.dat=merge(stg.dat,stg.dbkeys,"DBKEY")

stg.dat=stg.dat[,c("SITE","Date","Data.Value")]
stg.dat=rbind(stg.dat,R127)
stg.dat=merge(stg.dat,stage.sites[,c("SITE","Region","min.GndElev.NGVD29")],"SITE")
stg.dat$WL.GT.Gnd=with(stg.dat,ifelse(Data.Value>min.GndElev.NGVD29,1,0))
stg.dat$Depth.cm=with(stg.dat,ft.to.m(Data.Value-min.GndElev.NGVD29)*100)
stg.dat$WY=WY(stg.dat$Date)

stg.dat.WY.HP=ddply(stg.dat,c("SITE","Region","WY"),summarise,HP.days=sum(WL.GT.Gnd,na.rm=T),N.val=N.obs(SITE),meanDepth=mean(Depth.cm,na.rm=T))
stg.dat.WY.HP=subset(stg.dat.WY.HP,WY%in%seq(1996,2019,1))
stg.dat.WY.HP$HP.freq=with(stg.dat.WY.HP,HP.days/N.val)
subset(stg.dat.WY.HP,N.val!=365)

plot(HP.freq~WY,subset(stg.dat.WY.HP,SITE=="NP-112"))
plot(meanDepth~WY,subset(stg.dat.WY.HP,SITE=="NP-112"))
plot(HP.freq~WY,subset(stg.dat.WY.HP,SITE=="R127"))
plot(meanDepth~WY,subset(stg.dat.WY.HP,SITE=="R127"))
stg.dat.WY.HP.region=ddply(stg.dat.WY.HP,c("WY","Region"),summarise,meanHP.f=mean(HP.freq,na.rm=T),meanZ=mean(meanDepth,na.rm=T))

ddply(stg.dat.WY.HP.region,c("Region"),summarise,
      est=as.numeric(cor.test(WY,meanZ,method="kendall")$estimate),
      pval=cor.test(WY,meanZ,method="kendall")$p.value)



# Discharge Data ----------------------------------------------------------
SRS.dbkeys=data.frame(SITE=c(rep("S12A",3),rep("S12B",3),rep("S12C",3),rep("S12D",3),rep("S333",2),rep("S334",2),"S355A","S355B_P","S355B","S356","S335"),
                      DBKEY=c(c("FE771","P0796","03620"),c("FE772","P0950","03626"),c("FE773","P0951","03632"),c("FE774","P0952","03638"),c("15042","91487"),c("FB752","91488"),"MQ895","AM173","MQ896","64136","91489"),
                      Priority=paste0("P",c(rep(c(1,2,3),4),c(1,2),c(1,2),rep(1,5))),
                      WQSite=c(rep("S12A",3),rep("S12B",3),rep("S12C",3),rep("S12D",3),rep("S333",2),rep("S334",2),"S355A","S355B","S355B","S356-334","S356-334"))


srs.flow=data.frame()
for(i in 1:nrow(SRS.dbkeys)){
  tmp=DBHYDRO_daily(dates[1],dates[2],SRS.dbkeys$DBKEY[i])
  tmp$DBKEY=SRS.dbkeys$DBKEY[i]
  srs.flow=rbind(tmp,srs.flow)
  print(i)
}
srs.flow$WY=WY(srs.flow$Date)
srs.flow=merge(srs.flow,SRS.dbkeys,"DBKEY")
flow.xtab=data.frame(reshape::cast(srs.flow,Date+WY+SITE+WQSite~Priority,value="Data.Value",fun.aggregate=function(x) ifelse(sum(x,na.rm=T)==0,NA,sum(x,na.rm=T))))
flow.xtab$fflow.cfs=with(flow.xtab,ifelse(is.na(P1)==T&is.na(P2)==T,P3,ifelse(is.na(P2)==T&is.na(P3)==T,P1,ifelse(is.na(P1)==T&is.na(P3)==T,P2,P1))));#final flow value for analysis

srs.flow.da.xtab=reshape::cast(flow.xtab,Date+WY~SITE,value="fflow.cfs",fun.aggregate=function(x) mean(cfs.to.km3d(ifelse(x<0,NA,x)),na.rm=T))
srs.flow.da.xtab$S355B=rowSums(srs.flow.da.xtab[,c("S355B","S355B_P")],na.rm=T)
srs.flow.da.xtab=srs.flow.da.xtab[,c("Date","WY","S12A","S12B", "S12C", "S12D", "S333", "S334", "S355A", "S355B","S356","S335")]

plot(S12A~Date,srs.flow.da.xtab)
srs.flow.da.xtab[is.na(srs.flow.da.xtab)]<-0
srs.flow.da.xtab$minS356_S335=apply(srs.flow.da.xtab[,c('S356','S335')],1,min); # with(srs.flow.da.xtab,ifelse(S356>S335,S335,S356));
srs.flow.da.xtab$m15.ESRS.Q=rowSums(srs.flow.da.xtab[,c("S333","S355A","S355B","minS356_S335")],na.rm=T)
# srs.flow.da.xtab$S333R.m15=with(srs.flow.da.xtab,S333*(1-ifelse(m15.ESRS.Q>0,ifelse(S334/m15.ESRS.Q>1,1,S334/m15.ESRS.Q),0)))
# srs.flow.da.xtab$S355A.adj.m15=with(srs.flow.da.xtab,S355A*(1-ifelse(m15.ESRS.Q>0,ifelse(S334/m15.ESRS.Q>1,1,S334/m15.ESRS.Q),0)))
# srs.flow.da.xtab$S355B.adj.m15=with(srs.flow.da.xtab,S355B*(1-ifelse(m15.ESRS.Q>0,ifelse(S334/m15.ESRS.Q>1,1,S334/m15.ESRS.Q),0)))
# srs.flow.da.xtab$S356.adj.m15=with(srs.flow.da.xtab,minS356_S335*(1-ifelse(m15.ESRS.Q>0,ifelse(S334/m15.ESRS.Q>1,1,S334/m15.ESRS.Q),0)))
srs.flow.da.xtab$TFlow.m15=rowSums(srs.flow.da.xtab[,c("S12A","S12B","S12C","S12D","S333","S355A","S355B","minS356_S335")],na.rm=T)
srs.flow.da.xtab.WY=ddply(subset(srs.flow.da.xtab,WY%in%seq(1996,2019,1)),"WY",summarise,TFlow.km3=sum(TFlow.m15,na.rm=T))

plot(TFlow.km3~WY,srs.flow.da.xtab.WY,type="b")
with(srs.flow.da.xtab.WY,cor.test(WY,TFlow.km3,method="kendall"))
# write.csv(srs.flow.da.xtab.WY,paste0(export.path,"SRS_WYDischarge.csv"),row.names = F)
# https://www.sfwmd.gov/sites/default/files/documents/ag_item_4_TSCB_Method3.pdf
# https://www.sfwmd.gov/sites/default/files/documents/TS%26CB_TP_Compliance_1st_qtr_2020_Data_Report_method1%262%263_0.pdf

# TSCB
TS.dbkeys=data.frame(SITE=c("S332","S175","S332D","S332D","S332DX1","S328","G737","S18C"),
                     DBKEY=c("91486","15752","TA413","91485","91484","AN558","AN674","15760"),
                     Priority=c("P1","P1","P1","P2","P1","P1","P1","P1"))
TS.flow=data.frame()
for(i in 1:nrow(TS.dbkeys)){
  tmp=DBHYDRO_daily(dates[1],dates[2],TS.dbkeys$DBKEY[i])
  tmp$DBKEY=TS.dbkeys$DBKEY[i]
  TS.flow=rbind(tmp,TS.flow)
  print(i)
}
TS.flow$WY=WY(TS.flow$Date)
TS.flow=merge(TS.flow,TS.dbkeys,"DBKEY")
TS.flow$Date=date.fun(TS.flow$Date)

ts.flow.xtab=data.frame(reshape::cast(TS.flow,Date+WY+SITE~Priority,value="Data.Value",fun.aggregate=function(x) ifelse(sum(x,na.rm=T)==0,NA,sum(x,na.rm=T))))
ts.flow.xtab$P3=NA
ts.flow.xtab$fflow.cfs=with(ts.flow.xtab,ifelse(is.na(P1)==T&is.na(P2)==T,P3,ifelse(is.na(P2)==T&is.na(P3)==T,P1,ifelse(is.na(P1)==T&is.na(P3)==T,P2,P1))));#final flow value for analysis
ts.flow.xtab$fflow.cfs=with(ts.flow.xtab,ifelse(SITE%in%c("S332","S175")&Date>date.fun("1999-08-30"),NA,fflow.cfs))

ts.flow.da.xtab=reshape::cast(ts.flow.xtab,Date+WY~SITE,value="fflow.cfs",fun.aggregate=function(x) mean(cfs.to.km3d(ifelse(x<0,NA,x)),na.rm=T))
ts.flow.da.xtab[is.na(ts.flow.da.xtab)]<-0
ts.flow.da.xtab$DBasin.Flow=with(ts.flow.da.xtab,(S332D-S332DX1-S328))
ts.flow.da.xtab$DBasin.Flow=with(ts.flow.da.xtab,ifelse(DBasin.Flow<0,0,DBasin.Flow))
ts.flow.da.xtab$TFlow=with(ts.flow.da.xtab,DBasin.Flow+S328+G737+S18C+S332+S175)

ts.flow.da.xtab.WY=ddply(subset(ts.flow.da.xtab,WY%in%seq(1996,2019,1)),"WY",summarise,TFlow.km3=sum(TFlow,na.rm=T))

plot(TFlow.km3~WY,ts.flow.da.xtab.WY,type="b")
with(ts.flow.da.xtab.WY,cor.test(WY,TFlow.km3,method="kendall"))

# write.csv(ts.flow.da.xtab.WY,paste0(export.path,"TSPh_WYDischarge.csv"),row.names = F)
##
SRS.hydro=merge(subset(stg.dat.WY.HP.region,Region%in%c("UpperSRS","MidSRS","LowSRS")),
                srs.flow.da.xtab.WY,"WY")
range(SRS.hydro$meanHP,na.rm=T)
ddply(SRS.hydro,"Region",summarise,mean.val=mean(meanHP.f))
range(SRS.hydro$TFlow.km3)

plot(meanHP.f~TFlow.km3,SRS.hydro,type="n")
with(subset(SRS.hydro,Region=="UpperSRS"),points(TFlow.km3,meanHP.f,pch=21,bg="indianred1"))
with(subset(SRS.hydro,Region=="MidSRS"),points(TFlow.km3,meanHP.f,pch=21,bg="dodgerblue1"))
with(subset(SRS.hydro,Region=="LowSRS"),points(TFlow.km3,meanHP.f,pch=21,bg="forestgreen"))

ddply(SRS.hydro,c("Region"),summarise,
      est=as.numeric(cor.test(TFlow.km3,meanHP.f,method="spearman")$estimate),
      pval=cor.test(TFlow.km3,meanHP.f,method="spearman")$p.value)

srs.HP.Q=lm(meanHP.f~TFlow.km3+WY+Region,SRS.hydro)
layout(matrix(1:4,2,2));plot(srs.HP.Q)
gvlma::gvlma(srs.HP.Q)

TS.hydro=merge(subset(stg.dat.WY.HP.region,Region%in%paste0(c("Upper","Mid","Low"),"TS")),
               ts.flow.da.xtab.WY,"WY")

range(TS.hydro$meanHP.f,na.rm=T)
ddply(TS.hydro,"Region",summarise,mean.val=mean(meanHP.f))
range(TS.hydro$TFlow.km3)

ddply(TS.hydro,c("Region"),summarise,
      est=as.numeric(cor.test(TFlow.km3,meanHP.f,method="spearman")$estimate),
      pval=cor.test(TFlow.km3,meanHP.f,method="spearman")$p.value)


### hydro figure
# png(filename=paste0(plot.path,"Hydro_Trend.png"),width=6.5,height=4,units="in",res=200,bg="white")
par(family="serif",oma=c(0.5,0.5,0.5,1),mar=c(0.1,0.1,0.1,0.1),xpd=F)
layout(matrix(c(1:6),2,3,byrow=T))

SRS.stg.site=c("NP-201","NESRS1","NESRS2","NP-203","NP-P33","NP-P36")
TS.stg.site=c("NP-112","NP-TSB","NP-TSH")

bbox.lims=bbox(subset(sloughs,NAME=="Shark River"))
plot(ENP.shore,col="white",border="grey",lwd=0.05,ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)])
plot(sloughs.clp,col="grey90",border=NA,add=T)
plot(canal,col="grey",add=T)
plot(ENP,bg=NA,lwd=0.5,add=T)
plot(subset(structures,NAME%in%c(paste0("S12",LETTERS[1:4]),"S333","S334","S356","S355A","S355B")),add=T,pch=21,bg="black")
tmp=subset(wmd.mon,STATION%in%SRS.stg.site&ACTIVITY_S=="Stage")
tmp=tmp[order(match(tmp$STATION,SRS.stg.site)),]
plot(tmp,pch=21,bg=c(rep("indianred1",3),rep("dodgerblue1",2),"forestgreen"),add=T,lwd=0.01)
text(tmp,"SITE",pos=1,cex=0.5,halo=T)
mapmisc::scaleBar(utm17,"bottomleft",bty="n",cex=1,seg.len=4)
legend("topleft",legend=c("Discharge","Upper Slough Stage","Mid Slough Stage","Lower Slough Stage"),
       pch=c(21),lty=c(NA),lwd=c(0.01),
       pt.bg=c("black","indianred1","dodgerblue1","forestgreen"),
       pt.cex=1,ncol=1,cex=0.7,bg="white",box.lty=0,y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=1,title.adj = 0,title=" Monitoring Locations")
box(lwd=1)

xlim.val=c(1996,2019);by.x=10;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/by.x)
ylim.val=c(0.5,1);by.y=0.1;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
par(family="serif",mar=c(3,4,0.5,0.5),xpd=F)
plot(meanHP.f~WY,stg.dat.WY.HP.region,type="n",ylim=ylim.val,xlim=xlim.val,ann=F,axes=F)
abline(h=ymaj,v=xmaj,lty=3,col="grey",lwd=0.5)
with(subset(stg.dat.WY.HP.region,Region=="UpperSRS"),pt_line(WY,meanHP.f,1,"indianred1",1,21,"indianred1"))
with(subset(stg.dat.WY.HP.region,Region=="MidSRS"),pt_line(WY,meanHP.f,1,"dodgerblue1",1,21,"dodgerblue1"))
with(subset(stg.dat.WY.HP.region,Region=="LowSRS"),pt_line(WY,meanHP.f,1,"forestgreen",1,21,"forestgreen"))
axis_fun(1,xmaj,xmin,xmaj,line=-0.5)
axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=2,line=2.5,"Hydroperiod (Prop of WY)",cex=0.8)

ylim.val=c(0,3);by.y=1;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(TFlow.km3~WY,srs.flow.da.xtab.WY,type="n",ylim=ylim.val,xlim=xlim.val,ann=F,axes=F)
abline(h=ymaj,v=xmaj,lty=3,col="grey",lwd=0.5)
with(srs.flow.da.xtab.WY,pt_line(WY,TFlow.km3,1,"black",1.25,21,"black"))
axis_fun(1,xmaj,xmin,xmaj,line=-0.5)
axis_fun(2,ymaj,ymin,format(ymaj,nsmall=1));box(lwd=1)
mtext(side=2,line=2.5,"Discharge (km\u00B3 WY\u207B\u00B9)",cex=0.8)
text(xlim.val[2]+2,ylim.val[2],"A",xpd=NA,cex=1.5)

par(mar=c(0.1,0.1,0.1,0.1))
bbox.lims=bbox(subset(sloughs,NAME=="Taylor Sloug"))
plot(ENP.shore,col="white",border="grey",lwd=0.05,ylim=c(bbox.lims[c(2)],2818533),xlim=bbox.lims[c(1,3)])
plot(sloughs.clp,col="grey90",border=NA,add=T)
plot(canal,col="grey",add=T)
plot(ENP,bg=NA,lwd=0.5,add=T)
plot(subset(structures,NAME%in%c("S332","S332D","S332DX1","S328","G737","S18C","S175")),add=T,pch=21,bg="black")
tmp=subset(wmd.mon,STATION%in%TS.stg.site&ACTIVITY_S=="Stage")
tmp=tmp[order(match(tmp$STATION,TS.stg.site)),]
plot(tmp,pch=21,bg=c(rep("indianred1",2),"forestgreen"),add=T,lwd=0.01)
text(tmp,"SITE",pos=2,cex=0.5,halo=T)
EDEN.R127=SpatialPointsDataFrame(data.frame(UTMX=c(539626.9),UTMY=c(2804113.2)),data=data.frame(SITE="R127"),proj4string=CRS(SRS_string="EPSG:26917"))
text(EDEN.R127,"SITE",pos=2,cex=0.5,halo=T)
plot(EDEN.R127,add=T,pch=21,bg="dodgerblue1",lwd=0.01)
mapmisc::scaleBar(utm17,"bottomleft",bty="n",cex=1,seg.len=4)
box(lwd=1)

xlim.val=c(1996,2019);by.x=10;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/by.x)
ylim.val=c(0.5,1);by.y=0.1;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
par(family="serif",mar=c(3,4,0.5,0.5),xpd=F)
plot(meanHP.f~WY,stg.dat.WY.HP.region,type="n",ylim=ylim.val,xlim=xlim.val,ann=F,axes=F)
abline(h=ymaj,v=xmaj,lty=3,col="grey",lwd=0.5)
with(subset(stg.dat.WY.HP.region,Region=="UpperTS"),pt_line(WY,meanHP.f,1,"indianred1",1,21,"indianred1"))
with(subset(stg.dat.WY.HP.region,Region=="MidTS"),pt_line(WY,meanHP.f,1,"dodgerblue1",1,21,"dodgerblue1"))
with(subset(stg.dat.WY.HP.region,Region=="LowTS"),pt_line(WY,meanHP.f,1,"forestgreen",1,21,"forestgreen"))
axis_fun(1,xmaj,xmin,xmaj,line=-0.5)
axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=2,line=2.5,"Hydroperiod (Prop of WY)",cex=0.8)
mtext(side=1,line=1.75,"Water Year",cex=0.8)

ylim.val=c(0,1);by.y=0.2;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(TFlow.km3~WY,ts.flow.da.xtab.WY,type="n",ylim=ylim.val,xlim=xlim.val,ann=F,axes=F)
abline(h=ymaj,v=xmaj,lty=3,col="grey",lwd=0.5)
with(ts.flow.da.xtab.WY,pt_line(WY,TFlow.km3,1,"black",1.25,21,"black"))
axis_fun(1,xmaj,xmin,xmaj,line=-0.5)
axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=2,line=2.5,"Discharge (km\u00B3 WY\u207B\u00B9)",cex=0.8)
mtext(side=1,line=1.75,"Water Year",cex=0.8)
text(xlim.val[2]+2,ylim.val[2],"B",xpd=NA,cex=1.5)
dev.off()