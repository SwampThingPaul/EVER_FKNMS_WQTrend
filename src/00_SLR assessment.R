## 
## Sea level stress 
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

library(zoo)


# -------------------------------------------------------------------------

vaca.key.link="https://tidesandcurrents.noaa.gov/sltrends/data/8723970_meantrend.txt"
vaca.sl.dat=read.table(vaca.key.link,header=T)
date.fil=seq(date.fun("1971-01-01"),date.fun("2022-12-01"),"1 month")
vaca.sl.dat=merge(vaca.sl.dat,
                  data.frame(Year=as.numeric(format(date.fil,"%Y")),Month=as.numeric(format(date.fil,"%m"))),
                  c("Year","Month"),all.y=T)

vaca.sl.dat$Monthly_MSL=vaca.sl.dat$Monthly_MSL*1000
vaca.sl.dat$MA_MSL_6mon=with(vaca.sl.dat,c(rep(NA,5),rollapply(Monthly_MSL,width=6,FUN=function(x) mean(x,na.rm=T))))
vaca.sl.dat$MA_MSL_5YR=with(vaca.sl.dat,c(rep(NA,59),rollapply(Monthly_MSL,width=60,FUN=function(x) mean(x,na.rm=T))))
vaca.sl.dat$monCY_date=with(vaca.sl.dat,date.fun(paste(Year,Month,"01",sep="-")))


plot(Monthly_MSL~monCY_date,vaca.sl.dat)
lines(MA_MSL_5YR~monCY_date,vaca.sl.dat,col="red",lwd=2)
lines(MA_MSL_6mon~monCY_date,vaca.sl.dat,col="dodgerblue",lwd=2)


vaca.sl.dat$SLS=with(vaca.sl.dat,MA_MSL_6mon-MA_MSL_5YR)

plot(SLS~monCY_date,vaca.sl.dat)
