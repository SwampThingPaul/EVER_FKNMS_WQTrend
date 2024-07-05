## 
## Everglades/South Florida Nitrogen Synthesis
##   Exploring random forest
## 
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
library(openxlsx)
library(classInt)
library(zoo)

#GIS Libraries
# library(sp)
# library(rgdal)
# library(PROJ)
# library(rgeos)
# library(tmap)
# library(raster)
# 
# library(fields)
# library(spdep)
# 
# # GAM
# library(ggplot2)
# library(viridis)
# library(mgcv)
# library(gratia)

#Paths
wd="C:/Julian_LaCie/_GitHub/EVER_FKNMS_WQTrend"

paths=paste0(wd,c("/Plots/","/Export/","/Data/","/GIS","/src/"))
# Folder.Maker(paths);#One and done. Creates folders in working directory.
plot.path=paths[1]
export.path=paths[2]
data.path=paths[3]
GIS.path=paste0(wd,"/GIS")

gen.GIS="C:/Julian_LaCie/_GISData"
db.path=paste(gen.GIS,"/SFER_GIS_Geodatabase.gdb",sep=""); 

#Helper variables 
N.mw=14.0067
P.mw=30.973762
C.mw=12.0107

wgs84=CRS("+init=epsg:4326")
utm17=CRS("+init=epsg:26917")

tmap_mode("view")


# -------------------------------------------------------------------------

# load(paste0(export.path,"regiontrend_reanalysis.RData"))

dat.all.GM
dat.xtab=dcast(dat.all.GM,STATION+Region+WY~variable,value.var="GM",mean)


# Random Forest -----------------------------------------------------------
### following this 
## https://www.r-bloggers.com/2021/04/random-forest-in-r/

library(randomForest)
library(caret)


## Random Forest
head(dat.xtab)
# vars=c('Region',"TN","TP","TOC","WY")
vars=c('Region',"TN","TP","TOC",'Chla',"SALINITY")

dat.xtab.cl=na.omit(dat.xtab[,vars])
dat.xtab.cl$Region=as.factor(dat.xtab.cl$Region)

# data partition
set.seed(222)
ind <- sample(2, nrow(dat.xtab.cl), replace = TRUE, prob = c(0.7, 0.3))
train <- dat.xtab.cl[ind==1,]
test <- dat.xtab.cl[ind==2,]

## random forest
# rf=randomForest(Region~.,data=dat.xtab.cl,proximity=T)
rf=randomForest(Region~.,data=train,proximity=T)
print(rf)

# prediction & confusion matrix - train
p1=predict(rf,train)
caret::confusionMatrix(p1,train$Region)

# prediction & confusion matrix - test
p2=predict(rf,test)
caret::confusionMatrix(p2,test$Region)


## Error rate 
head(rf$err.rate)

plot(rf)
lines(1:500,rf$err.rate[,1],col="purple",lwd=2); #OOB
lines(1:500,rf$err.rate[,2],col="blue",lwd=2); # Coastal Mangroves
lines(1:500,rf$err.rate[,3],col="green",lwd=2); # FLBay
lines(1:500,rf$err.rate[,4],col="red",lwd=2); # Keys
lines(1:500,rf$err.rate[,5],col="grey",lwd=2); # Shelf

## No. of nodes
hist(treesize(rf))
# Variable Importance
varImpPlot(rf,
           sort = T,
           n.var = 5)
importance(rf)
## TOC the most important


partialPlot(rf, train, TOC, "FLBay")
partialPlot(rf, train, TOC, "Shelf")


## MDS plot
MDSplot(rf, train$Region)
