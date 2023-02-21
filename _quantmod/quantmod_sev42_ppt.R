######## Pulse metrics from quantmod for LTER data ########
######## LTER pulse dynamics working group, July 2022, Jennifer Rudgers ########

## This tutorial uses SEV climate data from one SEV ecosystem/met station as an example ##

######## Data steps  ########
rm(list=ls(all=TRUE)) #give R a blank slate

#set the working directory (wd) on your computer - optional time saver
#will save any files you create into that folder
setwd("C:/Users/jrudgers/Desktop/SEV Data/Pulse dynamics/SEV quantmod & fourier 2.0/")

# https://urldefense.com/v3/__http://www.quantmod.com/examples/intro/__;!!FjuHKAHQs5udqho!MsmKtf5lX-ZPH2JAyjRWFewbNfbOHPm-74v0qYzU6Hw3RDE9pRXSz2idBhwkaaPhFp8yWiRqX9FDvDjfhA$  

#Read data from your working directory
#or use read.csv(file.choose()) to import from browser
sev42ppt<-read.csv("sev42ppt_raw.csv")

#call libraries (make sure you have installed packages first)
library(tidyverse)
library(quantmod)
library(xts)
library(emmeans)
library(ggthemes)

#QUANTMOD requires a single time series ordered by the order of observations (time)
#!!! WARNING: xts will not accept any NAs, NaN, or Inf !!!
# for temperature datasets, you may want to fill missing values 
# in xts using last observation
# xts_last <- na.locf(xtsdata) 
# or interpolate NAs na.approx(xtsdata)  
# note: date is being read as character format, we need to format date
#make date data into ISO format
sev42ppt$date<-as.Date(sev42ppt$date)
summary(sev42ppt) 

# filter to location called 42 and remove NAs
sev42ppt<-data.frame(sev42ppt %>% filter (date>'1989-12-31',!is.na(date),!is.na(ppt)))
summary(sev42ppt) 
# check for and remove max T outliers
hist(sev42ppt$ppt)




# STEP 1: Convert to xts (time series) object -----------
#!!! WARNING assumes records of observations are ordered by time !!!
# this function sorts by date so the records are in the correct order 
# prior to analysis
sev42xtsppt <- xts(sev42ppt[,-1], order.by=sev42ppt[,1])


######## Apply quantmod to raw data ########
#### make a quantmod graphic - time series of pulses
chartSeries(sev42xtsppt,theme="white")

# STEP 2: set a threshold -------
threshold.peak<-20
# rationale: 20 mm rain event activates plants 

# STEP 3: Function to identify peaks[valleys] --------
sev42peaksppt<-findPeaks(sev42xtsppt, thresh=threshold.peak)
plot(sev42xtsppt[sev42peaksppt-1])

# GRAPH PEAKS
peaksppt<-as.data.frame(sev42xtsppt[sev42peaksppt-1])

#plot peaks onto raw data 
peaks_graphicppt<-ggplot(sev42ppt, aes(x = date, y = ppt)) +
  geom_line(colour='deepskyblue1') +
  labs(x = "Date",
       y = "Precipitation (mm)")+
  geom_point(data=peaksppt,aes(x = as.Date(row.names(peaksppt)), y = V1), 
             colour='black', size=2)+
  theme_tufte(base_size = 14, base_family="sans")
peaks_graphicppt
#save graphic - revise file name for your dataset
ggsave("sev42_ppt_peaks.jpg",peaks_graphicppt,dpi=300,width=10, height=4)

# STEP 4: Create a data frame to add to the group project -----------
#name your LTER
lter<-"SEV"
#name your site
site<-"sev42"
#name your driver
driver<-"ppt"
#units = units of driver (e.g., mm)
units<-"mm"
#are you analyzing peaks or valleys? [options: peak or valley]
pv="peak"
#save number of months and years
nmonths<-nmonths(sev42xtsppt)
nyears<-nyears(sev42xtsppt)

# STEP 5: Calculate metrics of peaks ----------
#how many peaks per total obs, which are in units of days (here, obs = 365 days*33 years) ?
peak_per_d<-length(peaksppt$V1)/length(sev42ppt$ppt)
peak_per_d

#how many peaks per year?
peaks_per_y<-length(peaksppt$V1)/nyears(sev42xtsppt)
peaks_per_y

#average peak magnitude (in mm for precip)
peak_mean<-mean(peaksppt$V1)
peak_mean

#peak standard deviation
peak_sd<-sd(peaksppt$V1)

#peak CV
peak_CV<-sd(peaksppt$V1)/mean(peaksppt$V1)
peak_CV

# STEP 6: Standardized regression models for temporal change--------
# peak number vs. time: Has the number of peaks increased/decreased over time?
# get slope of (number of peaks per year for each year) vs. year (and p-value) to look for temporal change in number of peaks

# add year and time columns to peaks dataset
peaksppt$time<-as.POSIXct(row.names(peaksppt))
peaksppt$year<-as.numeric(format(as.POSIXct(row.names(peaksppt)),"%Y"))
summary(peaksppt)

# first, add any missing years that had no peaks (add zeros) - probably a more efficient way to do this...
year<-min(as.numeric(format(as.POSIXct(sev42ppt$date),"%Y"))):max(as.numeric(format(as.POSIXct(sev42ppt$date),"%Y")))
years<-as.data.frame(year)
years
peak.sum<-peaksppt %>% group_by(year) %>% summarise(mean.peak=mean(V1), count=n())
peak.sum
peak.number<-merge(years,peak.sum,by.x="year",by.y="year",all.x=TRUE)
peak.number[is.na(peak.number)] <- 0 
peak.number

# second, build the stats models and save the slope and p as output
peak.number.lm<-lm(scale(count)~scale(year),data=peak.number)
summary(peak.number.lm)


# centering is done by subtracting the column means (omitting NAs) of V1 or year from their corresponding column
# scaling is done by dividing the (centered) column of V1 or year by their standard deviations 
# this returns a standardized regression coefficient that makes it possible to compare across drivers that differ in measured units and/or timescales

#plot(peak.number.lm) # turn this on to check model statistical assumptions
lmsum.number<-summary(peak.number.lm)
peak.number.slope<-peak.number.lm$coefficients[2]
peak.number.slope
peak.number.p<-lmsum.number$coefficients[2,4]
peak.number.p

peak.number.graph<-ggplot(data=peak.number, aes(x = year, y = count)) + 
  geom_point(color='darkslategrey',size=4) +
  geom_smooth(method = "lm",  se = FALSE, color="darkslategrey")+
  ylab("Number of peaks")+
  xlab("Year")+
  theme_tufte(base_size = 14, base_family="sans")
peak.number.graph
#save graphic - revise file name for your dataset
ggsave("sev42_ppt_peak.number.jpg",peak.number.graph,dpi=300,width=10, height=4)

# peak magnitude (all peaks) vs. time: Has the magnitude of peaks increased/decreased over time?
# get slope of magnitude of peaks vs. time (and p-value)
peak.magnitude.lm<-lm(scale(V1)~scale(time),data=peaksppt)

# plot(peak.magnitude.lm) # turn this on to check model statistical assumptions
lmsum.mag<-summary(peak.magnitude.lm)
peak.magnitude.slope<-peak.magnitude.lm$coefficients[2]
peak.magnitude.slope
peak.magnitude.p<-lmsum.mag$coefficients[2,4]
peak.magnitude.p

peak.mag.graph<-ggplot(data=peaksppt, aes(x = time, y = V1)) + 
  geom_point(color='darkslategrey',size=4) +
  geom_smooth(method = "lm",  se = FALSE, color="darkslategrey")+
  ylab("Average peak magnitude (mm)")+
  xlab("Year")+
  theme_tufte(base_size = 14, base_family="sans")
peak.mag.graph
#save graphic - revise file name for your dataset
ggsave("sev42_ppt_peak.magnitude.jpg",peak.mag.graph,dpi=300,width=10, height=4)

# STEP 7:  Save metrics into your site location data frame ------------
pulse_metrics_sev42ppt<-data.frame(lter,site,driver,units,pv,nyears,nmonths,peak_mean,peak_sd,peak_CV,peaks_per_y,peak_per_d,
                                peak.number.slope,peak.number.p,peak.magnitude.slope,peak.magnitude.p,threshold.peak)
write.csv(pulse_metrics_sev42ppt,"pulse_metrics_sev42ppt_quantmod.csv")
#then you will paste your data into our cross-site dataframe










