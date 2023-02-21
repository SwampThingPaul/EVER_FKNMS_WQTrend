#to install color packages
install.packages("scales")
install.packages("directlabels")
install.packages("tidyverse")
install.packages("dplyr")
install.packages("patchwork")
install.packages("ggplot2")
install.packages("reshape2")
install.packages("segmented")

#to open tidyverse and dplyr and patchwork
library(scales)
library(directlabels)
library(tidyverse)
library(dplyr)
library(patchwork)
library(ggplot2)
library(reshape2) 
library(segmented)
library(lubridate)

#Convert Date
df <- FCEP
FCEP2 = FCEP %>% 
  mutate(Date = mdy(Date))

#Fit breakpoint model TS Marsh
df <- FCEP2
fitline = lm(TSMarshP ~ Date, data = df)
segmented.fit = segmented(fitline, seg.Z = ~Date)
fit=numeric(length(df$Date))*NA
fit[complete.cases(rowSums(cbind(df$TSMarshP,df$Date)))] = broken.line(segmented.fit)$fit
#Make the plot
ggplot(df, aes(x=Date, y=TSMarshP))+
  geom_point(size=2,aes()) +
  geom_line (aes (x = Date, y = fit)) +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16), legend.text=element_text(size=14), strip.text.x = element_text(size = 14)) +
  xlab("Date\n") + ylab("\nTS/Ph Marsh Water TP (ug L-1)") +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  scale_x_date(breaks = "1 year", labels=date_format("%b-%Y"))
  
#Fit breakpoint model SRS Marsh
df <- FCEP2
fitline = lm(SRSMarshP ~ Date, data = df)
segmented.fit = segmented(fitline, seg.Z = ~Date)
fit=numeric(length(df$Date))*NA
fit[complete.cases(rowSums(cbind(df$SRSMarshP,df$Date)))] = broken.line(segmented.fit)$fit
#Make the plot
ggplot(df, aes(x=Date, y=SRSMarshP))+
  geom_point(size=2,aes()) +
  geom_line (aes (x = Date, y = fit)) +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16), legend.text=element_text(size=14), strip.text.x = element_text(size = 14)) +
  xlab("Date\n") + ylab("\nSRS Marsh Water TP (ug L-1)") +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  scale_x_date(breaks = "1 year", labels=date_format("%b-%Y"))

#Fit breakpoint model TS/Ph Seagrass
df <- FCEP2
fitline = lm(TSSeagrassP ~ Date, data = df)
segmented.fit = segmented(fitline, seg.Z = ~Date)
fit=numeric(length(df$Date))*NA
fit[complete.cases(rowSums(cbind(df$TSSeagrassP,df$Date)))] = broken.line(segmented.fit)$fit
#Make the plot
ggplot(df, aes(x=Date, y=TSSeagrassP))+
  geom_point(size=2,aes()) +
  geom_line (aes (x = Date, y = fit)) +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16), legend.text=element_text(size=14), strip.text.x = element_text(size = 14)) +
  xlab("Date\n") + ylab("\nTS/Ph Seagrass Water TP (ug L-1)") +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  scale_x_date(breaks = "1 year", labels=date_format("%b-%Y"))

#Fit breakpoint model TS/Ph Seagrass
df <- FCEP2
fitline = lm(TSMangroveP ~ Date, data = df)
segmented.fit = segmented(fitline, seg.Z = ~Date)
fit=numeric(length(df$Date))*NA
fit[complete.cases(rowSums(cbind(df$TSSeagrassP,df$Date)))] = broken.line(segmented.fit)$fit
#Make the plot
ggplot(df, aes(x=Date, y=TSMangroveP))+
  geom_point(size=2,aes()) +
  geom_line (aes (x = Date, y = fit)) +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16), legend.text=element_text(size=14), strip.text.x = element_text(size = 14)) +
  xlab("Date\n") + ylab("\nTS/Ph Mangrove Water TP (ug L-1)") +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  scale_x_date(breaks = "1 year", labels=date_format("%b-%Y"))

#Fit breakpoint model TS/Ph Seagrass
df <- FCEP2
fitline = lm(SRSMangroveP ~ Date, data = df)
segmented.fit = segmented(fitline, seg.Z = ~Date)
fit=numeric(length(df$Date))*NA
fit[complete.cases(rowSums(cbind(df$SRSMangroveP,df$Date)))] = broken.line(segmented.fit)$fit
#Make the plot
ggplot(df, aes(x=Date, y=SRSMangroveP))+
  geom_point(size=2,aes()) +
  geom_line (aes (x = Date, y = fit)) +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16), legend.text=element_text(size=14), strip.text.x = element_text(size = 14)) +
  xlab("Date\n") + ylab("\nSRS Mangrove Water TP (ug L-1)") +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  scale_x_date(breaks = "1 year", labels=date_format("%b-%Y"))


df <- fceperi

df$SITE=as.factor(df$SITE)

str(df)

df$SITE <- factor(df$SITE, levels=c("SRS1", "SRS2", "SRS3", "TS/Ph1", "TS/Ph2", "TS/Ph3", "TS/Ph9", "TS/Ph10", "TS/Ph11"))

df()
ggplot(df) + 
  geom_boxplot(aes (SITE, TOC), width = 0.5, position=position_dodge(0.8)) + 
  theme_bw() +
  xlab("Site") + ylab("TOC Accumulation (g m-2 d-1)") +
  ylim (0, 0.5) 

df()
ggplot(df) + 
  geom_boxplot(aes (SITE, TIC), width = 0.5, position=position_dodge(0.8)) + 
  theme_bw() +
  xlab("Site") + ylab("TIC Accumulation (g m-2 d-1)") +
  ylim (0, 0.5) 

df()
ggplot(df) + 
  geom_boxplot(aes (SITE, TP), width = 0.5, position=position_dodge(0.8)) + 
  theme_bw() +
  xlab("Site") + ylab("TP Accumulation (ug m-2 d-1)") +
  ylim (0, 60) 

df <- fceperi2

df$SITE=as.factor(df$SITE)
df$VAR=as.factor(df$VAR)

str(df)

df$SITE <- factor(df$SITE, levels=c("SRS1", "SRS2", "SRS3", "TS/Ph1", "TS/Ph2", "TS/Ph3", "TS/Ph9", "TS/Ph10", "TS/Ph11"))

df()
ggplot(df, aes (SITE, VALUE)) + 
  geom_boxplot(aes (fill=VAR), width = 0.5, position=position_dodge(0.8)) +
  theme_bw() +
  scale_fill_manual(values = c("grey", "white")) +
  labs(fill = "") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16,face="bold"), legend.text=element_text(size=14)) +
  xlab("Site") + ylab("Accumulation Rate (g/m2/d)") +
  ylim (0, 0.5) 

df <- fceperi2

df$SITE=as.factor(df$SITE)
df$VAR=as.factor(df$VAR)
df$TYPE=as.factor(d$TYPE)

str(df)

df$SITE <- factor(df$SITE, levels=c("SRS1", "SRS2", "SRS3", "TS/Ph1", "TS/Ph2", "TS/Ph3", "TS/Ph9", "TS/Ph10", "TS/Ph11"))


df()
ggplot(df, aes (TYPE, VALUE)) + 
  geom_boxplot(aes (fill=TYPE), width = 0.5, position=position_dodge(0.8)) +
  facet_wrap(~ VAR, nrow = 1) +
  scale_fill_manual(values = c("blue", "turquoise2", "green")) +
  theme_bw() +
  labs(fill = "") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16,face="bold"), legend.text=element_text(size=14), strip.text.x = element_text(size = 14)) +
  xlab("Site") + ylab("Accumulation Rate (g/m2/d)") + 
  ylim (0, 0.3) 

