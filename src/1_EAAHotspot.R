## Title:      EAA Hotspot Data analysis 
## Created by: Paul Julian (pjulian@evergladesfoundation.org)
## Created on: 02/13/2023

## BAD ## https://www.tidyverse.org/articles/2017/12/workflow-vs-script/
#Clears Everything...start fresh.
rm(list=ls(all=T));cat("\014");dev.off()

# Libraries
# devtools::install_github("SwampThingPaul/AnalystHelper")
library(AnalystHelper);
library(plyr)
library(reshape2)
library(openxlsx)

library(rgdal)
library(rgeos)
library(raster)

#Paths
wd="C:/Julian_LaCie/_GitHub/EAAHotSpot"
paths=paste0(wd,c("/Plots/","/Export/","/Data/","/GIS"))
plot.path=paths[1]
export.path=paths[2]
data.path=paths[3]
GIS.path=paths[4]

GIS.path.gen="C:/Julian_LaCie/_GISData"

# Helper variables
nad83.pro=CRS("+init=epsg:4269")
utm17=CRS("+proj=utm +zone=17 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

## GIS Data ---------------------------------------------------------------
wmd.mon=spTransform(readOGR(paste0(GIS.path.gen,"/SFWMD_Monitoring_20200221"),"Environmental_Monitoring_Stations"),utm17)
# wmd.mon=subset(wmd.mon,ACTIVITY_S=="Surface Water Grab")

subset(wmd.mon,SITE=="WP04.1TN01")@data

# ogrListLayers(paste(GIS.path.gen,"/AHED_release/AHED_20171102.gdb",sep=""))
basins=spTransform(readOGR(paste0(GIS.path.gen,"/AHED_release/AHED_20171102.gdb"),"WATERSHED"),utm17)

# -------------------------------------------------------------------------
# Structure level data provided by SFWMD. 

## farmID -----------------------------------------------------------------
farm.xwalk=read.xlsx(paste0(data.path,"Farm and structure.xlsx"))
farm.vars=c("X1", "Permit", "New.Permit", "Farm.ID", "Acres", "Old.acres", 
            "Basin", "Structures1", "Basin", "Structures2", "Comments", "Structures3", 
            "Structures4", "Structures5", "Structures6", "Structures7", "Structures8",
            "Structures9", "Structures10", "Structures11","Structures12", "Structures13",
            "Land.use", "Discharging", "Annual.allocation.(MG/Y)")
colnames(farm.xwalk)=farm.vars
farm.struct.xwalk=melt(farm.xwalk[,c("Farm.ID",paste0("Structures",1:13))],id.vars = "Farm.ID",value.name="structure")
farm.struct.xwalk=subset(farm.struct.xwalk,is.na(structure)==F)
unique(farm.struct.xwalk$structure)

str.vals=farm.struct.xwalk$structure[is.na(sapply(strsplit(farm.struct.xwalk$structure,"\\."),"[",2))]
str.tmp=with(farm.struct.xwalk,ifelse(structure%in%str.vals,NA,strsplit(structure,"-")))
farm.struct.xwalk$StructureName=sapply(str.tmp,"[",1);#trim any extra suffixes

ddply(farm.struct.xwalk,"Farm.ID",summarise,N.val=N.obs(StructureName))
farm.struct.xwalk=farm.struct.xwalk[,c("Farm.ID","structure","StructureName")]

## WQ Data ----------------------------------------------------------------
head.val=c("StructureName","date","TP.mgL","Q.milgal","RF.in","WQSampType")
dat=read.table(paste0(data.path,"eaawqdwn.dat"),col.names=head.val)
dat$date2=with(dat,date.fun(paste(substr(date,1,4),substr(date,5,6),substr(date,7,8),sep="-")))

## Data cleaning
dat[dat==-9]=NA #-9 is NA values

# check zero TP values
dates.vals=date.fun(c("1993-05-01","2020-02-29"))
sites=unique(subset(dat,TP.mgL==0)$StructureName)

tmp=DBHYDRO_WQ(dates.vals[1],dates.vals[2],sites[2],25)
unique(subset(dat,StructureName==sites[2])$TP.mgL)

subset(tmp,Date.EST==date.fun("2000-08-23"))
subset(dat,StructureName==sites[2]&date2==date.fun("2000-08-23"))

tmp2=subset(dat,StructureName==sites[2])
tmp2=merge(tmp2,tmp[,c("Date.EST","HalfMDL")],by.x="date2",by.y="Date.EST")

plot(TP.mgL~HalfMDL,tmp2,ylab="Data for paper",xlab="DBHYDRO")

# dat$TP.mgL[dat$TP.mgL==0]=NA

subset(dat,TP.mgL==0)
range(subset(dat,TP.mgL==0)$date2)
nrow(subset(dat,TP.mgL==0))
summary(dat)

dat$Q.cfs=dat$Q.milgal*1.54722865;# convert million gallons per day to cfs
dat$CY=as.numeric(format(dat$date2,"%Y"))
dat$WY=WY(dat$date2)

range(dat$date2)
range(dat$CY)
range(dat$WY)

dat2=subset(dat,WY%in%seq(2000,2019,1))

# data checking
subset(dat2,TP.mgL==0)
range(subset(dat2,TP.mgL==0)$date2)
nrow(subset(dat2,TP.mgL==0))
nrow(subset(dat,TP.mgL==0&Q.cfs!=0))

## 
unique(dat2$StructureName)
unique(substr(dat2$StructureName,1,2))

plot(subset(wmd.mon,SITE%in%unique(dat2$StructureName)))
plot(basins,add=T)

dat2=merge(dat2,farm.struct.xwalk[,c("StructureName","Farm.ID")],"StructureName",all.x=T)
unique(subset(dat2,is.na(Farm.ID))$StructureName)
length(unique(subset(dat2,is.na(Farm.ID))$StructureName))
paste(unique(subset(dat2,is.na(Farm.ID))$StructureName),collapse=", ")

subset(farm.struct.xwalk,substr(StructureName,1,4)=="BC19")
subset(farm.struct.xwalk,substr(StructureName,1,4)=="BC19")
subset(farm.struct.xwalk,substr(StructureName,1,4)=="OC09")

subset(wmd.mon,SITE=="OC09.5TN13")@data


### FWM values
cfs.to.Ld=function(x) x*(28.3168/(1/86400))
af.to.liters=function(x) x* 1233481.85532

head(dat2)
dat2=dat2[order(dat2$StructureName,dat2$date2),]
subset(dat2,StructureName=="BC00.1TN")
plot(TP.mgL~date2,subset(dat2,StructureName=="BC00.1TN"))

dat2$load.mt=with(dat2,Load.Calc.kg(Q.cfs,TP.mgL))

farm.level.load.WY=ddply(dat2,c("Farm.ID","WY"),summarise,TLoad.mt=sum(load.mt,na.rm=T),TFlow.AcFtd=sum(cfs.to.acftd(Q.cfs),na.rm=T))
farm.level.load.WY$FWM=with(farm.level.load.WY,TLoad.mt/af.to.liters(TFlow.AcFtd))*1e9

subset(farm.level.load.WY,FWM>1500)

plot(FWM~WY,farm.level.load.WY)
plot(FWM~WY,subset(farm.level.load.WY,Farm.ID=="26-006-01"))
plot(FWM~WY,subset(farm.level.load.WY,Farm.ID=="50-048-01"),log="y")


farm.level.load=ddply(dat2,c("Farm.ID","date2"),summarise,meanTP=mean(TP.mgL,na.rm=T),TFlow=sum(Q.cfs,na.rm=T))
farm.level.load$WY=WY(farm.level.load$date2)
farm.level.load$load.mt=with(farm.level.load,Load.Calc.kg(TFlow,meanTP))

farm.level.load.WY2=ddply(farm.level.load,c("Farm.ID","WY"),summarise,TLoad.mt=sum(load.mt,na.rm=T),TFlow.AcFtd=sum(cfs.to.acftd(TFlow),na.rm=T))
farm.level.load.WY2$FWM=with(farm.level.load.WY2,TLoad.mt/af.to.liters(TFlow.AcFtd))*1e9

vars=c('Farm.ID',"WY","FWM")
FWM.compare=merge(farm.level.load.WY[,vars],farm.level.load.WY2[,vars],c('Farm.ID',"WY"))
FWM.compare=merge(FWM.compare,ddply(farm.struct.xwalk,"Farm.ID",summarise,N.val=N.obs(StructureName)),"Farm.ID",all.x=T)

vars=c('Farm.ID',"WY","TLoad.mt")
load.compare=merge(farm.level.load.WY[,vars],farm.level.load.WY2[,vars],c('Farm.ID',"WY"))
load.compare=merge(load.compare,ddply(farm.struct.xwalk,"Farm.ID",summarise,N.val=N.obs(StructureName)),"Farm.ID",all.x=T)


plot(FWM.x~FWM.y,subset(FWM.compare,N.val>1));abline(0,1)


plot(TLoad.mt.x~TLoad.mt.y,subset(load.compare,N.val>1));abline(0,1)


# png(filename=paste0(plot.path,"LoadFWM_Compare.png"),width=6.5,height=3.5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,3,0.5,1),oma=c(2,1,1,0.25));
layout(matrix(1:2,1,2))

ylim.val=c(0,2e5);by.y=0.5e5;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=ylim.val;by.x=by.y;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
plot(TLoad.mt.y~TLoad.mt.x,subset(load.compare,N.val>1),ylim=ylim.val,xlim=ylim.val,ann=F,axes=F,type="n");
abline(h=ymaj,v=xmaj,lty=3,col="grey",lwd=0.5)
points(TLoad.mt.y~TLoad.mt.x,subset(load.compare,N.val>1),pch=21,bg=adjustcolor("indianred1",0.5),col="grey",lwd=0.1,cex=0.8)
abline(0,1)
axis_fun(1,xmaj,xmin,format(xmaj/1e5),line=-0.5)     
axis_fun(2,ymaj,ymin,format(ymaj/1e5));box(lwd=1)
mtext(side=3,adj=0,"Annual TP Load",font=2)
mtext(side=2,line=2.25,"Method 1: Avg Conc \u00D7 Q (x10\u2074 mt WY\u207B\u00B9)")
mtext(side=1,line=1.75,"Method 2: \u2211 Conc \u00D7 Q (x10\u2074 mt WY\u207B\u00B9)")

ylim.val=c(0,2000);by.y=500;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=ylim.val;by.x=by.y;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
plot(FWM.y~FWM.x,subset(FWM.compare,N.val>1),ylim=ylim.val,xlim=ylim.val,ann=F,axes=F,type="n");
abline(h=ymaj,v=xmaj,lty=3,col="grey",lwd=0.5)
points(FWM.y~FWM.x,subset(FWM.compare,N.val>1),pch=21,bg=adjustcolor("dodgerblue1",0.5),col="grey",lwd=0.1,cex=0.8)
abline(0,1)
axis_fun(1,xmaj,xmin,format(xmaj),line=-0.5)     
axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=3,adj=0,"Annual TP FWM",font=2)
mtext(side=2,line=2.5,"Method 1: Avg Conc \u00D7 Q (\u03BCg L\u207B\u00B9)")
mtext(side=1,line=1.75,"Method 2: \u2211 Conc \u00D7 Q (\u03BCg L\u207B\u00B9)")
dev.off()