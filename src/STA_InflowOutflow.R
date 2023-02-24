## Title:      EAA Hotspot Data analysis - S5,S6,S7,S8 structures
##             STA inflows
## Created by: Paul Julian (pjulian@evergladesfoundation.org)
## Created on: 02/16/2023

## BAD ## https://www.tidyverse.org/articles/2017/12/workflow-vs-script/
#Clears Everything...start fresh.
rm(list=ls(all=T));cat("\014");dev.off()

# Libraries
# devtools::install_github("SwampThingPaul/AnalystHelper")
library(AnalystHelper);
library(plyr)
library(reshape2)
library(zoo)

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
utm17=CRS("+init=epsg:26917")

# -------------------------------------------------------------------------
dates=date.fun(c("1999-05-01","2021-04-30"))

## STA Flow and WQ Site tables --------------------------------------------
STA1E.Q.sites=data.frame(
  DBKEY=c("TP366","91476","TP367","90974","TP368","91516","TP369","91517"),
  Site=c(rep('S319',2),rep('G311',2),rep('S361',2),rep('S362',2)),
  region=c(rep("inflow",6),rep("outflow",2)),
  STA="STA1E",
  Priority=c(rep(c("P1","P2"),4))
  )
STA1E.Q.sites$WQ=STA1E.Q.sites$Site
STA1W.Q.sites=data.frame(
  DBKEY=c("JW221","90941","M2901","90973","JW222","90934"),
  Site=c(rep("G302",2),rep("G310",2),rep("G251",2)),
  region=c(rep("inflow",2),rep("outflow",4)),
  STA="STA1W",
  Priority=c(rep(c("P1","P2"),3)),
  WQ=c(rep("G302",2),rep("G310",2),rep("ENR012",2))
)
ENR.Q.sites=data.frame(
  DBKEY=c("16222","90933","W3883","90932","JW222","90934"),
  Site=c(rep("G250",2),rep("G250S",2),rep("G251",2)),
  region=c(rep("inflow",4),rep("outflow",2)),
  STA="ENR",
  Priority=c(rep(c("P1","P2"),3))
)
STA2.Q.sites=data.frame(
  DBKEY=c("15034","91663","90327","91202","91201","90328","91208","N0659","91008","90329","91209"),
  Site=c(rep("S6",2),rep("G434",2),"G434S",rep("G435",2),rep("G335",2),rep("G436",2)),
  region=c(rep("inflow",7),rep("outflow",4)),
  STA="STA2",
  Priority=c(rep(c("P1","P2"),2),"P1",rep(c("P1","P2"),3)),
  WQ=c(rep("S6",2),rep("G434",2),"G434",rep("G435",2),rep("G335",2),rep("G436",2)))
STA2.Q.sites=subset(STA2.Q.sites,Site!="G434S");# removed G434 seepage pump
STA34.Q.sites=data.frame(
  DBKEY=c("TA438","91094","TA437","91105",
          "90348","91119","TA582","91120","90349","91121","90359","91122","TA583","91123","90360","91124",
          "90361","91140","TA584","91141","90362","91142","TA585","91143","90363","91144","W3981","91174",
          "90364","91151","TA586","91152","90365","91153","90366","91154","TA587","91155","90367","91156"),
  Site=c(rep("G370",2),rep("G372",2),
         rep("G376A",2),rep("G376B",2),rep("G376C",2),rep("G376D",2),rep("G376E",2),rep("G376F",2),
         rep("G379A",2),rep("G379B",2),rep("G379C",2),rep("G379D",2),rep("G379E",2),rep("G388",2),
         rep("G381A",2),rep("G381B",2),rep("G381C",2),rep("G381D",2),rep("G381E",2),rep("G381F",2)),
  region=c(rep("inflow",4),rep("outflow",36)),
  STA="STA34",
  Priority=c(rep(c("P1","P2"),20)),
  WQ=c(rep("G370",2),rep("G372",2),
       rep("G376B",6),rep("G376E",6),
       rep("G379B",6),rep("G379D",4),rep("G388",2),
       rep("G381B",6),rep("G381E",6))
)
STA34.Q.sites2=data.frame(
  DBKEY=c("W3964","91107","W3965","91108","W3966","91109","W3967","91110","W3968","91111","W3969","91112",
          "W3970","91125","W3971","91126","W3972","91127","W3973","91128","W3974","91129",
          "W3975","91145","W3976","91146","W3977","91147","W3978","91148","W3979","91149","W3980","91150"),
  Site=c(rep("G374A",2),rep("G374B",2),rep("G374C",2),rep("G374D",2),rep("G374E",2),rep("G374F",2),
         rep("G377A",2),rep("G377B",2),rep("G377C",2),rep("G377D",2),rep("G377E",2),
         rep("G380A",2),rep("G380B",2),rep("G380C",2),rep("G380D",2),rep("G380E",2),rep("G380F",2)),
  region=c("inflow"),
  STA="STA34",
  Priority=c(rep(c("P1","P2"),17)),
  WQ=c(rep("G374B",6),rep("G374E",6),
       rep("G377B",6),rep("G377D",4),
       rep("G380B",6),rep("G380E",6))
)
STA34.Q.sites=rbind(STA34.Q.sites2,STA34.Q.sites)#rbind(STA34.Q.sites2,subset(STA34.Q.sites,region!="inflow"))

A1FEB.Q.sites=data.frame(
  DBKEY=c("64172","P8674",
          "92164","92162","92166","92169","92170","92172","92174","92176","92178","92180"),
  Site=c("G721","G720",paste0("G724",LETTERS[1:10])),
  region=c(rep("inflow",2),rep("outflow",10)),
  STA="A1FEB",
  Priority=c("P1"),
  WQ=c("G370","G372",paste0("G724",LETTERS[1:10]))
)

L8FEB.Q.sites=data.frame(
  DBKEY=c("AM245","AO210"),
  Site=c("G538","G539"),
  region=c("inflow","outflow"),
  STA="L8FEB",
  Priority=c("P1")
)
L8FEB.Q.sites$WQ=L8FEB.Q.sites$Site

## To do: Add STA 5/6 Q and WQ site list



## Discharge Data ---------------------------------------------------------
Q.sites=rbind(
  STA1E.Q.sites,
  STA1W.Q.sites,
  STA2.Q.sites,
  STA34.Q.sites,
  A1FEB.Q.sites,
  L8FEB.Q.sites
)

Q.dat=data.frame()
# start_time=Sys.time()
for(i in 1:length(Q.sites$DBKEY)){
  tmp=DBHYDRO_daily(dates[1],dates[2],Q.sites$DBKEY[i])
  tmp$DBKEY=as.character(Q.sites$DBKEY[i])
  Q.dat=rbind(Q.dat,tmp)
  print(i)
}
# end_time=Sys.time()
# end_time - start_time

# length(unique(Q.dat$DBKEY))
# length(unique(Q.sites$DBKEY))

Q.dat=merge(Q.dat,Q.sites,"DBKEY")
Q.dat$Date.EST=date.fun(Q.dat$Date)

Q.dat.xtab=dcast(Q.dat,STA+region+Site+WQ+Date.EST~Priority,value.var="Data.Value",mean)
Q.dat.xtab$fflow.cfs=with(Q.dat.xtab,ifelse(is.na(P1)|is.nan(P1),P2,P1))

range(Q.dat.xtab$fflow.cfs,na.rm=T)
Q.dat.xtab$fflow.cfs=with(Q.dat.xtab,ifelse(fflow.cfs<0,0,fflow.cfs));# positive flow only
Q.dat.xtab$WY=WY(Q.dat.xtab$Date.EST)

#Fix STA 3/4 Inflow volumes post FEB
unique(subset(Q.dat.xtab,STA=="STA34"&Site%in%c("G370","G372"))&WY>=2005)$Site)

Q.dat.xtab$fflow.cfs=with(Q.dat.xtab,ifelse(STA=="STA34"&region=="inflow"&Site%in%c("G370","G372")&WY>=2015,0,fflow.cfs))
Q.dat.xtab$fflow.cfs=with(Q.dat.xtab,ifelse(STA=="STA34"&region=="inflow"&!(Site%in%c("G370","G372"))&WY<2015,0,fflow.cfs))

STA.WY.Q.tot=ddply(Q.dat.xtab,c("STA","region","WY"),summarise,
                   TFlow.AcFt=sum(cfs.to.acftd(fflow.cfs),na.rm=T))

subset(STA.WY.Q.tot,STA=="STA34")
subset(STA.WY.Q.tot,STA=="A1FEB")



## WQ Data ----------------------------------------------------------------
WQ.sites=unique(Q.sites$WQ)

#TP only
wq.dat=data.frame()
for(i in 1:length(WQ.sites)){
  tmp=DBHYDRO_WQ(dates[1],dates[2],WQ.sites[i],25)
  wq.dat=rbind(wq.dat,tmp)
  print(i)
}

## verify we got all WQ sites
# length(WQ.sites)
# length(unique(wq.dat$Station.ID))
# did we miss any?
# WQ.sites[!(WQ.sites%in%unique(wq.dat$Station.ID))]

wq.dat.xtab=dcast(wq.dat,Station.ID+Date.EST~Collection.Method,value.var="HalfMDL",mean)
head(wq.dat.xtab)
wq.dat.xtab[,c("Station.ID","Date.EST","G","ACF")]

plot(ACT~G,wq.dat.xtab);abline(0,1)

WQ.Q=merge(Q.dat.xtab,
           wq.dat.xtab[,c("Station.ID","Date.EST","G","ACF","ACT")],
           by.x=c("WQ","Date.EST"),by.y=c("Station.ID","Date.EST"),all.x=T)
head(WQ.Q)
# If flow within the last 7-days then use ACF otherwise G/ACT
WQ.Q$Q.window=with(WQ.Q,ave(fflow.cfs,Site,FUN=function(x) c(rep(NA,13),rollapply(x>0,width=14,sum,na.rm=T))))
WQ.Q$TP=with(WQ.Q,ifelse(fflow.cfs>0&is.na(ACF)==F&Q.window>0,ACF,ifelse(is.na(ACT)==T,G,ACT)))
# WQ.Q$TP=rowMeans(WQ.Q[,c("G","ACF")],na.rm=T)
# WQ.Q$TP=with(WQ.Q,ifelse(Q.cfs>0&is.na(ACF)==F,ACF,ifelse(is.na(ACT)==T,G,ACT)))
# subset(WQ.Q,Struct=="S5A"&Date.EST==date.fun("1979-11-05")) #spot check

WQ.Q$TP.inter=with(WQ.Q,ave(TP,Site,FUN=function(x)dat.interp(x)))
WQ.Q$TP.load=with(WQ.Q,Load.Calc.kg(fflow.cfs,TP.inter))

### FWM values
cfs.to.Ld=function(x) x*(28.3168/(1/86400))
af.to.liters=function(x) x* 1233481.85532

WQ.Q.WY=ddply(WQ.Q,c("STA","region","WY"),summarise,TLoad=sum(TP.load,na.rm=T),TFlow.AcFtd=sum(cfs.to.acftd(fflow.cfs),na.rm=T))
WQ.Q.WY$FWM=with(WQ.Q.WY,TLoad/af.to.liters(TFlow.AcFtd))*1e9


plot(FWM~WY,subset(WQ.Q.WY,STA=="L8FEB"&region=="inflow"))
plot(FWM~WY,subset(WQ.Q.WY,STA=="L8FEB"&region=="outflow"))

plot(FWM~WY,subset(WQ.Q.WY,STA=="A1FEB"&region=="inflow"))
plot(FWM~WY,subset(WQ.Q.WY,STA=="A1FEB"&region=="outflow"))

plot(FWM~WY,subset(WQ.Q.WY,STA=="STA2"&region=="inflow"))
plot(FWM~WY,subset(WQ.Q.WY,STA=="STA2"&region=="outflow"))

## Compare
STA.sfer=openxlsx::read.xlsx(paste0(data.path,"STA_data.xlsx"),sheet=1)

head(WQ.Q.WY)
head(STA.sfer)


STA.compare.Q=merge(dcast(WQ.Q.WY,WY+STA~region,value.var = "TFlow.AcFtd",mean),
                  STA.sfer[,c("WY","STA","InflowQ.acft","OutflowQ.acft")],c("WY","STA"))


plot(outflow~OutflowQ.acft,STA.compare.Q);abline(0,1)

STA.compare.load=merge(dcast(WQ.Q.WY,WY+STA~region,value.var = "TLoad",mean),
                    STA.sfer[,c("WY","STA","InflowLoad.mt","OutflowLoad.mt")],c("WY","STA"))
STA.compare.load$inflow=STA.compare.load$inflow/1000
STA.compare.load$outflow=STA.compare.load$outflow/1000

plot(inflow~InflowLoad.mt,STA.compare.load);abline(0,1)
plot(outflow~OutflowLoad.mt,STA.compare.load);abline(0,1)



# png(filename=paste0(plot.path,"STA_FlowLoad_Compare.png"),width=6.5,height=5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(2,2,1,0.5),oma=c(2,3,1,0.25));
layout(matrix(1:4,2,2,byrow=T))

ylim.val=c(0,8e5);by.y=2e5;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=ylim.val;by.x=by.y;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)

plot(inflow~InflowQ.acft,STA.compare.Q,ylim=ylim.val,xlim=ylim.val,ann=F,axes=F,type="n");
abline(h=ymaj,v=xmaj,lty=3,col="grey",lwd=0.5)
points(inflow~InflowQ.acft,STA.compare.Q,pch=21,bg=adjustcolor("indianred1",0.5),col="grey",lwd=0.1,cex=1.25)
abline(0,1)
axis_fun(1,xmaj,xmin,format(xmaj/1e5),line=-0.5)     
axis_fun(2,ymaj,ymin,format(ymaj/1e5));box(lwd=1)
mtext(side=3,"Inflow",font=2)
mtext(side=3,adj=0," Annual Flow",line=-1.25,font=3)
mtext(side=2,line=2.25,"Recaculated Q\n(x10\u2074 Ac-Ft WY\u207B\u00B9)")
mtext(side=1,line=1.5,"SFER Reported Q (x10\u2074 Ac-Ft WY\u207B\u00B9)")

plot(outflow~OutflowQ.acft,STA.compare.Q,ylim=ylim.val,xlim=ylim.val,ann=F,axes=F,type="n");
abline(h=ymaj,v=xmaj,lty=3,col="grey",lwd=0.5)
points(outflow~OutflowQ.acft,STA.compare.Q,pch=21,bg=adjustcolor("indianred1",0.5),col="grey",lwd=0.1,cex=1.25)
abline(0,1)
axis_fun(1,xmaj,xmin,format(xmaj/1e5),line=-0.5)     
axis_fun(2,ymaj,ymin,format(ymaj/1e5));box(lwd=1)
mtext(side=3,"Outflow",font=2)
mtext(side=3,adj=0," Annual Flow",line=-1.25,font=3)
mtext(side=1,line=1.5,"SFER Reported Q (x10\u2074 Ac-Ft WY\u207B\u00B9)")

ylim.val=c(0,100);by.y=20;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=ylim.val;by.x=by.y;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
plot(inflow~InflowLoad.mt,STA.compare.load,ylim=ylim.val,xlim=ylim.val,ann=F,axes=F,type="n");
abline(h=ymaj,v=xmaj,lty=3,col="grey",lwd=0.5)
points(inflow~InflowLoad.mt,STA.compare.load,pch=21,bg=adjustcolor("indianred1",0.5),col="grey",lwd=0.1,cex=1.25)
abline(0,1)
axis_fun(1,xmaj,xmin,format(xmaj),line=-0.5)     
axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=3,adj=0," Annual Load",line=-1.25,font=3)
mtext(side=2,line=2.25,"Recaculated TP Load\n(tons WY\u207B\u00B9)")
mtext(side=1,line=1.5,"SFER Reported Load (tons WY\u207B\u00B9)")

ylim.val=c(0,30);by.y=10;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=ylim.val;by.x=by.y;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
plot(outflow~OutflowLoad.mt,STA.compare.load,ylim=ylim.val,xlim=ylim.val,ann=F,axes=F,type="n");
abline(h=ymaj,v=xmaj,lty=3,col="grey",lwd=0.5)
points(outflow~OutflowLoad.mt,STA.compare.load,pch=21,bg=adjustcolor("indianred1",0.5),col="grey",lwd=0.1,cex=1.25)
abline(0,1)
axis_fun(1,xmaj,xmin,format(xmaj),line=-0.5)     
axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=3,adj=0," Annual Load",line=-1.25,font=3)
mtext(side=1,line=1.5,"SFER Reported Load (tons WY\u207B\u00B9)")
dev.off()