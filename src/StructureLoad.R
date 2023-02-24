## Title:      EAA Hotspot Data analysis - S5,S6,S7,S8 structures
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
utm17=CRS("+init=epsg:26917")

# -------------------------------------------------------------------------
dates=date.fun(c("1978-05-01","2022-04-30"))

Q.sites=data.frame(DBKEY=c("15040","AN602","15037","15034","TA868","00319","AN791"),
                   Site=c("S8","S7Temp","S7","S6","S5A_S5AS","S5A_S5AS","S5ATemp"),
                   Struct=c("S8","S7","S7","S6","S5A","S5A",'S5A'),
                   Priority=c("P1",'P1',"P1","P1","P2","P1","P1"))

Q.dat=data.frame()
for(i in 1:length(Q.sites$DBKEY)){
  tmp=DBHYDRO_daily(dates[1],dates[2],Q.sites$DBKEY[i])
  tmp$DBKEY=as.character(Q.sites$DBKEY[i])
  Q.dat=rbind(Q.dat,tmp)
  print(i)
}
  
Q.dat=merge(Q.dat,Q.sites,"DBKEY")
Q.dat$Date.EST=date.fun(Q.dat$Date)

Q.dat.xtab=dcast(Q.dat,Site+Struct+Date.EST~Priority,value.var="Data.Value",mean)
Q.dat.xtab$Q.cfs=with(Q.dat.xtab,ifelse(is.na(P1),P2,P1))
Q.dat.xtab$Inflow=with(Q.dat.xtab,ifelse(Q.cfs<0,0,Q.cfs))

Q.dat.da=ddply(Q.dat.xtab,c("Struct","Date.EST"),summarise,Q.cfs=sum(Inflow,na.rm=T))

plot(Q.cfs~Date.EST,subset(Q.dat.da,Struct=='S5A'))
plot(Q.cfs~Date.EST,subset(Q.dat.da,Struct=='S6'))
plot(Q.cfs~Date.EST,subset(Q.dat.da,Struct=='S7'))
plot(Q.cfs~Date.EST,subset(Q.dat.da,Struct=='S8'))


# 
WQ.sites=unique(Q.sites$Struct)
#TP only
wq.dat=data.frame()
for(i in 1:length(WQ.sites)){
  tmp=DBHYDRO_WQ(dates[1],dates[2],WQ.sites[i],25)
  wq.dat=rbind(wq.dat,tmp)
  print(i)
}


head(wq.dat)

wq.dat.xtab=dcast(wq.dat,Station.ID+Date.EST~Collection.Method,value.var="HalfMDL",mean)

wq.dat.xtab[,c("Station.ID","Date.EST","G","ACF","ACT")]

plot(ACT~G,wq.dat.xtab);abline(0,1)

WQ.Q=merge(Q.dat.da,
           wq.dat.xtab[,c("Station.ID","Date.EST","G","ACF","ACT")],
           by.x=c("Struct","Date.EST"),by.y=c("Station.ID","Date.EST"),all.x=T)
# If flow within the last 7-days then use ACF otherwise G/ACT
WQ.Q$Q.window=with(WQ.Q,ave(Q.cfs,Struct,FUN=function(x) c(rep(NA,6),rollapply(x>0,width=7,sum,na.rm=T))))
WQ.Q$TP=with(WQ.Q,ifelse(Q.cfs>0&is.na(ACF)==F&Q.window>0,ACF,ifelse(is.na(ACT)==T,G,ACT)))
# WQ.Q$TP=with(WQ.Q,ifelse(Q.cfs>0&is.na(ACF)==F,ACF,ifelse(is.na(ACT)==T,G,ACT)))
# subset(WQ.Q,Struct=="S5A"&Date.EST==date.fun("1979-11-05")) #spot check

WQ.Q$TP.inter=with(WQ.Q,ave(TP,Struct,FUN=function(x)dat.interp(x)))
WQ.Q$TP.load=with(WQ.Q,Load.Calc.kg(Q.cfs,TP.inter))
WQ.Q$WY=WY(WQ.Q$Date.EST)

### FWM values
cfs.to.Ld=function(x) x*(28.3168/(1/86400))
af.to.liters=function(x) x* 1233481.85532

WQ.Q.WY=ddply(WQ.Q,c("Struct","WY"),summarise,TLoad=sum(TP.load,na.rm=T),TFlow.AcFtd=sum(cfs.to.acftd(Q.cfs),na.rm=T))
WQ.Q.WY$FWM=with(WQ.Q.WY,TLoad/af.to.liters(TFlow.AcFtd))*1e9

plot(TFlow.AcFtd~WY,subset(WQ.Q.WY,Struct=="S5A"),type="l")
plot(TFlow.AcFtd~WY,subset(WQ.Q.WY,Struct=="S6"),type="l")
plot(TFlow.AcFtd~WY,subset(WQ.Q.WY,Struct=="S7"),type="l")
plot(TFlow.AcFtd~WY,subset(WQ.Q.WY,Struct=="S8"),type="l")



# png(filename=paste0(plot.path,"EAAStruct_FWM.png"),width=4.5,height=5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,3,0.5,1),oma=c(2,1.5,0.5,0.25));
layout(matrix(1:4,4,1))

# ylim.val=c(10,550);by.y=100;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
ylim.val=c(10,600);ymaj=log.scale.fun(ylim.val,"major");ymin=log.scale.fun(ylim.val,"minor")
xlim.val=c(1979,2022);by.x=10;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/by.x)
struct.vals=c("S5A","S6","S7","S8")
struct.labs=c("S-5A/S-5AS","S-6","S-7","S-8")
for(i in 1:4){
plot(FWM~WY,WQ.Q.WY,ylim=ylim.val,xlim=xlim.val,ann=F,axes=F,type="n",log="y");
abline(h=ymaj,v=xmaj,lty=3,col="grey",lwd=0.5)
with(subset(WQ.Q.WY,Struct==struct.vals[i]),pt_line(WY,FWM,2,"dodgerblue2",1,21,"dodgerblue2",pt.lwd=0.01,cex=1.25))
axis_fun(1,xmaj,xmin,xmaj,line=-0.5)     
axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
mtext(side=3,adj=0,line=-1.25,paste0(" ",struct.labs[i]),font=3,cex=0.75)
}
# mtext(side=2,line=2.5,"Annual TP FWM (\u03BCg L\u207B\u00B9)")
mtext(side=2,line=-0.5,outer=T,"Annual TP FWM (\u03BCg L\u207B\u00B9)")
mtext(side=1,line=1.5,"Water Year")
dev.off()


library(changepoint)
tmp.dat=subset(WQ.Q.WY,Struct==struct.vals[1])
ansmean=cpt.mean(tmp.dat$FWM,method="AMOC")
tmp.dat[ansmean@cpts[1],]
plot(ansmean)



WQ.Q.WY.2000=subset(WQ.Q.WY,WY>1999)
FWM.trend=ddply(WQ.Q.WY.2000,c("Struct"),summarise,
      est=cor.test(FWM,WY,method="kendall")$estimate,
      pval=cor.test(FWM,WY,method="kendall")$p.value)

struct.vals=c("S5A","S6","S7","S8")
pt.change=data.frame()
for(i in 1:length(struct.vals)){
  tmp.dat=subset(WQ.Q.WY.2000,Struct==struct.vals[i])
  pt.rslt=trend::pettitt.test(tmp.dat$FWM)
  
  rslt.sum=data.frame(Struct=struct.vals[i],
                      pt.stat=pt.rslt$statistic,
                      pt.pval=pt.rslt$p.value,
                      pt.change=tmp.dat$WY[pt.rslt$estimate])
  pt.change=rbind(pt.change,rslt.sum)
}

# png(filename=paste0(plot.path,"EAAStruct_FWM_2000.png"),width=4.5,height=5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,3,0.5,1),oma=c(2,1.5,0.5,0.25));
layout(matrix(1:4,4,1))

ylim.val=c(10,300);ymaj=log.scale.fun(ylim.val,"major");ymin=log.scale.fun(ylim.val,"minor")
xlim.val=c(2000,2022);by.x=5;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/by.x)
struct.vals=c("S5A","S6","S7","S8")
struct.labs=c("S-5A/S-5AS","S-6","S-7","S-8")
for(i in 1:4){
  plot(FWM~WY,WQ.Q.WY,ylim=ylim.val,xlim=xlim.val,ann=F,axes=F,type="n",log="y");
  abline(h=ymaj,v=xmaj,lty=3,col="grey",lwd=0.5)
  with(subset(WQ.Q.WY,Struct==struct.vals[i]),pt_line(WY,FWM,2,"dodgerblue2",1,21,"dodgerblue2",pt.lwd=0.01,cex=1.25))
  trend.dat=subset(FWM.trend,Struct==struct.vals[i])
  pt.dat=subset(pt.change,Struct==struct.vals[i])
  if(trend.dat$pval<0.05){
    mod=mblm::mblm(FWM~WY,subset(WQ.Q.WY.2000,Struct==struct.vals[i]),repeated=F)
    x.val=seq(2000,2022,1)
    pred.mod=predict(mod,data.frame(WY=x.val),interval="confidence")
    pred.mod2=ifelse(pred.mod<10,10,pred.mod)
    shaded.range(x.val,pred.mod2[,2],pred.mod2[,3],bg="grey",col.adj=0.5,lty=0)
    lines(x.val,pred.mod[,1])
    lines(x.val,pred.mod[,2],lty=2,lwd=0.5)
    lines(x.val,pred.mod[,3],lty=2,lwd=0.5)
  }
  if(pt.dat$pt.pval<0.05){abline(v=pt.dat$pt.change,col="indianred1")}
  
  axis_fun(1,xmaj,xmin,xmaj,line=-0.5)     
  axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
  mtext(side=3,adj=0,line=-1.25,paste0(" ",struct.labs[i]),font=3,cex=0.75)
}
# mtext(side=2,line=2.5,"Annual TP FWM (\u03BCg L\u207B\u00B9)")
mtext(side=2,line=-0.5,outer=T,"Annual TP FWM (\u03BCg L\u207B\u00B9)")
mtext(side=1,line=1.5,"Water Year")
dev.off()