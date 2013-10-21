################
#Analyze and Summarize the station climate projections from CSAG
#

library(ncdf4)

setwd("/media/Data/Work/Regional/CFR/Weather/ClimateInterpolation")

### open file 
nc=nc_open("StationData_raw.nc",write=F)

### import metadata
time=ncvar_get(nc,"time")
lat=as.vector(ncvar_get(nc,"lat"))
lon=as.vector(ncvar_get(nc,"lon"))

### import data
tmax=ncvar_get(nc,"tmax")
tmin=ncvar_get(nc,"tmin")
ppt=ncvar_get(nc,"ppt")

### Number of stations through time
st_time=data.frame(time=time,date=(as.Date("1800-01-01")+time),
  tmax=apply(tmax,1,function(x) sum(!is.na(x))),
  tmin=apply(tmin,1,function(x) sum(!is.na(x))),
  ppt=apply(ppt,1,function(x) sum(!is.na(x)))
  )

### mean values through time
time_mean=data.frame(time=time,date=(as.Date("1800-01-01")+time),
  tmax=apply(tmax,1,mean, na.rm=T),
  tmin=apply(tmin,1, mean, na.rm=T),
  ppt=apply(ppt,1, mean, na.rm=T)
  )


###### Draw some plots

pdf("output/stationTimeline.pdf",width=11,height=8.5,useDingbats=F)

#### Some summary plots
plot(ppt~date,data=st_time,type="l",col="blue",ylab="Observations per day",las=1,xlab="Date",xlim=c(as.Date("1980-01-01"),as.Date("2010-12-31")))
lines(tmin~date,data=st_time,type="l",col="orange")
lines(tmax~date,data=st_time,type="l",col="red")

plot(ppt~date,data=time_mean,type="l",col="blue",ylab="Mean",las=1,xlab="Date")
lines(tmin~date,data=time_mean,type="l",col="orange")
lines(tmax~date,data=time_mean,type="l",col="red")


spplot(obsbins2,zcol=2:7,asp=1,col.regions=c("blue","purple","purple","red","red"),pch=c(1,13,13,4,4))

dev.off()
system("ps2pdf output/stationTimeline.pdf output/stationTimeline_small.pdf")



