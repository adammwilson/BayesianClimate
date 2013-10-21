################
#Analyze and Summarize the station climate projections from CSAG
#
system("ssh -X protea")
R
library(ncdf4);library(raster);library(rgdal);library(reshape);library(geoR);library(sp)

setwd("/media/Data/Work/Regional/CFR/Weather/ClimateInterpolation")

### open file 
nc=nc_open("StationData.nc",write=F)
cl=nc_open("Climate.nc",write=F) 

### import metadata
time=ncvar_get(nc,"time")
lat=as.vector(ncvar_get(nc,"lat"))
lon=as.vector(ncvar_get(nc,"lon"))

### load CFR Boundary
cfr=readOGR("/media/Data/Work/Regional/CFR/BaseGISData/CFR.shp","CFR")

### Dates
### Subset to period of interest
start=as.Date("1990-01-01")
stop=as.Date("2009-12-31")
dates=seq(start,stop,1)
origin=as.Date("1800-01-01")  # used as the reference date in netcdfs

### Make variable table
vars=data.frame(
  name=c("ppt","tmax","tmin"),
  longname=c("Total Precipitation","Maxiumum Temperature","Miniumum Temperature"),
  unit=c("mm","C","C"),stringsAsFactors=F)
vars  

### Jobs object
jobs=expand.grid(date=dates,var=vars$name)
save(jobs,file="jobs.Rdata")

## update notdone to include all jobs
notdone=jobs
save(notdone,file="notdone.Rdata")

### Make prediction grid
bbox(cfr)
res=0.25 #125 #grid resolution
pgrid=expand.grid(lon=seq(17.5,26,by=res),lat=seq(-35.25,-30.5,by=res))
pgrid=SpatialPixelsDataFrame(pgrid,data=data.frame(id=1:nrow(pgrid),lon=pgrid$lon,lat=pgrid$lat),proj4string = CRS("+proj=longlat +datum=WGS84"))
pgrid$cfr=ifelse(pgrid$id%in%extract(raster(pgrid),cfr,weights=T,small=T)[[1]][,1],1,NA) #overlay CFR
## get distance to cfr to keep buffer of pixels around region - important to prevent extrapolation later
pdist=as(distance(raster(pgrid,layer="cfr")),"SpatialPixelsDataFrame")
pgrid$dcfr=pdist$values[overlay(pdist,pgrid)]/1000  #add distances in km
pgrid$cfr[pgrid$dcfr<=res*300]=1  #keep cells within X km of CFR
pgrid=pgrid[!is.na(pgrid$cfr),] #drop pixels that don't touch CFR
## update the index
pgrid$id=1:nrow(pgrid)

## check out the region
plot(raster(pgrid,"dcfr"),asp=1)
plot(cfr,add=T)

## number of prediction grid cells
pn <- nrow(coordinates(pgrid)) ;print(paste(pn," prediction cells"))

## save objects to send to cluster
save(pgrid,cfr,vars,dates,origin,res,file="Interp.Rdata")
writeGDAL(pgrid,"pgrid.tif")
### Fit variogram to get overall parameters


########################################################################
########################################################################
### Copy files to Davis
### copy directory of files to farm
system("scp StationData.nc zzwilson@davisfarm:/group/latimergrp/adam/cfr/interp")
system("scp Climate.nc zzwilson@davisfarm:/group/latimergrp/adam/cfr/interp")
system("scp CFR.nc zzwilson@davisfarm:/group/latimergrp/adam/cfr/interp")
system("scp Interp.Rdata zzwilson@davisfarm:/group/latimergrp/adam/cfr/interp/Interp.Rdata")
system("scp jobs.Rdata zzwilson@davisfarm:/group/latimergrp/adam/cfr/interp/jobs.Rdata")
#system("scp 4_KrigScript.R zzwilson@davisfarm:/group/latimergrp/adam/cfr/interp/4_KrigScript.R")
#system("scp /media/Data/Work/Regional/CFR/Weather/ClimateInterpolation/4_KrigScript.R zzwilson@davisfarm:/group/latimergrp/adam/cfr/interp/4_KrigScript.R")
### Now open 5_RunCluster to actually run the script, or...
### connect to the farm for maintenence things
system("git archive /media/Data/Work/Projects/WeatherInterpolation/.git HEAD:analysis 4_KrigScript.R | tar -x")

system("ssh davisfarm")
cd /group/latimergrp/adam/cfr/interp
module load gcc openmpi R szip hdf netcdf

###  Start R to install any needed packages
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/export/1/atm-apps/netcdf/lib:/share/atm-apps/netcdf/lib  # path for compiling netcdf


R
#library(Rsge)
### Install any needed packages
packages=c("geoR","sp","reshape","raster","plyr")
rlib="/group/latimergrp/R/x86_64-redhat-linux-gnu-library/2.13/"
for(i in packages[3]) install.packages(i,lib=rlib)
install.packages("ncdf4_1.2.tar.gz",lib=rlib)

install.packages("ncdf4",lib=rlib)
install.packages("ncdf",lib=rlib)

