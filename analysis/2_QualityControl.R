#############################################
#############################################
### Data Quality Control
### Following Schulze 2007 data quality control measures.  Quoted lines are from Shulze, 2007.

library(ncdf4)

setwd("/media/Data/Work/Regional/CFR/Weather/ClimateInterpolation")

### open file 
nc=nc_open("StationData_raw.nc",write=T)

### import metadata
lat=ncvar_get(nc,"lat")
lon=ncvar_get(nc,"lon")

### import data
tmax=ncvar_get(nc,"tmax")
tmin=ncvar_get(nc,"tmin")
ppt=ncvar_get(nc,"ppt")


###########################################
### Temperatures in degrees F - none found
### "Certain stations were found to still have their temperature data given in degree Fahrenheit, either for part of, or for the entire record, without this having been indicated. These records were identified by logical checks and then converted to °C."
max(tmax,na.rm=T)
max(tmin,na.rm=T)
###

##########################################
###  Tmax < Tmin
### Tmax should not be greater than Tmin

e1=!is.na(tmin)&!is.na(tmax)&tmin>tmax
paste(sum(e1,na.rm=T)," errors found")

tmax[e1]=NA
tmin[e1]=NA

##############################
### Tmax>45, neighboring mean <=40
### "Where daily maximum temperatures exceeded 45°C, these values were flagged and checked against same-day temperatures of surrounding stations. Where found to be anomalous, these records were deleted"

e2tmax=!is.na(tmax)&tmax>=45  # "error index" cases with temps < 20
e2dates=which(apply(e2tmax,1,any))  #get days with >0 errors
dists=spDists(cbind(lon,lat),cbind(lon,lat),longlat=T) #get distance matrix of all stations

chkneigh_TMAX<-function(date,tmax,lat,lon,nstation=5,upperthresh=45,neighborthresh=40){
  ttmax=tmax[date,] #get one days temp data
  tnotnull=!is.na(ttmax)  #identify nonnull stations
  distmat=apply(dists[tnotnull,tnotnull],1,function(x) rank(x)>1&rank(x)<=(nstation+1))  #get matrix of nstation closest stations for each non-null station
  tmeans=apply(distmat,2,function(x) mean(ttmax[tnotnull][x],na.rm=T))  #get means of 5 closest stations
  e2=rep(F,ncol(tmax)) #vector of falses
  e2[tnotnull][ttmax[tnotnull]>=upperthresh&tmeans<neighborthresh]=T  #set T for errors
  return(e2)
}

e2=matrix(F,ncol=ncol(tmax),nrow=nrow(tmax))  #matrix of Falses
e2[e2dates,]=do.call(rbind,lapply(e2dates,chkneigh_TMAX,tmax=tmax))  #update the errors with a true
attr(e2,"type")="Tmax>45, 5 neighbors<40"
paste(sum(e2,na.rm=T),"errors found")

tmax[e2]=NA

########################################################
### SmallDifference
### In the initial quality control the range between Tmxd and Tmnd on a given day was checked to be ≥ 1°C. In over 13 000 cases the daily range, Tra, was found to be < 1°C. On the assumption that a low temperature range was indicative of moist atmospheric conditions at low temperatures, days with Tra < 1°C were to be flagged as errors if both Tmxd and Tmnd were ≥ 20°C. When both Tmxd and Tmnd were <20°C, surrounding stations were checked for low Tra before deciding whether or not the data were to be flagged as errors."

table((tmax-tmin)<1)
table(tmax>=20&tmin>=20&tmax-tmin<1)  # 132 cases with temps > 20

### delete measurements with temps >= 20
e3=!is.na(tmax)&!is.na(tmin)&tmax>=20&tmin>=20&tmax-tmin<1
attr(e3,"type")="tra<1 tmax>20"

### Check surrounding stations for temps < 20
e3tra=!is.na(tmax)&!is.na(tmin)&tmax<20&tmin<20&tmax-tmin<1  # "error index" cases with temps < 20
e3dates=which(apply(e3tra,1,any))  #get days with >0 errors
dists=spDists(cbind(lon,lat),cbind(lon,lat),longlat=T) #get distance matrix of all stations

chkneigh_TRA<-function(date,tmax,tmin,lat,lon,nstation=5,neighborthresh=5){
  ttra=tmax[date,]-tmin[date,] #get one days temp data
  tnotnull=!is.na(ttra)  #identify nonnull stations
  distmat=apply(dists[tnotnull,tnotnull],1,function(x) rank(x)>1&rank(x)<=(nstation+1))  #get matrix of nstation closest stations for each non-null station
  tmeans=apply(distmat,2,function(x) min(ttra[tnotnull][x],na.rm=T))  #get means of 5 closest stations
  e3=rep(F,ncol(tmax)) #vector of falses
  e3[tnotnull][ttra[tnotnull]<1&tmeans>neighborthresh]=T  #set T for errors
  return(e3)
}

## udpate e3 created above to add errors when neighbors have big tra
e3[e3dates,]=do.call(rbind,lapply(e3dates,chkneigh_TRA,tmax=tmax,tmin=tmin))  #update the errors with a true
paste(sum(e3,na.rm=T),"errors found")

tmax[e3]=NA
tmin[e3]=NA

#########################################################
### RepetitiveStations - tmax, tmin, ppt
### "It is highly unlikely that identical temperatures (i.e. to 0.1°C) occur on more than 2 consecutive days. When identical Tmxd or Tmnd values were recorded for 3 or more consecutive days, those records were flagged as errors and later infilled synthetically from surrounding station values."

findreps=function(x){
  l2=c(NA,NA,x[1:(length(x)-2)])  #lag 2 days
  l1=c(NA,x[1:(length(x)-1)]) #lag 1 day
  tdup=x!=0 & x==l1 & x==l2
  return(tdup)
}

e4=apply(tmax,2,findreps)  #update the errors with a true
paste(sum(e4,na.rm=T),"errors found")
attr(e4,"type")="Repeated value"
tmax[e4]=NA

e5=apply(tmin,2,findreps)  #update the errors with a true
paste(sum(e5,na.rm=T),"errors found")
attr(e5,"type")="Repeated value"
tmin[e5]=NA

e6=apply(ppt,2,findreps)  #update the errors with a true
paste(sum(e6,na.rm=T),"errors found")
attr(e6,"type")="Repeated value"
ppt[e6]=NA



### make a new dataset
file.copy("StationData_raw.nc","StationData_QC.nc",overwrite=T)
nc_QC=nc_open("StationData_QC.nc",write=T)

ncvar_put(nc_QC,"tmax",tmax)
ncvar_put(nc_QC,"tmin",tmin)
ncvar_put(nc_QC,"ppt",ppt)

nc_close(nc_QC)
