################
#Analyze and Summarize the station climate projections from CSAG
#
library(reshape);library(rgdal);library(multicore);library(ncdf4);library(raster)

setwd("/media/Data/Work/Regional/CFR/Weather/ClimateInterpolation")

### load some functions
source("/media/Data/Work/projects/WeatherInterpolation/analysis/WeatherData_functions.R")

## get file names
files=list.files("../data/obs/20110712_StationData",recursive=T,full=T,pattern="[.]txt")
files=files[!grepl("station|error|nc|[.]R|grid|[.]pdf|[.]gz",files)]

## Build station list
meta=rbind.fill(mclapply(files,importmeta,verbose=F))

## get time series
origin=as.Date("1800-01-01")
time=as.integer(seq(min(meta$datestart),max(meta$datestop),1)-origin+1)
meta$startid=match(as.numeric(meta$datestart-origin+1),time)
meta$count=as.numeric(meta$datestop-meta$datestart+1)	

### construct table of unique stations
umeta=unique(meta[,c("station","lon","lat","elevation")])


#### make a netcdf file
## define variables
d_station=ncdim_def("station",units="",vals=1:nrow(umeta),unlim=T,create_dimvar=T,longname="measurement locations")
d_time=ncdim_def("time",units=paste("days since",origin,"00:00:00"),vals=time,unlim=F,create_dimvar=T,longname="time of measurement")
v_lat=ncvar_def("lat",units="degrees_north",dim=d_station,missval=-999,longname="station latitude")
v_lon=ncvar_def("lon",units="degrees_east",dim=d_station,missval=-999,longname="station longitude")
## Station metadata
v_stationID=ncvar_def("stationID",units="Station ID",dim=list(d_station),missval="NA",longname="Station ID",prec="char")
v_stationName=ncvar_def("stationName",units="Station Name",dim=list(d_station),missval="NA",longname="Station Name",prec="char")
v_stationElev=ncvar_def("elevation",units="m",dim=list(d_station),missval=-999,longname="Station Elevation")
v_stationStart=ncvar_def("timestart",units=paste("days since",origin,"00:00:00"),dim=list(d_station),missval=-999,
  longname="Time of first record")
v_stationStop=ncvar_def("timestop",units=paste("days since",origin,"00:00:00"),dim=list(d_station),missval=-999,
  longname="Time of last record")
### Weather data
v_tmax=ncvar_def("tmax",units="degrees_c",dim=list(d_time,d_station),missval=-999,longname="Maximum Temperature",compression=9)
v_tmin=ncvar_def("tmin",units="degrees_c",dim=list(d_time,d_station),missval=-999,longname="Minimum Temperature",compression=9)
v_ppt=ncvar_def("ppt",units="mm",dim=list(d_time,d_station),missval=-999,longname="Total Daily Precipitation",compression=9)

### Write the file
file.remove("StationData.nc")
nc_create("StationData.nc",list(v_tmax,v_tmin,v_ppt,v_lon,v_lat,v_stationName,v_stationID,v_stationElev,v_stationStart,v_stationStop))


### open file for writing
nc=nc_open("StationData.nc",write=T)
## Global Attributes
ncatt_put(nc,varid=0, "title","Daily weather station observations for the Cape Floristic Region of South Africa")
ncatt_put(nc,varid=0, "institution","Silander Lab, Ecology and Evolutionary Biology, University of Connecticut, Storrs, CT, USA")
ncatt_put(nc,varid=0, "source","Blurb, B. 2011. Daily dataset of 20th-century surface air temperature and precipitation series etc. 
Data and metadata available at http://www.csag.uct.ac.za")
ncatt_put(nc,varid=0, "references","http://www.csag.uct.ac.za/")
ncatt_put(nc,varid=0, "comment","All observations that failed quality control were removed from the dataset.  Data converted from station-by-station files to netCDF by Adam Wilson (adam.wilson@uconn.edu)")
## add meta data
ncvar_put(nc,"stationID",vals=umeta$station)
ncvar_put(nc,"stationName",vals=meta$stationname[match(umeta$station,meta$station)])
ncvar_put(nc,"elevation",vals=as.numeric(umeta$elevation))
## coordinates
ncvar_put(nc,"lat",vals=as.numeric(umeta$lat))
ncvar_put(nc,"lon",vals=as.numeric(umeta$lon))


## Fill the file with missing values (instead of default 0s) - takes a few minutes
lapply(c("tmax","tmin","ppt"),function(x) ncvar_put(nc,x,vals=rep(NA,prod(nc[["var"]][[x]][["size"]]))))

## add the data file by file from the meta table
lapply(1:nrow(meta),function(i){
  ncvar_put(nc,as.character(meta$var[i]),
            vals=importdata(meta$path[i]),
            start=c(meta$startid[i],match(as.character(meta$station[i]),as.character(umeta$station))),
            count=c(meta$count[i],1)
            )
  print(paste("Finished ",i, " out of ",nrow(meta)))
})

nc_close(nc)
#  system("ncpdq StationData_raw.nc StationData_raw_packed.nc")  #makes the file bigger!

#################################################################
#################################################################
### Build Climate NetCDF

### load CFR Boundary
cfr=readOGR("/media/Data/Work/Regional/CFR/BaseGISData/CFR.shp","CFR")

### Copy Schulze datasets to folder and reproject to WGS84
spath="/media/Data/Work/Regional/ZA/Schulze2007/GISData/grids/"
f=list.files(spath,pattern="tminave|tmaxave|gmednrfl",full=F)
f=data.frame(file=f[!grepl("aux",f)])
f$var=ifelse(grepl("tmin",f$file),"tmin",ifelse(grepl("tmax",f$file),"tmax",ifelse(grepl("gmednrfl",f$file),"ppt",NA)))
f$month=as.numeric(gsub("[a-z]","",f$file))
  
### Crop to CFR
cfrbbox=round(as.vector(bbox(cfr))+c(-.1,-.1,.1,.1),4)

mclapply(f$file,function(i)
  system(paste("gdalwarp -overwrite -ot Float32 -s_srs EPSG:4222 -tr .016667 .016667 -t_srs EPSG:4326 -dstnodata -3.40282e+38 -multi -r cubic -te ",
               paste(cfrbbox,collapse=" ")," ",spath,i," climate/",i,".tif",sep="")))
  
### Connect to monthly datasets
library(raster)

tr=readGDAL(paste("climate/",f$file[1],".tif",sep=""))

### take monthly climate surfaces and build a subsetted dataset for prediction and calculation of anomolies 
d_time=ncdim_def("time",units=paste("months since",origin-31,"00:00:00"),vals=1:12,unlim=F,calendar="standard",create_dimvar=T,longname="time")
d_lat=ncdim_def("lat",units="degrees_north",vals=sort(unique(coordinates(tr)[,2]),decreasing=F),longname="latitude")
d_lon=ncdim_def("lon",units="degrees_east",vals=sort(unique(coordinates(tr)[,1])),longname="longitude")

## Station metadata
v_tmax=ncvar_def("tmax",units="degrees_c",dim=list(d_lon,d_lat,d_time),missval=-999,longname="Mean Maximum Temperature",compression=9)
v_tmin=ncvar_def("tmin",units="degrees_c",dim=list(d_lon,d_lat,d_time),missval=-999,longname="Mean Minimum Temperature",compression=9)
v_ppt=ncvar_def("ppt",units="mm",dim=list(d_lon,d_lat,d_time),missval=-999,longname="Mean Monthly Precipitation",compression=9)

file.remove("Climate.nc")
nc_create("Climate.nc",list(v_tmax,v_tmin,v_ppt))
cl=nc_open("Climate.nc",write=T)

### Put in the data
  for (i in 1:nrow(f)){
    tr=readGDAL(paste("climate/",f$file[i],".tif",sep=""))
    ncvar_put(cl,as.character(f$var[i]),
            vals=as.matrix(tr)[,ncol(as.matrix(tr)):1],  #nrow(as.matrix(tr)):1
            start=c(1,1,f$month[i]),
            count=c(-1,-1,1)
            )
  }

  ## Global Attributes
ncatt_put(cl,varid=0, "title","Mean monthly weather variables for the Greater Cape Floristic Region of South Africa")
ncatt_put(cl,varid=0, "institution","Silander Lab, Ecology and Evolutionary Biology, University of Connecticut, Storrs, CT, USA")
ncatt_put(cl,varid=0, "source","Original Data from Schulze, 2007")
ncatt_put(cl,varid=0, "comment","Generated by Adam Wilson (adam.wilson@uconn.edu)")
ncatt_put(cl,varid=0, "Convention","CF-1.5")

### Variable attributes
ncatt_put(cl,varid="lon", "_CoordinateAxisType","lon")
ncatt_put(cl,varid="lon", "axis","X")
ncatt_put(cl,varid="lon", "standard_name","longitude")
ncatt_put(cl,varid="lat", "_CoordinateAxisType","lat")
ncatt_put(cl,varid="lat", "axis","Y")
ncatt_put(cl,varid="lat", "standard_name","latitude")
ncatt_put(cl,varid="time", "_CoordinateAxisType","time")
ncatt_put(cl,varid="time", "axis","T")
ncatt_put(cl,varid="tmax", "_CoordinateAxis","time lat lon")
ncatt_put(cl,varid="tmax", "coordinates","time lat lon")
nc_sync(cl)

system("ncdump -h Climate.nc")

nc_close(cl)

## Create CFR mask
cgrid=as(raster("Climate.nc",varname="tmax",layer=1),"SpatialGridDataFrame")
proj4string(cgrid)=CRS("+proj=longlat +datum=WGS84")
cgrid$id=1:nrow(cgrid)
cgrid$cfr=ifelse(cgrid$id%in%extract(raster(cgrid,layer="id"),cfr,weights=T,small=T)[[1]][,1],1,NA) #overlay CFR
## get distance to cfr to keep buffer of pixels around region - important to prevent extrapolation later
cdist=as(distance(raster(cgrid,layer="cfr")),"SpatialGridDataFrame")
cdist@data=data.frame(cfrdist_km=cdist$values/1000) #convert to km
writeRaster(raster(cdist),"CFR.nc",over=T,varname="CFR_dist",varunit="km",longname="Distance to Cape Floristic Region",yname="lat",xname="lon")

###########################################################################################
###########################################################################################
## Create lower resolution (0.05 degree ~ 5km) climate dataset for interpolations

### load CFR Boundary
cfr=readOGR("/media/Data/Work/Regional/CFR/BaseGISData/CFR.shp","CFR")

### Connect to monthly datasets
library(raster)
f$coarsefile=sub("c$","",sub("ave","",sub("gmednrfl","ppt",f$file)))

## create 0.05 files in g.climate script using grass
tr=readGDAL(paste("climate_0.05/",f$coarsefile[1],".tif",sep=""))

### take monthly climate surfaces and build a subsetted dataset for prediction and calculation of anomolies 
d_time=ncdim_def("time",units=paste("months since",origin-31,"00:00:00"),vals=1:12,unlim=F,calendar="standard",create_dimvar=T,longname="time")
d_lat=ncdim_def("lat",units="degrees_north",vals=sort(unique(coordinates(tr)[,2]),decreasing=F),longname="latitude")
d_lon=ncdim_def("lon",units="degrees_east",vals=sort(unique(coordinates(tr)[,1])),longname="longitude")

## Station metadata
v_tmax=ncvar_def("tmax",units="degrees_c",dim=list(d_lon,d_lat,d_time),missval=-999,longname="Mean Maximum Temperature",compression=9)
v_tmin=ncvar_def("tmin",units="degrees_c",dim=list(d_lon,d_lat,d_time),missval=-999,longname="Mean Minimum Temperature",compression=9)
v_ppt=ncvar_def("ppt",units="mm",dim=list(d_lon,d_lat,d_time),missval=-999,longname="Mean Monthly Precipitation",compression=9)

file.remove("Climate_0.05.nc")
nc_create("Climate_0.05.nc",list(v_tmax,v_tmin,v_ppt))
cl=nc_open("Climate_0.05.nc",write=T)

### Put in the data
  for (i in 1:nrow(f)){
    tr=readGDAL(paste("climate_0.05/",f$coarsefile[i],".tif",sep=""))
    ncvar_put(cl,as.character(f$var[i]),
            vals=as.matrix(tr)[,ncol(as.matrix(tr)):1],  #nrow(as.matrix(tr)):1
            start=c(1,1,f$month[i]),
            count=c(-1,-1,1)
            )
  }

  ## Global Attributes
ncatt_put(cl,varid=0, "title","Mean monthly weather variables for the Greater Cape Floristic Region of South Africa")
ncatt_put(cl,varid=0, "institution","Silander Lab, Ecology and Evolutionary Biology, University of Connecticut, Storrs, CT, USA")
ncatt_put(cl,varid=0, "source","Original Data from Schulze, 2007")
ncatt_put(cl,varid=0, "comment","Generated by Adam Wilson (adam.wilson@uconn.edu)")
ncatt_put(cl,varid=0, "Convention","CF-1.5")

### Variable attributes
ncatt_put(cl,varid="lon", "_CoordinateAxisType","lon")
ncatt_put(cl,varid="lon", "axis","X")
ncatt_put(cl,varid="lon", "standard_name","longitude")
ncatt_put(cl,varid="lat", "_CoordinateAxisType","lat")
ncatt_put(cl,varid="lat", "axis","Y")
ncatt_put(cl,varid="lat", "standard_name","latitude")
ncatt_put(cl,varid="time", "_CoordinateAxisType","time")
ncatt_put(cl,varid="time", "axis","T")
ncatt_put(cl,varid="tmax", "_CoordinateAxis","time lat lon")
ncatt_put(cl,varid="tmax", "coordinates","time lat lon")
nc_sync(cl)

system("ncdump -h Climate_0.05.nc")

nc_close(cl)


## Create CFR mask
cgrid=as(raster("Climate_0.05.nc",varname="tmax",layer=1),"SpatialGridDataFrame")
proj4string(cgrid)=CRS("+proj=longlat +datum=WGS84")
cgrid$id=1:nrow(cgrid)
cgrid$cfr=ifelse(cgrid$id%in%extract(raster(cgrid,layer="id"),cfr,weights=T,small=T)[[1]][,1],1,NA) #overlay CFR
## get distance to cfr to keep buffer of pixels around region - important to prevent extrapolation later
cdist=as(distance(raster(cgrid,layer="cfr")),"SpatialGridDataFrame")
cdist@data=data.frame(cfrdist_km=cdist$V1/1000) #convert to km
writeRaster(raster(cdist),"CFR_0.05.nc",over=T,varname="CFR_dist",varunit="km",longname="Distance to Cape Floristic Region",yname="lat",xname="lon")
