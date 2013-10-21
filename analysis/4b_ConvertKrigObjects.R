
##########################################################
## Convert output from GeoR's krige.bayes to netcdf format
## objects must have additional metadata appended in the attr(obj,"meta") location or the script will fail.  See 4_KrigScript.R for more details.

##First read in the arguments listed at the command line
args=(commandArgs(TRUE))
##args is now a list of character vectors
## First check to see if arguments are passed.
## if no args were given, print a warning and stop
if(length(args)==0) {stop("No parameters supplied, you must pass parameters")}

## Then cycle through each element of the list and evaluate the expressions.
eval(parse(text=args))
## now there is an i that corresponds to the row in the notdone object that will be processed.

### Keep all iterations of the predictions?
### F will only keep the mean and sd of the predictions
keepall=F

## load libraries
## running on farm?
farm=!grepl("protea",system("uname -a",intern=T))
if(!farm){
	library(sp);library(geoR);library(ncdf4);library(reshape);library(raster)
	setwd("/media/Data/Work/Regional/CFR/Weather/ClimateInterpolation")
        ##  Setup output directories
        version="v1.3b"
        if(keepall) {
          odir=paste("output/",version,"/nc_all/",sep="")
          if(!file.exists(odir)) dir.create(odir)}
        if(!keepall) {
          odir=paste("output/",version,"/nc/",sep="")
          if(!file.exists(odir)) dir.create(odir)}
      }
## running on farm
if(farm) {
	rlib="/group/latimergrp/R/x86_64-redhat-linux-gnu-library/2.13"
	require(plyr,lib=rlib)
	require(sp,lib.loc=rlib)
	require(geoR,lib.loc=rlib)
	require(ncdf4,lib.loc=rlib)
	require(reshape,lib.loc=rlib)
	require(raster,lib.loc=rlib)
        ##  Setup output directories
        if(keepall) {
          odir="output/nc/"
          if(!file.exists(odir)) dir.create(odir)}
        if(!keepall) {
          odir="output/nc_mean/"
          if(!file.exists(odir)) dir.create(odir)}
      }


## load the list of 'krige objects' to be processed
load("KO.Rdata")

## load the prediction grid, start dates, and cfr border
load("Interp.Rdata")

## load the krige object to be processed (it's called "kc")
load(KO[i])

idate=attr(kc,"meta")$date
vari=attr(kc,"meta")$variable

## reshape function
reshapekrige<-function(obj,grid){
    tobj=obj$predictive$simulations
    dimnames(tobj)=list(id=grid$id,iter=1:ncol(tobj))
    tobj2=melt(tobj)
    tobj2[c("lon","lat")]=grid@data[match(tobj2$id,grid$id),c("lon","lat")]
    tobj3=cast(tobj2,iter~lon~lat)
    return(tobj3)
  }

var_iter=reshapekrige(kc,pgrid) #var_pred
var_mean=as.matrix(SpatialPixelsDataFrame(cbind(lon=pgrid$lon,lat=pgrid$lat),data=data.frame(kc$predictive$mean)))
var_sd=sqrt(as.matrix(SpatialPixelsDataFrame(cbind(lon=pgrid$lon,lat=pgrid$lat),data=data.frame(kc$predictive$variance))))
var_q2.5=apply(var_iter,2:3,quantile,0.025,na.rm=T);var_q2.5=var_q2.5[,ncol(var_q2.5):1]
var_q97.5=apply(var_iter,2:3,quantile,0.975,na.rm=T);var_q2.5=var_q97.5[,ncol(var_q97.5):1]

		## rescale ppt 
if(vari=="ppt") {
  var_mean=(var_mean-1)  #subtract 1 added before kriging (still have to add 1 to climate later)
  var_iter=var_iter-1
  var_q2.5=var_q2.5-1
  var_q97.5=var_q97.5-1
}


## Set dimentions
d_time=ncdim_def("time",units=paste("days since",origin,"00:00:00"),vals=as.integer(idate-origin),unlim=T,calendar="standard",longname="time")
d_lat=ncdim_def("lat",units="degrees_north",vals=sort(unique(pgrid$lat),decreasing=F),longname="latitude")
d_lon=ncdim_def("lon",units="degrees_east",vals=sort(unique(pgrid$lon)),longname="longitude")
d_iter=ncdim_def("iter",units="iterations",vals=1:nrow(kc$posterior$sample),unlim=FALSE)
d_params=ncdim_def("parameters",units="",create_dimvar=F,vals=1:ncol(kc$posterior$sample),unlim=FALSE)  #posterior krige parameter summaries
#d_params2=ncdim_def("parametermoments",units="",create_dimvar=F,vals=1:2,unlim=FALSE)  #posterior krige parameter summaries
d_valid=ncdim_def("validation",units="",create_dimvar=F,vals=1:nrow(attr(kc,"meta")$valid),unlim=FALSE)
d_valid2=ncdim_def("validation2",units="",create_dimvar=F,vals=1:ncol(attr(kc,"meta")$valid),unlim=FALSE)

## Define variables
comp=9  #define compression level: 9 is the highest
if(keepall){
  v_var=ncvar_def(vari,units=vars$unit[vars$name==vari],dim=list(d_iter,d_lon,d_lat,d_time),missval=-999,
    longname=vars$longname[vars$name==vari],compress=comp)
}
  v_var_mean=ncvar_def(paste(vari,"_mean",sep=""),units=vars$unit[vars$name==vari],dim=list(d_lon,d_lat,d_time),missval=-999,
    longname=vars$longname[vars$name==vari],compress=comp)
  v_var_q2.5=ncvar_def(paste(vari,"_q2.5",sep=""),units=vars$unit[vars$name==vari],dim=list(d_lon,d_lat,d_time),missval=-999,
    longname=vars$longname[vars$name==vari],compress=comp)
  v_var_q97.5=ncvar_def(paste(vari,"_q97.5",sep=""),units=vars$unit[vars$name==vari],dim=list(d_lon,d_lat,d_time),missval=-999,
    longname=vars$longname[vars$name==vari],compress=comp)
  v_var_sd=ncvar_def(paste(vari,"_sd",sep=""),units=vars$unit[vars$name==vari],dim=list(d_lon,d_lat,d_time),missval=-999,
    longname="Standard Deviation of prediction",compress=comp)
## Add q2.5 and q97.5 here?
  v_var_valid=ncvar_def(paste(vari,"_valid",sep=""),units="values",dim=list(d_valid,d_valid2,d_time),missval=-999,
    longname="Stations held out for validation",compress=comp)
  v_krige=ncvar_def(paste(vari,"_params",sep=""),units="values",dim=list(d_iter,d_params,d_time),missval=-999,
    longname="posterior kriged parameters",compress=comp)

## set up nc file
   ncfile=paste(odir,vari,"_",idate,".nc",sep="")
   if(file.exists(ncfile)) file.remove(ncfile)
   if(keepall) nc_create(ncfile,vars=list(v_var,v_var_mean,v_var_sd,v_krige,v_var_valid),verbose=F)   #save every iteration
   if(!keepall) nc_create(ncfile,vars=list(v_var_mean,v_var_sd,v_var_q2.5,v_var_q97.5,v_krige,v_var_valid),verbose=F)   #save every iteration

   nc=nc_open(ncfile,write=T)
   print("NetCDF file created, adding data")

  ncatt_put(nc,paste(vari,"_params",sep=""), "names", paste(names(kc$posterior$sample),collapse=","))
  ncatt_put(nc,paste(vari,"_valid",sep=""),"columns",paste(colnames(attr(kc,"meta")$valid),collapse=","),prec="char")

  ## Add data
   if(keepall) ncvar_put(nc,vari,vals=var_iter,start=c(1,1,1,1),c(-1,-1,-1,1),verb=F)
   ncvar_put(nc,paste(vari,"_mean",sep=""),vals=var_mean[,ncol(var_mean):1],start=c(1,1,1),c(-1,-1,1),verb=F)
   ncvar_put(nc,paste(vari,"_sd",sep=""),vals=var_sd[,ncol(var_sd):1],start=c(1,1,1),c(-1,-1,1),verb=F)
   ncvar_put(nc,paste(vari,"_q2.5",sep=""),vals=var_q2.5[,ncol(var_q2.5):1],start=c(1,1,1),c(-1,-1,1),verb=F)
   ncvar_put(nc,paste(vari,"_q97.5",sep=""),vals=var_q97.5[,ncol(var_q97.5):1],start=c(1,1,1),c(-1,-1,1),verb=F)
   ncvar_put(nc,paste(vari,"_valid",sep=""),vals=t(as.matrix(attr(kc,"meta")$valid)),start=c(1,1,1),c(-1,-1,1),verb=F)
   #ncvar_put(nc,paste(vari,"_params",sep=""),vals=apply(kc$posterior$sample,2,function(x) c(mean=mean(x),var=var(x))),start=c(1,1,1),c(-1,-1,1),verb=F)
   ncvar_put(nc,paste(vari,"_params",sep=""),vals=t(as.matrix(kc$posterior$sample)),start=c(1,1,1),c(-1,-1,1),verb=F)
  
print("Data added, updating attributes")
  ################################
  ## Attributes
    ncatt_put(nc,"iter", "description","Posterior samples from predictive distribution",prec="character")
  ## Global Attributes
    ncatt_put(nc,varid=0, "title","Interpolated observed weather data for the Cape Floristic Region",prec="character")
    ncatt_put(nc,varid=0, "institution","Silander Lab, EEB, UCONN, Storrs, CT",prec="character")
    ncatt_put(nc,varid=0, "source","Interpolated Daily Station Data",prec="character")
    ncatt_put(nc,varid=0, "comment","Kriged by Adam Wilson (adam.wilson@uconn.edu) using the krige.bayes in library GeoR in R",prec="character")
    ncatt_put(nc,varid=0, "nstation",as.integer(attr(kc,"meta")$nStations),prec="short")
    ncatt_put(nc,varid=0, "RunningTime",paste(round(attr(kc,"meta")$duration,2),attr(attr(kc,"meta")$duration,"units")),prec="character")
 
## Close the file
    nc_sync(nc)
    nc_close(nc)

#system(paste("ncpdq -O -4 -L 9 ",ncfile,ncfile))


### Quit R
print("Finished, quitting R")
q("no")
