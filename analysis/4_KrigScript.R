### script to run on each node that loads one day's data,
### kriges it and converts it to a netcdf file

##First read in the arguments listed at the command line
args=(commandArgs(TRUE))
##args is now a list of character vectors
## First check to see if arguments are passed.
## if no args were given, print a warning and stop
if(length(args)==0) {stop("No parameters supplied, you must pass parameters")}

## Then cycle through each element of the list and evaluate the expressions.
eval(parse(text=args))
## now there is an i that corresponds to the row in the notdone object that will be processed.

odir="output/daybyday/"
 
## load libraries
## running on farm?
farm=!grepl("protea",system("uname -a",intern=T))
if(!farm){
	library(sp);library(geoR);library(ncdf4);library(reshape);library(raster)
	setwd("/media/Data/Work/Regional/CFR/Weather/ClimateInterpolation")}
## running on farm
if(farm) {
	rlib="/group/latimergrp/R/x86_64-redhat-linux-gnu-library/2.13"
	require(plyr,lib=rlib)
	require(sp,lib.loc=rlib)
        require(geoR,lib.loc=rlib)
	require(ncdf4,lib.loc=rlib)
	require(reshape,lib.loc=rlib)
	require(raster,lib.loc=rlib)
}


## load the prediction grid, start dates, and cfr border
load("Interp.Rdata")

## load the list of 'not done' dates
load("notdone.Rdata")

## identify date in progress
idate=as.Date(notdone$date[i])
vari=as.character(notdone$var[i])

## Start timer to compute how long it takes
starttime=Sys.time()

### Print some summary information (to log)
print(paste("Process running on ",system("uname -a",intern=T)))
print(paste("Process started at ",starttime))
print(paste("Processing variable:",vari,", array element:", i,"  Date:",idate,sep=""))
print(paste("Temporary Directory",tempdir()))


### open station data 
nc=nc_open("StationData.nc",write=F)
### import metadata
time=ncvar_get(nc,"time")
station_origin=as.Date(strsplit(ncatt_get(nc,"time")$units," ")[[1]][3])
itime=which(station_origin+time==idate)  #time index of idate

### hold out data for validation
vn=3  # of stations to hold out

### import Observed data for this day
data=data.frame(
  station=as.vector(ncvar_get(nc,"station")),
  lat=as.vector(ncvar_get(nc,"lat")),
  lon=as.vector(ncvar_get(nc,"lon")),
  value=as.numeric(ncvar_get(nc,vari,start=c(itime,1),count=c(1,-1))))
  data=data[!is.na(data$value)&data$value!=-999,]

  ## print some progress
  print(paste("Now working on ",vari))
  print(head(data))
  print(summary(data))

  ## correct column names on table if incorrect
  latc=grep("lat",colnames(data))
  lonc=grep("lon",colnames(data))
  if(data[1,"lat"]>0) {colnames(data)[latc]="lon";colnames(data)[lonc]="lat"}
  
	## Drop stations that are Xkm from CFR
	data$dcfr=as.vector(extract(raster("CFR.nc"),data[,c("lon","lat")]))
	data=data[data$dcfr<=75&!is.na(data$value)&!is.na(data$dcfr),]

	### drop duplicate coordinates
        dupes=dup.coords(data[,c('lon','lat')])
        if(!is.null(dupes)) data=data[!rownames(data)%in%dupes[2:nrow(dupes),],]

  ###########################
  ## apply climate correction
  clim=brick("Climate.nc",varname=vari);clim
  clim=subset(clim,subset=as.numeric(format(idate,"%m")))
  ## drop any stations outside the climate region
  data=data[data$lat>=bbox(clim)[2,1]&data$lat<=bbox(clim)[2,2]&data$lon>=bbox(clim)[1,1]&data$lon<=bbox(clim)[1,2],]
  data$clim=as.numeric(extract(clim,data[,c("lon","lat")],method="bilinear"))
	
  data=data[!is.na(data$clim),]
  if(vari%in%c("tmax","tmin"))   data$climanomaly=data$clim-data$value

    climateoffset=1
    climanomalyoffset=1  #set to 1 for log transform
if(vari=="ppt")  {
  ## Scale for precipitation (add 1 to climate to avoid dividing by 0, will be subtracted off later)
  ## Add 1 to final value to avoid transformation troubles (geoR will only transform if >0
    data$clim=ifelse(data$clim<0,0,data$clim)
    data$climanomaly=(data$value/(data$clim+climateoffset))+climanomalyoffset
  }

  print(summary(data$climanomaly))

  ## extract data for validation
  valid=data[data$station%in%sample(data$station,vn),]
  ## delete the validation data from the fitting data
  data=data[!data$station%in%valid$station,] 

  ## if all values are 0 (happens in precip), it fails at prediction, so add 1mm to one station
  if(all(data$climanomaly==1)) data$climanomaly[sample(1:nrow(data),1)]=2
  
  ## Set up data for geoR
  border=cfr@polygons[[1]]@Polygons[[1]]@coords
  dat=as.geodata(data,coords.col=c("lon","lat"),data.col="climanomaly",borders=border)  #fitting data
  
  ## identify prediction locations
   pdat=as.geodata(pgrid)
#  pdat=as.geodata(rbind(data,valid),coords.col=c("lon","lat"),data.col="climanomaly",borders=border)  #fitting data
#  pdat=as.geodata(cbind(coordinates(raster("CFR_0.05.nc")),
#    data=1:nrow(coordinates(raster("CFR_0.05.nc")))),coords.col=c("x","y"),data.col="data",borders=border)
  
  ## Run the model
  ##  Set up Model
  mod=model.control(trend.d=trend.spatial("1st",dat),trend.l=trend.spatial("1st",pdat),lambda=ifelse(vari=="ppt",0,1))
  out=output.control(n.posterior=1000,n.predictive=1000,signal=T)

##################
  ##  Set up Priors
if(vari%in%c("tmax","tmin"))  #For temperature
  pri=prior.control(beta.prior="normal",beta=rep(0,3),beta.var.std=diag(10,3),
    sigmasq.prior="reciprocal",
    phi.prior="exponential",phi=4,phi.discrete=seq(0,20,length.out=101), # plot(pri$priors.info$phi$probs~names(pri$priors.info$phi$probs),ylim=c(0,.05))
    tausq.rel.prior="reciprocal",tausq.rel.discrete=seq(0,1,length.out=51)) #tau goes 0-1

if(vari=="ppt")  #for Precipitation
  pri=prior.control(beta.prior="normal",beta=rep(0,3),beta.var.std=diag(10,3),
    sigmasq.prior="reciprocal",
    phi.prior="exponential",phi=0.5,phi.discrete=seq(0,3,length.out=81), #2
    tausq.rel.prior="reciprocal",tausq.rel.discrete=seq(0,1,length.out=21)) #tau goes 0-1

  ## Run final version and make the predictions
  gc(); kc=krige.bayes(dat,loc=pdat[[1]],model=mod,prior=pri,output=out); gc()

  ## stop timer
  stoptime=Sys.time()
  diftime=stoptime-starttime

  ## update metadata by adding non-standard attributes to kc (krige.bayes object)
  attr(kc,"meta")=list(variable=vari,date=idate,nStations=nrow(data),
        starttime=starttime,stoptime=stoptime,duration=diftime,
        valid=valid,climateoffset=climateoffset,climanomalyoffset=climanomalyoffset)

save(kc,file=paste(odir,vari,"_",idate,".Rdata",sep=""),compress="xz",compression_level=9)


  ## Some plots for validation
saveplot=F
if(saveplot){
  pdf(paste("output/plots/summary",i,vari,idate,".pdf",sep="_"))
  par(mfrow=c(2,5));hist(kc)
  plot(kc)
  plot(variog(dat, max.dist = 6,lambda=ifelse(vari=="ppt",0,1)))#, ylim=c(0, .75)) 
  lines(kc, summ = mean, lwd=1, lty=2);lines(kc, summary = "median", post = "parameters")
  my.summary <- function(x){quantile(x, prob = c(0.05, 0.5, 0.95))}
  lines(kc, summ = my.summary, ty="l", lty=c(2,1,2), col=1)
  #par(mfrow=c(1,1));
  ## extract validation data
  obsindex=overlay(SpatialPixelsDataFrame(pdat[[1]],data=data.frame(mean=kc$predictive$mean,row.names=rownames(pdat[[1]]))),
               SpatialPointsDataFrame(rbind(dat[[1]],valid[,c('lon','lat')]),
                                      data=data.frame(id=c(rownames(dat[[1]]),rownames(valid)),row.names=c(rownames(dat[[1]]),rownames(valid)))))
  pred=cbind.data.frame(mean=kc$predictive$mean[obsindex],sd=sqrt(kc$predictive$variance[obsindex]),obs=c(dat$data,valid$climanomaly))
  
  plot(mean~obs,data=pred,col=c(rep("black",length(dat$data)),rep("green",vn)),pch=16) #,ylim=c(1,3),xlim=c(1,3)
  abline(0,1,col="red")
  arrows(pred$obs,pred$mean-pred$sd,pred$obs,pred$mean+pred$sd,col="grey",length=0)
  dev.off()
}
  
	

### Quit R
print("Finished, quitting R")
q("no")

