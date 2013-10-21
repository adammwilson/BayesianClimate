### Import station metadata from text file 
importmeta<-function(f,verbose=T){
	if(verbose) print(paste("Starting file:",f))
	data.frame(
			"station"=scan(f,what="char",skip=54,nlines=1,quiet=T)[3],
			"lon"=scan(f,what="char",skip=57,nlines=1,quiet=T)[3],
			"lat"=scan(f,what="char",skip=56,nlines=1,quiet=T)[3],
			"elevation"=scan(f,what="char",skip=58,nlines=1,quiet=T)[3],
			"var"=tolower(scan(f,what="char",skip=52,nlines=1,quiet=T)[3]),
			"datestart"=as.Date(scan(f,what="char",skip=59,nlines=1,quiet=T)[3],"%Y%m%d"),
			"datestop"=as.Date(scan(f,what="char",skip=60,nlines=1,quiet=T)[3],"%Y%m%d"),
			"stationname"=scan(f,what="char",skip=55,nlines=1,quiet=T)[3],
			"path"=f,stringsAsFactors=F)
}
 
## put data in the file
importdata<-function(f){
	d=read.csv(f,h=T,skip=62)[,]
	#d$date=as.Date(as.character(d$date),"%Y%m%d")
	d$VAR[!grepl("_",d$QC)]=NA  #drop all measurements with problematic data
	return(as.numeric(d$VAR))
}




### Function to define the nc files
makenc<-function(o) {
  ## subset meta
  filename=paste("CFR_Africa.2.2_",o,".nc",sep="")
  subset=meta$output==o
  stations=unique(meta$station[subset])
  models=unique(meta$model[subset])
  ## Extract Dates
  refdate=as.Date("1960-01-01")
  dtemp=importdata(meta$path[subset][1])
  dates=seq(attr(dtemp,"datestart"),attr(dtemp,"datestop"),1)
  datecode=as.integer(dates-refdate)  #get date with respect to reference
  ## Get string lengths
  nchar=max(c(nchar(models),nchar(stations)))
  ## Set dimentions
  dimnchar=ncdim_def("nchar",create_dimvar=F,units="",vals=1:nchar, unlim=FALSE)
  dm=ncdim_def("model",create_dimvar=F,units="",vals=1:length(models), unlim=FALSE)
  ds=ncdim_def("station",create_dimvar=F,units="",vals=1:length(stations), unlim=FALSE)
  dt=ncdim_def("time",units=paste("days since",refdate),vals=datecode,calendar="standard",unlim=FALSE)
  ## Set dimension variables
  vm=ncvar_def("modelname",units="name",dim=list(dimnchar,dm),missval=NA,longname="Model Name",prec="char") #following example in manual - ugh.
  vs=ncvar_def("stationcode",units="name",dim=list(dimnchar,ds),missval=NA,longname="Station Code",prec="char")
  vlat=ncvar_def("latitude",units="degrees_south",dim=ds,missval=NA,longname="Station Latitude",prec="double")
  vlon=ncvar_def("longitude",units="degrees_east",dim=ds,missval=NA,longname="Station Longitude",prec="double")
   ## Set weather variables
  vppt=ncvar_def("ppt",units="mm",dim=list(dm,ds,dt),missval=-999,longname="Total Precipitation",prec="double")
  vtmax=ncvar_def("tmax",units="C",dim=list(dm,ds,dt),missval=-999,longname="Maxiumum Temperature",prec="double")
  vtmin=ncvar_def("tmin",units="C",dim=list(dm,ds,dt),missval=-999,longname="Minimum Temperature",prec="double")
  ## Make the netcdf file with all variables above
  nc_create(filename,vars=list(vm,vlat,vlon,vs,vppt,vtmax,vtmin),verbose=F)
  ## open file for writing and update attributes
  nc=nc_open(filename,write=T)
  ## Global Attributes
  ncatt_put(nc,varid=0, "title","Downscaled climate projections for the Cape Floristic Region")
  ncatt_put(nc,varid=0, "institution","Climate Systems Analysis Group, UCT, Cape Town, South Africa")
  ncatt_put(nc,varid=0, "source","Downscaled GCM models")
  ncatt_put(nc,varid=0, "references","http://www.csag.uct.ac.za/")
  ncatt_put(nc,varid=0, "comment","Converted to NetCDF by Adam Wilson (adam.wilson@uconn.edu)")
  ## Add the dimensional data
  ncvar_put(nc,"latitude",st$lat[st$code%in%stations])
  ncvar_put(nc,"longitude",st$lon[st$code%in%stations])
  ncvar_put(nc,"modelname",models)
  ncvar_put(nc,"stationcode",st$code[st$code%in%stations])
  ## Close the file
  nc_sync(nc)
  nc_close(nc)
  print(paste("Finished output:",o))
}

putdatanc<-function(o,vars=c("ppt","tmax","tmin")) {
  ## Put data into the netcdf file
  ## subset meta
  filename=paste("CFR_Africa.2.2_",o,".nc",sep="")
  ## Open the netcdf file
  nc=nc_open(filename,write=T)
  subset=meta$output==o
  models=ncvar_get(nc,varid="modelname")
  stations=ncvar_get(nc,varid="stationcode")
  tmeta=meta[subset,]
  ## add the data
  lapply(1:nrow(tmeta),function(i) {  #loop through files and add them to the correct place in model-station-time array
    dtemp=read.table(paste("data/",tmeta$path[i],sep=""),h=T)                            # read in data 
    tstation=which(stations==tmeta$station[i])                                           # get station id
    tmodel=which(models==tmeta$model[i])                                                 # get model id
    tvars=colnames(dtemp)[!grepl("date",colnames(dtemp))]                                # get variabiles in file
    lapply(tvars[tvars%in%vars],function(v)                                              # loop though variables and put them in the nc file
           ncvar_put(nc,v,vals=dtemp[,v],start=c(tmodel,tstation,1),c(1,1,nrow(dtemp))))
    if(i%in%pretty(1:nrow(tmeta),10)) print(paste("Finished ",i," out of ",nrow(tmeta)," files from output ",o))
  })
  nc_close(nc)
  print(paste("Finished output:",o))
}


## check if data was put in the netcdf file correctly   
checkdata<-function(o,vars=c("ppt_tot","ppt_med","ppt_days2","ppt_dryspell","tmax_avg","tmin_avg","tmax_max","pdsi")){
  filename=paste("CFR_Africa.2.2_",o,".nc",sep="")
  ## Open the netcdf file
  nc=nc_open(filename,write=F)
  models=ncvar_get(nc,varid="modelname")
  stations=ncvar_get(nc,varid="stationcode")
  tmeta=meta[meta$output==o,][sort(unlist(tapply(1:nrow(meta[meta$output==o,]),list(meta$output[meta$output==o],meta$model[meta$output==o]),sample,10))),]
  #tmeta=meta[meta$output==o,]
  same=lapply(1:nrow(tmeta),function(i) {  #loop through files and add them to the correct place in model-station-time array
    dtemp=read.table(paste("data/",tmeta$path[i],sep=""),h=T)                            # read in data 
    tstation=which(stations==tmeta$station[i])                                           # get station id
    tmodel=which(models==tmeta$model[i])                                                 # get model id
    tvars=colnames(dtemp)[!grepl("date",colnames(dtemp))]                                # get variabiles in file
    same=do.call(rbind,lapply(tvars[tvars%in%vars],function(v)                           # loop though variables and see if they're the same
      identical(dtemp[,v],as.numeric(ncvar_get(nc,v,start=c(tmodel,tstation,1),c(1,1,nrow(dtemp)))))))
    if(i%in%pretty(1:nrow(tmeta),10)) print(paste("Finished ",i," out of ",nrow(tmeta)," files from output ",o))
    return(same)
  })
  do.call(rbind,same)
}
