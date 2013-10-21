###########################################################
system("ssh -X protea")
R

###########################################################
#{{{ Libraries, working directories, etc.
library(reshape);library(lattice);library(rgdal);library(latticeExtra);library(multicore);library(ncdf4)
library(sp);library(geoR);library(raster);library(abind)

setwd("/media/Data/Work/Regional/CFR/Weather/ClimateInterpolation")

### load some functions
source("/media/Data/Work/Stat/Scripts/Functions.r")

### Set version number for file locations, keep this up to date with git repository to link code and output!
version="v1.3b"

### Mount ramdisk?
if(!file.exists("/tmp/ramdisk")) system("sudo mkdir /tmp/ramdisk; sudo chmod 777 /tmp/ramdisk")
system("sudo mount -t tmpfs -o size=20G tmpfs /tmp/ramdisk/")
scratch="/tmp/ramdisk/"
#system("sudo umount  /tmp/ramdisk/")

#}}}

###########################################################
#{{{ Download and process data to netcdf
## Copy data back here

if(!file.exists(paste("output/",version,"/nc",sep=""))) dir.create(paste("output/",version,"/nc",sep=""),recursive=T)
system(paste("rsync -uave ssh zzwilson@davisfarm:/group/latimergrp/adam/cfr/interp/output/nc_mean/ output/",version,"/nc",sep=""),wait=F)

## merge them by time
outputfile=paste("output/",version,"/CFR_Weather_Krige.nc",sep="")

## variables to process
## these variables are used to construct the netcdf files and everything below.  
vars=c("tmax","tmin","ppt")

## Concatenate the day-by-day files to a single timeseries for each variable
mclapply(vars,function(v)
         system(paste("ncrcat -O -4 -L 9 -x -v ",v," output/",version,"/nc/",v,"* output/",version,"/",v,"_krig.nc",sep=""))) 

## get time_origin
time_origin=as.Date(strsplit(ncatt_get(nc_open(paste("output/",version,"/ppt_krig.nc",sep="")),"time")$units," ")[[1]][3])


## Identical dates must be extracted from files before merging if not all days are complete.
## Get list of common dates to be extracted
ds=mclapply(vars,function(v)  {
  system(paste("cdo showdate output/",version,"/",v,"_krig.nc > output/",version,"/",v,"_dates.txt",sep=""))
  dates=as.Date(scan(paste("output/",version,"/",v,"_dates.txt",sep=""),what="char"))
  file.remove(paste("output/",version,"/",v,"_dates.txt",sep=""))
  return(data.frame(date=dates))
})
datekeep=merge(merge(ds[[1]],ds[[2]]),ds[[3]])$date

d2=lapply(vars,function(v){
  vd=time_origin+as.vector(ncvar_get(nc_open(paste("output/",version,"/",v,"_krig.nc",sep="")),"time"))
  which(vd%in%datekeep)
});names(d2)=vars


## merge all the variables into a single file taking care to line up timesteps if not all files have been processed
#system(paste("cdo -O -z zip merge ",
#             ifelse("ppt"%in%vars,paste(" -seltimestep,",paste(d2[["ppt"]],collapse=",")," output/",version,"/ppt_krig.nc ",sep=""),""),
#             ifelse("tmin"%in%vars,paste(" -seltimestep,",paste(d2[["tmin"]],collapse=",")," output/",version,"/tmin_krig.nc ",sep=""),""),
#             ifelse("tmax"%in%vars,paste(" -seltimestep,",paste(d2[["tmax"]],collapse=",")," output/",version,"/tmax_krig.nc ",sep=""),""),
#             "output/",version,"/weather_coarse_anom.nc",sep=""))
system(paste("cdo -O -z zip merge ",
             ifelse("ppt"%in%vars,paste(" output/",version,"/ppt_krig.nc ",sep=""),""),
             ifelse("tmin"%in%vars,paste(" output/",version,"/tmin_krig.nc ",sep=""),""),
             ifelse("tmax"%in%vars,paste(" output/",version,"/tmax_krig.nc ",sep=""),""),
             "output/",version,"/weather_coarse_anom.nc",sep=""))

## remove the temporary files
file.remove(c(paste("output/",version,"/ppt_krig.nc",sep=""),
              paste("output/",version,"/tmax_krig.nc",sep=""),
              paste("output/",version,"/tmin_krig.nc",sep="")))


## interpolate to Climate Grid
## get schulze grid
system("cdo griddes Climate.nc > climategrid.txt")
#system("cdo griddes -selvar,tmax output/daybyday/tmax_1990-01-01.nc") 


## remap it to climate grid and mask out pixels > 5km from CFR
system(paste("cdo -z zip ifthen -invertlat -ltc,5 -selvar,CFR_dist CFR.nc -remapbic,climategrid.txt ",
             "-selvar,",paste(rep(vars,each=2),rep(c("_mean","_sd"),times=length(vars)),sep="",collapse=",")," output/",
             version,"/weather_coarse_anom.nc output/",version,"/weather_fine_anom.nc",sep="")) 
	
## transform back to real units using the climate data
  system(paste("cdo -O -z zip merge ",
               ## tmax
                ifelse("tmax"%in%vars,paste("-ymonadd -mulc,-1 -selvar,tmax_mean output/",version,"/weather_fine_anom.nc -selvar,tmax Climate.nc ",
               ## convert SD							 
               "-selvar,tmax_sd output/",version,"/weather_fine_anom.nc ",sep=""),""),
                ifelse("tmin"%in%vars,paste("-ymonadd -mulc,-1 -selvar,tmin_mean output/",version,"/weather_fine_anom.nc -selvar,tmin Climate.nc ",
               ## convert SD
               "-selvar,tmin_sd output/",version,"/weather_fine_anom.nc ",sep=""),""),
               ## take maximum of 0 or ppt to drop negative values using gtc,9999.
               ## Then multiply by mean ppt+1 (1 was added before kriging to remove 0s from climate.)
                ifelse("ppt"%in%vars,paste(" -ymonmul -max -gtc,99999 -selvar,ppt_mean output/",version,
               "/weather_fine_anom.nc -selvar,ppt_mean output/",version,
               "/weather_fine_anom.nc -addc,1 -selvar,ppt Climate.nc ",
               ## convert SD
               "-ymonmul -selvar,ppt_sd output/",version,
               "/weather_fine_anom.nc -selvar,ppt Climate.nc ",sep=""),""),
               " output/",version,"/weather_fine_unpacked.nc",
               sep=""))

## TODO: Why does missing data region around coast change?!!

## Pack it.
system(paste("ncpdq -O -4 -L 9 output/",version,"/weather_fine_unpacked.nc output/",version,"/weather_fine.nc",sep=""),wait=F)
#}}}

## compare nc_short and nc_byte precision
nc1=nc_open("output/v1.3b/weather_fine.nc")
nc2=nc_open("output/v1.3b/weather_fine_byte.nc")
t1=as.vector(ncvar_get(nc1,"tmax_mean",c(1,1,1),c(-1,-1,5)))
t2=as.vector(ncvar_get(nc2,"tmax_mean",c(1,1,1),c(-1,-1,5)))
plot(t1,t2)

###########################################################
#{{{ extract validation data to compare with remapped data
inc=nc_open(paste("output/",version,"/weather_coarse_anom.nc",sep=""))
time=as.vector(ncvar_get(inc,"time"))
time_origin=as.Date(strsplit(ncatt_get(inc,"time")$units," ")[[1]][3])
dates=time_origin+time  #time index of idate

## import station data to use for validation
stnc=nc_open("StationData.nc",write=F)
### import metadata
time=ncvar_get(stnc,"time")
lat=as.vector(ncvar_get(stnc,"lat"))
lon=as.vector(ncvar_get(stnc,"lon"))

## variables to process
## these variables are used to construct the netcdf files and everything below.  
vars=c("tmax","tmin","ppt")


reshape_valid<-function(x,var){
  tx=do.call(rbind.data.frame,mclapply(1:dim(x)[3],function(p) cbind.data.frame(timeid=p,matrix(as.vector(x[,,p]),ncol=7,byrow=T))))
  colnames(tx)=c("timeid","station","lat","lon","obs","dcfr","clim","climanom")
  tx$var=var
  return(tx)}

valid=rbind.data.frame(
  reshape_valid(ncvar_get(inc,"tmax_valid"),"tmax"),
  reshape_valid(ncvar_get(inc,"tmin_valid"),"tmin"),
  reshape_valid(ncvar_get(inc,"ppt_valid"),"ppt"))


## table of validation sources
valid_sources=data.frame(
		vid=paste("V",1:6,sep=""),
		scale=c("coarseanom","coarseanom","fineanom","fineanom","finepred","finepred"),
		type=rep(c("mean","sd"),3),
		name=c("krig_anom_mean","krige_anom_sd","fine_anom_mean","fine_anom_sd","fine_pred_mean","fine_pred_sd"),
		file=c(
                  paste("output/",version,"/weather_coarse_anom.nc",sep=""),
                  paste("output/",version,"/weather_coarse_anom.nc",sep=""),
                  paste("output/",version,"/weather_fine_anom.nc",sep=""),
                  paste("output/",version,"/weather_fine_anom.nc",sep=""),
                  paste("output/",version,"/weather_fine_unpacked.nc",sep=""),
                  paste("output/",version,"/weather_fine_unpacked.nc",sep="")),stringsAsFactors=F)

pextract<-function(brick,points,pointid="station",...){
  ## brick is the raster brick you want to extract from
  ## points are the points you want to extract data for
  ## clusters is the number of clusters you want to use
  require(cluster)
  ## Define clusters
  clusters=max(24,round(nrow(points)/5))
  ## group stations into cluster to speed up extraction
  points$cluster=as.factor(clara(coordinates(points),clusters)$clustering)
  t1=do.call(rbind,mclapply(1:clusters,function(ki){
    print(paste("Processing Cluster:",ki," out of ",clusters," clusters"))
    tpoints=points[points$cluster==ki,]
    tpoints@data[c("lon","lat")]=coordinates(tpoints)
    ext=extent(min(tpoints$lon*.999),max(tpoints$lon*1.001),
      ifelse(min(tpoints$lat)<0,min(tpoints$lat*1.001),min(tpoints$lat*.999)),
      ifelse(max(tpoints$lat)<0,max(tpoints$lat*0.999),max(tpoints$lat*1.001)))
    ## lapply through validatation types and update the valid table with the predictions
    #print(paste("Cropping Cluster:",ki))
    cbrick=crop(brick,ext)
    #print(paste("Extracting Cluster:",ki," out of ",clusters," clusters"))
    tt1=extract(cbrick,tpoints)
    if(is.null(dim(tt1))) tt1=t(as.matrix(tt1))
    dimnames(tt1)[1]=list(station=tpoints@data[,pointid])
    return(tt1)
  }))
	gc()
	return(t1)}

getvalid<-function(v){
  vt=valid[valid$var==v,]
  st_u=unique(vt[,c("station","lat","lon")])
  coordinates(st_u)=c("lon","lat")
  for(r in 5:6){#1:nrow(valid_sources)){
    print(paste("####################   Now processing ",valid_sources$name[r],"for ",v))
    t1=pextract(brick=brick(valid_sources$file[r],varname=paste(v,valid_sources$type[r],sep="_")),points=st_u)
    tmatch=match(vt$station,rownames(t1))
    #timeids=match(vt$date,ncvar_get(inc,"time")+station_origin)
    vt[,valid_sources$name[r]]=as.numeric(do.call(c,lapply(1:nrow(vt),function(i) t1[tmatch[i],vt$timeid[i]])))
  }
  return(vt)
}

valid=do.call(rbind.data.frame,lapply(vars,getvalid))

                                       #
## drop NAs
#valid=valid[!is.na(valid$krig_anom),]
valid$date=(ncvar_get(inc,"time")+time_origin)[valid$timeid]

#TODO: figure out why some values are >200?!?!  Missing underlying data?  Problem with packing?
table(valid$fine_pred_mean[valid$var=="ppt"]<200)
table(valid$fine_pred_mean[valid$var=="ppt"]<1)

###########################
### Extract krig parameters
getparms<-function(vari="tmax",inc=nc_open(paste("output/",version,"/weather_coarse_anom.nc",sep=""))) {
	kparms=ncvar_get(inc,paste(vari,"_params",sep=""))
        kparms2=abind(mclapply(1:dim(kparms)[3],function(p) matrix(as.vector(kparms[,,p]),ncol=dim(kparms)[2],byrow=T)),along=3)
	dimnames(kparms2)=list(
                  iteration=1:dim(kparms2)[1],
                  parameter=c("beta0","beta1","beta2","sigmasq","phi","tausq.rel"),
                  date=as.character(dates[1:dim(kparms2)[3]]))
	return(kparms2)}

kparms=lapply(vars,getparms); names(kparms)=vars

## Summarize them
skparms=mclapply(kparms,function(x) cast(melt(apply(x,c(2,3),quantile,c(0,0.025,.5,.975,1))),parameter+date~Var.1))
names(skparms)=vars
skparms=cast(melt.list(skparms),L1+parameter+date~Var.1)
colnames(skparms)[1]="var"
skparms$date=as.Date(skparms$date)

#}}}

###########################################################
#{{{ Extract stations used in interpolations
stnc=nc_open("StationData.nc",write=F)
st_time=ncvar_get(stnc,"time")
station_origin=as.Date(strsplit(ncatt_get(stnc,"time")$units," ")[[1]][3])
st_date=which(station_origin+st_time>=as.Date("1990-01-01")&station_origin+st_time<=as.Date("2010-12-31"))  #time index of idate

st=data.frame(station=as.vector(ncvar_get(stnc,"station")),lat=as.vector(ncvar_get(stnc,"lat")),lon=as.vector(ncvar_get(stnc,"lon")))
st$temp=apply(t(ncvar_get(stnc,"tmax",start=c(min(st_date),1),count=c(length(st_date),-1))),1,function(x) any(!is.na(x)))
st$ppt=apply(t(ncvar_get(stnc,"ppt",start=c(min(st_date),1),count=c(length(st_date),-1))),1,function(x) any(!is.na(x)))
## limit to within 75km of CFR as was done in kriging
st$dcfr=as.vector(extract(raster("CFR.nc"),st[,c("lon","lat")]))
st=st[st$dcfr<=75&!is.na(st$dcfr),]
coordinates(st)=c("lon","lat")
st$type=ifelse(st$temp&st$ppt,"Temperature and Precipitation",ifelse(st$temp&!st$ppt,"Temperature","Precipitation"))
#}}}

###########################################################
#{{{ Save it
writeOGR(st,dsn="StationLocations.shp",layer="StationLocations",driver="ESRI Shapefile")
## Save validation data
save(valid,kparms,skparms,st,file=paste("output/",version,"/valid_0.25deg.Rdata",sep=""))
#}}}

###############################################################################################
#{{{ download climate metrics
if(!file.exists(paste("output/",version,"/metrics_fromfarm",sep=""))) dir.create(paste("output/",version,"/metrics_fromfarm",sep=""),recursive=T)
system(paste("rsync -uave ssh zzwilson@davisfarm:/group/latimergrp/adam/cfr/interp/output/metrics/ output/",version,"/metrics_fromfarm",sep=""),wait=F)
#}}}

###############################################################################################
#{{{ Process climate metrics

## get list of metrics to process
f=data.frame(file=list.files(paste("output/",version,"/metrics_fromfarm",sep="")),stringsAsFactors=F)
f[,c("metric","year","iter")]=do.call(rbind.data.frame,strsplit(sub(".nc","",f$file),"_"))
f$iter=as.numeric(as.character(f$iter))

cast(f,metric~year,fun.aggregate=length)

fu=unique(cbind.data.frame(metric=as.character(f$metric),year=as.character(f$year)));rownames(fu)=1:nrow(fu)
fu$done=fu$metric%in%sub("[.]nc*","",list.files(paste("output/",version,"/metrics_summary",sep="")))
unique(cbind(as.character(fu$metric),fu$done))

## find iterations not yet done
is=(0:999)[!(0:999%in%as.numeric(do.call(rbind,strsplit(sub(".nc","",list.files(scratch)),"_"))[,2]))]

for(m in unique(fu$metric[!fu$done])){
  mclapply(is,function(i) {
    if(i%in%pretty(0:999,5)) print(paste("Starting iteration",i,"of metric",m))
    system(paste("ncrcat -O  output/",version,"/metrics_fromfarm/",m,"_*_",i,".nc ",scratch,m,"_",i,".nc",sep=""))
  })
  system(paste("ncecat -u iteration -O ",scratch,m,"_* output/",version,"/metrics_summary/",m,".nc",sep=""))
  file.remove(list.files(scratch,pattern="[.]tmp",full=T))
  file.remove(list.files(scratch,pattern=m,full=T))
  print(paste("Finished metric",m))
}


## confirm they are all finished and update fu$done
file.remove(list.files(scratch,pattern="[.]tmp",full=T))  #remove temporary files
file.remove(list.files(paste("output/",version,"/metrics_summary",sep=""),full=T,pattern="[.]tmp"))  #remove temporary files
fu$done=fu$metric%in%sub("[.]nc*","",list.files(paste("output/",version,"/metrics_summary",sep="")))
table(fu$metric,fu$done)


#}}}

### Generate single summary file for all metrics from summaries
f=file=list.files(paste("output/",version,"/metrics_summary",sep=""),pattern="nc")
f=f[!grepl("mean",f)]

mclapply(f[-1],function(fi) { gc()
                           system(paste("ncpdq -O -4 -L 9 output/",version,"/metrics_summary/",fi," output/",version,"/metrics_summary/packed_",fi,sep=""),wait=T)
                           system(paste("ncwa -y avg -a iteration output/",version,"/metrics_summary/packed_",fi," output/",version,"/metrics_summary/mean_",fi,sep=""))} )


                           system(paste("ncwa -3 -y avg -a iteration output/",version,"/metrics_summary/",fi," output/",version,"/metrics_summary/mean_",fi,sep=""))

y=1990
m="CDD"

file.remove("output/v1.3b/metrics_fromfarm/CDD_1992_103.nc")

### remove files with structural problems 
#for(i in list.files(paste("output/",version,"/metrics_fromfarm",sep=""),full=T)) {
#   print(paste("Checking",i))
#   if(system(paste("cdo -s griddes ",i),intern=T)=="") {  #returns "" if there was an error
#     file.remove(i)  #remove it
#     print(paste("###########################    Removing ",i))
#   }
# }
   


system(paste("cdo merge ",paste("output/",version,"/metrics_summary/",f,sep="",collapse=",") ,"output/",version,"/metrics_summary/CFRMetricSummary.nc",sep=""))



################################################################################
#{{{ assess the normality of posterior predictions
### if they are normal, we could just use the mean/sd to generate posteriors rather than keeping them.  Could be muuch faster.
if(!file.exists(paste("output/",version,"/nc_all",sep=""))) dir.create(paste("output/",version,"/nc_all",sep=""),recursive=T)
system(paste("rsync -uave ssh --include \"*/\" --include \"*_1990-01*\" --exclude \"*\" zzwilson@davisfarm:/group/latimergrp/adam/cfr/interp/output/nc/ output/",version,"/nc_all",sep=""),wait=F)


i="ppt_1990-01-01.nc"
nc=nc_open(paste("output/",version,"/nc_all/",i,sep=""))
d=ncvar_get(nc,"ppt")

ntest=function(x) {
  if(all(is.na(x))) return(F)
  shapiro.test(x)$p.value<0.01
}

image(apply(d,c(2,3),ntest))

x=d[,1,13]
qqnorm(x);qqline(x)
shapiro.test(x)
#}}}


