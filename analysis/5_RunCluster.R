## Script to organize commands to run the cluster

#{{{ Load data, set directories, etc.
## If running on protea:
system("ssh -X protea")
R
setwd("/media/Data/Work/Regional/CFR/Weather/ClimateInterpolation")

### connect to the farm
system("ssh -X davisfarm")
module load gcc openmpi R szip hdf/5.1.8 netcdf cdo nco
cd /group/latimergrp/adam/cfr/interp

R

rlib="/group/latimergrp/R/x86_64-redhat-linux-gnu-library/2.13/"
library(ncdf4,lib=rlib);library(sp,lib=rlib);library(raster,lib=rlib)

#####################################################################
## make some directories needed for output
if(!file.exists("log")) dir.create("log")
if(!file.exists("output")) dir.create("output")
if(!file.exists("output/daybyday")) dir.create("output/daybyday")
if(!file.exists("output/nc")) dir.create("output/nc")
if(!file.exists("output/plots")) dir.create("output/plots")


#}}}
#{{{ Run the interpolations
## Load the interpolation data

load("Interp.Rdata")

## Load the jobs data
load("jobs.Rdata")
#tjobs=range(which(jobs$var%in%c("tmax","tmin")))
#pjobs=range(which(jobs$var=="ppt"))

#jobs=jobs[jobs$var%in%c("tmax","tmin"),];rownames(jobs)=1:nrow(jobs)
#jobs=jobs[jobs$var%in%c("ppt"),];rownames(jobs)=1:nrow(jobs)#"ppt",]

### update jobs doc to remove completed jobs
done=as.data.frame(do.call(rbind,strsplit(sub(".Rdata","",list.files("output/daybyday",pattern="Rdata")),"_")));colnames(done)=c("var","date")
notdone=jobs[!paste(jobs$var,jobs$date,sep="_")%in%paste(done$var,done$date,sep="_"),]; print(paste(nrow(notdone)," jobs remaining"))
if(nrow(done)==0) notdone=jobs

save(notdone,file="notdone.Rdata")

#
#notdone=expand.grid(date=dates[1:3],var=vars$name)

## get schulze grid
system("cdo griddes Climate.nc > climategrid.txt")

###################
### Generate script to submit with qsub
####^&*(####$ -tc 128 #max threads

cat(paste("
#!/bin/bash
#$ -t 1-",nrow(notdone),"  #These are the elements it loops through
#$ -cwd
#$ -S /bin/bash
#$ -r y
#$ -j y  #merge error and output in log
#$ -o log/
#$ -l mem_free=1.5G  #choose nodes with this much free
#$ -l h_vmem=1.8G  #kill jobs if they go over this
#$ -pe threaded 1  #how many threads to use (set to int(jobsize/1.4GB)+1)
#$ -m abe
#$ -M adam.wilson@uconn.edu

  ulimit -d 1500000 -m 1500000 -v 1500000  #limit memory usage
  export R_LIBS=\"/usr/lib64/R/library:/group/latimergrp/R/x86_64-redhat-linux-gnu-library/2.13/\"
  module load gcc openmpi R szip hdf netcdf
  ## check to check space in /tmp
       mkdir -p /scratch/$USER/$JOB_ID/$SGE_TASK_ID
       export TMPDIR=\"/scratch/$USER/$JOB_ID/$SGE_TASK_ID\"
       export TMP=\"/scratch/$USER/$JOB_ID/$SGE_TASK_ID\"
       export TEMP=\"/scratch/$USER/$JOB_ID/$SGE_TASK_ID\"
       echo \"Temporary Directory set to: $TMPDIR\"
  Rscript --verbose --vanilla 4_KrigScript.R i=\"$SGE_TASK_ID\" #> \"log/log_$SGE_TASK_ID.log\"
  rm -r $TMPDIR
",sep=""),file="InterpScript")

### Check the files
system("cat InterpScript")
#system("cat 4_KrigScript.R")

## Submit it!
system("qsub InterpScript")

## check progress
system("qstat -f") ; print(paste(max(0,length(system("qstat",intern=T))-2)," processes running"))
# system("ssh c0-8.farm.caes.ucdavis.edu")
# system("qalter -p +1024 25964")  #decrease priority of job to run extraction below.
system("cat log/InterpScript.o55934.2")

## check log
system(paste("cat",list.files("log",pattern="InterpScript",full=T)[100]))
#system(paste("cat",list.files("log",pattern="InterpScript",full=T)[13]," | grep \"Temporary Directory\""))

### Delete the logs
#system("rm log/Interp*")

### Delete the krige output - careful!
#system("rm -r output/daybyday/")

### kill jobs
#system("qdel -u zzwilson")

#############  Or run it from here:
library(multicore,lib=rlib)
head(jobs$var[order(jobs$var,decreasing=T)])

mclapply(1,function(i)
         system(paste("Rscript --verbose --vanilla /media/Data/Work/git/PlanetFlux/WeatherInterpolation/4_KrigScript.R i=",i,sep="")))

#{{{ figure out how long it took to run interpolations
system("qrsh")
cd /group/latimergrp/adam/cfr/interp
R
rlib="/group/latimergrp/R/x86_64-redhat-linux-gnu-library/2.13/"
library(sp,geoR,lib=rlib)
library(multicore,lib=rlib)
f=list.files("output/daybyday",full=T)

times=do.call(rbind.data.frame,mclapply(f,function(i){
  print(i)
  load(i)
  kcm=attr(kc,"meta")
  res=data.frame(
    starttime=kcm$starttime,stoptime=kcm$stoptime,
    duration=as.numeric(kcm$duration),units=units(kcm$duration),
    variable=kcm$variable,nstations=kcm$nStations)
  rm(kc)
#  if(i==f[1]) write.table(res,file="output/times.csv",row.names=F,sep=",")
#  if(i!=f[1]) write.table(res,file="output/times.csv",row.names=F,append=T,col.names=F,sep=",")
  return(res)
}))

write.csv(times,file="times.csv")
#}}}

#}}}

#################################################################
#{{{ Process the krige.bayes objects into netCDF files
if(!file.exists("log/convertlog")) dir.create("log/convertlog")
if(!file.exists("output/nc")) dir.create("output/nc")

## get list of files to process
KO=list.files("output/daybyday",full=T)

## subset to only those not yet converted 
KO=KO[!sub(".Rdata","",sub("output/daybyday/","",KO))%in%sub("[.]nc","",list.files("output/nc_mean"))];length(KO)

## save list to be imported by slaves
save(KO,file="KO.Rdata")

###################
### Generate script to submit with qsub
cat(paste("
#!/bin/bash
#$ -t 1-",length(KO),"  #These are the elements it loops through
#$ -cwd
#$ -S /bin/bash
#$ -r y
#$ -j y  #merge error and output in log
#$ -o log/convertlog/
module load gcc openmpi R szip hdf netcdf nco
ulimit -d 1500000 -m 1500000 -v 1500000
  Rscript --verbose --vanilla 4b_ConvertKrigObjects.R ",
          " keepall=F i=\"$SGE_TASK_ID\" 
",sep=""),file="ConvertScript")

### Check the files
#system("cat 4b_ConvertKrigObjects.R")
system("cat ConvertScript")

## Submit it!
system("qsub ConvertScript")

## check progress
system("qstat -f") ; print(paste(max(0,length(system("qstat",intern=T))-2)," processes running"))

## connect to a node
system("ssh c0-8.farm.caes.ucdavis.edu")
exit

## check a log
system(paste("cat",list.files("log/convertlog/",full=T)[6]))

### Delete the krige output - careful!
#system("rm -r output/nc/")

#}}}

#################################################################
#{{{ Process the netcdf files into climate metric objects
#if(!file.exists("log/convertlog")) dir.create("log/convertlog")
#if(!file.exists("output/iter")) dir.create("output/iter")

## get schulze grid
system("cdo griddes Climate.nc > climategrid.txt")

## get list of files to process
niter=999

## Variables to process
vars=c("tmax","tmin","ppt")
version="v1.3b"

## list of metric filenames and their variable
met=rbind.data.frame(
  c("GDD","tmin"),
  c("MATmin","tmin"),
  c("MMTmin","tmin"),
  c("MTmin","tmin"),
  c("SU","tmax"),
  c("CSU","tmax"),
  c("MATmax","tmax"),
  c("MMTmax","tmax"),
  c("MTmax","tmax"),
  c("CDD","ppt"),
  c("rr2","ppt"),
  c("sdii","ppt"),
  c("ECAr20mm","ppt"),
  c("MAP","ppt"),
  c("MMP","ppt"));colnames(met)=c("metric","var")

## full list of iteration-variable-metric-year to be processed
metrics=expand.grid(iter=0:niter,year=1990:2009,metric=met$metric)
metrics$var=met$var[match(metrics$metric,met$metric)]
metrics$name=paste(metrics$metric,metrics$year,metrics$iter,sep="_")

## add list of already finished values
farm=!grepl("protea",system("uname -a",intern=T))
if(!farm) metrics$done=metrics$name%in%sub(".nc","",list.files(paste("output/",version,"/metrics",sep=""),pattern="nc"))
if(farm) metrics$done=metrics$name%in%sub(".nc","",list.files(paste("output/metrics",sep=""),pattern="nc"))

## table of finished metrics
table(metrics$metric,metrics$done)

## create list of variables to process
iter=unique(metrics[!metrics$done,c("iter","year","var")])

## print status
print(paste("Finished ",sum(metrics$done),"out of",nrow(metrics),"(",
            round(sum(metrics$done)/nrow(metrics),2)*100,
            "%), will run another",nrow(iter),"jobs"))

save(iter,file="iternotdone.Rdata")

###################
### Generate script to submit with qsub
cat(paste("
#!/bin/bash
#$ -t 1:",nrow(iter),"  #These are the elements it loops through
#$ -cwd
#$ -S /bin/bash
#$ -r y
#$ -j y  #merge error and output in log
#$ -o log/convertlog/
#$ -pe threaded 16  #how many threads to use (set to int(jobsize/1.4GB)+1)
## export R_LIBS=\"/usr/lib64/R/library:/group/latimergrp/R/x86_64-redhat-linux-gnu-library/2.13/\"
module load gcc openmpi R szip/2.1 hdf/5.1.8 netcdf cdo nco
#ulimit -d 1500000 -m 1500000 -v 1500000  #one core
#ulimit -d 3000000 -m 3000000 -v 3000000 #two cores
ulimit -d 24000000 -m 24000000 -v 24000000 #16 cores
Rscript --verbose --vanilla 7_ClimateMetrics.R i=\"$SGE_TASK_ID\" 
",sep=""),file="ClimateScript")

### Check the files
#system("cat 4b_ConvertKrigObjects.R")
system("cat ClimateScript")

## Submit it!
system("qsub ClimateScript")

## check progress
system("qstat -f") ; print(paste(max(0,length(system("qstat",intern=T))-2)," processes running"))
system("cat log/convertlog/ClimateScript.o55871.3011")
list.files("output/metrics")
system("df -h ")
system(" qstat -j 56185 | grep vmem")
system("qstat -g t")
#system("qdel 55871.3561")
## disconnect from farm
q("no")
exit


### Climate metrics
## Mount ramdisk for processing http://www.linuxscrew.com/2010/03/24/fastest-way-to-create-ramdisk-in-ubuntulinux
if(!file.exists("/tmp/ramdisk")) system("sudo mkdir /tmp/ramdisk; sudo chmod 777 /tmp/ramdisk")
system("sudo mount -t tmpfs -o size=10G tmpfs /tmp/ramdisk/")

system("cdo griddes Climate_0.05.nc > climategrid_0.05.txt")

system.time(mclapply(1:nrow(iter),function(i) # takes 17 hours
         system(paste("Rscript --vanilla /media/Data/Work/projects/WeatherInterpolation/analysis/7_ClimateMetrics.R ",
                      "i=",i," >> metriclog.txt",sep=""))))

#system("sudo umount /tmp/ramdisk/")

# 165s for 24 jobs at 1min resolution; takes ~7s per iteration-year. Will take ~5.1 days to complete
# 43s  for 24 jobs at 2min resolution; takes ~2s per iteration-year. Will take 31 hours to complete
# 26s for 24 jobs at 0.05degree res;  18 hours for all files 
#}}}

#{{{ Some random junk
## distance of various resolutions
res=.05 
lat=-33
lon=18

round(distm(data.frame(lon=c(lon,lon,lon+res),lat=c(lat,lat-res,lat)))/1000,3)


#}}}
