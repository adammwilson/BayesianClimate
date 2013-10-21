
############################################################################################
############################################################################################
### Process the interpolated datasets to calculate indices using nco
#system("ssh davisfarm")
#qrsh #-pe threaded 2 #-now n
#module load gcc openmpi R szip/2.1 hdf/5.1.8 netcdf cdo nco

#R
#i=1396

##First read in the arguments listed at the command line
args=(commandArgs(TRUE))
##args is now a list of character vectors
## First check to see if arguments are passed.
## if no args were given, print a warning and stop
if(length(args)==0) {stop("No parameters supplied, you must pass parameters")}

## Then cycle through each element of the list and evaluate the expressions.
eval(parse(text=args))
## now there is an i that corresponds to the row in the notdone object that will be processed.

## load the list of 'not done' dates
load("iternotdone.Rdata")

## identify bin in progress
it=iter$iter[i]
vari=iter$var[i]
yeari=iter$year[i]


## load libraries
## running on farm?
farm=!grepl("protea",system("uname -a",intern=T))
if(!farm){
	library(sp);library(geoR);library(ncdf4)
	setwd("/media/Data/Work/Regional/CFR/Weather/ClimateInterpolation")
        ## Working directories
        scratch=paste("/tmp/ramdisk/",i,sep="")
        if(!file.exists(scratch)) dir.create(scratch,recursive=T)
        ## destination directory
        version="v1.3b"
        odir=paste("output/",version,"/metrics",sep="")
        verbose=T
        suppressMessages(require(Hmisc))
        suppressMessages(require(gsubfn))
        if(!file.exists(odir)) dir.create(odir,recursive=T)
        ## use coarse climate data data?
        coarse=T
        res="0.05"
      }
## running on farm
if(farm) {
	rlib="/group/latimergrp/R/x86_64-redhat-linux-gnu-library/2.13"
        setwd("/group/latimergrp/adam/cfr/interp")
#	require(plyr,lib=rlib)
#	require(sp,lib.loc=rlib)
#	require(geoR,lib.loc=rlib)
	require(ncdf4,lib.loc=rlib)
#	require(reshape,lib.loc=rlib)
#	require(raster,lib.loc=rlib)
        ## Working directories
        scratch=paste("/scratch/zzwilson/",i,sep="")
        if(!file.exists(scratch)) dir.create(scratch,recursive=T)
        ## destination directory
        version=""  #leave blank if not keeping different versions
        odir="output/metrics"
        if(!file.exists(odir)) dir.create(odir,recursive=T)
        verbose=T
        ## use coarse climate data data?
        coarse=F
        res=""
      }


### Print some summary information (to log)
if(verbose){
  print(paste("Process running on ",system("uname -a",intern=T)))
  print(paste("Processing array element: ", i,sep=""))
  print(paste("Processing variable: ", vari,sep=""))
  print(paste("Processing iteration: ", it,sep=""))
  print(paste("Processing year: ", yeari,sep=""))
  print(paste("Temporary Directory",tempdir()))
  print(system("df -h",intern=T))
}

#######################################################################
#######################################################################
### Combine daily files create a single timeseries for each iteration
if(verbose) print("pulling out timeseries for one iteration")

while(!file.exists(paste(scratch,"/weather_fine_unpacked.nc",sep=""))) {
  ## use while because sometimes cdo merge fails...  ???
  
if(!farm) system(paste("ncrcat -O ",paste("-d iter,",it,",",it," ",sep="")," output/",version,"/nc_all/",vari,"_",yeari,"* ",scratch,"/",vari,".nc",sep=""))
if(farm) system(paste("ncrcat -O ",paste("-d iter,",it,",",it," ",sep="")," output/nc/",vari,"_",yeari,"* ",scratch,"/",vari,".nc",sep=""))

### Delete the iteration dimension to make cdo tools happy later
system(paste("ncwa -O -a iter ",scratch,"/",vari,".nc ", scratch,"/",vari,".nc",sep=""))

## merge all variables 
time_origin=as.Date(strsplit(ncatt_get(nc_open(paste(scratch,"/",vari,".nc",sep="")),"time")$units," ")[[1]][3])


## remap it to climate grid and mask out pixels > 5km from CFR
if(verbose) print("Remapping to climate grid")

if(!coarse) system(paste("cdo -s ifthen -invertlat -ltc,5 -selvar,CFR_dist CFR.nc -remapbic,climategrid.txt ",
             " -selvar,",vari," ",scratch,"/",vari,".nc ",scratch,"/fine_anom.nc",sep="")) 
if(coarse) system(paste("cdo -s ifthen -invertlat -ltc,5 -selvar,CFR_dist CFR_",res,".nc -remapbic,climategrid_",res,".txt ",
             " -selvar,",vari," ",scratch,"/",vari,".nc ",scratch,"/fine_anom.nc",sep="")) 


## remove temporary files
file.remove(paste(scratch,"/",vari,".nc",sep=""))
	
## transform back to real units using the climate data
if(verbose) print("Transforming back to real units")
  system(paste("cdo -s -O merge ",
               ## tmax
                ifelse(vari=="tmax",
                       paste("-ymonadd -mulc,-1 -selvar,tmax ",scratch,"/fine_anom.nc -selvar,tmax ",
                             ifelse(coarse,paste(" Climate_",res,".nc ",sep="")," Climate.nc "),sep=""),""),
                ifelse(vari=="tmin",
                       paste("-ymonadd -mulc,-1 -selvar,tmin ",scratch,"/fine_anom.nc -selvar,tmin ",
                             ifelse(coarse,paste(" Climate_",res,".nc ",sep="")," Climate.nc "),sep=""),""),
               ## take maximum of 0 or ppt to drop negative values using gtc,9999.
               ## Then multiply by mean ppt+1 (1 was added before kriging to remove 0s.)
                ifelse(vari=="ppt",
                       paste(" -ymonmul -max -gtc,99999 -selvar,ppt ",
                             scratch,"/fine_anom.nc -selvar,ppt ",
                             scratch,"/fine_anom.nc -addc,1 -selvar,ppt",
                             ifelse(coarse,paste(" Climate_",res,".nc ",sep="")," Climate.nc "),
                             sep=""),
                       ""),
               " ",scratch,"/weather_fine_unpacked.nc",
               sep=""))
} #close while loop

## remove temporary files
file.remove(paste(scratch,"/fine_anom.nc",sep=""))

## remove temporary files
#file.remove(paste(scratch,"/weather_fine_unpacked.nc",sep=""))
            
##########################################################################
##########################################################################
if(verbose) print("Starting processing of Climate Metrics")

## input file for climate metrics                                            
ifile=paste(scratch,"/weather_fine_unpacked.nc",sep="")


if(vari=="tmin"){
##########################################
### Frost Days - FD
### The number of days where TN < 0 ℃ per time period
### https://code.zmaw.de/embedded/cdo/1.4.7/cdo.html#x1-6440002.16.2
  system(paste("cdo -s -eca_fd -addc,273.15 -selname,tmin ",ifile," ",odir,"/FD_",yeari,"_",it,".nc",sep=""))

  ##########################################
### Consecutive Frost Days - CFD
### The largest number of consecutive days where TN < 0 ℃ per time period
### https://code.zmaw.de/embedded/cdo/1.4.7/cdo.html#x1-6440002.16.2
  system(paste("cdo -s -eca_cfd -addc,273.15 -selname,tmin ",ifile," ",odir,"/CFD_",yeari,"_",it,".nc",sep=""))


  ##########################################
  ## Growing Degree Days - GDD
if(verbose)   print("Growing Degree Days")
  ## Growing degree days per time period
  ## use heating degree days function but invert by multiplying by -1 and using -10 threshold
  system(paste("cdo -s -eca_hd,-10 -mulc,-1 -addc,273.15 -selname,tmin ",ifile," ",odir,"/GDD_",yeari,"_",it,".nc",sep=""))
  ### Update variable name and description
#  system(paste("ncdump -h ",odir,"/GDD_",yeari,"_",it,".nc",sep=""))
  system(paste("ncrename  -v heating_degree_days_per_time_period,growing_degree_days ",odir,"/GDD_",yeari,"_",it,".nc",sep=""))
  system(paste("ncatted  -a long_name,growing_degree_days,o,c,",
               "\"Growing degree days are the sum of the difference between 10C and daily min temperature Y\" ",
               odir,"/GDD_",yeari,"_",it,".nc",sep=""))

#########################
## Mean Minimum Temperature (MATmin) 
if(verbose) print("Mean Annual Minimum Temperature")
  system(paste("cdo -s timmean -selname,tmin ",ifile," ",odir,"/MATmin_",yeari,"_",it,".nc",sep=""))

    #########################
  ## Mean Monthly Minimum Temperature (MMTmin)
if(verbose)  print("Mean Annual Minimum Temperature")
  system(paste("cdo -s monmean -selname,tmin ",ifile," ",odir,"/MMTmin_",yeari,"_",it,".nc",sep=""))

  #########################
  ## Monthly Min Temperature (MTmin)
if(verbose)  print("Monthly Minimum Temperature")
  system(paste("cdo -s monmin -selname,tmin ",ifile," ",odir,"/MTmin_",yeari,"_",it,".nc",sep=""))

  
} ## End Tmin

if(vari=="tmax"){
  #########################
  ## Summer Heat Waves (SU)
if(verbose)  print("Summer Heat Waves")
  ## ECASU - Summer days index per time period
  ## The  number of days where daily temps are greater than 25 ℃.
  ## https://code.zmaw.de/embedded/cdo/1.4.7/cdo.html#x1-6480002.16.3
  system(paste("cdo -s -eca_su  -addc,273.15 -selname,tmax ",ifile," ",odir,"/SU_",yeari,"_",it,".nc",sep=""))

  #########################
  ## Consecutive Summer Heat Waves (CSU)
if(verbose)  print("Summer Heat Waves")
  ## ECACSU - Consecutive summer days index per time period
  ## The largest number of consecutive days where daily temps
  ## are greater than 35 ℃.
  ## https://code.zmaw.de/embedded/cdo/1.4.7/cdo.html#x1-6480002.16.3
  system(paste("cdo -s -eca_csu,35  -addc,273.15 -selname,tmax ",ifile," ",odir,"/CSU_",yeari,"_",it,".nc",sep=""))

  #########################
  ## Mean Maximum Temperature (MATmax)
if(verbose)  print("Mean Annual Maximum Temperature")
  system(paste("cdo -s timmean -selname,tmax ",ifile," ",odir,"/MATmax_",yeari,"_",it,".nc",sep=""))

  #########################
  ## Mean Monthly Maximum Temperature (MMTmax)
if(verbose)  print("Mean Annual Maximum Temperature")
  system(paste("cdo -s monmean -selname,tmax ",ifile," ",odir,"/MMTmax_",yeari,"_",it,".nc",sep=""))

 #########################
  ## Monthly Max Temperature (MTmax)
if(verbose)  print("Monthly Maximum Temperature")
  system(paste("cdo -s monmax -selname,tmax ",ifile," ",odir,"/MTmax_",yeari,"_",it,".nc",sep=""))

  
} ## End tmax

if(vari=="ppt"){
#eca pd Precipitation days index per time period
  
#########################
  ## Consecutive dry days (CDD)
if(verbose)  print("Consecutive dry days")
  ## The largest number of consecutive days where RR is < 2 mm per year (by subtracting 1 before calculating)
  ## https://code.zmaw.de/embedded/cdo/1.4.7/cdo.html#x1-6400002.16.1
  system(paste("cdo -s setname,max_consecutive_dry_days -selvar,consecutive_dry_days_index_per_time_period -eca_cdd -subc,1 -selname,ppt ",ifile,
               " ",odir,"/CDD_",yeari,"_",it,".nc",sep=""))

  #########################
  ## Total Wet days (rr1)
if(verbose)  print("Total wet days")
  ## The number of days where RR is > 2 mm per year
  ## https://code.zmaw.de/embedded/cdo/1.4.7/cdo.html#x1-6400002.16.1
  system(paste("cdo -s -eca_rr1 -subc,1 -selname,ppt ",ifile," ",odir,"/rr2_",yeari,"_",it,".nc",sep=""))

  #########################
  ## Simple daily Preciptition intensity (sdii)
if(verbose)  print("Precipitation intensity")
  ## The mean precipititation for days where RR is > 1 mm
  ## https://code.zmaw.de/embedded/cdo/1.4.7/cdo.html#x1-6400002.16.1
  system(paste("cdo -s eca_sdii -selname,ppt ",ifile," ",odir,"/sdii_",yeari,"_",it,".nc",sep=""))

  #########################
  ## Precipitation Days Index (r20mm)
  ## The number of days per year where daily precipitation is at least 20 mm
  ## https://code.zmaw.de/embedded/cdo/1.4.7/cdo.html#x1-6980002.16.14
   system(paste("cdo -s -eca_r20mm -selname,ppt ",ifile," ",odir,"/ECAr20mm_",yeari,"_",it,".nc",sep=""))

  #########################
  ## Annual Precipitation (MAP)
if(verbose)  print("Annual Precipitation")
  system(paste("cdo -s yearsum -selname,ppt ",ifile," ",odir,"/MAP_",yeari,"_",it,".nc",sep=""))

  #########################
  ## Monthly Precipitation (MMP)
if(verbose)  print("Monthly Precipitation")
  system(paste("cdo -s monsum -selname,ppt ",ifile," ",odir,"/MMP_",yeari,"_",it,".nc",sep=""))
}

### Remove all temporary files so we don't blow anything up!
if(verbose)  print("Removing all temporary files")
  system(paste("rm -R ",scratch))
  
## Finished                                            
print(paste("################################################################    Finished job",i,"out of",nrow(iter),", quitting R"))
q("no")
