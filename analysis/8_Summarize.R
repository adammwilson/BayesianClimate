#system("ssh -X protea")
#R

setwd("/media/Data/Work/Regional/CFR/Weather/ClimateInterpolation")
library(reshape);library(lattice);library(rgdal);library(latticeExtra);library(multicore);library(ncdf4)
library(sp);library(geoR);library(raster);library(abind);library(gridBase)

version="v1.3b"
## get time_origin
time_origin=as.Date(strsplit(ncatt_get(nc_open(paste("output/",version,"/weather_fine.nc",sep="")),"time")$units," ")[[1]][3])

##############################################################
##############################################################

st=readOGR("StationLocations.shp","StationLocations")
load(paste("output/",version,"/valid_0.25deg.Rdata",sep=""))
vars=c("tmax","tmin","ppt")
## import station data to use for validation
stnc=nc_open("StationData.nc",write=F)

## Set precip that's <1mm to 0 to eliminate tiny precips

valid$fine_pred_mean[valid$var=="ppt"&valid$fine_pred_mean<1]=0
#valid=valid[!(valid$var=="ppt"&valid$fine_pred_mean>200),]

### import metadata
time=ncvar_get(stnc,"time")
lat=as.vector(ncvar_get(stnc,"lat"))
lon=as.vector(ncvar_get(stnc,"lon"))

table(st$type)

#{{{ How long did it take to run the interpolation?
times=read.csv(paste("output/",version,"/times.csv",sep=""))
times$duration[times$units=="mins"]=times$duration[times$units=="mins"]/60
times$units[times$units=="mins"]="hours"
bwplot(duration~variable,data=times,horizontal=F,scales=list(y=list(log=T)))
## total computation time in days
round(sum(times$duration)/24,1)

                                        #}}}

#{{{#############  Do sample krige to show parameters
vari="tmax"
idate=as.Date("2009-01-03")
itime=which(time_origin+time==idate)

## Example krige parameters
data=data.frame(
  station=as.vector(ncvar_get(stnc,"station")),
  lat=as.vector(ncvar_get(stnc,"lat")),
  lon=as.vector(ncvar_get(stnc,"lon")),
  value=as.numeric(ncvar_get(stnc,vari,start=c(itime,1),count=c(1,-1))))
data$dcfr=as.vector(extract(raster("CFR.nc"),data[,c("lon","lat")]))
data=data[data$dcfr<=75&!is.na(data$value)&!is.na(data$dcfr)&data$value!=-999,]

  ###########################
  ## apply climate correction
  clim=brick("Climate.nc",varname=vari)
  clim=subset(clim,subset=as.numeric(format(idate,"%m")))
  ## drop any stations outside the climate region
  data=data[data$lat>=bbox(clim)[2,1]&data$lat<=bbox(clim)[2,2]&data$lon>=bbox(clim)[1,1]&data$lon<=bbox(clim)[1,2],]
  data$clim=as.numeric(extract(clim,data[,c("lon","lat")],method="bilinear"))
  if(vari%in%c("tmax","tmin"))   data$climanomaly=data$clim-data$value


dat=as.geodata(data,  coords.col=c("lon","lat"),data.col="climanomaly")
dat=jitterDupCoords(dat,max=.01)
vg=variog(dat, max.dist = 12,trend="cte")
vf=variofit(vg,cov.model="exponential")
vf

mod=model.control(trend.d=trend.spatial("cte",dat),cov.model="exponential")
  out=output.control(n.posterior=100,n.predictive=100,signal=T)
  pri=prior.control(#beta.prior="normal",#beta=rep(0,3),beta.var.std=diag(10,3),
    sigmasq.prior="reciprocal",
    phi.prior="exponential",phi=8,phi.discrete=seq(0,15,length.out=201), # plot(pri$priors.info$phi$probs~names(pri$priors.info$phi$probs),ylim=c(0,.02))
    tausq.rel.prior="reciprocal",tausq.rel.discrete=seq(0,1,length.out=51)) #tau goes 0-1

  ## Run final version and make the predictions
  gc(); kc=krige.bayes(dat,model=mod,prior=pri,output=out); gc()
save(vg,vf,kc,file="ExampleKrige.Rdata")

## now there are vg,vf,kc objects for plotting later
#}}}

#{{{  Summarize Validation
## accuracy assessment function
accuracy<-function(obs,pred,ppt=F,predthresh=1,logoffset=1,verbose=T) {
	## Calculate confusion matrix of rain/no rain
	tt=addmargins(table(Observed=obs>2,Predicted=pred>=predthresh))
	## compute logs if desired for calculating R2
	if(ppt) {
		if(any(obs<=0|pred<=0)){	lobs=log(obs+logoffset);lpred=log(pred+logoffset)}
		if(!any(obs<=0|pred<=0)){   lobs=log(obs);lpred=log(pred)}
               }

	## generate some summary stats for PPT
	if(ppt) tm=round(matrix(c(
                  tt[2,2]/(tt[2,2]+tt[2,1]), #TruePositiveRate
                  tt[1,1]/(tt[1,1]+tt[1,2]), #TrueNegativeRate
                  tt[2,2]/(tt[2,2]+tt[1,2]), #Positive Predicted value
                  tt[1,1]/(tt[1,1]+tt[2,1]), #Negative Predicted value
                  sqrt(sum((obs-pred)^2,na.rm=T)/length(obs)), #RMSE 
                  (1/length(obs))*sum(abs(pred-obs),na.rm=T), #mae
                  (1/length(obs))*sum(pred-obs,na.rm=T), #mer
                  length(obs), #n
                  summary(lm(obs~pred))$r.squared, #R2)
                  summary(lm(lobs~lpred))$r.squared), # log-log R2)
                  dimnames=list(c("Sensitivity (True Positive Rate)",
                    "Specificity (True Negative Rate)",
                    "Positive Predicted Value (Predicted Wet that were Wet)",
                    "Negative Predicted Value (Predicted Dry that were Dry)",
                    "RMSE (Root Mean Square Errors)",
                    "MAE (Mean Absolute Error)",
                    "MER (Mean Error)",
                    "n",
                    "R^2",
                    "log-log R^2"), "Value")),3)
	## generate some summary stats for temp
	if(!ppt) tm=round(matrix(c(
                   sqrt(sum((obs-pred)^2,na.rm=T)/length(obs)), #RMSE 
                   (1/length(obs))*sum(abs(pred-obs),na.rm=T), #mae
                   (1/length(obs))*sum(pred-obs,na.rm=T), #mer
                   length(obs), #n
                   summary(lm(obs~pred))$r.squared), #R2)
                   dimnames=list(c("RMSE (Root Mean Square Errors)",
                     "MAE (Mean Absolute Error)",
                     "MER (Mean Error)",
                     "n",
                     "R^2"), "Value")),3)
        
	if(verbose){
		if(ppt) {print("Validation data (numbers indicate numbers of observations in each category)")
                          print(tt)}
		print("Summary of Confusion matrix:")
		print(tm)
	}
	return(list(confmat=tt,summary=tm))}

### Run validation function for daily hold-out data
ac=lapply(vars,function(vari) {
  tv=valid[valid$var==vari,]
  return(accuracy(obs=tv$obs,pred=tv$fine_pred_mean,ppt=vari=="ppt",predthresh=2))
})


#{{{ Validation of monthly data


## Get station data
vd=read.csv("/media/Data/Work/Regional/CFR/Weather/data/obs/IndependantWeatherStations/Swartberg_Weather_Data.csv")
vdl=melt(vd,id.vars=c("Station","Year","Month"),measure.vars=c("MaxTemp","MinTemp","Rainfall"))
vdl$date=as.Date(paste(vdl$Year,vd$Month,15),format="%Y %B %d")
vdl=vdl[!is.na(vdl$value),c("Station","variable","value","date")]

## identify unique stations
vds=unique(vd[!is.na(vd$Longitude),c("Station","Latitude","Longitude")])
coordinates(vds)=c("Longitude","Latitude")
vds@data[,c("lon","lat")]=coordinates(vds)

vdm=do.call(rbind.data.frame,mclapply(c("MMTmax","MMTmin","MMP"),function(m) {
  nc=nc_open(paste("output/",version,"/metrics_summary_coarse/",m,".nc",sep=""))
  vds$col=apply(matrix(vds$lon),1,function(x) which.min(abs(x-ncvar_get(nc,"lon"))))
  vds$row=apply(matrix(vds$lat),1,function(x) which.min(abs(x-ncvar_get(nc,"lat"))))
 ## loop through metrics
  vdm=do.call(rbind.data.frame,lapply(1:nrow(vds),function(i) {
    tvdm=ncvar_get(nc,NA,start=c(vds$col[i],vds$row[i],1,1),count=c(1,1,-1,-1))
    all=F
    if(all){  #return all iterations?
      time=ncvar_get(nc,"time")
      dimnames(tvdm)=list(time=time,iteration=1:1000)
      tvdml=melt(tvdm)
      tvdml$value[tvdml$value==-999]=NA
      tvdml$location=vds$Station[i]
      return(tvdml)
    }
    if(!all){
      tvdm=melt(apply(tvdm,1,quantile,c(0.025,0.5,0.975)))
      colnames(tvdm)=c("quantile","time","value")
      tvdm$time=ncvar_get(nc,"time")[tvdm$time]
      tvdm$value[tvdm$value==-999]=NA
      tvdm$location=vds$Station[i]
      return(tvdm)
    }
  }))
  vdm$metric=m
  print(paste("Finished",m))
  return(vdm)
}))

## add dates
time_origin=as.Date(strsplit(ncatt_get(nc_open(paste("output/",version,"/metrics_summary/MATmax.nc",sep="")),"time")$units," ")[[1]][3])
vdm$date=as.Date(format(time_origin+vdm$time,"%Y-%m-15"))
vdm$metric=ifelse(vdm$metric=="MTmax","tmax",ifelse(vdm$metric=="MTmin","tmin","ppt"))


## merge with vd
vdl$variable=ifelse(vdl$variable=="MaxTemp","tmax",ifelse(vdl$variable=="MinTemp","tmin","ppt"))
vdc=merge(
  vdm[,c("location","metric","quantile","value","date")],
  vdl,by.x=c("date","metric","location"),by.y=c("date","variable","Station"))
colnames(vdc)[grep("value.x",colnames(vdc))]="Predicted"
colnames(vdc)[grep("value.y",colnames(vdc))]="Observed"
vdc=cast(vdc,date+metric+location+Observed~quantile,value="Predicted")
vdc$month=factor(format(vdc$date,"%B"),levels=months(as.Date(paste(2000,1:12,15),"%Y %m %d")),ordered=T)
vdc$season=ifelse(format(vdc$date,"%m")%in%c("12","01","02"),"Summer",
  ifelse(format(vdc$date,"%m")%in%c("03","04","05"),"Fall",
         ifelse(format(vdc$date,"%m")%in%c("06","07","08"),"Winter",
                ifelse(format(vdc$date,"%m")%in%c("09","10","11"),"Spring",NA))))
vdc$season=factor(vdc$season,levels=c("Spring","Summer","Fall","Winter"),ordered=T)
colnames(vdc)[grep("%",colnames(vdc))]=c("q2.5","q50","q97.5")

### Drop bad data
vdc2=vdc[!grepl("Blesberg",vdc$location),]
vdc2=vdc2[!(grepl("Top",vdc2$location)&vdc2$metric%in%c("tmax","tmin")),]
vdc2=vdc2[!(grepl("Albertsberg",vdc2$location)&vdc2$metric=="tmin"),]
vdc2=vdc2[!(grepl("De Wetsvlei",vdc2$location)&vdc2$metric=="tmin"),]
vdc2=vdc2[!(grepl("De Wetsvlei",vdc2$location)&vdc2$metric=="tmin"),]
vdc2=vdc2[!(grepl("phuisv",vdc2$location)&vdc2$metric=="tmin"),]
vdc2=vdc2[!grepl("Gamkapoort",vdc2$location),]
vdc2=vdc2[!grepl("Dwyka hek",vdc2$location),]
vdc2=vdc2[!grepl("Auto Besemfontein",vdc2$location),]
vdc2=vdc2[!(vdc2$date>as.Date("1998-01-01")&grepl("^Auto$",vdc2$location)&vdc2$metric!="ppt"),]
vdc2=vdc2[!(vdc2$date>as.Date("1995-01-01")&grepl("Wetsvlei",vdc2$location)&vdc2$metric=="tmax"),]
vdc2=vdc2[!(vdc2$date>as.Date("2000-01-01")&vdc2$metric!="ppt"),]



### Run validation function for daily hold-out data
acm=lapply(vars,function(vari) {
  tv=vdc2[!is.na(vdc2$Observed)&vdc2$metric==vari,]
  return(accuracy(obs=tv$Observed,pred=tv$q50,ppt=vari=="ppt",predthresh=2))
})


## print pdf
pdf(paste("output/",version,"/monthlyvalidation.pdf",sep=""),width=11,height=8.5)


### show time series
vdcl=melt(as.data.frame(vdc2),measure.vars=c("Observed","q2.5","q50","q97.5"))
useOuterStrips(xyplot(value~date|metric+location,data=vdcl,panel=function(x,y,subscripts){
  td=vdcl[subscripts,]
  if(nrow(td)>1){
    panel.arrows(td$date[td$variable=="q2.5"],td$value[td$variable=="q2.5"],td$date[td$variable=="q2.5"],td$value[td$variable=="q97.5"],len=0,col="blue")
    panel.xyplot(td$date[td$variable=="q50"],y[td$variable=="q50"],pch=16,cex=.5,type="l",col="blue",lwd=.5)
    panel.xyplot(td$date[td$variable=="Observed"],y[td$variable=="Observed"],pch=16,cex=.5,type="l",col="red",lwd=1.5)
  }
},scales=list(y=list(relation="free"),log=F),
             main="Model Validation using independent monthly data \n from the Swartberg Region",sub="Grey bars are 95% credible intervals",
             ylab="Predicted",xlab="Observed",layout=c(length(unique(vdcl$location)),length(unique(format(vdcl$date,"%y"))),length(unique(vdcl$metric)))))

### Show scatterplots
for(m in unique(vdc$metric))
  print(xyplot(q50~Observed|location,data=vdc2[!is.na(vdc2$Observed)&vdc2$metric==m,],panel=function(x,y,subscripts){
  td=vdc2[!is.na(vdc2$Observed)&vdc2$metric==m,][subscripts,]
  panel.arrows(td$Observed,td$q2.5,td$Observed,td$q97.5,len=0,col="grey")
#  panel.arrows(log(td$Observed),log(td$q2.5),log(td$Observed),log(td$q97.5),len=0,col="grey")
  panel.xyplot(x,y,pch=16,cex=.5)
  panel.abline(0,1,col="red")
  if(td$metric[1]=="ppt") {
#    tlm=lm(I(log(td[,"50%"]+1))~I(log(td$Observed+1)))
    tlm=lm(td$q50~td$Observed)
    panel.key(text = c(paste("R^2=",round(summary(tlm)$r.squared,2))),lines=F,points=F,corner = c(0,.98))
    panel.abline(tlm,col="blue")
  }
  if(td$metric[1]!="ppt") {
    tlm=lm(td$q50~td$Observed)
    panel.key(text = c(paste("R^2=",round(summary(tlm)$r.squared,2))),
              lines=F,points=F,corner = c(0,.98))
    panel.abline(tlm,col="blue")
  }
},scales=list(log=F),
             main=paste(m,"Model Validation using independent monthly data \n from the Swartberg Region"),sub="Grey bars are 95% credible intervals",
             ylab="Predicted",xlab="Observed",asp=1))

dev.off()
#}}}


## make table of validation results
library(xtable)
tab= cbind.data.frame(Variable=c("Maximum Temperature","Minimum Temperature","Precipitation"),
                   RMSE=do.call(c,lapply(ac,function(x) x$summary["RMSE (Root Mean Square Errors)",])),
                   MAE=do.call(c,lapply(ac,function(x) x$summary["MAE (Mean Absolute Error)",])),
                   MER=do.call(c,lapply(ac,function(x) x$summary["MER (Mean Error)",])),
                   r2=do.call(c,lapply(ac,function(x) x$summary["R^2",])),
                   RMSE=do.call(c,lapply(acm,function(x) x$summary["RMSE (Root Mean Square Errors)",])),
                   MAE=do.call(c,lapply(acm,function(x) x$summary["MAE (Mean Absolute Error)",])),
                   MER=do.call(c,lapply(acm,function(x) x$summary["MER (Mean Error)",])),
                   r2=do.call(c,lapply(acm,function(x) x$summary["R^2",]))
  )

my.add.to.row <- function(x) {
    first <- "\\hline \\multicolumn{1}{c}{} & "
    middle <- paste(paste("\\multicolumn{2}{c}{(", seq(x), ")}", sep = ""), collapse = " & ")
    last <- paste("\\\\ \\cline {", 2, "-", 1 + 2 * x, "}", collapse = "")
    string <- paste(first, middle, last, collapse = "")
    list(pos = list(-1), command = string)
}

print(xtable(tab,label="tab:valid",
             caption="Validation results for each variable. The comparison of daily data are from the daily observations not included in model fitting. The monthly data compare predicted and observed total monthly precipitation and monthly maximum and minium temperature at a set of remote stations.  The poorer fit for the monthly temperature is expected as it represents the ability of the model to estimate the single daily maximum and minimum in each month, while the monthly precipitation comparison is the aggregated total for the month. RMSE: Root Mean Squared Errors, MAE: Mean Absolute Error, MER: Mean Error",
             align="llrrrr|rrrr"),
             add.to.row=list(pos=list(-1), command="\\hline & \\multicolumn{4}{c|}{Daily data} & \\multicolumn{4}{c}{Monthly data} \\\\"),
             include.rownames=F,file=paste("output/",version,"/ValidationTable.tex",sep=""))

                                        #}}}

table(vdl$variable)  #number of observations in monthly validation dataset


#{{{  Extract data from points A & B to illustrate uncertainty
expts=data.frame(lat=c(-32.278,-33.959),lon=c(19.112,18.4090),name=c("Cederberg","Table Mountain"))
coordinates(expts)=c("lon","lat")
expts@data[,c("lon","lat")]=coordinates(expts)
## summary of station distances for ms (distance (km) of 10 closest stations)
apply(spDists(expts,st,longlat=T),1,function(x) sort(x)[1:15])

exgrid=paste(expand.grid(vars,c("mean","sd"))[,1],expand.grid(vars,c("mean","sd"))[,2],sep="_")
expred=mclapply(exgrid,function(x){
  source=paste("output/",version,"/weather_fine_unpacked.nc",sep="")
  #source=paste("output/",version,"/weather_fine_unpacked.nc",sep="")
  t=extract(brick(source,varname=x),expts)
  dimnames(t)=list(name=expts$name,date=as.character(time_origin+ncvar_get(nc_open(source),"time")))
  return(t)
  })
names(expred)=exgrid
expredl=melt.list(expred);expredl$date=as.Date(expredl$date)
expredl[,c("var","type")]=do.call(rbind,strsplit(expredl$L1,split="_"))
expredl=cast(expredl,date+name+var~type,value="value")

                                        #}}}

#{{{ Explore climate metrics

##test one
##m="CDD"
##nc=nc_open(paste("output/",version,"/metrics_summary/",m,".nc",sep=""))

## extract all values of all metrics for locations in expts
exmet=do.call(rbind.data.frame,lapply(sub(".nc","",list.files(paste("output/",version,"/metrics_summary",sep=""),pattern="[.]nc$")),function(m) {
  nc=nc_open(paste("output/",version,"/metrics_summary/",m,".nc",sep=""))
  expts$col=apply(matrix(expts$lon),1,function(x) which.min(abs(x-ncvar_get(nc,"lon"))))
  expts$row=apply(matrix(expts$lat),1,function(x) which.min(abs(x-ncvar_get(nc,"lat"))))
  d2=do.call(rbind.data.frame,lapply(1:nrow(expts),function(i) {
    d2=ncvar_get(nc,NA,start=c(expts$col[i],expts$row[i],1,1),count=c(1,1,-1,-1))
    dimnames(d2)=list(year=1990:2009,iteration=1:ncol(d2))
    d2l=melt(d2)
    d2l$value[d2l$value==-999]=NA
    d2l$location=expts$name[i]
    d2l$year=as.factor(d2l$year)
    return(d2l)
  }))
  d2$metric=m
  print(paste("Finished",m))
  return(d2)
}))

## clean up values
exmet$value=as.numeric(exmet$value)
exmet=exmet[!is.na(exmet$value),]
exmet$year=as.numeric(as.character(exmet$year))

## save them
save(exmet,file=paste("output/",version,"exmet.Rdata"))

## load them
load(paste("output/",version,"exmet.Rdata"))

## toss weird values
exmet$value[exmet$metric=="CDD"&exmet$value>=365]=NA
exmet$value[exmet$metric=="MATmax"&exmet$value==0]=NA
exmet$value[exmet$metric=="MATmin"&exmet$value==0]=NA
exmet=exmet[exmet$metric!="CSU",]
      
#{{{ brick is confused by 4D arrays, couldnt' make it work correctly
#r=do.call(rbind.data.frame,lapply(1:20,function(i) {
#  b=readAll(brick(paste("output/",version,"/metrics_summary/GDD.nc",sep=""),level=i))
#  r=extract(b,expts,layer=1,nl=c(1000,1))
#  r=t(r)
#  dimnames(r)=list(iteration=1:1000,location=expts$name)
#  rl=melt(r)
#  rl$date=i
#  rl$value[rl$value==-999]=NA
#  return(rl)
#}))
#}}}

pdf(paste("output/",version,"/MetricEvaluationAtPoints.pdf",sep=""),width=11,height=8.5)
       
## xyplot(value~year|location,data=exmet,type="l",panel=function(x,y,subscripts){
##   dt=exmet[subscripts,]
##   dt=dt[order(dt$iteration,dt$year),]
##   for(i in 1:1000) panel.xyplot(dt$year[dt$iteration==i],dt$value[dt$iteration==i],type="l",col="#C1CDCD30")
## },subset=T,scales=list(relation="free"),layout=c(1,2),ylim=c(0,150))

#bwplot(value~as.factor(year)|location+metric,data=exmet,horizontal=F,layout=c(2,1),fill="grey",pch="|",
#       ylab="value",xlab="Year",las=1,scales=list(x=list(rot=45))) ###       ylab="Longest period of cumulative dry days",

length(unique(exmet$location))*length(unique(exmet$metric))

useOuterStrips(combineLimits(bwplot(value~as.factor(year)|location+metric,data=exmet,horizontal=F,layout=c(2,7),fill="grey",pch="|",
       ylab="Climate Metric Value",xlab="Year",las=1,
                                    scales = list(x=list(rot=45),y=list(rot=0,relation="free")),
                                    par.settings = list(plot.symbol = list(pch = 16,cex=.3,col="black"),
                                      box.rectangle=list(col="black"),box.umbrella=list(col="black"),strip.background=list(col="transparent"),
                                      par.strip.text = list(cex =1))),margin.x = c(1, 2)))

### summary stats for paper
texmet=exmet[exmet$metric=="CDD",]
aggregate(texmet$value,list(year=texmet$year,location=texmet$location),
       function(x) c(mean=mean(x,na.rm=T),quantile(x,0.025,na.rm=T),quantile(x,0.975,na.rm=T),IQW=quantile(x,0.975,na.rm=T)-quantile(x,0.025,na.rm=T)))

## compare standard deviations for two locations
#sds=melt(tapply(exmet$value,list(year=exmet$year,metric=exmet$metric,location=exmet$location),function(x) sd(x,na.rm=T)/mean(x,na.rm=T)))
#xyplot(value~year|metric,groups=location,data=sds,type="l",scales=list(relation="free"),auto.key=T,main="Comparison of the Coefficient of Variation between sites")


dev.off()

                                        #}}}

smallpdf<-function(x) system(paste("convert -density 1200 -resize 3000x2000 ",x," tiff:- | ",
             "convert tiff:- -quality 95 -compress jpeg ",sub(".pdf","",x),"_small.pdf",sep=""),wait=F)


########################################################################
#{{{ plot maps
## plot the station maps
library(grid);library(gridBase)
    cfr=readOGR("/media/Data/Work/Regional/CFR/BaseGISData/CFR.shp","CFR")
    africacountries=readOGR("/media/Data/Work/Regional/USGS_data/Africa/database/mapbase/f_pol/vpf/polbnda.shp","polbnda")
    africacountries=africacountries[africacountries@data$AREA>10000000000,]  #get rid of really small areas
    coast=readOGR("/media/Data/Work/Regional/CFR/BaseGISData/CFR_coastline/coastline.shp","coastline")
#load("/media/Data/Work/Regional/CFR/BiomassAccumulation/BiomassAccumulation_500m/analysis/data/covariates.Rdata")
#load("/media/Data/Work/Regional/CFR/BiomassAccumulation/BiomassAccumulation_500m/analysis/data/ndvi.Rdata")
#ndvil=do.call(rbind.data.frame,mclapply(ndvi,function(x) unique(cbind(lat=x$lat,lon=x$lon))))
#coordinates(ndvil)=c("lon","lat")
#ndvil=ndvil[!is.na(overlay(ndvil,cfr)),]

    ## stations
    st=readOGR("StationLocations.shp","StationLocations")
    st$pch=ifelse(st$temp&st$ppt,8,ifelse(st$temp&!st$ppt,3,4))
    st$col=ifelse(st$temp&st$ppt,"black",ifelse(st$temp&!st$ppt,"red","blue"))

pdf(paste("output/",version,"/maps.pdf",sep=""),width=11,height=8.5)
### plot with two example points 
plot(cfr,col=grey(.95),new=T,ylab="test");box()
axis(1,las=1,xlab="Longitude");axis(2,las=1,ylab="Latitude")
mtext("Longitude",1,padj=4);mtext("Latitude",2,padj=-4)
#lines(coast)
points(st,col=st$col,pch=st$pch)
points(expts,col="black",cex=2,pch=16)
points(vds,col="black",cex=2)
#text(x=coordinates(expts),labels=expts$name,col="black",cex=2,pos=4)
text(20.4,-31.6,labels="Cederberg",col="black",cex=2,pos=4)
arrows(20.4,-31.6,19.112,-32.278,length=.4,lwd=3)
text(18.5,-35,labels="Table\nMountain",col="black",cex=2,pos=4)
arrows(18.75,-34.5,18.409,-33.959,length=.4,lwd=3)
legend("bottomright",legend=
     c("CFR","Precipitation","Temperature","Temperature and Precipitation","Validation Stations"),
       pch=c(15,4,3,8,1),col=c("grey","blue","red","black","black"),bty="n",cex=1.5)
pushViewport(viewport(x=.8,y=.72,width=.25,height=.25,name="d1"))
par(new=T,fig=gridFIG(),mar=c(.5,.1,.4,0)) #,mgp=c(1,3,2))
plot(africacountries,col=grey(.9),lty=1,lwd=.1,new=F,xlim=c(-20,60))
box()
plot(cfr,col="black",add=T,lty=0,bty="o")
arrows(0,-34,17,-34,length=.1)
popViewport();par(new=T,opar)  #reset the plotting parameters
dev.off()

     smallpdf(paste("output/",version,"/maps.pdf",sep=""))
#}}}

#########################################################################
#{{{ plot kriging parameters
pdf(paste("output/",version,"/kriging.pdf",sep=""),width=11,height=8.5)

for(i in 1:4){
   plot(vg,ylim=c(0, 30),xlim=c(0,11),pch=16,col="blue",
       main="Example variogram of single day's temperature data",xlab="Distance (degrees)")
  lines(vg,col="blue")
  legend("topright",legend=c("Empirical","Conventional","Bayesian"),col=c("blue","red","purple"),lty=1,bty="n",cex=1.5,lwd=2)
  if(i==2){
    lines(vf,col="red",lwd=2)
    arrows(vf$cov.pars[2],12,vf$cov.pars[2],0,col="red",length=0.1)
    text(vf$cov.pars[2],0,label=expression(paste("Range (",phi,")")),pos=4,col="red")
    arrows(20,vf$nugget+vf$cov.pars[1],-.3,vf$nugget+vf$cov.pars[1],col="red",length=.1)
    text(0,vf$nugget+vf$cov.pars[1]+.5,label=expression(paste("Sill (",tau+sigma^2,")")),pos=4,col="red")
    arrows(0,vf$nugget,-.3,vf$nugget,col="red",length=.1)
    text(0,vf$nugget,label=expression(paste("Nugget (",tau,")")),pos=1,col="red")
  
  }
  if(i%in%c(3)){
    # sigmasq
    kc.sig=density(kc$posterior$sample$sigmasq)
    polygon(c(kc.sig$y*100,rep(0,length(kc.sig$y))),c(kc.sig$x,rev(kc.sig$x)),col="#8B008B20")
    text(3,25,label=expression(paste("Sill (",tau+sigma^2,")")),pos=4,col="#8B008B")
    ## tausq
    kc.tau=density(kc$posterior$sample$tausq.rel*kc$posterior$sample$sigmasq)
    polygon(c(kc.tau$y*2,rep(0,length(kc.tau$y))),c(kc.tau$x,rev(kc.tau$x)),col="#8B008B20")
    text(1.3,3.8,label=expression(paste("Nugget (",tau,")")),pos=1,col="#8B008B")
    ## phi
    kc.phi=density(kc$posterior$sample$phi)
    polygon(c(kc.phi$x,rev(kc.phi$x)),c(kc.phi$y*100,rep(0,length(kc.phi$y))),col="#8B008B20")
    text(7,11,label=expression(paste("Range (",phi,")")),pos=1,col="#8B008B")
    ## lines
    my.summary <- function(x){quantile(x, prob = c(0.025, 0.5, 0.975))}
    lines(kc, summ = my.summary, ty="l", lty=c(2,1,2), col="#8B008B",lwd=2)
  }
}

## Single Figure for paper
par(mfrow=c(2,1),mar=c(0,4,0,2)+0.1, oma=c(5,0,3,0)+0.1)
ccol=grey(.4)
plot(vg,ylim=c(0, 20),xlim=c(-.25,8),pch=16,col="black",type="o",lwd=2,xlab="",cex.lab=1.5,las=1,xaxt="n",ylab="")
  legend("topleft",legend=c("Empirical","Conventional"),col=c("black",ccol),lty=1,pch=c(16,NA),bty="n",cex=1.5,lwd=2)
axis(4,las=1)
      lines(vf,col=ccol,lwd=2)
      arrows(vf$cov.pars[2],7.7,vf$cov.pars[2],0,col=ccol,length=0.1) #vf$nugget+vf$cov.pars[1]
      text(vf$cov.pars[2],0,label=expression(paste("Range (",phi,")")),pos=4,col=ccol)
      arrows(7.2,vf$nugget+vf$cov.pars[1]+.1,8.3,vf$nugget+vf$cov.pars[1],col=ccol,length=.1)
    text(7.1,vf$nugget+vf$cov.pars[1]+.5,label=expression(paste("Sill (",tau+sigma^2,")")),pos=4,col=ccol)
    arrows(0,vf$nugget,-.5,vf$nugget,col=ccol,length=.1)
    text(-0.2,3.2,label=expression(paste("Nugget (",tau,")")),pos=1,col=ccol)
    lines(vf,col=ccol,lwd=2)
bcol="#8B008B"
bcol=grey(.4)
bt=20
    plot(vg,ylim=c(0, 22),xlim=c(-.25,8),pch=16,col="black",type="o",lwd=2,xlab="Distance (arc-degrees)",cex.lab=1.5,las=1,ylab="")
    legend("topleft",legend=c("Empirical","Bayesian","95% Credible Intervals"),col=c("black",bcol,bcol),lty=c(1,1,20),pch=c(16,NA,NA),bty="n",cex=1.5,lwd=2)
    axis(4,las=1)
    text(4,1.1,label=expression(paste("Range (",phi,")")),pos=4,col=bcol) #phi
    boxplot(kc$posterior$sample$phi,horizontal=T,add=T,at=0,col=paste(bcol,bt,sep=""),boxwex=2,las=1) #phi
    text(7.1,vf$nugget+vf$cov.pars[1]+.5,label=expression(paste("Sill (",tau+sigma^2,")",sep="")),pos=4,col=bcol) #Sill
    boxplot(kc$posterior$sample$tausq.rel+kc$posterior$sample$sigmasq,horizontal=F,add=T,boxwex=0.5,at=8.1,col=paste(bcol,bt,sep=""),las=1) #sill
    text(-0.2,5.2,label=expression(paste("Nugget (",tau,")")),pos=1,col=bcol) #nugget
    boxplot(kc$posterior$sample$tausq.rel*kc$posterior$sample$sigmasq,horizontal=F,add=T,boxwex=0.5,at=-.3,col=paste(bcol,bt,sep=""),las=1) #nugget
    lines(kc, summ = my.summary, ty="l", lty=c(2,1,2), col=bcol,lwd=2) #lines
mtext("Semivariance",2,cex=1.5,adj=1.2,padj=-3)    
mtext("Distance (arc-degrees)",1,cex=1.5,padj=3)    

dev.off()
smallpdf(paste("output/",version,"/kriging.pdf",sep=""))
#}}}

#{{{ plot predictions for example locations
pdf(paste("output/",version,"/exloc.pdf",sep=""),width=11,height=8.5)

### Compare example locations
datelim=c(as.Date(c("2001-01-01","2002-01-01")))
tsubset=expredl$var%in%c("tmax","tmin")&expredl$date>=datelim[1]&expredl$date<=datelim[2]
useOuterStrips(xyplot(mean~date|name+var,data=expredl[tsubset,],type="l",
       panel=function(x,y,subscripts=subscripts){
         td=expredl[tsubset,][subscripts,]
           for(i in 1:1000) panel.xyplot(td$date,rnorm(nrow(td),td$mean,td$sd),col="#00000020",type="l",pch=16,cex=.3)
           #panel.arrows(x0=x,x1=x,y0=y-td$sd,y1=y+td$sd,col="grey",length=0)
           panel.xyplot(td$date,td$mean,type="l")
       },subscripts=T,
       main=paste("Time series of posterior temperature predictions"),
       scales=list(relation="sliced",log=F),layout=c(2,2),
                      ylab="Daily Temperature (C)",xlab="Date (2000)",sub="Grey lines are 1000 posterior time series"))


### Compare example locations
for(v in vars){
tsubset=expredl$var==v&expredl$date>=datelim[1]&expredl$date<=datelim[2]
print(xyplot(mean~date|name,data=expredl[tsubset,],type="l",
       panel=function(x,y,subscripts=subscripts){
         td=expredl[tsubset,][subscripts,]
         if(td$var[1]=="ppt") {
           panel.arrows(x0=x,x1=x,y0=ifelse(y-td$sd<0,0,y-td$sd),y1=y+td$sd,col="grey",length=0)
           panel.xyplot(td$date,td$mean,type="p",pch=16,cex=.3)
         }
           if(td$var[1]!="ppt") {
           panel.arrows(x0=x,x1=x,y0=y-td$sd,y1=y+td$sd,col="grey",length=0)
           panel.xyplot(td$date,td$mean,type="l")
         }
       },subscripts=T,
       main=paste("Comparison of predictions for locations near and far from the nearest station"),
       sub=paste("Grey bars are +/- 1 SD"),
                      layout=c(1,2)))
}
dev.off()
smallpdf(paste("output/",version,"/exloc.pdf",sep=""))
#}}}

#{{{ plot validation data
pdf(paste("output/",version,"/validation.pdf",sep=""),width=11,height=8.5)


### Predicted vs. Observed
xyplot(fine_pred_mean~obs|var,data=valid,asp=1,ylab="Predicted",xlab="Observed",
       panel=function(x,y,subscripts=subscripts){
         td=valid[subscripts,]
         #print(head(td))
         if(td$var[1]=="ppt") panel.arrows(x0=x,x1=x,y0=ifelse(y-td$fine_pred_sd<0,0,y-td$fine_pred_sd),y1=y+td$fine_pred_sd,col="grey",length=0)
         if(td$var[1]!="ppt") panel.arrows(x0=x,x1=x,y0=y-td$fine_pred_sd,y1=y+td$fine_pred_sd,col="grey",length=0)
         panel.xyplot(td$obs,td$fine_pred_mean,pch=16,cex=.3,col="black")
         panel.abline(0,1,col="red",lty="dashed")
         t=accuracy(td$obs,td$fine_pred_mean)
#         panel.key(text = c(
#                     bquote( paste(RMSE== .(t$summary[1]))),
#                     bquote( paste(MAE== .(t$summary[2]))),
#                     bquote( paste(MER==.(t$summary[3]))),
#                     bquote( paste(R^2==.(t$summary[4])))),
         panel.key(text = c(
                     paste("RMSE=",round(t$summary[1],1)),
                     paste("MAE=",round(t$summary[2],1)),
                     paste("MER=",round(t$summary[3],1)),
                     paste("R^2=",round(t$summary[5],1))),
                   lines=F,points=F,corner = c(0,.98))
       },subscripts=T,
       main=paste("Predicted vs. Observed Daily Weather for",nrow(valid)," hold out observations\n(Kriged at 1/4 degree resolution)"),
       sub=paste("Grey bars are +/- 1 SD \n RMSE: Root Mean Squared Errors, MAE: Mean Absolute Error, MER: Mean Error"),
       strip=strip.custom(factor.levels=c("Precipitation (mm)","Maximum Temperature (C)","Minimum Temperature (C)")),
       par.settings = list(strip.background=list(col="transparent"), par.xlab.text=list(cex=1.5),par.ylab.text=list(cex=1.5)),
       scales=list(rot=0,cex=1.5,relation="free",log=F),par.strip.text = list(cex =1.2),layout=c(3,1))

### Predicted vs. Observed precip with log axis for precip
xyplot(I(fine_pred_mean+1)~I(obs+1)|var,data=valid[valid$var=="ppt",],asp=1,ylab="Predicted",xlab="Observed",
       panel=function(x,y,subscripts=subscripts){
         td=valid[valid$var=="ppt",][subscripts,]
         panel.xyplot(x,y,pch=16,cex=.4)
         panel.abline(0,1,col="red")
         t=accuracy(x,y,verbose=F)
         panel.key(text = c(
                     paste("RMSE=",round(t$summary[1],1)),
                     paste("MAE=",round(t$summary[2],1)),
                     paste("MER=",round(t$summary[3],1)),
                     paste("R^2=",round(t$summary[5],1))),
                   lines=F,points=F,corner = c(0,.98))
       },subscripts=T,
       main=paste("Predicted vs. Observed Daily Weather for",nrow(valid)," hold out observations\n(Kriged at 1/4 degree resolution)\nLog Axis"),
       sub=paste("Grey bars are +/- 1 SD \n RMSE: Root Mean Squared Errors, MAE: Mean Absolute Error, MER: Mean Error"),
       scales=list(relation="free",log=T))#,ylim=log(c(1,200)),xlim=log(c(1,200)))


### Predicted vs. Observed by month
useOuterStrips(xyplot(fine_pred_mean~obs|factor(format(valid$date,"%B"),ordered=T,levels=month.name)+var,data=valid,asp=1,ylab="Predicted",xlab="Observed",
       panel=function(x,y,subscripts=subscripts){
         td=valid[subscripts,]
         if(td$var[1]=="ppt") panel.arrows(x0=x,x1=x,y0=ifelse(y-td$fine_pred_sd<0,0,y-td$fine_pred_sd),y1=y+td$fine_pred_sd,col="grey",length=0)
         if(td$var[1]!="ppt") panel.arrows(x0=x,x1=x,y0=y-td$fine_pred_sd,y1=y+td$fine_pred_sd,col="grey",length=0)
         panel.xyplot(td$obs,td$fine_pred_mean,pch=16,cex=.2)
         panel.abline(0,1,col="red",lty="dashed")
       },subscripts=T,
main=paste("Predicted vs. Observed Daily Weather for",nrow(valid)," hold out observations\n(Kriged at 1/4 degree resolution)"),
sub=paste("Grey bars are +/- 1 SD \n RMSE: Root Mean Squared Errors, MAE: Mean Absolute Error, MER: Mean Error"),layout=c(12,3),
                      scales = list(
                        x=list(rot=45,at = rep(list(TRUE, NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL), times=3)),
                        y=list(rot=0,at = rep(list(TRUE, NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL), times=3)),
                        relation = "free", limits =
                        rep(tapply(valid$fine_pred_mean,list(valid$var),
                                   function(x) round(range(x,na.rm=T))),each=12)),
                      par.settings = list(plot.symbol = list(pch = 16,cex=.3,col="black"),
                        strip.background=list(col="transparent"),
                        layout.widths = list(ylab=10,left.padding=1,axis.panel = rep(c(1, 0), c(0, 1)))),
                      par.strip.text = list(cex =.75)),
               strip.left=strip.custom(par.strip.text=list(cex=.5),factor.levels=c("Precipitation","Max Temperature","Min Temperature")))

    
## Krig Parameters over time
for(parm in unique(skparms$parameter)){
  ylim=range(c(skparms[skparms$parameter==parm,"0%"],skparms[skparms$parameter==parm,"100%"]))
  print(xyplot(skparms[skparms$parameter==parm,"50%"]~date|var,data=skparms[skparms$parameter==parm,],
         panel=function(x,y,subscripts){
           td=skparms[skparms$parameter==parm,][subscripts,]
           panel.arrows(td$date,td[,"0%"],td$date,td[,"100%"],col=grey(.9),length=0)
           panel.arrows(td$date,td[,"2.5%"],td$date,td[,"97.5%"],col=grey(.7),length=0)
           panel.xyplot(td$date,td[,"50%"],pch=16,cex=.2)
         },subscripts=T,scales=list(relation="free"),layout=c(1,3),ylab="Value",xlab="Date",main=parm,ylim=ylim))
}  

### Summarize temporally
valid$pred_error=valid$fine_pred_mean-valid$obs
bwplot(pred_error~format(date,"%m")|var,data=valid,scales=list(relation="free"),layout=c(3,1),main="Residuals by month",xlab="Month",ylab="Residuals")
bwplot(skparms[,"50%"]~format(date,"%m")|parameter+var,data=skparms,layout=c(6,3),
   scales=list(relation="free"),main="Parameters by month",xlab="Month",ylab="Median posterior parameter value")

svalid=aggregate(valid$pred_error,
  list(month=format(valid$date,"%m"),lat=valid$lat,lon=valid$lon,var=valid$var),mean,na.rm=T)
#ccut=as.character(cut(svalid$x,seq(-5,35,length=11),lab=rainbow(10)))
levelplot(x~lon*lat|month+var,data=svalid,contour=T,col.regions=heat.colors(100))

## scatterplot matrix of mean posterior parameter values
for(vari in unique(skparms$var))
print(splom(cast(date~parameter,value="50%",data=skparms[skparms$var==vari,])[,-1],main=vari))

dev.off()
smallpdf(paste("output/",version,"/validation.pdf",sep=""))

                                        #shrinkpdf(paste("output/",version,"/validation.pdf",sep=""),maxsize=1)
#}}}

#{{{ plot mean and sd as a map

mets=sub("[.]nc*","",list.files(paste("output/",version,"/metrics_summary_coarse",sep=""),pattern=".*nc$"))

iyear=2000

getpred=function(metric,year){
  nc=nc_open(paste("output/",version,"/metrics_summary_coarse/",metric,".nc",sep=""))
  d=ncvar_get(nc,NA,start=c(1,1,which(year==1990:2009),1),count=c(-1,-1,1,-1))
  d[d==-999]=NA
  ds=apply(d,1:2,function(i) c(Mean=mean(i,na.rm=T),SD=sd(i,na.rm=T)))
  dimnames(ds)=list(type=c("Mean","SD"),lon=ncvar_get(nc,"lon"),lat=ncvar_get(nc,"lat"))
  dl=as.data.frame(cast(melt(ds),lon+lat~type,value="value"))
  dl2=melt(dl,measure.vars=c("Mean","SD"))
  dl2$metric=metric
  dl2$year=year
  dl2=dl2[!is.na(dl2$value),]
  print(paste("Finished processing",metric,"for year",year))
  return(dl2)
}

smet=do.call(rbind.data.frame,mclapply(mets,getpred,year=2000))
save(smet,file=paste("output/",version,"/metrics_summary/smet.Rdata",sep=""))

load(paste("output/",version,"/metrics_summary/smet.Rdata",sep=""))
smet$variable=as.character(smet$variable)
smet$variable[smet$variable=="mean"]="Mean"
smet$variable[smet$variable=="sd"]="SD"
smet$metric[smet$metric=="ECAr20mm"]="ECA"

## add coefficient of variation
smet2=as.data.frame(cast(smet,lon+lat+metric+year~variable,value="value"))
smet2$CV=smet2$SD/smet2$Mean
smet2$CV[is.nan(smet2$CV)]=0
smet2l=melt(smet2,measure.vars=c("Mean","SD","CV"))
smet2l2=smet2l[smet2l$variable!="SD",]
smet2l2=smet2l2[smet2l2$metric%in%c("CDD","CFD","CSU","ECA"),]

### summary stats for paper
tsmet=smet[smet$variable=="CV"&smet$year==2000&smet$metric=="CDD",]
tsmet=cast(smet[smet$year==2000&smet$metric=="CDD",],lon+lat~variable,value="value")
tsmet$CV=tsmet$SD/tsmet$Mean
hist(tsmet$CV)
aggregate(tsmet$Mean,list(bin=cut(tsmet$lat,c(-36,-34,-30),c("South","North"))),
       function(x) c(mean=mean(x,na.rm=T),sd=sd(x,na.rm=T),CV=sd(x,na.rm=T)/mean(x,na.rm=T),quantile(x,0.025,na.rm=T),quantile(x,0.975,na.rm=T),IQW=quantile(x,0.975,na.rm=T)-quantile(x,0.025,na.rm=T)))


g.col=colorRampPalette(c("blue","green","yellow","red","red","darkred","darkred"));ncol=100
#g.col=colorRampPalette(c(grey(.9),grey(.05)));ncol=100
coast= readOGR("/media/Data/Work/Regional/CFR/BaseGISData/CFR_coastline/coastline.shp","coastline")
png(paste("output/",version,"/Metrics.png",sep=""),width=1780,height=1080)
plot(1:10)  #random plot workaround because of new=T bug in panel function
useOuterStrips(levelplot(value~lon*lat|metric+variable,data=smet[smet$metric%in%c("CDD","CFD","CSU","ECA"),],panel=function(x,y,z,subscripts){
#  at=unique(quantile(z[subscripts],seq(0,1,len=ncol),na.rm=T))
  at=seq(range(z[subscripts])[1],range(z[subscripts])[2],len=ncol) # if((max(at)-min(at))<.5)
  cols=level.colors(z[subscripts], at, col.regions=g.col(length(at)))
  panel.levelplot(x,y,z,subscripts=subscripts,col.regions=g.col(length(at)),at=at)
#  panel.xyplot(coordinates(st[st$ppt==1,])[,1],coordinates(st[st$ppt==1,])[,2],pch=4,col="black")
  ## Draw histogram
  pushViewport(viewport(x=.69,y=.75,width=.55,height=.5))
  first=current.row()==1&current.column()==1
  par(fig=gridFIG(),mar=c(2,.5,.5,.5),new=T)
  par(fig=gridFIG(),mar=c(2,.5,.5,.5))
  pat=pretty(z[subscripts],3)
  hist(z[subscripts],breaks=at,col=g.col(length(at)),border=NA,yaxt="n",ylab="",xlab="",cex.axis=1,main="",new=T,axes=F)
  axis(1,at=pat,cex.axis=2)
  popViewport()
},subscripts=T,colorkey=F,scales=list(draw=T),xlim=c(17.8,26),ylim=c(-34.85,-30.8),
       par.settings=list(strip.background=list(col="transparent"),axis.text=list(cex=2),
       par.xlab.text=list(cex=2),par.ylab.text=list(cex=2)),par.strip.text=list(cex=2),as.table=T,asp=1,ylab="Latitude",xlab="Longitude"))
dev.off()

### version with contours 
pdf(paste("output/",version,"/CDD.pdf",sep=""),width=11,height=8.5)
plot(1:10)  #random plot workaround because of new=T bug in panel function
#g.col=colorRampPalette(c(grey(.95),grey(.3)));ncol=100
#g.col=colorRampPalette(c(rgb(153,213,148,m=255),rgb(255,255,191,m=255),rgb(252,141,89,m=255)));ncol=100
g.col=colorRampPalette(c(rgb(37,132,182,m=255),rgb(168,221,181,m=255),rgb(204,223,199,m=255)));ncol=100
levelplot(Mean~lon*lat,data=smet2[smet2$metric=="CDD",],panel=function(x,y,z,subscripts){
  td=smet2[smet2$metric=="CDD",][subscripts,]
  at=unique(quantile(z[subscripts],seq(0,1,len=ncol),na.rm=T))
  ##at=seq(range(z[subscripts])[1],range(z[subscripts])[2],len=ncol) # if((max(at)-min(at))<.5)
  cols=level.colors(z[subscripts], at, col.regions=g.col(length(at)))
  panel.levelplot(td$lon,td$lat,td$Mean,subscripts=1:nrow(td),col.regions=g.col(length(at)),at=at)
#  cat=round(unique(quantile(td$CV,seq(0,1,len=10),na.rm=T)),3)
  cat2=pretty(td$CV,10)#round(seq(range(td$CV)[1],range(td$CV)[2],len=5),3) # if((max(at)-min(at))<.5)
  cat=pretty(td$CV,3)
  catl=cat
#  catl[seq(1:length(cat),by=2)]=NA
  panel.contourplot(td$lon,td$lat,td$CV,subscripts=1:nrow(td),region=F,contour=T,at=cat2,labels=F,lwd=.5)
  panel.contourplot(td$lon,td$lat,td$CV,subscripts=1:nrow(td),region=F,contour=T,at=cat,labels=list(labels=catl,cex=.8,col="darkblue"),label.style="flat",shrink=c(10,10))
  panel.xyplot(coordinates(st)[st$dcfr==0&grepl("Precipitation",st$type),1],coordinates(st)[st$dcfr==0&grepl("Precipitation",st$type),],cex=.75,pch=3,col="black",lwd=2)
  sp.lines(coast)
  ## Draw histogram
  pushViewport(viewport(x=.6,y=.75,width=.6,height=.45))
  first=current.row()==1&current.column()==1
  par(fig=gridFIG(),mar=c(2,.5,.5,.5),new=T)
  par(fig=gridFIG(),mar=c(5,5,1,1))
  pat=pretty(z[subscripts],3)
  plot(td$Mean,td$CV,xlab=expression(paste("Mean",(mu))),ylab=expression(paste("CV",(sigma/mu))),
       col=as.character(cut(td$Mean,breaks=c(0,at,Inf),labels=g.col(length(at)+1))),
       yaxp=c(0,.6,3),pch=16,cex=.75,new=F,bty="n",las=1,cex.lab=1.5,cex.axis=1.5)
#  hist(z[subscripts],breaks=at,col=g.col(length(at)),border=NA,yaxt="n",ylab="",xlab="",cex.axis=1,main="",new=T,axes=F)
#  axis(1,at=pat,cex.axis=2)
  popViewport()
},subscripts=T,colorkey=F,scales=list(draw=T),xlim=c(17.8,26),ylim=c(-34.85,-30.8),
       par.settings=list(strip.background=list(col="transparent"),axis.text=list(cex=2),
       par.xlab.text=list(cex=2),par.ylab.text=list(cex=2)),par.strip.text=list(cex=2),as.table=T,asp=.7,ylab="Latitude",xlab="Longitude")
dev.off()

### version with mean and CV
#g.col=colorRampPalette(c(grey(.9),grey(.05)));ncol=100
#png(paste("output/",version,"/Metrics3.png",sep=""),width=1780,height=1080)
ncol=25
pdf(paste("output/",version,"/Metrics3.pdf",sep=""),width=11,height=8.5)
plot(1:10)  #random plot workaround because of new=T bug in panel function
useOuterStrips(levelplot(value~lon*lat|metric+variable,data=smet2l2[smet2l2$variable%in%c("Mean","CV"),],panel=function(x,y,z,subscripts){
  #at=unique(quantile(z[subscripts],seq(0,1,len=ncol),na.rm=T))
  at=seq(range(z[subscripts])[1],range(z[subscripts])[2],len=ncol) # if((max(at)-min(at))<.5)
  cols=level.colors(z[subscripts], at, col.regions=g.col(length(at)))
  panel.levelplot(x,y,z,subscripts=subscripts,col.regions=g.col(length(at)),at=at)
  sp.lines(coast)
  ## Draw histogram
  pushViewport(viewport(x=.69,y=.75,width=.55,height=.5))
  first=current.row()==1&current.column()==1
  par(fig=gridFIG(),mar=c(2,.5,.5,.5),new=T)
  par(fig=gridFIG(),mar=c(2,.5,.5,.5))
#  at2=at[-c(1,(length(at)-10):(length(at)))]
#  at2=at[-c(1,length(at))]
  at2=at
  pat=pretty(z[subscripts][z[subscripts]>=min(at2)&z[subscripts]<=max(at2)],4)
  hist(z[subscripts][z[subscripts]>=min(at2)&z[subscripts]<=max(at2)],breaks=at2,col=g.col(length(at2)),border=NA,yaxt="n",ylab="",xlab="",cex.axis=1,main="",new=T,axes=F)
  axis(1,at=pat,cex.axis=1.2)
  popViewport()
},subscripts=T,colorkey=F,scales=list(draw=T),xlim=c(17.8,26),ylim=c(-34.85,-30.8),
       par.settings=list(strip.background=list(col="transparent"),axis.text=list(cex=2),
       par.xlab.text=list(cex=1.5),par.ylab.text=list(cex=1.5)),par.strip.text=list(cex=2),as.table=T,asp=1,ylab="Latitude",xlab="Longitude"))
dev.off()

## scatterplot matrix
splom(cast(smet[smet$variable=="Mean",],lon+lat~metric)[,c("MATmax","MATmin","MAP","CDD","CFD","CSU","ECA")],pch=16,col="black",cex=.25)

smetw=cast(smet[smet$variable=="Mean",],lon+lat+year~metric)
splom(~smetw[c("MATmax","MATmin","MAP","CDD","CFD","CSU","ECA")],groups=year,smetw,pch=16,cex=.25)
splom(~smetw[4:11],groups=year,smetw,pch=16,cex=.25,panel=panel.superpose)


## With SD Hillshade
dr=brick(dl)
projection(dr)="+proj=lonlat +ellps=WGS84"
fg=dr;res(fg)=0.008
fsd=resample(subset(dr,"sd"),fg,method="bilinear")
hs=hillShade(terrain(fsd,opt="slope"),terrain(fsd,opt="aspect"))
pptcol<-colorRampPalette(c("blue","green","yellow","red","darkred","purple"))

#pdf(paste("output/",version,"/CDD.pdf",sep=""),width=11,height=8.5)
png(paste("output/",version,"/CDD.png",sep=""),width=1280,height=1080)
par(fig=c(0,1,0,1))
brk=quantile(subset(dr,"mean"),seq(0,1,len=ncol+1))
plot(hs, col=grey(0:100/100), legend=FALSE,las=1,cex.lab=1.5,xlab="Longitude",ylab="Latitude",cex.axis=1.5)
plot(subset(dr,"mean"), breaks=brk,lab.breaks=ncol,col=adjustcolor(pptcol(100),alpha.f=0.5), add=TRUE,legend=F)
points(st[st$dcfr==0&grepl("Precipitation",st$type),],pch=16)
par(fig=c(.3,.9,.5,.9), new=TRUE)#mar=c(1,3,1,1),
plot(dl$mean,dl$sd/dl$mean,xlab="Mean Value",ylab="Coefficient of Variation",
     col=adjustcolor(as.character(cut(dl$mean,breaks=brk,labels=pptcol(ncol))),alpha.f=0.5),
     yaxp=c(0,.6,3),pch=16,cex=.75,new=F,bty="n",las=1,cex.lab=1.5,cex.axis=1.5)
dev.off()

smallpdf(paste("output/",version,"/CDD.pdf",sep=""))
         
## Densities
par(fig=c(.05,1,.05,1))
brk=quantile(subset(dr,"mean"),seq(0,1,len=ncol+1))
plot(subset(dr,"mean"),breaks=brk,col=pptcol(ncol),lab.breaks=ncol,legend=F,las=1)
par(fig=c(.3,.9,.6,.8), new=TRUE)#mar=c(1,3,1,1),
plot(dl$mean,dl$sd/dl$mean,xlab="Mean Value",ylab="Standard Deviation",col=as.character(cut(dl$mean,breaks=brk,labels=pptcol(ncol))),pch=16,cex=.5,new=F,bty="n",las=1)



### make color bar
#x=dl$mean
#y=dl$sd/dl$mean
#ncol=500
#    at=seq(min(x,na.rm=T),max(x,na.rm=T),len=ncol)
#    col=data.frame(t(col2rgb(as.character(cut(x,breaks=ncol,labels=pptcol(ncol))))))
#    col=as.data.frame(t(rgb2hsv(r=col$red,g=col$green,b=col$blue,maxColorValue=255)))
#   col$s=as.numeric(as.character(cut(y,breaks=quantile(y,seq(0,1,len=ncol),na.rm=T),labels=seq(1,.3,len=ncol-1)))) #saturation
#   col$v=as.numeric(as.character(cut(y,breaks=quantile(y,seq(0,1,len=ncol),na.rm=T),labels=seq(1,0.4,len=ncol-1))))  #value
#   col[is.na(col)]=0
#   col.val=hsv(h=col$h,s=col$s,v=col$v)
#plot(dl$mean,dl$sd,xlab="Mean Value",ylab="Standard Deviation",col=col.val,pch=16,cex=1,new=F,bty="n")
#plot(dl$mean,dl$sd/dl$mean,xlab="",ylab="",col=col.val,pch=16,cex=1,new=F,bty="n")
#fcex=apply(as.matrix((dl$mean/dl$sd)*.1),1,function(i) min(1,i,na.rm=T)); hist(fcex)
#plot(coordinates(dl)[,1],coordinates(dl)[,2],xlab="",ylab="",col=col.val,pch=16,cex=fcex,new=T)


### Plot the map
#spplot(dl,zcol="mean",colours = col.val, panel = function(x, y, colours) {
#  sp.grid(dl,col = colours,pch=16)
#  panel.xyplot(x,y,col =adjustcolor("white",alpha.f=.8),pch=16,cex=fcex)
                                      #sp.lines(as(coast,"SpatialLines"),lty=1,col="black")
#},colorkey=F,scales=list(draw=T),par.settings=list(axis.text=list(cex=1.8)),xlim=c(17.8,26),ylim=c(-34.85,-30.8))


#}}}

#####  Make animation of weather data
tmax_mean=brick(paste("output/",version,"/weather_fine.nc",sep=""),varname="tmax_mean")
tmax_sd=brick(paste("output/",version,"/weather_fine.nc",sep=""),varname="tmax_sd")

spplot(as(subset(tmax_mean,1),"SpatialGridDataFrame"))
