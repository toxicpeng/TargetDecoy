######### Head #########

###Reset all variables and clear all objects from R memory###
rm(list=ls())

###Load packages for code to run. Only need to run at the start of your R session.###
library(xcms)
library(MassSpecWavelet)
library(Rcpp)
library(RcppArmadillo)
library(isopat)
library(gtools)
library(caTools)
data(iso_list)

###Set file pathways. If you are not currently in your root directory or have differently named folders, you will get an error###
###You can use setwd() to change the current directory to the root.###
###If your folder structure is different, you can change the folder names here.###
rootdir<-getwd()
path<-rootdir
path.source<-paste(path,"/data/20190426", sep="")
path.destination<-paste(path,"/analysis/products analysis/lockmass", sep="")
path.lock<-paste(path,"/analysis/products analysis/lockmass", sep="")

###Load function files and set universal variables.###
###If you make changes to function files, you will need to run these 2 lines again###
setwd(path)
source("Nontargeted_fun.r")
#LockMass<-c(131.9614,165.0188,172.956991,227.201104,255.232405,312.972263,411.12912)##Lock Mass, C8H5O4, C14H27O2, C16H31O2, humic acid


###Pre-load mzXML raw data into a list. This saves time when running other functions which require multiple xcmsRaw objects.###
xrawdata<-NULL
xrawdata<-list()
setwd(path.source)
msfiles<-list.files(pattern="\\.mzXML$", ignore.case=TRUE)
msfiles<-msfiles[145:147]
for (i in 1:length(msfiles)){
  xdata<-xcmsRaw(msfiles[i],includeMSn=TRUE)
  xrawdata[i]<-xdata
}
#for Jeans samples, NEG=xrawdata[1:9], POS=xrawdata[10:20]
######### Body #########

polarity<-1##if neg -1, if pos 1

######Lockmass Detection######

###Identify lockmasses in raw data and create a csv file###
setwd(path.source)
#xset<-xcmsSet(msfiles,method='centWave',ppm=2.5,peakwidth=c(5,20),snthresh=10,polarity="positive")##data files,peak detection algorithm,mass error,peak width (min and max width in seconds),signal/noise threshold,polarity) 
xset<-xcmsSet(msfiles,method='centWave',ppm=2.5,peakwidth=c(5,20),snthresh=10,polarity="positive")
result<-findlock(xrawdata,xset,2000,0.001)##list of xcmsRaw objects, xcmsSet object, intensity threshold, mzstep. See Nontargeted_fun.R for "findlock" code.
setwd(path.destination)
write.table(result, file="lockmassproductsPOS_sample43to45.csv", sep = ',',row.names=FALSE,col.names=c("mz","minintensity","sampleID"))
#Reference Note: This function found 802 masses when run on the 100_300_Neg sample set (no blanks)

######Extra Functions for Confirmation#####

###Plot lockmasses###
setwd(path.source)
msfiles<-list.files()
msfiles<-msfiles[3:10]
testMass<-255.2324
polarity<--1##if neg -1, if pos 1
LockMass.POS<-c(81.070425,93.070425,105.070425,139.11229,151.042199,171.138505,413.26696,445.120583)##Lock Mass for positive ion mode
setwd(path.lock)
LockMass.NEG<-read.table("lockmassProd.csv",header=TRUE,sep=',')
LockMass.NEG<-LockMass.NEG$Lock
if (polarity==1){
  LockMass<-LockMass.POS
}else{
  LockMass<-LockMass.NEG
}
setwd(path.source)
xrawdata<-xcmsRaw(msfiles[1])
ppmshift<-plotlock(xrawdata,LockMass,10)

##############initial curve fitting using overall lockmass############
ppmshift<-plotlock(xrawdata,LockMass,10)
lock.shift<-data.frame(lockmass=LockMass,shift=rep(0,length(LockMass)))
index.save<-NULL
for (i in 1:length(LockMass)){
  temp<-ppmshift[i,]
  index<-which(temp<15)
  if (length(index)<1){
    index.save<-c(index.save,i)
    next}
  lock.shift$shift[i]<-median(temp[index])
}
if (length(index.save)>0){
  lock.shift<-lock.shift[-index.save,]}##delete those lockmass not detected
lock.shift$shift<-lock.shift$shift*10^(-6)
LockMass.cal<-fitlock(lock.shift,lock.shift$lockmass,NULL)
LockMass.cal<-LockMass.cal[[1]]
lock.shift<-(lock.shift$lockmass-LockMass.cal)*10^6/LockMass.cal
plot(LockMass.cal,lock.shift)

