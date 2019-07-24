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
path.data<-paste(path,"/data", sep="")
path.dust<-paste(path,"/dust analysis", sep="")
path.lock<-paste(path,"/products analysis/lockmass", sep="")
path.prodref<-paste(path, "/products analysis/reference compounds analysis", sep="")
path.out<-paste(path,"/refCal", sep="")
path.db<-paste(path,"/SMILES_DATABASE", sep="")
path.lockdata<-paste(path.data,"/20190426", sep="")

###Load function files and set universal variables.###
###If you make changes to function files, you will need to run these 2 lines again###
setwd(path)
source("Nontargeted_fun.r")
polarity<--1##if neg -1, if pos 1
#LockMass<-c(131.9614,165.0188,172.956991,227.201104,255.232405,312.972263,411.12912)##Lock Mass, C8H5O4, C14H27O2, C16H31O2, humic acid
electron<-0.0005485799
ppm<-2
ppm.ms2<-3##ms2 spectra accuracy

###Pre-load mzXML raw data into a list. This saves time when running other functions which require multiple xcmsRaw objects.###
xrawdata<-NULL
xrawdata<-list()
setwd(path.lockdata)
msfiles<-list.files(pattern="\\.mzXML$", ignore.case=TRUE)
msfiles<-msfiles[96]
for (i in 1:length(msfiles)){
  xdata<-xcmsRaw(msfiles[i],includeMSn=TRUE)
  xrawdata[i]<-xdata
}
#for Jeans samples, NEG=xrawdata[1:9], POS=xrawdata[10:20]
######### Body #########

######Lockmass Code######

###Identify lockmasses in raw data and create a csv file###
setwd(path.lockdata)
xset<-xcmsSet(msfiles,method='centWave',ppm=2.5,peakwidth=c(5,20),snthresh=10,polarity="negative")##data files,peak detection algorithm,mass error,peak width (min and max width in seconds),signal/noise threshold,polarity) 
result<-findlock(xrawdata[10:20],xset,2000,0.001)##list of xcmsRaw objects, xcmsSet object, intensity threshold, mzstep. See Nontargeted_fun.R for "findlock" code.
setwd(path.jeans)
write.table(result, file="lockmassjeansPOS.csv", sep = ',',row.names=FALSE,col.names=c("mz","minintensity","sampleID"))
#Reference Note: This function found 802 masses when run on the 100_300_Neg sample set (no blanks)

###Plot lockmasses###
setwd(path.lockdata)
msfiles<-list.files()
msfiles<-msfiles[1:96]
testMass<-255.2324
polarity<--1##if neg -1, if pos 1
LockMass.POS<-c(81.070425,93.070425,105.070425,139.11229,151.042199,171.138505,413.26696,445.120583)##Lock Mass for positive ion mode
setwd(path.prod)
LockMass.NEG<-read.table("lockmassProd.csv",header=TRUE,sep=',')
LockMass.NEG<-LockMass.NEG$Lock
if (polarity==1){
  LockMass<-LockMass.POS
}else{
  LockMass<-LockMass.NEG
}
setwd(path.lockdata)
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


####################Mass Calibration###############
setwd(path.lock)
LockMass.NEG<-read.table("lockmassProd.csv",header=TRUE,sep=',')
LockMass.NEG<-LockMass.NEG$Lock
path.rawsource<-paste(path, "/data/20190426", sep="")
path.caldata<-paste(path,"/caldata/20190426 products",sep="")

setwd(path.rawsource)
for (i in 1:length(msfiles)){
  #xrawdata<-xcmsRaw(msfiles[i],includeMSn=TRUE)
  if (polarity==1){
    xrawcaldata<-MassCal(xrawdata[[i]],LockMass.POS,'POS',3)}###lockshift is used to preliminarly adjust searching space
  if (polarity==-1){
    xrawcaldata<-MassCal(xrawdata[[i]],LockMass.NEG,'NEG',3)}
  setwd(path.caldata)
  name<-msfiles[i]
  write.mzdata(xrawcaldata,name)####save the calibrated data to new files, msndata cannot be stored
  print(c(i,"out of",length(msfiles)))
}


##############Data extraction for quick peak searching##############
path.rawsource<-paste(path, "/data/20190426", sep="")
path.calsource<-paste(path, "/caldata/20190426 products/caldata_July18_singlepluslinear", sep="") #Change source folder as needed
path.destination<-paste(path, "/products analysis/reference compounds analysis/uncalibrated peaks", sep="") #Change destination folder as needed

setwd(path.calsource)#Change data set as required
msfiles<-list.files()
msfiles<-msfiles[6] #Change file as needed (should be a single file)
#xset<-xcmsSet(msfiles,method='centWave',ppm=2.5,peakwidth=c(5,20),snthresh=3)##data files,peak detection algorithm,mass error,peak width (min and max width in seconds),signal/noise threshold,polarity) 
xraw<-xcmsRaw(msfiles)
#setwd(path.destination)
#write.table(xset@peaks[,c(1,4,9)], file="raw xcmsRaw sample100test.csv", sep = ',',row.names=FALSE) #Change file name as needed
RawFrame<-data.frame(mz=xraw@env$mz,intensity=xraw@env$intensity)
#write.table(RawFrame, file="raw xcmsRaw sample99.csv", sep = ',',row.names=FALSE) #Change file name as needed

#####Scan through extracted xcmsRaw data frame to pick out best match based on highest intensity#####
RawFrame<-RawFrame[order(RawFrame$mz),]

setwd(path.prodref)
refmz<-read.table("CPref.csv",header=TRUE,sep=',')
refmz<-refmz$mz_e
refmz<-sort(refmz)
minus3ppm<-refmz-(refmz*0.000003)
plus3ppm<-refmz+(refmz*0.000003)
refmatch<-data.frame(mzref=NA,mzmatch=NA,mzshift=NA)
for (i in 1:length(refmz)){
  indexfind<-which(RawFrame$mz>minus3ppm[i] & RawFrame$mz<plus3ppm[i])
  if(length(indexfind)==0){next}
  ppmarray<-((RawFrame$mz[indexfind]-refmz[i])/refmz[i])*10^6
  lowshift<-which.min(abs(ppmarray))
  mzmatch<-RawFrame$mz[indexfind[lowshift]]
  refmatch[i,]<-rbind(c(refmz[i],mzmatch,ppmarray[lowshift]))
}
nomatch<-which(is.na(refmatch$mzref))
refmatch<-refmatch[-nomatch,]
write.table(refmatch, file="cal CP reference matches sample99_singlepluslinear.csv", sep = ',',row.names=FALSE)


