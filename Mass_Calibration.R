
####################for drinking water, cpp funtion, the oxygen number could be 1.2*carbon+3, cutoff e5, S/N>5, change is.BrCl
library(xcms)
library(MassSpecWavelet)
library(Rcpp)
library(RcppArmadillo)
library(isopat)
library(gtools)
library(caTools)
data(iso_list)
rootdir<-getwd()##If you opened the current file via the R project, the following directory assignments will be correct
path<-rootdir
path.data<-paste(path,"/data", sep="")
path.dust<-paste(path,"/dust analysis", sep="")
path.prod<-paste(path,"/products analysis", sep="")
path.out<-paste(path,"/refCal", sep="")
path.db<-paste(path,"/SMILES_DATABASE", sep="")

source("Nontargeted_fun.r")
polarity<--1##if neg -1, if pos 1
LockMass<-c(131.9614,165.0188,172.956991,227.201104,255.232405,312.972263,411.12912)##Lock Mass, C8H5O4, C14H27O2, C16H31O2, humic acid
electron<-0.0005485799
ppm<-2
ppm.ms2<-3##ms2 spectra accuracy

path.lockdata<-paste(path.data,"/20190426", sep="")
xrawdata<-NULL
xrawdata<-list()
setwd(path.lockdata)
msfiles<-list.files()
msfiles<-msfiles[199:294]
for (i in 1:length(msfiles)){
  xdata<-xcmsRaw(msfiles[i],includeMSn=TRUE)
  xrawdata[i]<-xdata
}

##########identify lockmass#############
mzwin<-5###2.5ppm for mz cutoff
timewin<-0.5###30 sec for rt cutoff
timewin2<-2####60 sec for rt cutoff, since library was established for a long time
xset<-xcmsSet(msfiles,method='centWave',ppm=2.5,peakwidth=c(5,20),snthresh=10,nSlaves=1,polarity="negative")##peak width, the min and max range of chromatographic peaks in seconds
result<-findlock(xrawdata,xset,2000,0.002)##xset, intensity threshold, mzstep
setwd(path.prod)
write.table(result, file="lockmassProd2000_neg300_500_noblanks.csv", sep = ',',row.names=FALSE,col.names=c("mz","minintensity","sampleID"))
##Reference Note: This function found 802 masses when run on the 100_300_Neg sample set (no blanks)

#################plot LockMass######################
setwd(path.lockdata)
msfiles<-list.files()
msfiles<-msfiles[1:102]
testMass<-255.2324
polarity<--1##if neg -1, if pos 1
LockMass.POS<-c(81.070425,93.070425,105.070425,139.11229,151.042199,171.138505,413.26696,445.120583)##Lock Mass for positive
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
  lock.shift$shift[i]<-mean(temp[index])
}
lock.shift<-lock.shift[-index.save,]##delete those lockmass not detected
lock.shift$shift<-lock.shift$shift*10^(-6)
LockMass.cal<-fitlock(lock.shift,lock.shift$lockmass,NULL)
LockMass.cal<-LockMass.cal[[1]]
lock.shift<-(lock.shift$lockmass-LockMass.cal)*10^6/LockMass.cal


####################Mass Calibration###############


setwd(path.prod)
LockMass.NEG<-read.table("lockmassProd.csv",header=TRUE,sep=',')
LockMass.NEG<-LockMass.NEG$Lock
path.caldata<-paste(path.data,"/testCal",sep="")

setwd(path.lockdata)
for (i in 1:length(msfiles)){
  setwd(path.lockdata)
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


