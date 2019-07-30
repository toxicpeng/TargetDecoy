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
path.proddata<-paste(path.data,"/20190426", sep="")
path.jeansdata<-paste(path.data,"/20190614",sep="")

###Load function files and set universal variables.###
###If you make changes to function files, you will need to run these 2 lines again###
setwd(path)
source("Nontargeted_fun.r")
polarity<--1##if neg -1, if pos 1
#LockMass<-c(131.9614,165.0188,172.956991,227.201104,255.232405,312.972263,411.12912)##Lock Mass, C8H5O4, C14H27O2, C16H31O2, humic acid


###Pre-load mzXML raw data into a list. This saves time when running other functions which require multiple xcmsRaw objects.###
xrawdata<-NULL
xrawdata<-list()
setwd(path.jeansdata)
msfiles<-list.files(pattern="\\.mzXML$", ignore.case=TRUE)
msfiles<-msfiles[96]
for (i in 1:length(msfiles)){
  xdata<-xcmsRaw(msfiles[i],includeMSn=TRUE)
  xrawdata[i]<-xdata
}
#for Jeans samples, NEG=xrawdata[1:9], POS=xrawdata[10:20]
######### Body #########






##############Data extraction for quick peak searching##############
path.rawsource<-paste(path, "/data/20190614", sep="")
path.calsource<-paste(path, "/caldata/20190614 jeans", sep="") #Change source folder as needed
path.destination<-paste(path, "/jeans analysis/reference compounds", sep="") #Change destination folder as needed

datatype<-"CAL" #Change data set as required

if(datatype=="RAW"){setwd(path.rawsource)}
if(datatype=="CAL"){setwd(path.calsource)}
msfiles<-list.files()
msfiles<-msfiles[1] #Change file as needed (should be a single file)
xset<-xcmsSet(msfiles,method='centWave',ppm=2.5,peakwidth=c(5,20),snthresh=3)##data files,peak detection algorithm,mass error,peak width (min and max width in seconds),signal/noise threshold,polarity) 
xraw<-xcmsRaw(msfiles)
#setwd(path.destination)
#write.table(xset@peaks[,c(1,4,9)], file="raw xcmsRaw sample100test.csv", sep = ',',row.names=FALSE) #Change file name as needed
RawFrame<-data.frame(mz=xraw@env$mz,intensity=xraw@env$intensity)
#write.table(RawFrame, file="cal xcmsRaw July22 singlepluslinearto500 sample100.csv", sep = ',',row.names=FALSE) #Change file name as needed

#####Scan through extracted xcmsRaw data frame to pick out best match based on highest intensity#####
RawFrame<-RawFrame[order(RawFrame$mz),]

if(datatype=="RAW"){
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
  setwd(path.destination)
  write.table(refmatch, file="raw CP reference matches jeanssampleA.csv", sep = ',',row.names=FALSE)
}

if (datatype=="CAL"){
  setwd(path.prodref)
  refmz<-read.table("CPref.csv",header=TRUE,sep=',')
  refmz<-refmz$mz_e
  refmz<-sort(refmz)
  refmatch<-data.frame(matrix(data=NA,nrow=121,ncol=3))
  for (i in 1:length(refmz)){
    if(refmz[i]<500){
      linfitppm<-(refmz[i]*6.777697e-09)+(-3.3888485e-06)
      shiftmass<-refmz[i]*(1-linfitppm)
      minus3ppm<-shiftmass-(shiftmass*0.000003)
      plus3ppm<-shiftmass+(shiftmass*0.000003)
      indexfind<-which(RawFrame$mz>minus3ppm & RawFrame$mz<plus3ppm)
      if(length(indexfind)==0){next}
      ppmarray<-((RawFrame$mz[indexfind]-shiftmass)/shiftmass)*10^6
    }else{
      minus3ppm<-refmz[i]-(refmz[i]*0.000003)
      plus3ppm<-refmz[i]+(refmz[i]*0.000003)
      indexfind<-which(RawFrame$mz>minus3ppm & RawFrame$mz<plus3ppm)
      if(length(indexfind)==0){next}
      ppmarray<-((RawFrame$mz[indexfind]-refmz[i])/refmz[i])*10^6
    }
    if(length(ppmarray)==0){next}
    lowshift<-which.min(abs(ppmarray))
    mzmatch<-RawFrame$mz[indexfind[lowshift]]
    refmatch[i,]<-rbind(c(refmz[i],mzmatch,ppmarray[lowshift]))
  }
  names(refmatch)[1:3]<-c("mzref","mzmatch","ppm")
  nomatch<-which(is.na(refmatch$mzref))
  refmatch<-refmatch[-nomatch,]
  setwd(path.destination)
  write.table(refmatch, file="cal CP reference matches jeanssampleA.csv", sep = ',',row.names=FALSE)
}



