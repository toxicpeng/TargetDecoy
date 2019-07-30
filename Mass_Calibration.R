###Step 0: Reset R Environment-----------------------------------------------------------------------------------------

###Reset all variables and clear all objects from R memory###
rm(list=ls())

###Step 1: Load packages and set root directory (only run at start of R session)---------------------------------------

###Load packages for code to run###
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
path.lock<-paste(path,"/analysis/products analysis/lockmass", sep="")
###Load function files and set universal variables.###
###If you make changes to function files, you will need to run these 2 lines again###
setwd(path)
source("Nontargeted_fun.r")
polarity<--1##if neg -1, if pos 1
#LockMass<-c(131.9614,165.0188,172.956991,227.201104,255.232405,312.972263,411.12912)##Lock Mass, C8H5O4, C14H27O2, C16H31O2, humic acid

###Step 2: Change source (raw data) and target (calibrated data) folders-------------------------------------------------------

path.source<-paste(path,"/data/20190614/Neg700900", sep="")
path.destination<-paste(path,"/caldata/20190614 jeans/Neg700900", sep="")

###Step 3: Pre-load raw data---------------------------------------------------------------------------------------------------

###Pre-load mzXML raw data into a list. This saves time when running other functions which require multiple xcmsRaw objects.###
xrawdata<-NULL
xrawdata<-list()
setwd(path.source)
msfiles<-list.files(pattern="\\.mzXML$", ignore.case=TRUE)
for (i in 1:length(msfiles)){
  xdata<-xcmsRaw(msfiles[i],includeMSn=TRUE)
  xrawdata[i]<-xdata
}

###Step 4: Mass calibration and writing to new mzXML files---------------------------------------------------------------------

###Mass Calibration###
setwd(path.lock)
LockMass.NEG<-read.table("lockmassProd.csv",header=TRUE,sep=',')
LockMass.NEG<-LockMass.NEG$Lock

for (i in 1:length(msfiles)){
  #xrawdata<-xcmsRaw(msfiles[i],includeMSn=TRUE)
  if (polarity==1){
    xrawcaldata<-MassCal(xrawdata[[i]],LockMass.POS,'POS',3)}###lockshift is used to preliminarly adjust searching space
  if (polarity==-1){
    xrawcaldata<-MassCal(xrawdata[[i]],LockMass.NEG,'NEG',3)}
  setwd(path.destination)
  name<-msfiles[i]
  write.mzdata(xrawcaldata,name)####save the calibrated data to new files, msndata cannot be stored
  print(c(i,"out of",length(msfiles)))
}

