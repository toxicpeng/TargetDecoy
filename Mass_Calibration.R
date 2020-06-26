###Step 0: Reset R Environment-----------------------------------------------------------------------------------------

###Reset all variables and clear all objects from R memory###
rm(list=ls())

###Step 1: Load packages and set root directory (only run at start of R session)---------------------------------------

###Load packages for code to run###
library(MassSpecWavelet)
library(Rcpp)
library(RcppArmadillo)
library(isopat)
library(gtools)
library(caTools)
library(here)
library(xcms)
data(iso_list)

###Set file pathways. If you are not currently in your root directory or have differently named folders, you will get an error###
###You can use setwd() to change the current directory to the root.###
###If your folder structure is different, you can change the folder names here.###
path<-here()
path.lock<-paste(path,"/analysis/dust analysis", sep="")
###Load function files and set universal variables.###
###If you make changes to function files, you will need to run these 2 lines again###
setwd(path)
source("Nontargeted_fun.r")
polarity<--1##if neg -1, if pos 1
#LockMass<-c(131.9614,165.0188,172.956991,227.201104,255.232405,312.972263,411.12912)##Lock Mass, C8H5O4, C14H27O2, C16H31O2, humic acid

###Step 2: Change source (raw data) and target (calibrated data) folders-------------------------------------------------------

path.source<-paste(path,"/data/Hui_20151207_dust/rawdata/400-505", sep="")
path.destination<-paste(path,"/data/Hui_20151207_dust/caldata/400-505", sep="")

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
if(polarity==1){
  setwd(path.lock)
  LockMass.POS<-read.table("lockmassdust.csv",header=TRUE,sep=',') ##Note! The lock masses should be the most abundant isotopic mass (no e- added)
  LockMass.POS<-LockMass.POS$Lock}
if(polarity==-1){
setwd(path.lock)
LockMass.NEG<-read.table("lockmassdust.csv",header=TRUE,sep=',') ##Note! The lock masses should be the most abundant isotopic mass (no e- added)
LockMass.NEG<-LockMass.NEG$Lock}

#allstats<-list()
for (i in 1:length(msfiles)){
  xrawdata<-xcmsRaw(msfiles[i],includeMSn=TRUE)
  if (polarity==1){
    xrawcaldata<-MassCal(xrawdata[[i]],LockMass.POS,'POS',3)}###lockshift is used to preliminarly adjust searching space
  if (polarity==-1){
    xrawcaldata<-MassCal(xrawdata[[i]],LockMass.NEG,'NEG',3)}
  setwd(path.destination)
  name<-msfiles[i]
  write.mzdata(xrawcaldata,name)####save the calibrated data to new files, msndata cannot be stored
#  allstats[[i]]<-xrawcaldata
    print(c(i,"out of",length(msfiles)))
}

allstats_table<-data.frame(matrix(data=NA,nrow=24,ncol=11))
for (m in 1:length(allstats)){
  allstats_table[m,]<-glance(allstats[[m]])
}
names(allstats_table)<-names(glance(allstats[[1]]))
setwd(path)
write.table(allstats_table,file="calstats.csv",sep=',',row.names = TRUE)
