###################
#this script is used for Target-Decoy Search and untargeted chemical analysis
###Hui, 20190307
#--------------------------------------------------------
#Caldata is the calibrated mass of 400-500
#only need the MS2 data from each window
#--------------------------------------------------------
###parameters#####
rm(list=ls())

####################for drinking water, cpp funtion, the oxygen number could be 1.2*carbon+3, cutoff e5, S/N>5, change is.BrCl
library(xcms)
library(MassSpecWavelet)
library(Rcpp)
library(RcppArmadillo)
library(isopat)
library(gtools)
library(caTools)
library(seqinr)
library(rcdk)
data(iso_list)
path<-getwd()##If you opened the current file via the R project, the following directory assignments will be correct
path.data<-paste(path,"/data", sep="")
path.out<-paste(path,"/refCal", sep="")
path.db<-paste(path,"/SMILES_DATABASE", sep="")
path.sirius<-paste(path,"/Sirius", sep="")
path.siriusresult<-paste(path,"/Sirius/results", sep="")
path.calms2<-paste(path,"/400_500", sep="")
source("Nontargeted_fun.r")
polarity<--1##if neg -1, if pos 1
adducts.input<-c('[M]-','[M-H]-','[M-Br+O]-','[M-H-H2O]-','[M+Cl]-','[M+CH2O2-H]-')
Ionsource<-'ESI'
Intensitycut<-10^5###intensity cutoff to pick the peaks for matching
electron<-0.0005485799
ppm<-2
ppm.ms2<-3##ms2 spectra accuracy

###############plot LockMass######################
setwd(path.data)
msfiles<-list.files()
testMass<-255.2324
polarity<--1##if neg -1, if pos 1
LockMass.POS<-c(81.070425,93.070425,105.070425,139.11229,151.042199,171.138505,413.26696,445.120583)##Lock Mass for positive
setwd(path)
LockMass.NEG<-read.table("lockmass.csv",header=TRUE,sep=',')
LockMass.NEG<-LockMass.NEG$Lock
if (polarity==1){
  LockMass<-LockMass.POS
}else{
  LockMass<-LockMass.NEG
}

#----------------------------------------------
#detect peaks from mass spec raw files
#---------------------------------------------
setwd(path)
Library<-read.table("dustlibrary.csv",header=TRUE,sep=',')###the library of peaks
setwd(path.out)
msfiles<-list.files()
xset.raw<-xcmsSet(msfiles,method='centWave',ppm=2.5,peakwidth=c(10,30),snthresh=10,nSlaves=1)##peak width, the min and max range of chromatographic peaks in seconds
xtest<-xcmsRaw(msfiles[1])
mzrange<-xtest@mzrange
####delete those peaks with extreme m/z values, if this step is not done, some issues may happen during getEIC step#######
xset<-xset.raw
index<-which(xset.raw@peaks[,1]<mzrange[1]+1)
if (length(index)>0){
  xset@peaks<-xset.raw@peaks[-index,]}
index<-which(xset@peaks[,1]>mzrange[2]-1)
if (length(index)>0){
  xset@peaks<-xset@peaks[-index,]}
peaklist<-xset@peaks#into is the peak area, maxo is the max intensity                                                                      
len<-length(xset@peaks[,1])

##################group peaks
xset1<-group(xset,bw=60,minsamp=1,minfrac=1/length(msfiles),mzwid=0.001)

######################retention time correction
#xset2<-retcor(xset1,family="symmetric")
#xset2<-group(xset2,bw=30,minsamp=1,minfrac=1/length(msfiles),mzwid=0.001)###re-group using
xset2<-xset1 
######################
setwd(path.out)
xset3<-fillpeak(xset2,10,20)
test<-unlist(xset3@groupidx)##peak ID for each group
#peak_target<-which(test>=min(range_target)&&test<=max(range_target2))
len<-length(xset3@groupidx)#group number
len2<-length(msfiles)#data files
final<-array(rep(0,len*(len2+2)),dim=c(len,(len2+2)))##columns are m/z, rt, averaged intensity, and then 8 ratios
for (i in 1:len){
  print(c('PeakID...',i,'of...',len))
  temp<-unlist(xset3@groupidx[i])
  len3<-length(temp)
  for (j in 1:len3){
    index1<-xset3@peaks[temp[j],11]
    final[i,1]<-xset3@peaks[temp[j],1]##mz
    final[i,2]<-xset3@peaks[temp[j],4]/60##rt 
    final[i,index1+2]<-max(xset3@peaks[temp[j],9],final[i,index1+2])##intensity,select the maximum one if there are multiple peaks for one sample
  }}
final[which(final==0)]<-100##replace 0 values

Library<-final
colnames(Library)<-c('mz','rt',msfiles)
Library<-data.frame(Library)
Library$SampleID<-rep(1,nrow(Library))
index.save<-NULL
for (i in 1:nrow(Library)){##determine the file with the highest abundance
  index<-which.max(Library[i,3:ncol(Library)])
  if (max(Library[i,3:ncol(Library)])>Intensitycut){index.save<-c(index.save,i)}
  Library$SampleID[i]<-index[1]
}
if(length(index.save)>0){###only those peaks with abundances higher than cutoff were selected
  Library<-Library[index.save,]
}


####import TSCA database###
setwd(path.db)
library(ChemmineR)
library(ChemmineOB)
library(RSQLite)
mydb<-dbConnect(RSQLite::SQLite(),'TSCADB_all.db')
Database<-dbReadTable(mydb,"all.TSCA")

##prepare Decoy Database by replacing all oxygen with sulfur
Decoydb<-Database
for (i in 1:nrow(Database)){
  formula<-Decoydb$formula[i]
  Decoydb$formula[i]<-gsub('O','S',formula)
  
  Smile<-Decoydb$smiles[i]
  Decoydb$smiles[i]<-gsub('O','S',Smile)
  
  Smile<-grep('O',unlist(strsplit(Smile,'')))
  mz<-Decoydb$mz[i]+length(Smile)*15.97716##plus the mz of sulfur
  Decoydb$mz[i]<-mz
}


############revaluate the exact mass and retention time using raw file#############
setwd(path.out)
msfiles<-list.files()
Library.new<-mzSmooth(Library,msfiles,120,10)
index<-which(is.na(Library.new$mz))
if (length(index)>0){
Library.new<-Library.new[-index,]}


###########add fragmentation#################
setwd(path)
Library<-Library.new
setwd(path.data)###msn data can not be stored by write.mzdata, so we need to get raw files
msfiles<-list.files()
#################get the fragments by correlation in DIA window####
FragList<-GetFrag(Library,path.data,msfiles,10^(-5),LockMass)

#############double check the results across samples######
Check.Frag<-CheckFrag(FragList,path.data,msfiles,10^(-5),45)

#############remove isotopic peaks of fragments#########
primary.ID<-unique(Check.Frag$libraryid)
index.save<-NULL
for (i in 1:length(primary.ID)){
  index<-which(Check.Frag$libraryid==primary.ID[i])
  if (length(index)<=1){next}
  for (j in 2:length(index)){
    iso.mz<-Check.Frag$mz[index[j]]
    check.mz<-Check.Frag$mz[1:(index[j]-1)]
    check.intensity<-Check.Frag$intensity[1:(index[j]-1)]
    check.Carbon<-sapply(check.mz,FindCarbon,iso.mz,1.00335)
    iso.intensity<-Check.Frag$intensity[index[j]]
    for (k in 1:length(check.mz)){
    if (check.Carbon[k]==1&&iso.intensity<(0.3*check.intensity[k])){
      index.save<-c(index.save,index[j])###if there is any isotope
      break
    }}}}
Check.Frag<-DelList(Check.Frag,index.save)###delete those isotopic fragments

#############mass calibration of fragments##############
setwd(path.calms2)###the datasets for MS1 calibration establishment
msfiles<-list.files()
Cal.Frag<-FragCal(Check.Frag,path.calms2,msfiles,LockMass)

############mass calibration of MS1#############
setwd(path.out)
msfiles<-list.files()
Library.new<-mzSmooth(Library,msfiles,120,10)
index<-which(is.na(Library.new$mz))
if (length(index)>0){
  Library.new$mz[index]<-1
}


###############find isotope#########################################
Isotope.Data<-IsotopeFind(Library.new,path.out,10^(-5))###Identify Isotope Peaks
IsotopeData<-IsotopeAcross(path.out,msfiles,Isotope.Data,10^(-5),30)#####remove noises based on correlations across samples

#############use the primary carbon isotope##########################
primary.ID<-unique(IsotopeData$ID)
for (i in 1:length(primary.ID)){
  index<-which(IsotopeData$ID==primary.ID[i])
  raw.mz<-IsotopeData$mz[index[1]]
  raw.intensity<-IsotopeData$intensity[index[1]]
  for (j in 1:length(index)){
      iso.mz<-IsotopeData$Isotope[index[j]]
      if (iso.mz>raw.mz){break}
      check.Carbon<-FindCarbon(raw.mz,iso.mz,1.00335)##carbon isotope
      check.halogen<-FindCarbon(raw.mz,iso.mz,1.9975)##halogen isotope
      iso.intensity<-IsotopeData$intensity.iso[index[j]]
      if (check.Carbon==1&&raw.intensity<(0.3*iso.intensity)){#carbon
        IsotopeData$mz[index]<-raw.mz-1.00335*check.Carbon
        ID.lib<-IsotopeData$ID[index[j]]
        Library.new$mz[ID.lib]<-iso.mz
      }
      if (check.halogen==1&&raw.intensity>(0.1*iso.intensity)){#carbon
        IsotopeData$mz[index]<-iso.mz
        ID.lib<-IsotopeData$ID[index[j]]
        Library.new$mz[ID.lib]<-iso.mz
      }
      }
  }

##################exclude adducts#######################################
if (polarity==1){
  Adducts<-c('[M+H]+','[M+H-H2O]+','[M+NH4]+','[M+Na]+','[M+CH3OH+H]+')
  MW.adducts<-c(1.007825,-17.00274,18.03437,22.98977,33.03404)
  Adducts.Find<-AdductsFind(Library.new,path.out,MW.adducts,10^(-5))
}

if (polarity==-1){
  Adducts<-c('[M]-','[M-H]-','[M-Br+O]-','[M-H-H2O]-','[M+Cl]-','[M+CH2O2-H]-')
  MW.adducts<-c(0,-1.007825,-62.92342,-19.01839,34.96885,44.99765)
  Adducts.Find<-AdductsFind(Library.new,path.out,MW.adducts,2*10^(-6),Adducts)
}

#-----------------------------------------------
#Target Search
#-----------------------------------------------
cutoff<-5000
mwoffset<-0
mylib.Target<-DatabaseSearching(cutoff,polarity,Database,mwoffset)
setwd(path)

#------------------------------------------------
#Decoy Results
#------------------------------------------------
cutoff<-5000
mwoffset<-0#Ar
mylib.Decoy<-DatabaseSearching(cutoff,polarity,Decoydb,mwoffset)
setwd(path)

#-----------------------
#prepare data for Sirius Searching
#-----------------------
setwd(path.data)
msfiles<-list.files()
xraw<-xcmsRaw(msfiles[1],includeMSn=TRUE)
precursor<-preclist(xraw)
weightK<-c(1,1,1,0.8,0.6,1)#weight for MS1, ms2, ionmode, neutral, characteristic, and adducts, rt
#build MS files
setwd(path.sirius)
Sirius.build(mylib.Target,Cal.Frag,IsotopeData)
setwd("C:/Rprogram/Target_Decoy/Decoy")
Sirius.build(mylib.Decoy,Cal.Frag,IsotopeData)

#------------------------------------
#import MS2 data
#------------------------------------
setwd("C:/Rprogram/Target_Decoy/Decoy/results")
MS2files<-list.files()
mylib.output<-Import.MS2(mylib.Decoy,MS2files,Cal.Frag)
mylib.output$rtscore<-rep(0,nrow(mylib.output))
mylib.output2<-Finalscore(mylib.output,weightK,precursor)
setwd(path)
write.table(mylib.output2,file='Decoy_200.csv',sep=',',row.names = FALSE)

setwd("C:/Rprogram/Target_Decoy/Sirius/results")
MS2files<-list.files()
mylib.output<-Import.MS2(mylib.Target,MS2files,Cal.Frag)
mylib.output$rtscore<-rep(0,nrow(mylib.output))
mylib.output3<-Finalscore(mylib.output,weightK,precursor)
setwd(path)
write.table(mylib.output3,file='Target_200.csv',sep=',',row.names = FALSE)

#------------------------------------------
#initial ID without rt
#------------------------------------------
setwd(path)
Target<-mylib.output3
Decoy<-mylib.output2
cutoff<-Find.cut(Target,Decoy)

#----------------------------------------
#final ID with rt
#----------------------------------------
RT.coeff<-Predict.RT(Target,Database,cutoff)#rt coefficients and standard deviation
if (RT.coeff[4]>0.85){##if correlation is not good use the old one
  RT.coeff<-c(7.6826,0.4694,0.000446,0.85)
}
Target.rt<-Score.RT(Target,Database,RT.coeff)
Decoy.rt<-Score.RT(Decoy,Decoydb,RT.coeff)
cutoff.rt<-Find.cut(Target.rt,Decoy.rt)#recalculate score cutoff
output<-Output(Target.rt,cutoff.rt)
write.table(output,file='FinalID_200.csv',sep=',',row.names = FALSE)

#--------------------------
#unique ID
#-------------------------
Allcpd<-read.table("Sulfur.csv",header=TRUE,sep=',',fill=TRUE)
Uniqueid<-UniqueID(Allcpd)
write.table(Uniqueid,file='UniqueID.csv',sep=',',row.names = FALSE)
write.table(mylib.Target,file='AllID_1ppm.csv',sep=',',row.names = FALSE)



