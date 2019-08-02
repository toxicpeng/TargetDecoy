####
#this script is used for Target-Decoy Search and untargeted chemical analysis
####Hui, 20190307; Steven, 20190729
#--------------------------------------------------------]
#Directions:
  #Check that the folder assignments/file names in Steps 1 and 2 are correct.
  #Delete/move any old SIRIUS files
  #Run Steps 1 through 7
  #Follow Step 8 in SIRIUS
  #Run Steps 9 through 11
#--------------------------------------------------------]

#Note to self: All NTA run up to July 30 inclusive had the following "deleted" variables still present in .Rdata workspace as:
    #ppm<-1.5e-05
    #ppm.ms2<-3
  #These variables have been reinitialized (and ppm correctly set to 2) within the DatabaseSearching function.

#Note to self 2: The combine.mz function (line 1297) seems incredibly inefficient (3 sets of nested if/else statements).
  #This is LOW priority, but the code here could probably be cleaned up a little.
  #Also should probably evaluate the C++ function at line 937 at some point

#Note to self 3: A possible reason for no positive matches in the 700-900 range is that block of code which
  #deletes any masses from the list which fall outside the current range for MS2 peaks. 
  #I don't remember which function that was in though...


###Step 1: Run Once per Sample Set----------------------------------------------------------------------------------------

#######for drinking water, cpp funtion, the oxygen number could be 1.2*carbon+3, cutoff e5, S/N>5, change is.BrCl
library(MassSpecWavelet)
library(Rcpp)
library(RcppArmadillo)
library(isopat)
library(gtools)
library(caTools)
library(seqinr)
library(rcdk)
library(here)
library(xcms)
data(iso_list)

###File paths which will likely remain constant###
path<-here()  ##You should update your folder structure so that the following assignments are correct
path.lock<-paste(path,"/analysis/products analysis/lockmass", sep="") #This assignment can be changed if necessary
path.db<-paste(path,"/SMILES_DATABASE", sep="")
path.siriusTarget<-paste(path,"/Sirius/Target", sep="")
path.siriusTargetresults<-paste(path,"/Sirius/Target/results", sep="")
path.siriusDecoy<-paste(path,"/Sirius/Decoy", sep="")
path.siriusDecoyresults<-paste(path,"/Sirius/Decoy/results", sep="")
###Import function files
setwd(path)
source("Nontargeted_fun.r")

adducts.input<-c('[M]-','[M-H]-','[M-Br+O]-','[M-H-H2O]-','[M+Cl]-','[M+CH2O2-H]-')
Ionsource<-'ESI'
Intensitycut<-10^5###intensity cutoff to pick the peaks for matching

###Step 2: Modify for each run as necessary:------------------------------------------------------------------------------------------

###File paths/names which must be changed with every sample set
path.rawdata<-paste(path,"/data/20190614/Neg700900", sep="")#raw data
path.caldata<-paste(path,"/caldata/20190614 jeans/Neg700900", sep="")#calibrated data
path.finaloutput<-paste(path,"/jeans analysis/Neg700900_ranked", sep="")#folder for ID csv files

target.file<-"Target_Neg700900_ranked.csv" #output for Sirius Target matches with weighted scores
#decoy.file<-"Decoy_Neg700900_nort.csv" #output for Sirius Decoy matches with weighted scores
#RT.ID.file<-"FinalID_Neg700900_nort.csv" #output for Target matches with final scores (including RT prediction)
#uniqueresults.file<-'UniqueID_Neg700900_nort.csv' #output for final match list (duplicates removed)
#allresults.file<-'AllID_Neg700900_nort.csv' #output for all Library data

polarity<--1 ##if neg -1, if pos 1

###Step 3: Import Lockmasses and xcmsRaw files into memory-----------------------------------------------------------------------------

if(polarity==-1){
  setwd(path.lock)
  LockMass.NEG<-read.table("lockmassProd.csv",header=TRUE,sep=',')
  LockMass<-LockMass.NEG$Lock
}
if(polarity==1){
  setwd(path.lock)
  LockMass.POS<-read.table("lockmassProdPOS.csv",header=TRUE,sep=',')
  LockMass<-LockMass.POS$Lock
}

xrawdata<-NULL
xrawdata<-list()
setwd(path.rawdata)
msfiles<-list.files(pattern="\\.mzXML$", ignore.case=TRUE)
for (i in 1:length(msfiles)){
  xdata<-xcmsRaw(msfiles[i],includeMSn=TRUE)
  xrawdata[i]<-xdata
}

xcaldata<-NULL
xcaldata<-list()
setwd(path.caldata)
msfiles<-list.files(pattern="\\.mzXML$", ignore.case=TRUE)
for (i in 1:length(msfiles)){
  xdata<-xcmsRaw(msfiles[i],includeMSn=TRUE)
  xcaldata[i]<-xdata
}

###Step 4: Detect and group calibrated peaks into a Library----------------------------------------------------------------
#---------------------------------------------]
#detect peaks from mass spec calibrated files
#---------------------------------------------]
setwd(path.caldata)
#Library<-read.table("dustlibrary.csv",header=TRUE,sep=',')###the library of peaks
msfiles<-list.files()
xset.raw<-xcmsSet(msfiles,method='centWave',ppm=2.5,peakwidth=c(5,20),snthresh=3)##peak width, the min and max range of chromatographic peaks in seconds
xtest<-xcaldata[[1]]
mzrange<-xtest@mzrange
###delete those peaks with extreme m/z values, if this step is not done, some issues may happen during getEIC step
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
####
setwd(path.caldata)
xset3<-fillpeak(xset2,10,20,xcaldata)
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

############revaluate the exact mass and retention time using calibrated data##
Library.new<-mzSmooth(Library,120,10,xcaldata)
index<-which(is.na(Library.new$mz))
if (length(index)>0){
Library.new<-Library.new[-index,]}
Library<-Library.new


###Step 5: Use MS2 (raw files) fragmentation to recalibrate MS1 (cal files)------------------------------------------------

###########add fragmentation###
#msn data can not be stored by write.mzdata, so we need to get raw files
#get the fragments by correlation in DIA window###
FragList<-GetFrag(Library,xrawdata,10^(-5),LockMass)

#############double check the results across samples##
Check.Frag<-CheckFrag(FragList,xrawdata,10^(-5),45)

#############remove isotopic peaks of fragments##
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

#############mass calibration of fragments using raw data files (same files as used for calibrating MS1)##
Cal.Frag<-FragCal(Check.Frag,xrawdata,LockMass)

############mass calibration of MS1###
Library.new<-mzSmooth(Library,120,10,xcaldata)
index<-which(is.na(Library.new$mz))
if (length(index)>0){
  Library.new$mz[index]<-1
}

###Step 6: Add isotope and adduct data to the Library-----------------------------------------------------------------

###############find isotope###
Isotope.Data<-IsotopeFind(Library.new,xcaldata,10^(-5))###Identify Isotope Peaks
IsotopeData<-IsotopeAcross(xcaldata,Isotope.Data,10^(-5),30)#####remove noises based on correlations across samples

#############use the primary carbon isotope###
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

##################exclude adducts###
if (polarity==1){
  Adducts<-c('[M+H]+','[M+H-H2O]+','[M+NH4]+','[M+Na]+','[M+CH3OH+H]+')
  MW.adducts<-c(1.007825,-17.00274,18.03437,22.98977,33.03404)
  Adducts.Find<-AdductsFind(Library.new,xcaldata,MW.adducts,10^(-5))
}

if (polarity==-1){
  Adducts<-c('[M]-','[M-H]-','[M-Br+O]-','[M-H-H2O]-','[M+Cl]-','[M+CH2O2-H]-')
  MW.adducts<-c(0,-1.007825,-62.92342,-19.01839,34.96885,44.99765)
  Adducts.Find<-AdductsFind(Library.new,xcaldata,MW.adducts,2*10^(-6),Adducts)
}

###Step 7: Use Library to find matches in Target and Decoy databases------------------------------------------------------------

####import TSCA database###
setwd(path.db)
library(ChemmineR)
library(ChemmineOB)
library(RSQLite)
mydb<-dbConnect(RSQLite::SQLite(),'TSCADB_all.db')
Database<-dbReadTable(mydb,"TSCA")

##prepare Decoy Database by replacing all oxygen with sulfur
#Decoydb<-Database
#for (i in 1:nrow(Database)){
#  formula<-Decoydb$formula[i]
#  Decoydb$formula[i]<-gsub('O','S',formula)
#  
#  Smile<-Decoydb$smiles[i]
#  Decoydb$smiles[i]<-gsub('O','S',Smile)
#  
#  Smile<-grep('O',unlist(strsplit(Smile,'')))
#  mz<-Decoydb$mz[i]+length(Smile)*15.97716##plus the mz of sulfur
#  Decoydb$mz[i]<-mz
#  if(i%%1000==0){print(c("converted",i,"out of",nrow(Database)))}
#}


#-----------------------------------------------]
#Target Search
#-----------------------------------------------]
cutoff<-5000
mwoffset<-0
mylib.Target<-DatabaseSearching(cutoff,polarity,Database,mwoffset,xcaldata)


#------------------------------------------------]
#Decoy Results
#------------------------------------------------]
#cutoff<-5000
#mwoffset<-0#Ar
#mylib.Decoy<-DatabaseSearching(cutoff,polarity,Decoydb,mwoffset,xcaldata)


#-----------------------]
#prepare data for Sirius Searching
#-----------------------]
xraw<-xrawdata[[3]]
precursor<-preclist(xraw)
weightK<-c(1,1,1,0.8,0.6,1)#weight for MS1, ms2, ionmode, neutral, characteristic, and adducts, rt
#build MS files
setwd(path.siriusTarget)
Sirius.build(mylib.Target,Cal.Frag,IsotopeData)
#setwd(path.siriusDecoy)
#Sirius.build(mylib.Decoy,Cal.Frag,IsotopeData)

###Step 8: Import target/decoy results into SIRIUS 4 and run Formula Identification (NOT IN R)----------------------------------
  #Batch Import ".ms" files
  #Compute All
    #Add Fluorine (inf)
    #Set "instrument" to Orbitrap (default settings)
    #Turn on CSI:FingerID
      #Search in TSCA or Decoy - Sulfur (based on dataset imported)
      #Make sure the following adducts are selected (negative mode): 
        #[M+Cl]-, [M-Br+O]-, [M-H]-, [M]-, [M-H2O-H]-, [M+CH2O2-H]-
  #Export Results to correct sub-folder (check assignment from step 2)
    #Check off "export tree results"
###Step 9: Use exported SIRIUS csv files to generate statistical match scores---------------------------------------------------

#------------------------------------]
#import MS2 data
#------------------------------------]
#setwd(path.siriusDecoyresults)
#MS2files<-list.files()
#mylib.output<-Import.MS2(mylib.Decoy,MS2files,Cal.Frag)
#mylib.output$rtscore<-rep(0,nrow(mylib.output))
#mylib.output2<-Finalscore(mylib.output,weightK,precursor)
#setwd(path.finaloutput)
#write.table(mylib.output2,file=decoy.file,sep=',',row.names = FALSE)

setwd(path.siriusTargetresults)
MS2files<-list.files()
mylib.output<-Import.MS2(mylib.Target,MS2files,Cal.Frag)
#mylib.output$rtscore<-rep(0,nrow(mylib.output))
#mylib.output3<-Finalscore(mylib.output,weightK,precursor)
mylib.output3<-Finalscorerank(mylib.output,weightK,precursor)
setwd(path.finaloutput)
write.table(mylib.output3,file=target.file,sep=',',row.names = FALSE)

#------------------------------------------]
#initial ID without rt
#------------------------------------------]
#Target<-mylib.output3
#Decoy<-mylib.output2
#cutoff<-Find.cut(Target,Decoy)

###Step 10: Add predicted retention time scores to the compound ID files and apply score cutoff----------------------------------

#----------------------------------------]
#final ID with rt
#----------------------------------------]

#RT.coeff<-Predict.RT(Target,Database,cutoff)#rt coefficients and standard deviation
#if (RT.coeff[4]>0.85){##if correlation is not good use the old one
#  RT.coeff<-c(7.6826,0.4694,0.000446,0.85)
#}
#Target.rt<-Score.RT(Target,Database,RT.coeff)
#Decoy.rt<-Score.RT(Decoy,Decoydb,RT.coeff)
#cutoff.rt<-Find.cut(Target.rt,Decoy.rt)#recalculate score cutoff
#output<-Output(Target,cutoff)
#setwd(path.finaloutput)
#write.table(output,file=RT.ID.file, sep=',',row.names = FALSE)

###Step 11: Eliminate duplicates and give a final compound match for each mass---------------------------------------------------

#-------------------------]
#unique ID
#-------------------------]
#setwd(path.finaloutput)
#Allcpd<-read.table(RT.ID.file,header=TRUE,sep=',',fill=TRUE)
#Uniqueid<-UniqueID(output)
#write.table(Uniqueid,file=uniqueresults.file,sep=',',row.names = FALSE)
#write.table(mylib.Target,file=allresults.file,sep=',',row.names = FALSE)


###Extra: List of all printed counters---------------------------------------------------------------------------------
  #Step 4
    #fillpeak = "filled peak id [i] out of [number of files]"
    #mzSmooth = "smoothed [i] out of [number of files]"
  #Step 5
    #GetFrag = "got fragments in [i] out of [number of files]"
    #CheckFrag = "checked fragments in [i] out of [number of files]"
    #FragCal = "calibrated fragments in [i] out of [number of files]"
    #mzSmooth = "smoothed [i] out of [number of files]"
  #Step 6
    #IsotopeFind = "found isotopes in [i] out of [number of files]"
    #IsotopeAcross = "checked isotope noise in [i] out of [number of files]"
    #AdductsFind = "found adducts in [i] out of [number of files]"
  #Step 7
    #Decoydb = "converted [i] out of [number of rows in the database]"
    #Database Searching (Target)
      # = "calculated isotopes for [i] out of [number of rows in the target database]"
      # = "calculated scores for [i] out of [number of files]"
      # = "scored ionmode in [i] out of [number of rows in the library]"
    #Database Searching (Decoy)
      # = "calculated isotopes for [i] out of [number of rows in the decoy database]"
      # = "calculated scores for [i] out of [number of files]"
      # = "scored ionmode in [i] out of [number of rows in the library]"
    #Sirius.build (Target) = "built [i] out of [number of rows in the library]"
    #Sirius.build (Decoy) = "built [i] out of [number of rows in the library]

  #Step 10 (SKIPPED for now)
    #Predict.RT = "current R^2 is [R coefficient]"
    #Score.RT (Target) = "generated RT score for [i] out of [number of rows in the Target match file]"
    #Score.RT (Decoy) = "generated RT score for [i] out of [number of rows in the Decoy match file]"
  #Step 11
    #UniqueID = "finished [i] out of [number of rows in the library]"



