######### Head #########

###Reset all variables and clear all objects from R memory###
rm(list=ls())

###Load packages for code to run. Only need to run at the start of your R session.###
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


##Import NTA Ranked Target csv file, delete all columns except mass, Smiles, and score
##Put each match in its own cell, then delete any rows with a 0 or negative score

path<-here() 
path.analysis<-paste(path,"/jeans analysis/Neg700900_ranked", sep="")
setwd(path.analysis)
rankedtable<-read.csv("Target_Neg700900_ranked.csv",header=TRUE,sep=',',fill=TRUE, stringsAsFactors = FALSE)

expandtable<-as.data.frame(matrix(data=0,nrow=2500,ncol=4))
names(expandtable)<-c("mz","formula","SMILES","matchscore")
p<-1
for (i in 1:nrow(rankedtable)){
  scoreset<-strsplit(as.character(rankedtable$allscore[i]),';')
  scoreset<-as.numeric(unlist(scoreset))
  formset<-strsplit(as.character(rankedtable$formula.pred[i]),';')
  formset<-unlist(formset)
  smileset<-strsplit(as.character(rankedtable$SMILES[i]),';')
  smileset<-unlist(smileset)
  mz<-rankedtable$mz[i]
  
  if(length(scoreset)==0){next}
  if(length(scoreset)==1){
    expandtable$matchscore[p]<-scoreset
    expandtable$mz[p]<-mz
    expandtable$formula[p]<-formset
    expandtable$SMILES[p]<-smileset
    p<-p+1
    next
  }
  
  for(j in 1:length(scoreset)){
    expandtable$matchscore[p]<-scoreset[j]
    expandtable$mz[p]<-mz
    expandtable$formula[p]<-formset[j]
    expandtable$SMILES[p]<-smileset[j]
    p<-p+1
  }
}

lowscoreindex<-NULL
for(k in 1:nrow(expandtable)){
  if(expandtable$matchscore[k]<=0){
    lowscoreindex<-c(lowscoreindex,k)
  }
}

expandtable<-expandtable[-lowscoreindex,]
write.csv(expandtable,file="compoundmatches_700to900.csv",row.names = FALSE)



##Determining Formulae - for lock mass calibration figure##
formcalnum<-function(number_input,Peakinfo,index.input,iso_list,rawdata,ppm,ppm.ms2,IsotopeList,Database,Fragment,adducts,IsotopeDB,mwoffset){##number input is the limit of element composition number, isolist is the mw of element
  #########restrict the element number
  mz<-Peakinfo$mz[index.input]
  
  ###################constrain adducts###################
  mw.offset<-t(sapply(adducts,DefineAdducts))
  formula.index<-NULL
  Adducts.pred<-NULL
  mserror<-NULL
  for (k in 1:length(adducts)){
    mz.pred<-mz+0.0005485799*polarity-mw.offset[k,12]###the mz used for formula prediction, should consider electrons
    isotope.temp<-which(abs(mz.pred-IsotopeDB$mz)<mz.pred*ppm*10^(-6))##find the peaks
    if(length(isotope.temp)<1){next}
    formula.index.temp<-IsotopeDB$ID[isotope.temp]
    ###bromine###########
    if (adducts[k]=='[M-Br+O]-'){
      formula.temp<-Database$formula[formula.index.temp]
      Br.id<-grep('Br',formula.temp)
      if (length(Br.id)==0){next}##M-Br+O detected, but no bromine, deleted
      formula.index.temp<-formula.index.temp[Br.id]
    }
    mserror<-c(mserror,(mz.pred-IsotopeDB$mz[isotope.temp])*10^6/mz.pred)
    formula.index<-c(formula.index,formula.index.temp)
    Adducts.pred<-c(Adducts.pred,rep(adducts[k],length(formula.index.temp)))
  }
  if (length(formula.index)==0){return(NULL)}##no possible formulas
  formula.pred<-as.character(Database$formula[formula.index])
  SMILES.pred<-as.character(Database$smiles[formula.index])
  DatabaseID<-formula.index
  
  ###############delete those formulas with rare elements#############
  test.element<-'CHONPSClBrI0123456789'
  index<-NULL
  for (i in 1:length(formula.pred)){
    element.test<-strsplit(formula.pred[i],'')
    element.test<-unlist(element.test)
    match.element<-sapply(element.test,grep,test.element)
    if (length(match.element)==length(unlist(match.element))){##no rare elements
      index<-c(index,i)
    }}
  if(is.null(index)){return(NULL)}
  formula.pred<-formula.pred[index]
  Adducts.pred<-Adducts.pred[index]
  SMILES.pred<-SMILES.pred[index]
  mserror<-mserror[index]
  
  
  ####################isotope distribution##############################
  formula.split<-sapply(formula.pred,form.split)
  RMSE.iso1<-0
  if (length(ncol(formula.split))>0){
    form.offset<-sapply(Adducts.pred,DefineAdducts)
    formula.split<-formula.split+form.offset[1:11,]
    index.save<-NULL
    for (kk in 1:ncol(formula.split)){#Check if the number of element is less than -1
      if (min(formula.split[,kk])<0){
        index.save<-c(index.save,kk)
      }
    }
    if (length(index.save)>0){
      formula.split<-formula.split[,-index.save]
    }
    formula.split<-matrix(formula.split,nrow=11)
    if (length(formula.split)==0){return(NULL)}
    RMSE.iso1<-apply(formula.split,2,deiso.formula,iso_list,index.input,Peakinfo,IsotopeList)
  }
  if(RMSE.iso1[1]==0){return(NULL)}###error in isotope calculation because of rare element
  RMSE.iso<-exp((-0.5)*RMSE.iso1)
  
  ###############calculate the ms1 error
  MS1.score<-rep(0,length(mserror))
  index.small<-which(abs(mserror)<1)
  if (length(index.small)>0){MS1.score[index.small]<-1}##when smaller than 1ppm, no difference
  index.big<-which(abs(mserror)>=1)
  if (length(index.big)>0){
    MS1.score[index.big]<-exp((-0.5)*(mserror[index.big]/0.5)^2)}##2ppm for delta
  
  
  ###############calculate the subset of formulas and ms2 error##
  fragments<-which(Fragment$libraryid==index.input)
  index.fragments<-which(Fragment$intensity[fragments]>2000)##only abundant fragments are used
  ms2score<-rep(0,length(formula.pred))
  for (j in 1:length(formula.pred)){
    if (length(index.fragments)<1){break}##no fragments
    fragments<-fragments[index.fragments]
    number.limit<-form.split(formula.pred[j])##subset of the precusor formula
    form.offset<-sapply(Adducts.pred,DefineAdducts)
    number.limit<-number.limit+form.offset[1:11,]
    number.limit[8]<-number.limit[7]
    number.limit[10]<-number.limit[9]
    element.DAT<-c("C","H","N","O","P","S","35Cl","37Cl","79Br","81Br","I")
    mw.element<-iso.element(element.DAT,iso_list,1)##calculate mw for each element
    number.low<-rep(0,11)
    numberset<-cbind(number.low,number.limit,element.DAT)
    number_list<-c(1,3,4,3,2,3,2,2,2,2,2,2)#the ratio of element to carbon
    ###################fragments##################
    score.temp<-0
    for (k in 1:length(fragments)){
      mz.frag<-Fragment$mz[fragments[k]]
      mz.pred<-mz.frag+0.0005485799*polarity###the mz used for formula prediction, should consider electrons
      formula.frag<-formpred(mz.pred,ppm.ms2,numberset,mw.element[,3],number_list,mwoffset)##primary calculation
      if (max(formula.frag)==0){score.temp<-score.temp-1}else{
        score.temp<-score.temp+1
      }}
    ms2score[j]<-score.temp/length(fragments)}
  
  #################total score####################################
  MS1.score<-log(MS1.score+0.01)
  RMSE.iso<-log(RMSE.iso+0.01)
  all.score<-MS1.score+RMSE.iso+ms2score
  if (length(all.score)<1){
    return(NULL)}
  rank1<-order(-all.score)
  index<-rank1[1:length(rank1)]
  if (length(index)>1){
    formula.final<-cbind(formula.pred[index],Adducts.pred[index],SMILES.pred[index],DatabaseID[index],MS1.score[index],RMSE.iso[index],ms2score[index],all.score[index])
  }else{
    formula.final<-c(formula.pred[index],Adducts.pred[index],SMILES.pred[index],DatabaseID[index],MS1.score[index],RMSE.iso[index],ms2score[index],all.score[index])
  }
  formula.final<-matrix(formula.final,ncol=8)
  return(formula.final)}###just output the best one



##Gathering Formula prediction numbers for lock mass calibration figure##
DatabaseSearchingFormulae<-function(cutoff,polarity,Database,mwoffset,xcallist){#cutoff for intensity
    element<-c("C","H","N","O","P","S","35Cl","37Cl","79Br","81Br","I")
    number_input<-t(rbind(c(1,0,0,1,0,0,0,0,0,0,0),c(50,80,10,10,0,0,2,1,7,3,0)))
    rownames(number_input)<-element
    isotope<-1.9978##the value of isotope peaks
    inputvalue<-"yes"##for some compounds we don't need isotope information
    peakwidth<-10#peak width, seconds
    samplenum<-21##the number of sample
    cutratio<-5###cutoff for ratio difference
    timewin<-12###12 sec for rt cutoff
    mw.element<-iso.element(element,iso_list,1)
    
    ######calcualte isotope pattern for halogenated compounds###
    IsotopeDB<-NULL####save isotopic mz to IsotopeDB
    for (i in 1:nrow(Database)){
      formula<-Database$formula[i]
      carbon.index<-grep('C',formula)##must contain carbon
      if (length(carbon.index)==0){next}
      
      formula.split<-form.split(formula)
      formula.cal<-form.comb(formula.split)
      if (max(formula.split)==0) {next}
      iso.pattern<-isopattern_cal(iso_list,formula.cal,1e-4)
      iso.pattern[,2]<-iso.pattern[,2]/max(iso.pattern[,2])
      index<-which(iso.pattern[,2]>0.05)##only keep the isotopic pattern
      if (length(index)>1){
        IsotopeDB$mz<-c(IsotopeDB$mz,iso.pattern[index,1])
        IsotopeDB$ID<-c(IsotopeDB$ID,rep(i,length(index)))
      }
      if(i%%1000==0){print(c('calculated isotopes for',i,"out of",nrow(Database)))}
    }
    
    Database$mz<-Database$mz+mwoffset
    IsotopeDB$mz<-IsotopeDB$mz+mwoffset
    
    ppm<-1
    ppm.ms2<-3
    
    mylib<-Library.new
    mylib$formula.pred<-rep(0,nrow(Library.new))
    mylib$totalformulas<-rep(0,nrow(Library.new))
    mylib$isoscore<-rep(0,nrow(Library.new))
    mylib$mserror<-rep(0,nrow(Library.new))
    mylib$ms1score<-rep(0,nrow(Library.new))
    mylib$SMILES<-rep(0,nrow(Library.new))
    mylib$DbID<-rep(0,nrow(Library.new))
    mylib$Adducts<-rep(0,nrow(Library.new))
    for (i in 1:length(xcallist)){
      index<-which(Library.new$SampleID==i)
      if (length(index)==0){next}
      xrawdata<-xcallist[[i]]
      for (j in 1:length(index)){
        formula<-NULL
        formula<-formcalnum(number_input,Library.new,index[j],iso_list,xrawdata,ppm,ppm.ms2,IsotopeData,Database,Cal.Frag,adducts.input,IsotopeDB,mwoffset)
        if (length(formula)<2){next}
        #max.form<-nrow(formula) Not used???
        
        num_formulae<-length(formula[,1])
        mylib$totalformulas[index[j]]<-num_formulae
        
        mylib$formula.pred[index[j]]<-paste(formula[,1],collapse=';')
        mylib$isoscore[index[j]]<-paste(formula[,6],collapse=';')
        mylib$mserror[index[j]]<-paste(formula[,5],collapse=';')
        mylib$ms1score[index[j]]<-paste(formula[,8],collapse=';')
        mylib$SMILES[index[j]]<-paste(formula[,3],collapse=';')
        mylib$DbID[index[j]]<-paste(formula[,4],collapse=';')
        mylib$Adducts[index[j]]<-paste(formula[,2],collapse=';')
      }
      print(c('calculated scores for',i,"out of",length(xcallist)))
    }

setwd(path)
write.csv(mylib,file="housedustmatches_1ppm.csv",row.names = FALSE)

}


######Code for deleting zero values in a table, and reducing the number of data points for easier graphing######

tempmatches<-read.csv("formulamatchestemp.csv",header=TRUE,sep=',',fill=TRUE, stringsAsFactors = FALSE)

form10ppm<-tempmatches[,1:2]
form5ppm<-tempmatches[,3:4]
form1ppm<-tempmatches[,5:6]

currenttable<-form10ppm
indexf<-which(currenttable[,2]==0)
currenttable<-currenttable[-indexf,]


newtable<-NULL
for(i in seq(151,1001,by=1)){
  massrange<-c(i-1,i)
  rangeindex<-which(currenttable[,1]>=massrange[1] & currenttable[,1]<massrange[2])
  if(length(rangeindex)<1){next}
  newcount<-median(currenttable[rangeindex,2])
  newcount<-round(newcount)
  newtable<-rbind(newtable,c(massrange[1],newcount))
}
newtable<-as.data.frame(newtable)
names(newtable)<-c("mz10ppm","matches10ppm")
write.csv(newtable, file="10ppmmatches.csv",row.names=FALSE)

