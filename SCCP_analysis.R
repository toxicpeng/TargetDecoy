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
path_data<-'E:/Steven/Raw Data/20190918'
path_results<-'E:/Steven/Target Decoy/analysis/products analysis'
path_constants<-'E:/Steven/Target Decoy/analysis/products analysis'
setwd(path_data)
msfiles<-list.files()
msfiles<-mixedsort(msfiles)###sort the file name according to number

msfiles<-msfiles[164]

features.filename<-"SCCP feature list.csv"
xreport.filename<-"SCCPsPeakResults.csv"

#####Create list of SCCPs [M+Cl]- formulae, abundant isotopic mass, isotope %, [M] formula, average mass, average chlorine mass, chlorine mass%#####
CPlist<-data.frame(matrix(ncol=8,nrow=0))
Cnum<-1
Cmax<-25
Clnum<-1
CPformulas<-NULL
CPformulas_plusCl<-NULL
Clformulas<-NULL
while(Cnum<=25){

    Clmax<-(2*Cnum)+2
    for(j in 1:Clmax){
      Hnum<-(2*Cnum)+2-j
      if(Hnum<1){
        formula<-paste0("C",Cnum,"Cl",j)
        formula_plusCl<-paste0("C",Cnum,"Cl",j+1)
        Cls<-paste0("Cl",j)
      }else{
        formula<-paste0("C",Cnum,"H",Hnum,"Cl",j)
        formula_plusCl<-paste0("C",Cnum,"H",Hnum,"Cl",j+1)
        Cls<-paste0("Cl",j)
      }
      CPformulas<-c(CPformulas,formula)
      CPformulas_plusCl<-c(CPformulas_plusCl,formula_plusCl)
      Clformulas<-c(Clformulas,Cls)
    }
  
  Cnum<-Cnum+1
}

for(k in 1:length(CPformulas_plusCl)){
  alliso1<-isopattern(iso_list,CPformulas_plusCl[k],0.01)
  isomax<-which.max(alliso1[,2])
  CPiso<-c(CPformulas_plusCl[k],alliso1[isomax,1],alliso1[isomax,1]+0.00055,alliso1[isomax,2])
  CPlist[k,1:4]<-CPiso
  if(k%%100==0){print(k)}
}

for(m in 1:length(CPformulas)){
  alliso2<-isopattern(iso_list,CPformulas[m],0.01)
  isoweighted<-NULL
  isoweightedCl<-NULL
  for(n in 1:length(alliso2[,1])){
    weightedmass<-alliso2[n,1]*alliso2[n,2]
    isoweighted<-c(isoweighted,weightedmass)
  }
  avgmass<-sum(isoweighted)
  Cliso<-isopattern(iso_list,Clformulas[m],0.01)
  for(p in 1:length(Cliso[,1])){
    weightedmassCl<-Cliso[p,1]*Cliso[p,2]
    isoweightedCl<-c(isoweightedCl,weightedmassCl)
  }
  avgClmass<-sum(isoweightedCl)
  
  CPavg<-c(CPformulas[m],avgmass,avgClmass,(avgClmass/avgmass))
  CPlist[m,5:8]<-CPavg
  if(m%%100==0){print(m)}
}

#names(CPlist)<-c("Compound (+ Cl-)","Most Abundant Mass", "Mass + e-","% Abundance","Compound","Average Molecular Weight","Cl Average Weight","Cl Content (mass %)")
setwd(here())
write.table(CPlist, file="ALLCPsList.csv", sep = ',',row.names=FALSE,col.names=cbind("Compound (+ Cl-)","Most Abundant Mass", "Mass + e-","% Abundance","Compound","Average Molecular Weight","Cl Average Weight","Cl Content (mass %)"))


#####CP Concentration Calculations#####
setwd(path_constants)
CPref<-read.table(features.filename,header=TRUE,sep=',',fill=TRUE)

setwd(path_data)
CPareas<-read.table(xreport.filename,header=TRUE,sep=',',fill=TRUE)
CPareas[,4]<-as.numeric(CPareas[,4])
key<-which(CPareas[,2]==CPareas[1,2])
for(i in 1:length(key)){
  startline<-key[i]
  ISarea<-CPareas[startline,4]
  sampleCPs<-a
  
  
  
  
}



