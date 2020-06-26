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



#####RF-Cl linear estimation for all congeners (instead of mixture)#####

path_data<-"C:/Users/Steven Desktop/Documents/MSc2/TargetDecoy/analysis/products analysis/reference compounds analysis"
setwd(path_data)
CPknown<-read.table("SCCP linear estimation values.csv",header=TRUE,sep=',',fill=TRUE)
CongAreas<-CPknown[,2] #Congener Peak Areas (normalized to IS)
CongCl<-CPknown[,3] #Congener mass% Cl (based on formulas and average molecular mass)
Components<-nrow(CPknown)

MixConc<-0.1 #Concentration of the technical mixture used to collect this data (ppm)
MixCl<-0.515 #Manufacturer %Cl of the technical mixture
LS<-MixCl #Cl Concentration in the technical mixture (ppm)

CongConcLow<-0.00001 #Lowest expected concentration (in ppm) of individual congeners
CongConcHigh<-MixConc #Highest possible concentration (in ppm) of individual congeners (not possible to be higher than the concentration of the technical mixture)
xLow<-min(CongCl)
xHigh<-max(CongCl)
yLow<-(min(CongAreas))/CongConcHigh #Lowest possible y value, based on input data
yHigh<-(max(CongAreas))/CongConcLow #Highest possible y value based on input data
SlopeHigh<-(yHigh-yLow)/(xHigh-xLow) #delta y/delta x
SlopeLow<-(yLow-yHigh)/(xHigh-xLow) #denominator stays the same (otherwise the signs would cancel out and slopeLow would equal slopeHigh)
yIntLow<-yLow-(SlopeHigh*xHigh) #y=mx+b rearranges to b=y-mx. To get the lowest possible yInt, y must be small and mx must be large.
yIntHigh<-yHigh-(SlopeLow*xHigh) #To get the highest possible yInt, y must be large and mx must be small (or a large negative). xHigh is used here because Slope Low will either be negative or zero.

Rounding<-100 #Rounding off for the purposes of faster testing
#SlopeLow<-round(SlopeLow, digits=-log10(Rounding))
SlopeHigh<-round(SlopeHigh, digits=-log10(Rounding))
SlopeLow<-100000
#SlopeHigh<-100000
yIntLow<-round(yIntLow, digits=-log10(Rounding))
yIntHigh<-round(yIntHigh, digits=-log10(Rounding))

runcount<-(SlopeHigh-SlopeLow)
PassFactor<-1000 #A set of parameters will have passed the test if RS=LS within a factor of [PassFactor]
PassList<-matrix(data=NA,nrow=0,ncol=4)
digits<-floor(log10(runcount)) #Order of magnitude of the number of iterations of the main for-loop (for print counter purposes; see below)
printcounter<-10^(digits-2) #Order of magnitude for which the print counter will return <1000 values (see below)
print(c(runcount,printcounter))

start_time<-Sys.time()
for(Slope in seq(SlopeLow,SlopeHigh,by=Rounding)){
  for(yInt in seq(yIntLow,yIntHigh,by=Rounding)){
    if(yInt==0){next}
    RStop<-0
    RSbot<-0
    RS<-0
    CongenersTop<-NULL
    CongenersBot<-NULL
    for(k in 1:Components){
      if(CongAreas[k]==0){next}
      CongenersTop<-c(CongenersTop,((CongAreas[k]*CongCl[k])/(Slope*CongCl[k]+yInt)))
    }
    if(length(which(CongenersTop<0))>0){next}
    for(m in 1:Components){
      if(CongAreas[m]==0){next}
      CongenersBot<-c(CongenersBot,(CongAreas[m]/(Slope*CongCl[m]+yInt)))
    }
    if(length(which(CongenersBot<0))>0){next}
    if(RStop==0 | RSbot==0){next}
    RS<-RStop/RSbot
    if(RS>(LS/PassFactor) & RS<(LS*PassFactor)){
      ppmAcc<-abs(((RS-LS)/LS)*(10^6))
      PassList<-rbind(PassList,c(Slope,yInt,RS,ppmAcc))
    }
  }
  if(Slope==SlopeLow){print(c("First iteration complete"))}
  if(Slope%%printcounter==0 & Slope>SlopeLow){
    print(c("current slope =",Slope,"# of results so far=",nrow(PassList)))
    timetogo(start_time,(Slope-SlopeLow),runcount)
    }
}
print(c("number of results =",nrow(PassList)))
BestRow<-which.min(PassList[,4])
print(PassList[BestRow,])

#Function which uses the time passed between groups of for loop iterations to estimate the time remaining.
#***Requires the following code to be run immediately before the 1st iteration of the for loop: [starttime variable]<-Sys.time()
timetogo<-function(starttime,currentcount,totalcount){  #starttime = system time at the start of the for loop, currentcount = current for loop iteration, totalcount = total # of for loop iterations
  time_running<-difftime(Sys.time(),starttime,units="mins") #Calculate the time difference in minutes between now and the start of the for loop
  time_running<-as.numeric(time_running)
  fraction_complete<-round((currentcount/totalcount),digits=10) #Determine what fraction of the total # of iterations has been completed to this point
  total_time<-time_running/fraction_complete #Estimate the total runtime of the for loop based on how long it took to get to the current iteration
  minutes_to_go<-round((total_time-time_running),digits=1) #The remaining runtime will be the difference between total runtime and time passed thus far
  if(minutes_to_go>60){print(c("Time remaining =",round(minutes_to_go/60,digits=2),"hours"))}else{print(c("Time remaining =",minutes_to_go,"minutes"))}
}

  #For a single run of 51.5% Cl
#For Rounding==100 and PassFactor==10, nrow(PassList)==6
#For Rounding==10 and PassFactor==10, nrow(PassList)==227
#For Rounding==10 and PassFactor==5, nrow(PassList)==70
#For Rounding==10 and PassFactor==2, nrow(PassList)==15
#For Rounding==1 and PassFactor==2, nrow(PassList)==771
#For Rounding==1 and PassFactor==1.1, nrow(PassList)==166
  #At this point, upper bounds for Slope and y-Int were both set to 100, based on max values from earlier tests
#For Rounding==0.1 and PassFactor==1.1, nrow(PassList)==12921, min ppm = 31.6
#For Rounding==0.1 and PassFactor==1.01, nrow(PassList)==2280, min ppm = 12.3
#For Rounding==0.1 and PassFactor==1.001, nrow(PassList)==257, min ppm = 4.08
#For Rounding==0.01 and PassFactor==1.001, nrow(PassList)==24998, min ppm = 0.0650
#For Rounding==0.01 and PassFactor==1.0001, nrow(PassList)==2562, min ppm = 0.0650
#For Rounding==0.01 and PassFactor==1.00001, nrow(PassList)==285, min ppm = 0.0650
  #At this point, the rest of the SCCP values (10 runs total of 51.5% Cl) were added to the imported data list
  #Upper bounds for Slope and y-Int were reset to calculated values
#For Rounding==10 and PassFactor==2, nrow(PassList)==628, min ppm = 452
  #At this point, lower estimate of congener concentration was lowered from 0.1 ppb to 0.01 ppb
#For Rounding==10 and PassFactor==1.1, nrow(PassList)==130, min ppm = 452
#For Rounding==1 and PassFactor==1.1, nrow(PassList)==10517, min ppm = 12.9
#For Rounding==0.1 and PassFactor==1.1, nrow(PassList)==923125, min ppm = 0.133
#For Rounding==0.1 and PassFactor==1.001, nrow(PassList)==20176, min ppm = 0.0699
#For Rounding==0.1 and PassFactor==1.00001, nrow(PassList)==218, min ppm = 0.0699
  #New upper and lower bounds calculated for Slope and yInt