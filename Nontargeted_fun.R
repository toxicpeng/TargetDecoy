###########function to read formulae

##########FUNCTION to plot lockmass shift########
plotlock<-function(xraw,Lockmass,ppm){
  ppmmatrix<-matrix(rep(100,length(xraw@scanindex)*length(Lockmass)),nrow=length(Lockmass),ncol=length(xraw@scanindex))
  Lockmass<-Lockmass-0.0005485799*polarity
  for (i in 2:length(xraw@scanindex)){##mass calibration using the closest lockmass
    scanNum<-c(xraw@scanindex[i-1],xraw@scanindex[i])
    correctindex<-(scanNum[1]+1):scanNum[2]
    for (j in 1:length(Lockmass)){
      mz.index<-which(abs(xraw@env$mz[correctindex]-Lockmass[j])<ppm*10^(-6)*Lockmass[j])
      if (length(mz.index)==0){
        next}
      index2<-which.max(xraw@env$intensity[correctindex[mz.index]])##if multiple data exsit, pick out the most abundant one, questionable?
      mz.index<-mz.index[index2]
      mz.index<-mz.index[1]
      mz.index<-scanNum[1]+mz.index
      ppmmatrix[j,i]<-10^6*(xraw@env$mz[mz.index]-Lockmass[j])/Lockmass[j]}         
  }
  return(ppmmatrix)
}


##############function to identify lockmass, written by Steven 2017/12/20
findlock<-function(xrawAll,xset.input,intthresh,stepmz){
  newpeak<-NULL
  p<-0
  msfile<-filepaths(xset.input)
  for (k in 1:length(unlist(phenoData(xset.input)))){
    #xraw<-xcmsRaw(msfile[k])  Old code, kept just in case
    xraw<-xrawAll[[k]]
    for (mz.value in seq(from=300,to=500,by=stepmz)){
      mz.min<-mz.value-0.0005 
      mz.min<-max(mz.min, xraw@mzrange[1])
      mz.max<-mz.value+0.0005
      mz.max<-min(mz.max,xraw@mzrange[2])
      rt.min<-(min(xraw@scantime)-60) #range of analysis is the entire chromatograph
      rt.max<-(max(xraw@scantime)-60) #the first and last 1 minute are skipped to avoid sections with no analytes
      peaks<-rawEIC(xraw,mzrange=cbind(mz.min,mz.max),rtrange=cbind(rt.min, rt.max))
      if(mz.value==300 || mz.value==400){
        print(c("current mass=",mz.value,"sample=",k,"results so far=",p))
      }
      
      peaklist<-unlist(peaks$intensity)
      newlist<-NULL
      for(i in 1:length(peaklist)){
        if (peaklist[i]>0){newlist<-c(newlist,peaklist[i])}
        else if (is.na(peaklist[i+1])==TRUE||is.na(peaklist[i+2])==TRUE||is.na(peaklist[i+3])==TRUE||is.na(peaklist[i+4])==TRUE){next}
        else if (peaklist[i+1]>0||peaklist[i+2]>0||peaklist[i+3]>0||peaklist[i+4]>0){next} ##allows mz values through which do not have consecutive zero intensity values
        else {newlist<-c(newlist,peaklist[i])}}
      if (length(newlist)==0){next}
     
      minintensity<-min(newlist)
      
      if (minintensity>intthresh){
        newpeak<-rbind(newpeak,c(mz.value,minintensity,k))
        p<-p+1
        #print(c("running...",p,"sample =",k))
      }
    }
  }
  return (newpeak)
}

##############################Useful small functions#########################

##Function for manually calculating the AICc value of a curve fit
manualAICc<-function(lmobject){
  RSS<-sum(lmobject$residuals^2)
  nobs<-nrow(lmobject$model)
  K<-length(lmobject$coefficients)
  r1<-nobs*log(RSS/nobs)
  r2<-2*K
  r3<-(2*K*(K+1))/(nobs-K-1)
  result<-r1+r2+r3
  return(result)
}  
##Function for calculating the AICc of two combined curves
CalAIC<-function(RSS,K,N){
  r1<-N*log(RSS/N)
  r2<-2*K
  r3<-(2*K*(K+1))/(N-K-1)
  AICs<-r1+r2+r3
  return(AICs)
}

#########Calfunctions#######
CalFun<-function(method,x,y){
  if (method==1){
    fits<-lm(y~x)}
  if (method==2){
    fits<-lm(y~poly(x,2,raw=TRUE))}
  if (method==3){
    fits<-lm(y~poly(x,2,raw=TRUE))}
  if (method==4){
    fits<-lm(y~log(x))}
  if (method==5){
    fits<-lm(y~I(1/x))}
  return(fits)}

#########predicfunctions#######
PredictFun<-function(obj,x){
  library(stats)
  y<-predict(obj,data.frame(x=x))
  return(y)}     


###################curve fitting#####################
fitlock<-function(lockdf,mzlist,inputfun){
  ##If there is only 1 lockmass in the current column, create a flat curve and return it.
  table<-lockdf
  if(length(table[,1])<2){
    calmass<-mzlist/(1+table[1,2])
    return (calmass)
  }
  
  ##Attempts to first determine if a single curve can be used to fit all lock masses
  R.single<-0
  table.single<-table
  kk<-0
  while(R.single<0.8&&kk<3){##do not delete too many points
    kk<-kk+1
    fitlistall<-list()
    for (k in 1:5){
      results<-CalFun(k,table.single[,1],table.single[,2])
      fitlistall[[k]]<-results
    }
    AIC.single<-sapply(fitlistall, manualAICc)
    index.single<-which.min(AIC.single)
    AIC.single<-min(AIC.single)
    S.single<-summary(fitlistall[[index.single]])
    R.single<-S.single$r.squared
    res.single<-S.single$residuals
    if (R.single<0.8){
      index<-which.max(abs(res.single))##the maximal residual error
      table.single<-table.single[-index,]}
  }
  
#####save the curve fitting if it is pretty good########
savefun<-inputfun
if (nrow(table.single)>8&&R.single>0.9){
  savefun<-fitlistall[[index.single]]
}

  
  #####predict new mass###############
  table<-table.single
    fit.single<-fitlistall[[index.single]]
    calmass1<-NULL
    calmass2<-NULL
    calmass3<-NULL
    if (length(savefun)>0){
    if (R.single<0.8||nrow(table.single)<5){####if fitting is poor and inputfun exists, use the fitted curve saved before
      fit.single<-savefun
    }}
    index<-which(mzlist<max(table[,1]))
    index.rep<-which(mzlist[index]>min(table[,1]))
    index<-index[index.rep]
    if (length(index)>0){
      masserror<-PredictFun(fit.single,mzlist[index])
      calmass1<-mzlist[index]*(1-masserror)}
    index2<-which(mzlist>=max(table[,1]))
#    masserror2<-PredictFun(fit.single,269.1389)###predict at 269.1389 <<<<<<<<<<<<<<<<<<<<<<<<<<This was the breakpoint between polynomial and linear fitting.  Change this value to automatic detection.
    masserror2<-(mzlist[index2]*6.778e-09)+(-3.389e-06)
    if (length(index2)>0){
      calmass2<-mzlist[index2]*(1-masserror2)}
    index3<-which(mzlist<=min(table.single[,1]))###lesser than the minium one, use fitted function
    if (length(index3)>0){
      masserror3<-PredictFun(fit.single,mzlist[index3])
      calmass3<-mzlist[index3]*(1-masserror3)
    }
    calmass<-c(calmass3,calmass1,calmass2)
  return(list(calmass,savefun))}


##############formula prediction##############
form.parse<-function(rawformula){##formulae with two column, the second column is formulae
  formulae.save<-NULL
  for (i in 1:nrow(rawformula)){
    if (rawformula[2]=='0'){
      formulae.save<-rbind(formulae.save,rep(0,9))
      next}
    formula.input<-as.character(rawformula[2])
    charinput<-"C"
    Cnumber<-form.read(formula.input,charinput)
    charinput<-"H"
    Hnumber<-form.read(formula.input,charinput)
    charinput<-"N"
    Nnumber<-form.read(formula.input,charinput)
    charinput<-"O"
    Onumber<-form.read(formula.input,charinput)
    charinput<-"P"
    Pnumber<-form.read(formula.input,charinput)
    charinput<-"S"
    Snumber<-form.read(formula.input,charinput)
    charinput<-"Cl"
    Clnumber<-form.read(formula.input,charinput)
    charinput<-"Br"
    Brnumber<-form.read(formula.input,charinput)
    charinput<-"I"
    Inumber<-form.read(formula.input,charinput)
    formulae.cal<-c(Cnumber,Hnumber,Nnumber,Onumber,Pnumber,Snumber,Clnumber,Brnumber,Inumber)
    formulae.save<-rbind(formulae.save,formulae.cal)
  }
  return(formulae.save)}

form.ReadOne<-function(rawformula){##if it is CH, put it to C1H1 for isopattern function
  formula.save<-NULL
  formula.input<-as.character(rawformula)
  charinput<-"C"
  Cnumber<-form.read(formula.input,charinput)
  if (Cnumber>0){
    formula.save<-paste(formula.save,charinput,Cnumber,sep='')}
  charinput<-"H"
  Hnumber<-form.read(formula.input,charinput)
  if (Hnumber>0){
    formula.save<-paste(formula.save,charinput,Hnumber,sep='')}
  charinput<-"N"
  Nnumber<-form.read(formula.input,charinput)
  if (Nnumber>0){
    formula.save<-paste(formula.save,charinput,Nnumber,sep='')}
  charinput<-"O"
  Onumber<-form.read(formula.input,charinput)
  if (Onumber>0){
    formula.save<-paste(formula.save,charinput,Onumber,sep='')}
  charinput<-"P"
  Pnumber<-form.read(formula.input,charinput)
  if (Pnumber>0){
    formula.save<-paste(formula.save,charinput,Pnumber,sep='')}
  charinput<-"S"
  Snumber<-form.read(formula.input,charinput)
  if (Snumber>0){
    formula.save<-paste(formula.save,charinput,Snumber,sep='')}
  charinput<-"Cl"
  Clnumber<-form.read(formula.input,charinput)
  if (Clnumber>0){
    formula.save<-paste(formula.save,charinput,Clnumber,sep='')}
  charinput<-"Br"
  Brnumber<-form.read(formula.input,charinput)
  if (Brnumber>0){
    formula.save<-paste(formula.save,charinput,Brnumber,sep='')}
  charinput<-"I"
  Inumber<-form.read(formula.input,charinput)
  if (Inumber>0){
    formula.save<-paste(formula.save,charinput,Inumber,sep='')}
  return(formula.save)}

MzCal<-function(formula.input,cutoff){###find out all possible isotopic mz for a given formula
  iso.pattern<-isopattern_cal(iso_list,formula.input,1e-3)
  iso.pattern<-isopattern_Merge(iso.pattern,10^(-5))##merge the isotopic peaks within mass error
  index<-which(iso.pattern[,2]>cutoff)
  return(iso.pattern[index,])##return the mass
}

####This function is used to construct the TargetDatabase for a given formula, to calculate the major m/z of isotopic peaks
TargetConstruct<-function(Library){
  Database<-Library
  MZ.Matrix<-matrix(rep(0,nrow(Library)*6),nrow=nrow(Library),ncol=6)##6 maximal mz
  Database<-cbind(Library,MZ.Matrix)
  for (i in 1:nrow(Library)){
    mz.cal<-MzCal(form.ReadOne(Library[2]),0.1)
    mz.cal<-mz.cal[,1]
    for (j in 1:min(length(mz.cal),6)){
      Database[3+j]<-mz.cal[j]}
  }
  Database[Database==0]<-0.0005485799*polarity
  Database[,4:ncol(Database)]<-Database[,4:ncol(Database)]-0.0005485799*polarity##correct the m/z by loss of addition of electron
  return(Database)                            
}

####This function is used to construct the Decoy database
DecoyConstruct<-function(TargetDatabase,ImAdducts){
  Database<-TargetDatabase
  for (i in 1:nrow(TargetDatabase)){
    mz<-sample(ImAdducts,1,replace=F)##pick the random adduct from the pool
    Database[3]<-mz
    Database[4:ncol(Database)]<-Database[4:ncol(Database)]+mz
  }
  Database[TargetDatabase==0]<-0     
  return(Database)
}

##########DatabaseMatch###############################
DatabaseMatch<-function(Peaks,Database,ppm,msfiles){
  #######match the peak with exact mass########
  Database1<-as.matrix(Database[,4:ncol(Database)])
  score.save<-cbind(Peaks,Peaks[,1:5])
  score.save[,(ncol(score.save)-4):ncol(score.save)]<-0
  for (i in 1:length(msfiles)){
    xrawdata<-xcmsRaw(msfiles[i])
    index.sample<-which(Peaks[,3]==i)##the peak with maximal abundance detected in msfiles i
    if (length(index.sample)<1){next}
    for (k in 1:length(index.sample)){
      index<-which(abs(Peaks[index.sample[k],1]-Database1)<ppm*Peaks[index.sample[k],1])
      if (length(index)>0){
        score.init<--1000000000000
        for (j in 1:length(index)){###for each ID, just pick out the most possible ones
          #####score1, exact mass##############
          delta.mz<-Database1[index[j]]-Peaks[index.sample[k],1]
          score.mz<-exp((-0.5)*(delta.mz/(0.5*10^(-6)*Peaks[index.sample[k],1]))^2)##1ppm for delta
          
          ####score2, retention time##########
          index.data<-index[j]%%(nrow(Database1))
          if (index.data==0){index.data<-nrow(Database)}##if it is 0, it is the last row
          delta.rt<-Database[index.data,1]-Peaks[index.sample[k],2]
          score.rt<-exp((-0.5)*(delta.rt/1)^2)##1min for detal of rt
          
          ###isotopic peaks###################
          formula.mz<-Database[index.data,2]
          isotope<-MzCal(form.ReadOne(formula.mz),0.01)##the isotopic peak for a given formula
          if (Database[index.data,3]>2){##DecoyDatabase
            isotope[,1]<-Database[index.data,3]+isotope[,1]###addition of the adduct to the mz
          }
          RMSE.iso<-IsoscoreCal(xrawdata,Peaks[index.sample[k],2],isotope,Peaks[index.sample[k],1])
          score.iso<-RMSE.iso
          score.total<-score.mz+score.rt+score.iso
          if (score.total>score.init){###just put the best match to the database
            score.init<-score.total
            score.save[index.sample[k],ncol(score.save)]<-score.total
            score.save[index.sample[k],ncol(score.save)-1]<-score.iso
            score.save[index.sample[k],ncol(score.save)-2]<-score.rt
            score.save[index.sample[k],ncol(score.save)-3]<-score.mz
            score.save[index.sample[k],ncol(score.save)-4]<-index.data}}}
    }}
  return(score.save)
}

#############more efficiency to calculate isotope pattern by decreasing the dimension of peaks
isopattern_cal<-function(iso_list, compound, limit){
  isos <- iso_list
  isos <- isos[isos[, 4] != 0, ]
  getit <- seq(1, length(isos[, 1]), 1)
  getthat <- c()
  element <- c()
  number <- c()
  ende <- nchar(compound)
  i <- c(1)
  while (i <= ende) {
    warn <- TRUE
    if (substr(compound,  i) == c("[")) {
      k <- i
      while (any(substr(compound,  i) == c("]")) != TRUE) {
        i <- c(i + 1)}
      while (any(substr(compound,  i) == c("0", "1",
                                             "2", "3", "4", "5", "6", "7", "8", "9")) != TRUE) {
        i <- c(i + 1)}
      m <- c(i - 1)
      element <- c(element, substr(compound, k, m))
      getthat <- c(getthat, getit[isos[, 1] == substr(compound,
                                                      k, m)])
      warn <- FALSE}
    if (any(substr(compound,  i) == c("0", "1", "2", "3",
                                        "4", "5", "6", "7", "8", "9")) != TRUE) {
      k <- i
      while (any(substr(compound,  i) == c("0", "1",
                                             "2", "3", "4", "5", "6", "7", "8", "9")) != TRUE) {
        i <- c(i + 1)}
      m <- c(i - 1)
      i <- c(i - 1)
      element <- c(element, substr(compound, k, m))
      getthat <- c(getthat, getit[isos[, 1] == substr(compound,
                                                      k, m)])
      warn <- FALSE}
    if (any(substr(compound,  i) == c("0", "1", "2", "3",
                                        "4", "5", "6", "7", "8", "9")) == TRUE) {
      k <- i
      while (any(substr(compound,  i) == c("0", "1",
                                             "2", "3", "4", "5", "6", "7", "8", "9")) == TRUE) {
        i <- c(i + 1)}
      m <- c(i - 1)
      i <- i - 1
      number <- c(number, as.numeric(substr(compound, k,
                                            m)))
      warn <- FALSE}
    if (warn == TRUE) {
      stop("Calculation interrupted: compound chemical formula incorrect!")
    }
    i <- i + 1}
  for (i in 1:length(element)) {
    if (any(isos[, 1] == element[i]) == FALSE) {
      stop("Calculation interrupted: one element not found in iso_list")}}
  isos <- t(isos[getthat, ])
  monomass <- as.double(0)
  monoabund <- as.double()
  for (i in 1:length(element)) {
    getit <- isos[, isos[1, ] == element[i]]
    if (length(dim(getit)) > 0) {
      monomass <- c(monomass + (as.double(getit[3, as.numeric(getit[4,
                                                                    ]) == max(as.numeric(getit[4, ]))]) * as.numeric(number[i])))
      monoabund <- c(monoabund, (as.double(getit[4, as.numeric(getit[4,
                                                                     ]) == max(as.numeric(getit[4, ]))])^as.numeric(number[i])))}
    else {
      monomass <- c(monomass + (as.double(getit[3]) * as.numeric(number[i])))
      monoabund <- c(monoabund, (as.double(getit[4])^as.numeric(number[i])))}}
  if (monomass == 0) {
    stop("Calculation interrupted: monoisotopic mass could not be calculated!")}
  if (length(monoabund) == 0) {
    stop("Calculation interrupted: monoisotopic abundance incorrect!")}
  if (length(monoabund) != length(number)) {
    stop("Calculation interrupted: not all elements found in iso_list")}
  getit <- seq(1, length(isos[1, ]), 1)
  road <- matrix(nrow = length(element), ncol = 10, 0)
  rownames(road) = element
  colnames(road) = c("monoiso", rep("nonmono", 9))
  for (i in 1:length(element)) {
    getthat <- getit[isos[1, ] == element[i]]
    road[ 1] = getthat[as.numeric(isos[4, isos[1, ] ==
                                           element[i]]) == max(as.numeric(isos[4, isos[1, ] ==
                                                                                 element[i]]))]
    getthat <- getthat[as.numeric(isos[4, isos[1, ] == element[i]]) !=
                         max(as.numeric(isos[4, isos[1, ] == element[i]]))]
    if (length(getthat) > 0) {
      road[ c(2:(1 + length(getthat)))] = getthat}}
  leng1 <- length(isos[1, ])
  peaks <- matrix(ncol = (4 + leng1), nrow = 5e+03, as.double(0))
  leng2 <- length(peaks[1, ])
  leng3 <- length(peaks[, 1])
  peaks[1, 1] = monomass
  peaks[1, 2] = prod(monoabund)
  peaks[1, 3] = 0
  peaks[1, 4] = 0
  peaks[1, 4 + road[, 1]] = number
  colnames(peaks) = c("mass", "abundance", "generation", "stop?",
                      isos[2, ])
  if (length(monoabund) == 1 && monoabund == 1) {
    peaks[1, 4] = 1}
  k <- c(1)
  m <- c(2)
  b_start <- c(2)
  b_end <- c(1)
  g <- c(1)
  minlimit <- limit
  isonumber <- c()
  getit <- rep(1, length(road[1, ]))
  for (i in 1:length(road[, i])) {
    isonumber <- c(isonumber, length(getit[road[ ] != 0]))}
  for (i in 1:length(road[, 1])) {
    if (isonumber[i] > 1) {
      for (j in 2:isonumber[i]) {
        peaks[m, c(5:(4 + leng1))] = peaks[k, c(5:(4 +
                                                     leng1))]
        peaks[m, 3] = g
        peaks[m, 1] = c(peaks[k, 1] - (as.numeric(isos[3,
                                                       road[ 1]])) + (as.numeric(isos[3, road[
                                                                                                j]])))
        peaks[m, 2] = c(peaks[k, 2] * as.numeric(isos[4,
                                                      road[ j]])/(peaks[m, 4 + road[ j]] + 1) *
                          peaks[m, 4 + road[ 1]]/as.numeric(isos[4,
                                                                   road[ 1]]))
        peaks[m, 4 + road[ 1]] = c(peaks[m, 4 + road[
                                                       1]] - 1)
        peaks[m, 4 + road[ j]] = c(peaks[m, 4 + road[
                                                       j]] + 1)
        if (peaks[m, 2] < minlimit) {
          peaks[m, 4] = 1}
        m <- m + 1}}}
  g <- c(g + 1)
  k <- c(k + 1)
  b_end <- c(m - 1)
  while (any(peaks[peaks[, 3] == (g - 1), 4] == 0) && m < leng3 &&
         b_end >= b_start) {
    while (k <= b_end) {
      if (peaks[k, 2] >= limit) {
        for (i in 1:length(road[, 1])) {
          if (isonumber[i] > 1) {
            if (peaks[k, 4 + road[ 1]] != 0 && m <
                leng3) {
              for (j in 2:isonumber[i]) {
                peaks[m, c(5:(4 + leng1))] = peaks[k,
                                                   c(5:(4 + leng1))]
                peaks[m, 3] = g
                peaks[m, 2] = c(peaks[k, 2] * as.numeric(isos[4,
                                                              road[ j]])/(peaks[m, 4 + road[
                                                                                              j]] + 1) * peaks[m, 4 + road[ 1]]/as.numeric(isos[4,
                                                                                                                                                  road[ 1]]))
                peaks[m, 1] = round(peaks[k, 1] - (as.numeric(isos[3,
                                                                   road[ 1]])) + (as.numeric(isos[3,
                                                                                                    road[ j]])), digits = 9)
                peaks[m, 4 + road[ 1]] = c(peaks[m,
                                                   4 + road[ 1]] - 1)
                peaks[m, 4 + road[ j]] = c(peaks[m,
                                                   4 + road[ j]] + 1)
                if (peaks[m, 2] < minlimit) {
                  peaks[m, 4] = 1}
                m <- m + 1}}}}
        k <- c(k + 1)}
      else {
        k <- c(k + 1)}}
    b_end <- c(m - 1)
    b_start <- c(k)
    if (b_end > b_start) {
      getit <- seq(b_start, b_end, 1)
      getthat <- c()
      back <- getit[order(peaks[getit, 1], decreasing = FALSE)]
      peaksort <- peaks[back, ]
      for (i in 1:(length(peaksort[, 1]) - 1)) {
        if (peaksort[ 1] == peaksort[i + 1, 1]) {
          if (round(peaksort[ 2], digits = 3) == round(peaksort[i +
                                                                  1, 2], digits = 3)) {
            if (all(peaksort[ c(5:leng2)] == peaksort[i +
                                                        1, c(5:leng2)])) {
              getthat <- c(getthat, back[i + 1])}}}}
      leng4 <- length(getthat)
      if (leng4 > 0) {
        peaks <- peaks[-getthat, ]
        m <- c(m - leng4)
        b_end <- c(b_end - leng4)
        leng3 <- c(leng3 - leng4)
        rm(peaksort)}}
    g <- c(g + 1)}
  if (m >= leng3) {
    warning("Storage maximum for number of peaks reached!")}
  rm(road, isonumber, k, b_start, b_end, leng1, leng2, leng3)
  peaks <- peaks[peaks[, 1] != 0, ]
  if (m > 2) {
    peaks <- peaks[peaks[, 2] != 0, ]
    peaks <- peaks[order(peaks[, 1], decreasing = FALSE),
                   ]}
  return(peaks)}

#########Merge the isotope peaks with mass error###############
isopattern_Merge<-function(isotope,ppm){
  if (nrow(isotope)<2){return(isotope)}##only one isotope
  isotope<-isotope[order(isotope[,2],decreasing = TRUE),]
  index.save<-NULL
  for (i in 2:nrow(isotope)){
    index<-which(abs(isotope[1]-isotope[1:(i-1),1])<ppm*isotope[1])###overlapping with previous peaks, delete
    if (length(index)>0){
      index.save<-c(index.save,i)
      isotope[index[1],2]<-isotope[index[1],2]+isotope[2]##addition of the relative abundance
      isotope[1]<-0}}
  if (length(index.save)==0){return(isotope)}    
  return(isotope[-index.save,])}

##############mass calibration using lockmass        
MassCal<-function(xraw,LockMass,input,ppminput){
  ppmid<-ppminput*10^(-6)
  Lockmass<-sort(LockMass)
  Lockmass<-Lockmass-0.0005485799*polarity
  temp.fun<-NULL
  
  #########msn calibration using the average shift########
  ppmshift<-plotlock(xraw,LockMass,10)
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
  msnmzlist<-xraw@env$msnMz
  mzlist<-cbind(msnmzlist,1:length(msnmzlist))
  mzlist<-mzlist[order(mzlist[,1],decreasing=FALSE),]##sort the m/z for prediction
  output<-fitlock(lock.shift,mzlist[,1],temp.fun)
  calmz<-output[[1]]
  temp.fun<-output[[2]]
  mzlist[,1]<-calmz
  mzlist<-mzlist[order(mzlist[,2],decreasing=FALSE),]##sort the m/z for prediction
  xraw@env$msnMz<-mzlist[,1]
  
  ####calculate the primary shift for matching############
  LockMass.cal<-fitlock(lock.shift,lock.shift$lockmass,NULL)
  LockMass.cal<-LockMass.cal[[1]]
  lock.shift<-(lock.shift$lockmass-LockMass.cal)*10^6/LockMass.cal
  
  ##########calibrate ms1 data######################
  for (i in 2:length(xraw@scanindex)){##mass calibration using the closest lockmass
    scanNum<-c(xraw@scanindex[i-1],xraw@scanindex[i])
    correctindex<-(scanNum[1]+1):scanNum[2]
    mzlist<-xraw@env$mz[correctindex]
    
    ###################calculate the mass shift for each lockmass####################
    massshift<-rep(0,length(Lockmass))
    for (k in 1:length(Lockmass)){##define the shift for each lock mass
      mz.index<-which(abs(mzlist*(1-lock.shift*10^(-6))-Lockmass[k])<ppmid*Lockmass[k])
      if (length(mz.index)==0){
        massshift[k]<-100
        next}
      index2<-which.max(xraw@env$intensity[correctindex[mz.index]])##if multiple data exsit, pick out the most abundant one, questionable?
      mz.index<-mz.index[index2]
      mz.index<-mz.index[1]
      mz.index<-scanNum[1]+mz.index
      if (xraw@env$intensity[mz.index]<0){
        massshift[k]<-100
        next}###if the mass lock intensity too low, delete it
      massshift[k]<-(xraw@env$mz[mz.index]-Lockmass[k])/Lockmass[k]}##calculate the mass shift
    
    if (input=='POS'){###positive ion mode calibration
      
      shiftsave<-massshift[1:6]
      for (m in 1:100){
        if (length(shiftsave)<2){break}
        avg<-mean(shiftsave)
        SD<-sd(shiftsave)
        if (SD<1.5*10^(-6)){break}
        index<-which.max(abs(shiftsave))
        shiftsave<-shiftsave[-index]}
      print(c("cal",i))
      index.cor<-which(xraw@env$mz[correctindex]<200)####calibrate the mass below 250
      index.cor1<-which(xraw@env$mz[correctindex]>200)####calibrate the mass below 250
      index.save<-which(xraw@env$mz[correctindex[index.cor1]]<400)
      index.cor1<-index.cor1[index.save]
      index.cor2<-which(xraw@env$mz[correctindex]>400)####calibrate the mass below 250
      if (avg<10^(-5)){
        xraw@env$mz[correctindex[index.cor]]<-mzlist[index.cor]*(1-avg)}
      
      shiftsave<-massshift###for m/z between 200 and 400, use the averaged m/z 
      for (m in 1:100){
        if (length(shiftsave)<2){break}
        avg<-mean(shiftsave)
        SD<-sd(shiftsave)
        if (SD<2*10^(-6)){break}
        index<-which.max(abs(shiftsave))
        shiftsave<-shiftsave[-index]}
      if (avg<10^(-5)){
        xraw@env$mz[correctindex[index.cor1]]<-mzlist[index.cor1]*(1-avg)}
      
      
      shiftsave<-massshift[7:8]######for m/z above 400
      index.del<-which(shiftsave==100)
      if (length(index.del)>0){
        shiftsave<-shiftsave[-index.del]}
      if (length(shiftsave)>0&&length(index.cor2)>0){
        avg2<-mean(shiftsave)
        xraw@env$mz[correctindex[index.cor2]]<-mzlist[index.cor2]*(1-avg2)}
    }   
    
    if (input=='NEG'){
      mzlist<-sort(xraw@env$mz[correctindex])
      Lockmass1<-Lockmass                                     
      index1<-which(massshift==100)
      if (length(index1)>0){
        Lockmass1<-Lockmass1[-index1]
        massshift<-massshift[-index1]} ##Remove the massshift values of '100' from the list (and the corresponding Lockmasses)
      if (length(Lockmass1)==0){next}
      lockdf<-data.frame(lockmass=Lockmass1,ppm=massshift)
      output<-fitlock(lockdf,mzlist,temp.fun)
      calmz<-output[[1]]
      temp.fun<-output[[2]]
      xraw@env$mz[correctindex]<-calmz}
    print(i)
    }
  return(xraw)}

#############function to smooth the mz across three data points############
mzSmooth<-function(Library,msfiles,rtwin,ppmwin){
  ppmwin<-ppmwin*10^(-6)
Library.new<-Library
for (i in 1:length(msfiles)){
  index<-which(Library$SampleID==i)
  if (length(index)==0){next}
  xrawdata<-xcmsRaw(msfiles[i])
  for (j in 1:length(index)){
    temp.index<-index[j]
    mz<-Library$mz[temp.index]
    rt<-Library$rt[temp.index]*60
    
    #######smooth the exact mass and retention time############
    mz.min<-mz*(1-ppmwin)
    mz.min<-max(mz.min, xrawdata@mzrange[1])
    mz.max<-mz*(1+ppmwin)
    mz.max<-min(mz.max,xrawdata@mzrange[2])
    rt.min<-max(min(xrawdata@scantime),rt-rtwin) #range of analysis is the entire chromatograph
    rt.max<-min(max(xrawdata@scantime),rt+rtwin) #the first and last 3 minutes are skipped to avoid sections with no analytes
    peaks<-rawEIC(xrawdata,mzrange=cbind(mz.min,mz.max),rtrange=cbind(rt.min, rt.max))
    scan.max<-which.max(peaks$intensity)
    scan.max<-peaks$scan[scan.max[1]]
    if (max(peaks$intensity)<1000){
      Library.new$mz[temp.index]<-NA
      next}##delete the data if the compound was detected in the last scan
    if (scan.max>=(length(xrawdata@scanindex)-2)){
      Library.new$mz[temp.index]<-NA
      next}##delete the data if the compound was detected in the last scan
    scanNum<-c(xrawdata@scanindex[scan.max],xrawdata@scanindex[scan.max+1])
    correctindex<-(scanNum[1]+1):scanNum[2]
    mzlist<-xrawdata@env$mz[correctindex]
    intensity<-xrawdata@env$intensity[correctindex]
    index.mz<-which(abs(mzlist-mz)<ppmwin*mz)
    mz.scan1<-mzlist[index.mz[1]]###the mz of the maxiaml peak
    
    ####the mz for scaning point -1 to peak top######################
    if (scan.max==1){
      mz.scan2<-mz.scan1}else{
        scanNum<-c(xrawdata@scanindex[scan.max-1],xrawdata@scanindex[scan.max])
        correctindex<-(scanNum[1]+1):scanNum[2]
        mzlist<-xrawdata@env$mz[correctindex]
        intensity<-xrawdata@env$intensity[correctindex]
        index.mz<-which(abs(mzlist-mz)<ppmwin*mz)
        mz.scan2<-mzlist[index.mz[1]]}###the mz of the maxiaml peak
    
    ####the mz for scaning point 1 to peak top######################
    scanNum<-c(xrawdata@scanindex[scan.max+1],xrawdata@scanindex[scan.max+2])
    correctindex<-(scanNum[1]+1):scanNum[2]
    mzlist<-xrawdata@env$mz[correctindex]
    intensity<-xrawdata@env$intensity[correctindex]
    index.mz<-which(abs(mzlist-mz)<ppmwin*mz)
    mz.scan3<-mzlist[index.mz[1]]###the mz of the maxiaml peak
    
    #######replace the mz and rt with newly smoothed results    
    Library.new$mz[temp.index]<-sum(mz.scan1,mz.scan2,mz.scan3)/3
    Library.new$rt[temp.index]<-xrawdata@scantime[scan.max]/60
  }}
return(Library.new)}


####small function to find the isotope series#######
FindHomo<-function(mz1,mz2){
  delta.mz<-abs(mz1-mz2)
  for (x in 0:10){##chlorine
    for (y in 0:10){##bromine
      for (z in 0:2){###carbon
        delta<-x*1.99795+y*1.99705+1.00335*z
        if (abs(delta.mz-delta)<0.003){return(1)} ###<0.003,H and carbon     
        }}}
  return(0)
  }

FindCarbon<-function(mz1,mz2,deltamz){##check if two masses are differentiated by carbon isotope
  delta.mz<-mz1-mz2
  for (x in 1:5){##carbon
        delta<-deltamz*x
        if (abs(delta.mz-delta)<0.003){return(x)}      
      }
  return(0)
}


IsotopeFind<-function(Library.new,path.out,mz_tol){####delete the secondary isotopic peaks
  setwd(path.out)
  msfiles<-list.files()
  mz_tol<-3*10^(-6)
  index.save<-NULL
  IDsave<-0
  List.Isotope<-NULL
  for (i in 1:length(msfiles)){
    print(c('MS file...',i))
    xrawdata<-xcmsRaw(msfiles[i])
    index<-which(Library.new$SampleID==i)
    if (length(index)<1){next}
    for (j in 1:length(index)){
      temp<-index[j]
      mz<-Library.new$mz[temp]
      rt<-Library.new$rt[temp]*60
    mzrange<-xrawdata@mzrange
    minmz<-mzrange[1]
    maxmz<-mzrange[2]
    mzmin<-max(minmz,mz-mz*mz_tol)
    mzmax<-min(maxmz,mz+mz*mz_tol)
    rtrange<-xrawdata@scantime
    rtmin<-max(min(rtrange),rt-30)
    rtmax<-min(max(rtrange),rt+30)
    if (rtmax<0){next}
    peaks<-rawEIC(xrawdata,mzrange=cbind(mzmin,mzmax),rtrange=cbind(rtmin,rtmax))
    scan.max<-which.max(peaks$intensity)
    scan.max<-peaks$scan[scan.max[1]]
    scanNum<-c(xrawdata@scanindex[scan.max],xrawdata@scanindex[scan.max+1])
    correctindex<-(scanNum[1]+1):scanNum[2]
    mzlist<-xrawdata@env$mz[correctindex]
    mz.iso<-which(abs(mzlist-mz)<13)###range for isotope peak finding
    mz.iso<-mzlist[mz.iso]
    if (length(mz.iso)<1){next}
    
    
    ######only consider carbon, Br, and Cl isotope########
    save.iso<-NULL
    for (mm in 1:length(mz.iso)){
      match.iso<-FindHomo(mz.iso[mm],mz)
      if(match.iso==1){save.iso<-c(save.iso,mm)}##it is the isotope
    }
    if (length(save.iso)<1){next}
    mz.iso<-mz.iso[save.iso]
    
    native.peak<-peaks$intensity
    kk<-0
    for (k in 1:length(mz.iso)){
      mz0<-mz.iso[k]
      mzmin<-max(minmz,mz0-mz0*mz_tol)
      mzmax<-min(maxmz,mz0+mz0*mz_tol)
      isotope.peak<-rawEIC(xrawdata,mzrange=cbind(mzmin,mzmax),rtrange=cbind(rtmin,rtmax))
      isotope.peak<-isotope.peak$intensity
      if (sd(isotope.peak)==0||sd(native.peak)==0){next}###0 values
        corr<-cor(native.peak,isotope.peak)
        if (corr>0.64){
          List.Isotope$mz<-c(List.Isotope$mz,mz)
          List.Isotope$sampleID<-c(List.Isotope$sampleID,i)
          List.Isotope$Isotope<-c(List.Isotope$Isotope,mz0)
          List.Isotope$rt<-c(List.Isotope$rt,rt/60)
          List.Isotope$cor<-c(List.Isotope$cor,corr)
          List.Isotope$intensity.iso<-c(List.Isotope$intensity.iso,max(isotope.peak))
          List.Isotope$intensity<-c(List.Isotope$intensity,max(native.peak))
          List.Isotope$ID<-c(List.Isotope$ID,temp)
        }###indicate primary isotopic peaks detected
      }}}
  return(List.Isotope)}

###############function to remove isotope noises###################
IsotopeAcross<-function(path,msfiles,Isotope.Data,mz_tol,rtwin){
  setwd(path)
  if (length(msfiles)==1){return(Isotope.Data)}
  
  iso.save<-matrix(rep(0,length(Isotope.Data$Isotope)*length(msfiles)),ncol=length(msfiles))
  prec.save<-matrix(rep(0,length(Isotope.Data$mz)*length(msfiles)),ncol=length(msfiles))
  for (i in 1:length(msfiles)){
    print(c('checkisotope...',i))
    xraw<-xcmsRaw(msfiles[i])
    
    for (j in 1:length(Isotope.Data$mz)){
      mz.iso<-Isotope.Data$Isotope[j]###fragments
      mz.prec<-Isotope.Data$mz[j]###precursor
      
      #######setting up searching############
      mzrange<-xraw@mzrange
      minmz<-mzrange[1]
      maxmz<-mzrange[2]
      rt<-Isotope.Data$rt[j]*60
      mzmin<-max(minmz,mz.prec-mz.prec*mz_tol)
      mzmax<-min(maxmz,mz.prec+mz.prec*mz_tol)
      rtrange<-xraw@scantime
      rtmin<-max(min(rtrange),rt-rtwin)
      rtmax<-min(max(rtrange),rt+rtwin)
      
      #######intensity for precursor########
      peaks.prec<-rawEIC(xraw,mzrange=cbind(mzmin,mzmax),rtrange=cbind(rtmin,rtmax))
      scan.max<-which.max(peaks.prec$intensity)
      prec.save[j,i]<-max(peaks.prec$intensity)
      
      #######intensity for fragments########    
      mzmin<-max(minmz,mz.iso-mz.iso*mz_tol)
      mzmax<-min(maxmz,mz.iso+mz.iso*mz_tol)
      iso.peak<-rawEIC(xraw,mzrange=cbind(mzmin,mzmax),rtrange=cbind(rtmin,rtmax))
      iso.save[j,i]<-iso.peak$intensity[scan.max]
    }}
  
  
  #########correlations##################
  saveid<-NULL
  for (k in 1:nrow(iso.save)){
    iso.peak<-iso.save[k,]
    prec.peak<-prec.save[k,]
    if (sd(prec.peak)==0||sd(iso.peak)==0){next}###0 values
    corr<-cor(prec.peak,iso.peak)
    if (corr>0.85){saveid<-c(saveid,k)}
  }
  if (length(saveid)==0){return(NULL)}
  
  prec_list<-Isotope.Data
  Frag_list<-list()
  Frag_list$Isotope<-prec_list$Isotope[saveid]
  Frag_list$mz<-prec_list$mz[saveid]
  Frag_list$rt<-prec_list$rt[saveid]
  Frag_list$intensity<-prec_list$intensity[saveid]
  Frag_list$cor<-prec_list$cor[saveid]
  Frag_list$sampleID<-prec_list$sampleID[saveid]
  Frag_list$ID<-prec_list$ID[saveid]
  Frag_list$intensity.iso<-prec_list$intensity.iso[saveid]
  return (Frag_list)
  }

##############################find adducts###################################
AdductsFind<-function(Library.new,path.out,mw.adducts,mz_tol,Adducts){#find the adducts
  setwd(path.out)
  msfiles<-list.files()
  index.save<-NULL
  IDsave<-0
  List.Isotope<-NULL
  for (i in 1:length(msfiles)){
    print(c('MS file...',i))
    xrawdata<-xcmsRaw(msfiles[i])
    index<-which(Library.new$SampleID==i)
    if (length(index)<1){next}
    for (j in 1:length(index)){
      temp<-index[j]
      mz<-Library.new$mz[temp]
      rt<-Library.new$rt[temp]*60
      mzrange<-xrawdata@mzrange
      minmz<-mzrange[1]
      maxmz<-mzrange[2]
      mzmin<-max(minmz,mz-mz*mz_tol)
      mzmax<-min(maxmz,mz+mz*mz_tol)
      rtrange<-xrawdata@scantime
      rtmin<-max(min(rtrange),rt-30)
      rtmax<-min(max(rtrange),rt+30)
      if (rtmax<0){next}
      peaks<-rawEIC(xrawdata,mzrange=cbind(mzmin,mzmax),rtrange=cbind(rtmin,rtmax))
      scan.max<-which.max(peaks$intensity)
      scan.max<-peaks$scan[scan.max[1]]
      scanNum<-c(xrawdata@scanindex[scan.max],xrawdata@scanindex[scan.max+1])
      correctindex<-(scanNum[1]+1):scanNum[2]
      mz.iso<-xrawdata@env$mz[correctindex]
      if (length(mz.iso)<1){next}
      
      
      ######Find adducts########
      save.iso<-data.frame()
      mzadducts<-NULL
      adductssave<-NULL
      for (mm in 1:length(mz.iso)){
        for (kk in 2:length(mw.adducts)){##skip M and M-H
          mzdiff<-mw.adducts[kk]-mw.adducts[2]###difference between M-H or M+H
          if (abs(mz-mz.iso[mm]-mzdiff)<mz*mz_tol){
            mzadducts<-c(mzadducts,mz.iso[mm])
            adductssave<-c(adductssave,kk)
          }}}
      if (length(mzadducts)<1){next}
      save.iso<-data.frame(mz=mzadducts,adducts=Adducts[adductssave])
      mz.iso<-save.iso$mz
      
      native.peak<-peaks$intensity
      kk<-0
      for (k in 1:length(mz.iso)){
        mz0<-mz.iso[k]
        mzmin<-max(minmz,mz0-mz0*mz_tol)
        mzmax<-min(maxmz,mz0+mz0*mz_tol)
        isotope.peak<-rawEIC(xrawdata,mzrange=cbind(mzmin,mzmax),rtrange=cbind(rtmin,rtmax))
        isotope.peak<-isotope.peak$intensity
        if (sd(isotope.peak)==0||sd(native.peak)==0){next}###0 values
        corr<-cor(native.peak,isotope.peak)
        if (corr>0.64){
          List.Isotope$mzadducts<-c(List.Isotope$mzadducts,mz)
          List.Isotope$sampleID<-c(List.Isotope$sampleID,i)
          List.Isotope$mzMH<-c(List.Isotope$mzMH,mz0)
          List.Isotope$rt<-c(List.Isotope$rt,rt/60)
          List.Isotope$cor<-c(List.Isotope$cor,corr)
          List.Isotope$intensityMH<-c(List.Isotope$intensityMH,max(isotope.peak))
          List.Isotope$intensityadduct<-c(List.Isotope$intensityadduct,max(native.peak))
          List.Isotope$ID<-c(List.Isotope$ID,temp)
          List.Isotope$Adducts<-c(List.Isotope$Adducts,as.character(save.iso$adducts[k]))
        }###indicate primary isotopic peaks detected
      }}}
  return(List.Isotope)}


#################Rcpp function to enhance computation efficiency###############
cppFunction(
  'NumericMatrix itercal(NumericMatrix numberset, NumericVector mz_list,double ppm, double mz, double mwoffset){
  NumericMatrix output(500,12);
  int i1;
  int i2;
  int i3;
  int i4;
  int i5;
  int i6;
  int i7;
  int i8;
  int i9;
  int i10;
  int i11;
  int i12; 
  int kk=0;
  double temp1;
  double temp2;
  double temp3;
  double temp4;
  double temp;
  double value=100;
  double RDBE;
  for (i1=numberset(0,0);i1<=numberset(0,1);i1++){
  for (i2=numberset(1,0);i2<=(numberset(1,1)<(i1*2+3)?numberset(1,1):(i1*2+3));i2++){
  for (i3=numberset(2,0);i3<=(numberset(2,1)<(i1*1.3)?numberset(2,1):(i1*1.3));i3++){
  for (i4=numberset(3,0);i4<=(numberset(3,1)<(i1*1.2+2)?numberset(3,1):(i1*1.2+2));i4++){
  for (i5=numberset(4,0);i5<=(numberset(4,1)<(i1*0.3)?numberset(4,1):(i1*0.3));i5++){
  for (i6=numberset(5,0);i6<=numberset(5,1);i6++){
  for (i7=numberset(6,0);i7<=numberset(6,1);i7++){
  temp1=i1*mz_list[0]+i2*mz_list[1]+i3*mz_list[2]+i4*mz_list[3]+i5*mz_list[4]+i6*mz_list[5]+i7*mz_list[6];
  if (temp1>(mz+1)){
  break;}
  for (i8=numberset(7,0);i8<=numberset(7,1);i8++){
  temp2=temp1+i8*mz_list[7];
  if (temp2>(mz+1)){
  break;}
  for (i9=numberset(8,0);i9<=numberset(8,1);i9++){
  temp3=temp2+i9*mz_list[8];
  if (temp3>(mz+1)){
  break;}
  for (i10=numberset(9,0);i10<=numberset(9,1);i10++){
  temp4=temp3+i10*mz_list[9];
  if (temp4>(mz+1)){
  break;}
  for (i12=0;i12<=1;i12++){
  for (i11=numberset(10,0);i11<=numberset(10,1);i11++){
  temp=temp4+i11*mz_list[10]+i12*mwoffset;
  value=abs(1000000*(mz-temp)/temp);
  RDBE=i1-0.5/(i2+i7+i8+i9+i10+i11+0.0001)+0.5/(i3+i5+0.5)+1;
  if (value<ppm&&RDBE>=0&&RDBE<40){
  NumericVector out=NumericVector::create(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,value);
  output.row(kk)=out;
  kk=kk+1;
  if (kk>=500){return(output(Range(0,kk-1),Range(0,11)));}
  }
  if (temp>(mz+1)){break;}
  }}}}}}}}}}}}
  if (kk==0){
  kk=1;}
  return(output(Range(0,kk-1),Range(0,11)));
  }')

######################function to calculate formula##############
formpred<-function(mz,ppm,numberset,mz_list,number_list,mwoffset){
  numberset<-matrix(as.numeric(numberset[,1:2]),ncol=2,nrow=nrow(numberset))
  numberset[2,1]<-max(numberset[2,1],floor(numberset[1,1]*0.5-sum(numberset[7:11,2])))##the minimum number of bromine
  formulae<-rep(0,nrow(numberset))
  for.pred<-itercal(numberset,mz_list,ppm,mz,mwoffset)##mwoffset is for prediction of rare element for fragment
  if (length(for.pred)<2)
    for.pred<-NULL
  return(for.pred)
}


iso.otheratom<-rbind(1:6,1:6)###save the isotopic distribution of other atoms
iso.otheratom[1,1]<-iso_list[4,3]-iso_list[3,3]##mz difference of C
iso.otheratom[2,1]<-iso_list[4,4]/iso_list[3,4]##relative abundance of C
iso.otheratom[1,2]<-iso_list[2,3]-iso_list[1,3]##mz difference of H
iso.otheratom[2,2]<-iso_list[2,4]/iso_list[1,4]##relative abundance of H
iso.otheratom[1,3]<-iso_list[6,3]-iso_list[5,3]##mz difference of N
iso.otheratom[2,3]<-iso_list[6,4]/iso_list[5,4]##relative abundance of N
iso.otheratom[1,4]<-iso_list[9,3]-iso_list[7,3]##mz difference of O
iso.otheratom[2,4]<-iso_list[9,4]/iso_list[7,4]##relative abundance of O
iso.otheratom[1,5]<-0##mz difference of P
iso.otheratom[2,5]<-0##relative abundance of P
iso.otheratom[1,6]<-iso_list[12,3]-iso_list[10,3]##mz difference of S
iso.otheratom[2,6]<-iso_list[12,4]/iso_list[10,4]##relative abundance of S

#############more efficiency to calculate isotope pattern by decreasing the dimension of peaks
isopattern_cal<-function(iso_list, compound, limit){
  isos <- iso_list
  isos <- isos[isos[, 4] != 0, ]
  getit <- seq(1, length(isos[, 1]), 1)
  getthat <- c()
  element <- c()
  number <- c()
  ende <- nchar(compound)
  i <- c(1)
  while (i <= ende) {
    warn <- TRUE
    if (substr(compound, i, i) == c("[")) {
      k <- i
      while (any(substr(compound, i, i) == c("]")) != TRUE) {
        i <- c(i + 1)}
      while (any(substr(compound, i, i) == c("0", "1",
                                             "2", "3", "4", "5", "6", "7", "8", "9")) != TRUE) {
        i <- c(i + 1)}
      m <- c(i - 1)
      element <- c(element, substr(compound, k, m))
      getthat <- c(getthat, getit[isos[, 1] == substr(compound,
                                                      k, m)])
      warn <- FALSE}
    if (any(substr(compound, i, i) == c("0", "1", "2", "3",
                                        "4", "5", "6", "7", "8", "9")) != TRUE) {
      k <- i
      while (any(substr(compound, i, i) == c("0", "1",
                                             "2", "3", "4", "5", "6", "7", "8", "9")) != TRUE) {
        i <- c(i + 1)}
      m <- c(i - 1)
      i <- c(i - 1)
      element <- c(element, substr(compound, k, m))
      getthat <- c(getthat, getit[isos[, 1] == substr(compound,
                                                      k, m)])
      warn <- FALSE}
    if (any(substr(compound, i, i) == c("0", "1", "2", "3",
                                        "4", "5", "6", "7", "8", "9")) == TRUE) {
      k <- i
      while (any(substr(compound, i, i) == c("0", "1",
                                             "2", "3", "4", "5", "6", "7", "8", "9")) == TRUE) {
        i <- c(i + 1)}
      m <- c(i - 1)
      i <- i - 1
      number <- c(number, as.numeric(substr(compound, k,
                                            m)))
      warn <- FALSE}
    if (warn == TRUE) {
      stop("Calculation interrupted: compound chemical formula incorrect!")
    }
    i <- i + 1}
  for (i in 1:length(element)) {
    if (any(isos[, 1] == element[i]) == FALSE) {
      stop("Calculation interrupted: one element not found in iso_list")}}
  isos <- t(isos[getthat, ])
  monomass <- as.double(0)
  monoabund <- as.double()
  for (i in 1:length(element)) {
    getit <- isos[, isos[1, ] == element[i]]
    if (length(dim(getit)) > 0) {
      monomass <- c(monomass + (as.double(getit[3, as.numeric(getit[4,
                                                                    ]) == max(as.numeric(getit[4, ]))]) * as.numeric(number[i])))
      monoabund <- c(monoabund, (as.double(getit[4, as.numeric(getit[4,
                                                                     ]) == max(as.numeric(getit[4, ]))])^as.numeric(number[i])))}
    else {
      monomass <- c(monomass + (as.double(getit[3]) * as.numeric(number[i])))
      monoabund <- c(monoabund, (as.double(getit[4])^as.numeric(number[i])))}}
  if (monomass == 0) {
    stop("Calculation interrupted: monoisotopic mass could not be calculated!")}
  if (length(monoabund) == 0) {
    stop("Calculation interrupted: monoisotopic abundance incorrect!")}
  if (length(monoabund) != length(number)) {
    stop("Calculation interrupted: not all elements found in iso_list")}
  getit <- seq(1, length(isos[1, ]), 1)
  road <- matrix(nrow = length(element), ncol = 10, 0)
  rownames(road) = element
  colnames(road) = c("monoiso", rep("nonmono", 9))
  for (i in 1:length(element)) {
    getthat <- getit[isos[1, ] == element[i]]
    road[i, 1] = getthat[as.numeric(isos[4, isos[1, ] ==
                                           element[i]]) == max(as.numeric(isos[4, isos[1, ] ==
                                                                                 element[i]]))]
    getthat <- getthat[as.numeric(isos[4, isos[1, ] == element[i]]) !=
                         max(as.numeric(isos[4, isos[1, ] == element[i]]))]
    if (length(getthat) > 0) {
      road[i, c(2:(1 + length(getthat)))] = getthat}}
  leng1 <- length(isos[1, ])
  peaks <- matrix(ncol = (4 + leng1), nrow = 5e+03, as.double(0))
  leng2 <- length(peaks[1, ])
  leng3 <- length(peaks[, 1])
  peaks[1, 1] = monomass
  peaks[1, 2] = prod(monoabund)
  peaks[1, 3] = 0
  peaks[1, 4] = 0
  peaks[1, 4 + road[, 1]] = number
  colnames(peaks) = c("mass", "abundance", "generation", "stop?",
                      isos[2, ])
  if (length(monoabund) == 1 && monoabund == 1) {
    peaks[1, 4] = 1}
  k <- c(1)
  m <- c(2)
  b_start <- c(2)
  b_end <- c(1)
  g <- c(1)
  minlimit <- limit
  isonumber <- c()
  getit <- rep(1, length(road[1, ]))
  for (i in 1:length(road[, i])) {
    isonumber <- c(isonumber, length(getit[road[i, ] != 0]))}
  for (i in 1:length(road[, 1])) {
    if (isonumber[i] > 1) {
      for (j in 2:isonumber[i]) {
        peaks[m, c(5:(4 + leng1))] = peaks[k, c(5:(4 +
                                                     leng1))]
        peaks[m, 3] = g
        peaks[m, 1] = c(peaks[k, 1] - (as.numeric(isos[3,
                                                       road[i, 1]])) + (as.numeric(isos[3, road[i,
                                                                                                j]])))
        peaks[m, 2] = c(peaks[k, 2] * as.numeric(isos[4,
                                                      road[i, j]])/(peaks[m, 4 + road[i, j]] + 1) *
                          peaks[m, 4 + road[i, 1]]/as.numeric(isos[4,
                                                                   road[i, 1]]))
        peaks[m, 4 + road[i, 1]] = c(peaks[m, 4 + road[i,
                                                       1]] - 1)
        peaks[m, 4 + road[i, j]] = c(peaks[m, 4 + road[i,
                                                       j]] + 1)
        if (peaks[m, 2] < minlimit) {
          peaks[m, 4] = 1}
        m <- m + 1}}}
  g <- c(g + 1)
  k <- c(k + 1)
  b_end <- c(m - 1)
  while (any(peaks[peaks[, 3] == (g - 1), 4] == 0) && m < leng3 &&
         b_end >= b_start) {
    while (k <= b_end) {
      if (peaks[k, 2] >= limit) {
        for (i in 1:length(road[, 1])) {
          if (isonumber[i] > 1) {
            if (peaks[k, 4 + road[i, 1]] != 0 && m <
                leng3) {
              for (j in 2:isonumber[i]) {
                peaks[m, c(5:(4 + leng1))] = peaks[k,
                                                   c(5:(4 + leng1))]
                peaks[m, 3] = g
                peaks[m, 2] = c(peaks[k, 2] * as.numeric(isos[4,
                                                              road[i, j]])/(peaks[m, 4 + road[i,
                                                                                              j]] + 1) * peaks[m, 4 + road[i, 1]]/as.numeric(isos[4,
                                                                                                                                                  road[i, 1]]))
                peaks[m, 1] = round(peaks[k, 1] - (as.numeric(isos[3,
                                                                   road[i, 1]])) + (as.numeric(isos[3,
                                                                                                    road[i, j]])), digits = 9)
                peaks[m, 4 + road[i, 1]] = c(peaks[m,
                                                   4 + road[i, 1]] - 1)
                peaks[m, 4 + road[i, j]] = c(peaks[m,
                                                   4 + road[i, j]] + 1)
                if (peaks[m, 2] < minlimit) {
                  peaks[m, 4] = 1}
                m <- m + 1}}}}
        k <- c(k + 1)}
      else {
        k <- c(k + 1)}}
    b_end <- c(m - 1)
    b_start <- c(k)
    if (b_end > b_start) {
      getit <- seq(b_start, b_end, 1)
      getthat <- c()
      back <- getit[order(peaks[getit, 1], decreasing = FALSE)]
      peaksort <- peaks[back, ]
      for (i in 1:(length(peaksort[, 1]) - 1)) {
        if (peaksort[i, 1] == peaksort[i + 1, 1]) {
          if (round(peaksort[i, 2], digits = 3) == round(peaksort[i +
                                                                  1, 2], digits = 3)) {
            if (all(peaksort[i, c(5:leng2)] == peaksort[i +
                                                        1, c(5:leng2)])) {
              getthat <- c(getthat, back[i + 1])}}}}
      leng4 <- length(getthat)
      if (leng4 > 0) {
        peaks <- peaks[-getthat, ]
        m <- c(m - leng4)
        b_end <- c(b_end - leng4)
        leng3 <- c(leng3 - leng4)
        rm(peaksort)}}
    g <- c(g + 1)}
  if (m >= leng3) {
    warning("Storage maximum for number of peaks reached!")}
  rm(road, isonumber, k, b_start, b_end, leng1, leng2, leng3)
  peaks <- peaks[peaks[, 1] != 0, ]
  if (m > 2) {
    peaks <- peaks[peaks[, 2] != 0, ]
    peaks <- peaks[order(peaks[, 1], decreasing = FALSE),
                   ]}
  return(peaks)}


####function to group isotope value for selcted element######
iso.element<-function(element,iso_list,input){
  index<-NULL
  for (i in 1:length(element)){
    indextemp<-which(iso_list[,input]==element[i])
    if (length(indextemp)<1){
      indextemp<-which(iso_list[,2]==element[i])}##the isotope of Br and Cl
    indextemp<-min(indextemp)
    index<-c(index,indextemp)}
  return(iso_list[index,])}


###########function to read formulae
form.comb<-function(formulainput){
  formulainput[7]<-sum(formulainput[7:8])
  formulainput[8]<-sum(formulainput[9:10])
  formulainput[9]<-formulainput[11]
  charinput<-c("C","H","N","O","P","S","Cl","Br","I")
  formulainput<-as.character(formulainput)
  form.output<-NULL
  for (i in 1:length(charinput)){
    if (formulainput[i]!=0){
      form.output<-paste(form.output,charinput[i],formulainput[i],sep="")}}
  return(form.output)}

################funtion to calculate isotope distribution of proposed formulae########
deiso.formula<-function(vector1,iso_list,index.input,Peakinfo,Isotope.Data){
  mz_tol<-5*10^(-6)
  formula.cal<-form.comb(vector1[1:11])
  if (max(vector1)==0) {return(10)}
  iso.pattern<-isopattern_cal(iso_list,formula.cal,1e-4)
  iso.pattern<-iso.pattern[,1:2]
  iso.pattern[,2]<-iso.pattern[,2]/(max(iso.pattern[,2]))
  iso.pattern<-iso.pattern[which(iso.pattern[,2]>0.01),]
  iso.pattern<-combine.mz(iso.pattern,0.002)###combine those peaks could not be distinguished by mass spec
  iso.pattern[,1]<-iso.pattern[,1]-0.0005485799*polarity+mwoffset##plus the decoy adducts

  
  index<-which(Isotope.Data$ID==index.input)
  if(length(index)==0)(return(10))
  
  iso.intensity<-Isotope.Data$intensity.iso[index]###isotopic peaks from real spectra
  iso.mz<-Isotope.Data$Isotope[index]
  mz<-Peakinfo$mz[index.input]
  index<-which(abs(iso.pattern[,1]-mz)<mz_tol*mz)
  if (length(index)==0){return(10)}
  normal.ratio<-iso.pattern[index[1],2]###the ratio for mz
  index<-which(abs(iso.mz-mz)<mz_tol*mz)
  if (length(index)==0){return(10)}
  iso.intensity<-iso.intensity*normal.ratio/iso.intensity[index[1]]

  
  match.1<-data.frame(iso.pattern)###match expt spectra to predicted spectra
  match.1$match<-rep(0,nrow(iso.pattern))
  for (i in 1:nrow(iso.pattern)){
    index<-which(abs(match.1$mass[i]-iso.mz)<mz_tol*match.1$mass[i])
    if (length(index)>0){
      match.1$match[i]<-iso.intensity[index[1]]
    }
  }
  
  match.2<-data.frame(mz=iso.mz,intensity=iso.intensity)###match expt spectra to predicted spectra
  match.2$match<-rep(0,nrow(match.2))
  for (i in 1:nrow(match.2)){
    index<-which(abs(match.2$mz[i]-match.1$mass)<mz_tol*match.2$mz[i])
    if (length(index)>0){
      match.2$match[i]<-match.1$abundance[index[1]]
    }
  }
  
  RMSE.1<-sum((match.1[,2]-match.1[,3])^2)
  RMSE.2<-sum((match.2[,2]-match.2[,3])^2)
  return(RMSE.1+RMSE.2)
}


##############function to combine mz using mz_tol###########
combine.mz<-function(mz.list,mz.tol){
  mz.new<-NULL
  for (i in 1:nrow(mz.list)){
    if (i==1){
      mz.new<-mz.list[1,]
      }else{
        if (length(mz.new)==ncol(mz.list)){
          index<-which(abs(mz.new-mz.list[i,1])<mz.tol)
          }else{
            index<-which(abs(mz.new[,1]-mz.list[i,1])<mz.tol)}
        if (length(index)<1){
          mz.new<-rbind(mz.new,mz.list[i,])
          }else{
            if (length(mz.new)==ncol(mz.list)){
              mz.value<-sum(mz.new[1],mz.list[i,1])/2
              intensity<-sum(mz.new[2],mz.list[i,2])
              mz.new<-c(mz.value,intensity)
            }
            else{
              mz.value<-sum(mz.new[index,1],mz.list[i,1])/2
              intensity<-sum(mz.new[index,2],mz.list[i,2])
              mz.new[index[1],]<-c(mz.value,intensity)}
          }}}
  return(mz.new)}


#####################calculate the score for a given isotope distribution##############
#####################for actural library, adjusted by the predicted mz at first, and then normalized by summed value from the library##########
######for the first adjustment, if use maxmal peak or summed peak for adjustment, there may be big effects from noise########
######for the second adjustment, if use maximal peak for adjustment, there will be bias towards low isotopes###########
isotope.score<-function(iso.pattern,iso.info){##iso.pattern and iso.info is the library and small number of peaks, 
  index1<-ceiling(length(iso.info)/4)
  index<-which(abs(iso.pattern[,1]-iso.info[index1,1])<0.006)#index is the ID of mz for alignment
  factor<-iso.pattern[index,2]/sum(iso.pattern[,2])
  iso.pattern[,2]<-iso.pattern[,2]/sum(iso.pattern[,2])###this step is to use sum value for adjustment to avoid any bias
  iso.info[,2]<-iso.info[,2]*factor/(iso.info[index1,2])###normalized based on known mz
  RMSE<-0
  index.act<-rep(0,length(iso.info)/2)
  for (i in 1:(length(iso.pattern)/2)){##calculate the RMSE from library
    index.temp<-which(abs(iso.info[,1]-iso.pattern[1])<0.006)
    if (length(index.temp)>0){
      RMSE.temp<-(iso.pattern[2]-iso.info[index.temp,2])^2
      index.act[index.temp]<-i}else{##record the actual peaks have been calcualted
        RMSE.temp<-(iso.pattern[2])^2}
    RMSE<-RMSE+RMSE.temp}
  for (i in 1:(length(iso.info)/2)){##calculate the RMSE from unaligned actual peaks
    if (index.act[i]==0){##not aligned to library
      RMSE<-RMSE+(iso.info[2])^2}}
  return(RMSE)}


####################calculate probability of distribution of a given formulae based on prior information
density.form<-function(vector1){##vector1 is the number of element, C, H, O, N, P, S, Cl, 37Cl, Br, 81Br
  Hnumber<-vector1[2]+sum(vector1[7:11])
  Hnumber<-Hnumber/vector1[1]
  vector1<-vector1/vector1[1]
  prob<-exp((-0.5)*(Hnumber-1.6)^2/1.6^2)+exp((-0.5)*(vector[5]-0)^2/0.15^2)+exp((-0.5)*(vector[5]-0.1)^2/0.4^2)
}


#################formula calculation based on exact mass, isotope distribution#################
formcal<-function(number_input,Peakinfo,index.input,iso_list,rawdata,ppm,ppm.ms2,Isotope.Data,Database,Fragment,adducts,IsotopeDB,mwoffset){##number input is the limit of element composition number, isolist is the mw of element
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
    for (kk in 1:ncol(formula.split)){#Check if the number of elment is lesser than -1
      if (min(formula.split[,kk])<0){
        index.save<-c(index.save,kk)
      }
    }
    if (length(index.save)>0){
      formula.split<-formula.split[,-index.save]
    }
    formula.split<-matrix(formula.split,nrow=11)
    if (length(formula.split)==0){return(NULL)}
    RMSE.iso1<-apply(formula.split,2,deiso.formula,iso_list,index.input,Peakinfo,Isotope.Data)
  }
  if(RMSE.iso1[1]==0){return(NULL)}###error in isotope calculation because of rare element
  RMSE.iso<-exp((-0.5)*RMSE.iso1)
  
  ###############calculate the ms1 error
  MS1.score<-rep(0,length(mserror))
  index.small<-which(abs(mserror)<1)
  if (length(index.small)>0){MS1.score[index.small]<-1}##when smaller than 1ppm, no difference
  index.big<-which(abs(mserror)>1)
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
    formula.final<-cbind(formula.pred[index],Adducts.pred[index],SMILES.pred[index],DatabaseID[index],MS1.score[index],RMSE.iso[index],ms2score[index],all.score[index])}else{
      formula.final<-c(formula.pred[index],Adducts.pred[index],SMILES.pred[index],DatabaseID[index],MS1.score[index],RMSE.iso[index],ms2score[index],all.score[index])}
  formula.final<-matrix(formula.final,ncol=8)
  return(formula.final)}###just output the best one

################parse the formulas##################################
form.produce<-function(vector){
  formula.paste<-NULL
  if (vector[1]>0){
    if (vector[1]==1){formula.paste<-paste(formula.paste,"C",sep='')}else{
      formula.paste<-paste(formula.paste,"C",vector[1],sep='')}}
  if (vector[2]>0){
    if (vector[2]==1){formula.paste<-paste(formula.paste,"H",sep='')}else{
      formula.paste<-paste(formula.paste,"H",vector[2],sep='')}}
  if (vector[3]>0){
    if (vector[3]==1){formula.paste<-paste(formula.paste,"N",sep='')}else{
      formula.paste<-paste(formula.paste,"N",vector[3],sep='')}}
  if (vector[4]>0){
    if (vector[4]==1){formula.paste<-paste(formula.paste,"O",sep='')}else{
      formula.paste<-paste(formula.paste,"O",vector[4],sep='')}}
  if (vector[5]>0){
    if (vector[5]==1){formula.paste<-paste(formula.paste,"P",sep='')}else{
      formula.paste<-paste(formula.paste,"P",vector[5],sep='')}}
  if (vector[6]>0){
    if (vector[6]==1){formula.paste<-paste(formula.paste,"S",sep='')}else{
      formula.paste<-paste(formula.paste,"S",vector[6],sep='')}}
  if (vector[11]>0){
    if (vector[11]==1){formula.paste<-paste(formula.paste,"I",sep='')}else{
      formula.paste<-paste(formula.paste,"I",vector[11],sep='')}}
  if ((vector[7]+vector[8])>0){
    if ((vector[7]+vector[8])==1){formula.paste<-paste(formula.paste,"Cl",sep='')}else{
      formula.paste<-paste(formula.paste,"Cl",vector[7]+vector[8],sep='')}}
  if ((vector[9]+vector[10])>0){
    if ((vector[9]+vector[10])==1){formula.paste<-paste(formula.paste,"Br",sep='')}else{
      formula.paste<-paste(formula.paste,"Br",vector[9]+vector[10],sep='')}}
  if (length(formula.paste)<1){formula.paste<-0}
  return(formula.paste)}

#########function to extract fragment from DIA window#######
GetFrag<-function(Library,path,msfiles,mz_tol,LockMass){
  fragments<-NULL
  setwd(path)
  for (k in 1:length(msfiles)){
    print(c('getfragment...',k))
    xraw<-xcmsRaw(msfiles[k],includeMSn=TRUE)
    
  index<-which(Library$SampleID==k)
  if (length(index)==0){next}
  precursor<-preclist(xraw)
  len<-length(precursor)
  for (j in 1:length(index)){
  mz<-Library$mz[index[j]]
  DIAwin<-which(abs(mz-precursor)<=2.5)
  if (length(DIAwin)<1){next}
  DIAwin<-precursor[DIAwin[1]]
  ####get fragments from each window based correlation >0.85#####
  fragments<-getfrag(xraw,fragments,index[j],mz_tol,Library,DIAwin,30)
  }
  }
  return(fragments)}

########find precursor lists
preclist<-function (xmsn){
  x<-xmsn
  precmz<-xmsn@msnPrecursorMz
  len<-length(precmz)
  precur<-precmz[1]
  for (i in 2:len){
    if (length(which(precur==precmz[i]))==0){
      precur<-c(precur, precmz[i])
    }}
  return(precur)
}
###############################


#####################get fragments from DIA window only#############################
getfrag<-function(rawdata,prec_list,index2,mz_tol,Library,DIAmzwin,rtwindow){#get fragment spectra for each precursor ions
  Library$rt<-Library$rt*60
  index<-index2
  precurmz<-Library$mz[index2]
  
  DIAdata<-ms2copy(rawdata,DIAmzwin)##read rawdata for each DIA window
  mzrange<-DIAdata@mzrange
  minmz<-mzrange[1]
  maxmz<-mzrange[2]
  mzmin<-max(minmz,precurmz-precurmz*mz_tol)
  mzmax<-min(maxmz,precurmz+precurmz*mz_tol)
  rtrange<-DIAdata@scantime
  rtmin<-max(min(rtrange),Library$rt[index]-30)
  rtmax<-min(max(rtrange),Library$rt[index]+30)
  peaks<-rawEIC(DIAdata,mzrange=cbind(mzmin,mzmax),rtrange=cbind(rtmin,rtmax))
  scan.max<-which.max(peaks$intensity)
  scan.max<-peaks$scan[scan.max[1]]
  scanNum<-c(DIAdata@scanindex[scan.max],DIAdata@scanindex[scan.max+1])
  correctindex<-(scanNum[1]+1):scanNum[2]
  mz.frag<-DIAdata@env$mz[correctindex]
  if (length(mz.frag)<1){next}
  index<-which(mz.frag<precurmz-10)
  if (length(index)<1){next}
  mz.frag<-mz.frag[index]
  
  native.peak<-peaks$intensity
  kk<-0
  for (k in 1:length(mz.frag)){
    mz0<-mz.frag[k]
    mzmin<-max(minmz,mz0-mz0*mz_tol)
    mzmax<-min(maxmz,mz0+mz0*mz_tol)
    frag.peak<-rawEIC(DIAdata,mzrange=cbind(mzmin,mzmax),rtrange=cbind(rtmin,rtmax))
    frag.peak<-frag.peak$intensity
    if (sd(native.peak)==0||sd(frag.peak)==0){next}###0 values
    if (max(frag.peak)<2000){next}
    corr<-cor(native.peak,frag.peak)
    if (corr>0.8){#save the precursor ion to the list, if chromatographic is good
      list.frag<-list(precursor=precurmz,mz=mz.frag[k],rt=Library$rt[index2],intensity=max(frag.peak),score=corr,sampleid=Library$SampleID[index2],libraryid<-index2)
      if (length(prec_list)<1){
        prec_list<-list.frag}
      else{
        prec_list$precursor<-c(prec_list$precursor,list.frag$precursor)
        prec_list$mz<-c(prec_list$mz,list.frag$mz)
        prec_list$rt<-c(prec_list$rt,list.frag$rt)
        prec_list$intensity<-c(prec_list$intensity,list.frag$intensity)
        prec_list$score<-c(prec_list$score,list.frag$score)
        prec_list$sampleid<-c(prec_list$sampleid,Library$SampleID[index2])
        prec_list$libraryid<-c(prec_list$libraryid,index2)}}}
  return (prec_list)}

##############fragment calibration###############
FragCal<-function(Check.Frag,path,msfiles,LockMass){
  setwd(path)
  Save.Frag<-Check.Frag
  
  #########msn calibration using the average shift########
  for (i in 1:length(msfiles)){
    print(c('Calfragment...',i))
    index<-which(Check.Frag$sampleid==i)
    if(length(index)<1){next}
    
  xraw<-xcmsRaw(msfiles[i])
  ppmshift<-plotlock(xraw,LockMass,10)
  lock.shift<-data.frame(lockmass=LockMass,shift=rep(0,length(LockMass)))
  index.save<-NULL
  temp.fun<-NULL
  for (j in 1:length(LockMass)){
    temp<-ppmshift[j,]
    index.temp<-which(temp<15)
    if (length(index.temp)<1){
      index.save<-c(index.save,j)
      next}
    lock.shift$shift[j]<-mean(temp[index.temp])
  }
  lock.shift<-lock.shift[-index.save,]##delete those lockmass not detected
  lock.shift$shift<-lock.shift$shift*10^(-6)
  
  msnmzlist<-Check.Frag$mz[index]
  mzlist<-cbind(msnmzlist,1:length(msnmzlist))
  mzlist<-mzlist[order(mzlist[,1],decreasing=FALSE),]##sort the m/z for prediction
  mzlist<-matrix(mzlist,ncol=2)
  output<-fitlock(lock.shift,mzlist[,1],temp.fun)
  calmz<-output[[1]]
  temp.fun<-output[[2]]
  mzlist[,1]<-calmz
  mzlist<-mzlist[order(mzlist[,2],decreasing=FALSE),]##sort the m/z for prediction
  mzlist<-matrix(mzlist,ncol=2)
  Save.Frag$mz[index]<-mzlist[,1]}
  return(Save.Frag)
  }


#####################check fragments across samples, the correlation should be good among samples#############################
CheckFrag<-function(Fragments,path,msfiles,mz_tol,rtwin){#get fragment spectra for each precursor ions
  setwd(path)
  if (length(msfiles)==1){return(Fragments)}
  
  frag.save<-matrix(rep(0,length(Fragments$mz)*length(msfiles)),ncol=length(msfiles))
  prec.save<-matrix(rep(0,length(Fragments$mz)*length(msfiles)),ncol=length(msfiles))
  for (i in 1:length(msfiles)){
    print(c('checkfragment...',i))
    xraw<-xcmsRaw(msfiles[i],includeMSn=TRUE)
    precursor<-preclist(xraw)
    
    for (j in 1:length(Fragments$mz)){
      mz.frag<-Fragments$mz[j]###fragments
      mz.prec<-Fragments$precursor[j]###precursor
      
      #####precursor window and copy ms2 data######
        DIAwin<-which(abs(mz.prec-precursor)<=2.5)
        if (length(DIAwin)<1){next}
        DIAwin<-precursor[DIAwin[1]]
        DIAdata<-ms2copy(xraw,DIAwin)
        
      #######setting up searching############
        mzrange<-DIAdata@mzrange
        minmz<-mzrange[1]
        maxmz<-mzrange[2]
        rt<-Fragments$rt[j]
        mzmin<-max(minmz,mz.prec-mz.prec*mz_tol)
        mzmax<-min(maxmz,mz.prec+mz.prec*mz_tol)
        rtrange<-DIAdata@scantime
        rtmin<-max(min(rtrange),rt-rtwin)
        rtmax<-min(max(rtrange),rt+rtwin)
      
        #######intensity for precursor########
      peaks.prec<-rawEIC(DIAdata,mzrange=cbind(mzmin,mzmax),rtrange=cbind(rtmin,rtmax))
      scan.max<-which.max(peaks.prec$intensity)
      prec.save[j,i]<-max(peaks.prec$intensity)
      
        #######intensity for fragments########    
        mzmin<-max(minmz,mz.frag-mz.frag*mz_tol)
        mzmax<-min(maxmz,mz.frag+mz.frag*mz_tol)
        frag.peak<-rawEIC(DIAdata,mzrange=cbind(mzmin,mzmax),rtrange=cbind(rtmin,rtmax))
        frag.save[j,i]<-frag.peak$intensity[scan.max[1]]
    }}
        
  #########correlations##################
  saveid<-NULL
  for (k in 1:nrow(frag.save)){
    frag.peak<-frag.save[k,]
    prec.peak<-prec.save[k,]
    if (sd(prec.peak)==0||sd(frag.peak)==0){next}###0 values
    corr<-cor(prec.peak,frag.peak)
    if (corr>0.8){saveid<-c(saveid,k)}
  }
  if (length(saveid)==0){return(NULL)}
  
  prec_list<-Fragments
  Frag_list<-list()
  Frag_list$precursor<-prec_list$precursor[saveid]
  Frag_list$mz<-prec_list$mz[saveid]
  Frag_list$rt<-prec_list$rt[saveid]
  Frag_list$intensity<-prec_list$intensity[saveid]
  Frag_list$score<-prec_list$score[saveid]
  Frag_list$sampleid<-prec_list$sampleid[saveid]
  Frag_list$libraryid<-prec_list$libraryid[saveid]
  return (Frag_list)}

#################fucntion to copy MS2 data to MS1 matrix, but without precursor information
ms2copy <-function(xmsn,precursor) {
  x<-new("xcmsRaw")######if x<-xcms, the environment will be the same
  x@env<-new.env(parent = .GlobalEnv)
  index<-which(xmsn@msnPrecursorMz==precursor)
  x@tic <- xmsn@msnAcquisitionNum[index]
  x@scantime <- xmsn@msnRt[index]
  x@acquisitionNum <- xmsn@msnAcquisitionNum[index]
  x@polarity<-xmsn@polarity[1:length(index)]
  len2<-length(xmsn@msnPrecursorMz)
  index_total<-0
  index3<-0
  for (j in 1:length(index)){
    if (index[j]==len2){
      index2<-(xmsn@msnScanindex[index[j]]+1):length(xmsn@env$msnMz)########since the start of scan index is 0
    }
    else{
      index2<-(xmsn@msnScanindex[index[j]]+1):xmsn@msnScanindex[index[j]+1]
    }
    index3<-c(index3,index3[length(index3)]+length(index2))
    index_total<-c(index_total,index2)
  }
  index_total<-index_total[-1]
  index3<-index3[-length(index3)]
  x@env$mz <- xmsn@env$msnMz[index_total]
  x@env$intensity <- xmsn@env$msnIntensity[index_total]
  x@mzrange<-c(min(x@env$mz),max(x@env$mz))
  x@scanindex <-as.integer(index3)
  return(x)
}
#######################################

#-----------------------metfragscore---------------
#calculate score from data matrix
#--------------------------------------------------
metfragcal<-function(Library,Fragments,Database,ppm,charge){
  ms2score<-rep(0,nrow(Library))
  for (i in 1:nrow(Library)){
    mz<-Library$mz[i]
    Dbindex<-Library$DbID[i]
    if (Dbindex==0){next}##no searching results
    Dbindex<-strsplit(Dbindex,';')
    Dbindex<-Dbindex[[1]]
    index<-which(Fragments$libraryid==i)
    if (length(index)==0){next}##no fragments
    ms2spectra<-data.frame(mz=Fragments$mz[index],intensity=Fragments$intensity[index])##ms2 spectra
    mode<-charge
    smiles.match<-Database$SMILES[Dbindex]
    results<-metfrag_score(mz,ms2spectra,smiles.match,ppm,charge,mode)
    ms2score[i]<-paste(results$scores,sep=';')
  }
  return(ms2score)}


#------------------------------------
#use metfrag to calculate the ms2 score
#-----------------------------------
metfrag_score<-function(mz,ms2spectra,candidates,ppm,charge,mode){#ms2 spectra is the experimental results, candidates are the SMILES of compounds 
  cpd<-parse.smiles(candidates)
  mzs <-ms2spectra$mz
  ints <-ms2spectra$intensity
  exact.mass<-mz
  searching<-score.molecules.from.container(cpd, mzs, ints, exact.mass, 
                                        number.threads=1, mz.abs=0.003, 
                                        mz.ppm=ppm, pos.charge=charge, 
                                        mode=mode, tree.depth=2, score.names=c("FragmenterScore"), scoreWeights=c(1.0))
  if (length(smiles)==0){return(NULL)}
  results<-data.frame()
  scores<-getScores(score,scoreprop="Score")
  mols<-NULL
  for (i in 1:length(searching)){
  mols<-c(mols,get.smiles(searching[[i]]))
  }
  results<-data.frame(smiles=mols,score=scores)
  return(results)
  }

#------------------------
#read CAS
#------------------------
CAS.read<-function(CAS){
  CAS<-CAS[[1]]
  CAS<-as.character(CAS)
  temp<-strsplit(CAS,'')
  temp<-temp[[1]]
  temp1<-paste(c(temp[1:(length(temp)-3)]),collapse='')
  temp2<-paste(temp[(length(temp)-2):(length(temp)-1)],collapse='')
  CAS.new<-paste(c(temp1,temp2,temp[length(temp)]),collapse='-')
  return(CAS.new)
}

#---------------------
#delete List
#---------------------
DelList<-function(mylist,index){
  if (length(index)==0){return(mylist)}
  len<-length(mylist)
  for (i in 1:len){
    temp<-mylist[[i]]
    mylist[[i]]<-temp[-index]
  }
  return(mylist)}

#-----------------------
#split formula to matrix
#-----------------------
form.split<-function(formula){
  formula.paste<-rep(0,11)
  element<-c('C','H','N','O','P','S','Cl','37Cl','Br','81Br','I')
  for (i in 1:length(element)){
    if (length(grep(element[i],formula))<1){next}##no element
    temp<-strsplit(formula,element[i])
    temp<-temp[[1]]
    if (length(temp)==1){
      formula.paste[i]<-1##the end of the formula
      next}
    temp<-temp[2]
    temp<-strsplit(temp,'')
    temp<-temp[[1]]
    number.ele<-grep('[0-9]',temp)
    if (length(number.ele)<1){##all 1 atom
      formula.paste[i]<-1
      next
    }
    if (number.ele[1]>1){##1
      formula.paste[i]<-1
      next
    }
    temp.save<-number.ele[1]
    if (length(number.ele)==1){
      formula.paste[i]<-temp[temp.save]
      next}
    for (k in 2:length(number.ele)){
      if (number.ele[k]>(number.ele[k-1]+1)){break}
      temp.save<-c(temp.save,number.ele[k])}
    formula.paste[i]<-paste(temp[temp.save],collapse='')}
  return(as.numeric(formula.paste))}

#-----------------------build database-------------------
#-------------------------------------------------------
DB.build<-function(TSCA){
  xLogP<-rep(0,nrow(TSCA))
  for (i in 1:nrow(TSCA)){
    mysmiles<-smi.set@smilist[i]
    try(cpd<-parse.smiles(unlist(mysmiles)))
    try(xLogP[i]<-get.xlogp(cpd[[1]]))
  }
  
  smiles.list<-idsmiles(smi.set)
  molecule<-parse.smiles(smiles.list)
  formula<-lapply(molecule,getformula)
  
  formula.save<-NULL
  MASS<-NULL
  for (i in 1:length(formula)){
    formula.save<-c(formula.save,formula[[i]]@string)
    MASS<-c(MASS,formula[[i]]@mass)
  }
  
  results<-data.frame(smiles=idsmiles(smi.set),id=c(1:nrow(TSCA)),LogP=xLogP,formula=formula.save,mz=MASS)
  return(results)
}

getformula<-function(smile){
  formula<-(try(get.mol2formula(smile,charge=0)))
  if ('try-error' %in% class(formula)){
  molecule1 <- parse.smiles('N')[[1]]
  convert.implicit.to.explicit(molecule1)
  formula <- get.mol2formula(molecule1,charge=0)}
  return(formula)#return'NH3'
}

#------------------------------
#define adducts
#------------------------------
DefineAdducts<-function(adducts){
  element<-NULL
  if (adducts=='[M-H]-'){
    element<-c(0,-1,0,0,0,0,0,0,0,0,0,-1.007825)
  }
  if (adducts=='[M-H-H2O]-'){
    element<-c(0,-3,0,-1,0,0,0,0,0,0,0,-19.01839)
  }
  if (adducts=='[M+Cl]-'){
    element<-c(0,0,0,0,0,0,1,0,0,0,0,34.96885)
  }
  if (adducts=='[M+CH2O2-H]-'){
    element<-c(1,1,0,2,0,0,0,0,0,0,0,44.99765)
  }
  if (adducts=='[M]-'){
    element<-c(0,0,0,0,0,0,0,0,0,0,0,0)
  }
  if (adducts=='[M-Br+O]-'){
    element<-c(0,0,0,1,0,0,0,0,-1,0,0,-62.92342)
  }
  if (adducts=='[M+H]+'){
    element<-c(0,0,0,1,0,0,0,0,0,0,0,1.007825)
  }
  return(element)
}

#-------------------------------------
#neutral loss ID
#------------------------------------
LossCal<-function(Neutral,Fragment){
  loss.save<-NULL
  for (i in 1:nrow(Library)){
    index<-which(Fragment$libraryid==i)
    if (length(index)==0){next}##no fragment
    mz.prec<-Fragment$precursor[index[1]]
    mz.frag<-Fragment$mz[index]
    loss.mz<-mz.prec-mz.frag
    if (length(mz.frag)>1){
    for (j in 1:(length(mz.frag)-1)){
       for (m in (j+1):length(mz.frag)){
         loss.mz<-c(loss.mz,abs(mz.frag[m]-mz.frag[j]))
       }
    }}
    for (k in 1:length(loss.mz)){
    index2<-which(abs(loss.mz[k]-Neutral$mz)<0.002)##identify loss
    if (length(index2)>0){
      loss.save$mz<-c(loss.save$mz,Neutral$mz[index2])
      loss.save$libid<-c(loss.save$libid,i)
      loss.save$fragment<-c(loss.save$fragment,as.character(Neutral$frag[index2]))
    }}}
    return(loss.save)
  }
 
#------------------------------------------------
#Neutral Loss, and characteristic pattern seaching against SMARTS
#------------------------------------------------
pattern.neutral<-function(smiles,fragment){
  sdf<-smiles2sdf(smiles)
  sdf<-as(sdf,'SDFset')
  if (fragment=='CO2'){
    Match.loss= smartsSearchOB(sdf,"[CX3](=O)[OX1H0-,OX2H1]",uniqueMatches=FALSE)##carboxylic acid or conjugate salt
    if (max(Match.loss)>0){return(1)}##matched
  }
  if (fragment=='H2O'){
    Match.loss1= smartsSearchOB(sdf,"[CX3](=O)[OX1H0-,OX2H1]",uniqueMatches=FALSE)##carboxylic acid or conjugate salt
    Match.loss2= smartsSearchOB(sdf,"[#6][CX3](=O)[OX2H0][#6]",uniqueMatches=FALSE)##ester or aldehyde
    Match.loss<-c(Match.loss1,Match.loss2)
    if (max(Match.loss)>0){return(1)}##matched
  }
  if (fragment=='NH3'){
    Match.loss= smartsSearchOB(sdf,"[NX3;H2;!$(NC=[!#6]);!$(NC#[!#6])][#6]",uniqueMatches=FALSE)##primary amine
    if (max(Match.loss)>0){return(1)}##matched
  }
  if (fragment=='CH2O'){
    Match.loss= smartsSearchOB(sdf,"[CX3H1](=O)[#6]",uniqueMatches=FALSE)##Aldehyde
    if (max(Match.loss)>0){return(1)}##matched
  }
  if (fragment=='CH3OH'){
    Match.loss= smartsSearchOB(sdf,"[#6][CX3](=O)[OX2H0][#6]",uniqueMatches=FALSE)#ester or aldehyde
    if (max(Match.loss)>0){return(1)}##matched
  }
  if (fragment=='H2S'){
    Match.loss= smartsSearchOB(sdf,"[#16X2H]",uniqueMatches=FALSE)#THIOL
    if (max(Match.loss)>0){return(1)}##matched
  }
  if (fragment=='HCl'){
    Match.loss= smartsSearchOB(sdf,"[#6][Cl]",uniqueMatches=FALSE)#Cl
    if (max(Match.loss)>0){return(1)}##matched
  }
  if (fragment=='NO2'){
    Match.loss= smartsSearchOB(sdf,"[$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8]",uniqueMatches=FALSE)#NITRO GROUP
    if (max(Match.loss)>0){return(1)}##matched
  }
  if (fragment=='HCOOH'||fragment=='PO3'){
    Match.loss= smartsSearchOB(sdf,"[CX3](=O)[OX1H0-,OX2H1]",uniqueMatches=FALSE)#carboxylic acid or conjugate salt
    if (max(Match.loss)>0){return(1)}##matched
  }
  if (fragment=='H3PO4'){
    Match.loss= smartsSearchOB(sdf,"[$(P(=[OX1])([$([OX2H]),$([OX1-]),$([OX2]P)])([$([OX2H]),$([OX1-]),$([OX2]P)])[$([OX2H]),$([OX1-]),$([OX2]P)]),$([P+]([OX1-])([$([OX 2H]),$([OX1-]),$([OX2]P)])([$([OX2H]),$([OX1-]),$([OX2]P)])[$([OX2H]),$([OX1-]),$([OX2]P)])]",uniqueMatches=FALSE)#Phosphate
    if (max(Match.loss)>0){return(1)}##matched
  }
  if (fragment=='CO'){
    Match.loss= smartsSearchOB(sdf,"[CX3](=[OX1])C",uniqueMatches=FALSE)#Carboxilic acid, ester,anhydride and ketone
    if (max(Match.loss)>0){return(1)}##matched
  }
  if (fragment=='HBr'||fragment=='Br'){
    Match.loss= smartsSearchOB(sdf,"[#6][Br]",uniqueMatches=FALSE)#Cl
    if (max(Match.loss)>0){return(1)}##matched
  }
  if (fragment=='SO3'){
    Match.loss= smartsSearchOB(sdf,"[$([#16X4](=[OX1])(=[OX1])([#6])[OX2H,OX1H0-]),$([#16X4+2]([OX1-])([OX1-])([#6])[OX2H,OX1H0-])]",uniqueMatches=FALSE)#Sulfonic acid and conjugate
    if (max(Match.loss)>0){return(1)}##matched
  }
  return (0)
}

#------------------------------------------------
#Neutral Loss score
#------------------------------------------------
Score.neutral<-function(Library,fragment,Neutral){
  Neutralloss<-LossCal(Neutral,fragment)##determine which compounds show neutral loss
  if (length(Neutralloss)==0){return(NULL)}##No neutral loss
  Library$neutralscore<-rep(0,nrow(Library))
  for (i in 1:length(Neutralloss$mz)){
    smiles<-Library$SMILES[Neutralloss$libid[i]]##SMILES
    if (smiles==0){next}##no library compound
    smiles<-strsplit(smiles,';')
    smiles<-smiles[[1]]
    score.save<-NULL
    Loss.pattern<-Neutralloss$fragment[i]##the formula for neutral loss
    for (k in 1:length(smiles)){
    score.temp<-pattern.neutral(smiles[k],Loss.pattern)/2##0.5 for each neutral loss
    score.save<-paste(score.save,score.temp,sep=';')
    }
    Library$neutralscore[Neutralloss$libid[i]]<-score.save}
  return(Library)
}

#------------------------------------------------
#Characteristic ions
#------------------------------------------------
Score.charac<-function(Library,fragment,Character.ion){
  index.match<-NULL
  ions.match<-NULL
  for (i in 1:length(fragment$mz)){###determine if there is characteristic ions
    index<-which(abs(Character.ion$mz-fragment$mz[i])<0.002)
    if (length(index)>0){
      index.match<-c(index.match,i)
      ions.match<-c(ions.match,as.character(Character.ion$frag[index[1]]))}
  }
  if (length(index.match)==0){return(NULL)}##No neutral loss
  Library$chracaterscore<-rep(0,nrow(Library))
  for (i in 1:length(index.match)){
    index.lib<-fragment$libraryid[index.match[i]]
    smiles<-Library$SMILES[index.lib]##SMILES
    if (smiles==0){next}##no library compound
    smiles<-strsplit(smiles,';')
    smiles<-smiles[[1]]
    score.save<-NULL
    Loss.pattern<-ions.match[i]##the formula for neutral loss
    for (k in 1:length(smiles)){
      score.temp<-pattern.neutral(smiles[k],Loss.pattern)/2##0.5 for each neutral loss
      if(is.null(score.save)){
        score.save<-score.temp
        next
      }
      score.save<-paste(score.save,score.temp,sep=';')
    }
    Library$chracaterscore[index.lib]<-score.save}
  return(Library)
}

#------------------------------------------------
#ion mode pattern match
#------------------------------------------------
pattern.ionmode<-function(smiles,Ionmode,Ionsource,Adducts){
  sdf<-smiles2sdf(smiles)
  sdf<-as(sdf,'SDFset')
  if (Ionmode==-1){
    if (Adducts=='[M-H2O-H]-'){
      Match.loss1= smartsSearchOB(sdf,"[CX3](=O)[OX1H0-,OX2H1]",uniqueMatches=FALSE)##carboxylic acid or conjugate salt
      Match.loss2= smartsSearchOB(sdf,"[OX2H]",uniqueMatches=FALSE)##phenol and hydroxyl, acidic group including sulfonic and phosphoric
      Match.loss<-c(Match.loss1,Match.loss2)
      if (max(Match.loss)==0){return(-1)}##matched
      return (0)
    }
    if (Adducts=='[M+Cl]-'){
      Match.loss<-smartsSearchOB(sdf,"[#6][Cl]",uniqueMatches=FALSE)##chloride,M+Cl
      if (max(Match.loss)==0){return(-1)}##matched
      return (0)
    }
    if (Adducts=='[M-H]-'){
    Match.loss1= smartsSearchOB(sdf,"[CX3](=O)[OX1H0-,OX2H1]",uniqueMatches=FALSE)##carboxylic acid or conjugate salt
    Match.loss2= smartsSearchOB(sdf,"[CX3H1](=O)[#6]",uniqueMatches=FALSE)##aldehyde
    Match.loss3= smartsSearchOB(sdf,"[OX2H]",uniqueMatches=FALSE)##phenol and hydroxyl, acidic group including sulfonic and phosphoric
    Match.loss4= smartsSearchOB(sdf,"[NX3;H2;!$(NC=[!#6]);!$(NC#[!#6])][#6]",uniqueMatches=FALSE)##amine
    Match.loss5= smartsSearchOB(sdf,"[#16X2H]",uniqueMatches=FALSE)##thiol
    Match.loss6= smartsSearchOB(sdf,"[$([OH]-*=[!#6])]",uniqueMatches=FALSE)##hydroxy acidic
    Match.loss<-c(Match.loss1,Match.loss2,Match.loss3,Match.loss4,Match.loss5,Match.loss6)
    if (max(Match.loss)==0){return(-1)}##matched
    return (0)
    }
    if (Ionsource=='APCI'){
      Match.loss1<-smartsSearchOB(sdf,"[$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8]",uniqueMatches=FALSE)##nitro group
      Match.loss3<-smartsSearchOB(sdf,"[NX4](=O)",uniqueMatches=FALSE)##nitro group
      Match.loss2<-smartsSearchOB(sdf,"[#6][Br]",uniqueMatches=FALSE)##halogen M-X+O
      if (Adducts=='[M]-'){
        Match.loss<-c(Match.loss1,Match.loss3)
        if (max(Match.loss)==0){return(-1)}##matched
        return (0)}##NO2
      if (Adducts=='[M-Br+O]-'){
        Match.loss<-Match.loss2
        if (max(Match.loss)==0){return(-1)}##matched
        return (0)}##bromine
    }
  }
  if (Ionmode==1){
    Match.loss1= smartsSearchOB(sdf,"[NX3;H2,H1,H0;!$(NC=O)]",uniqueMatches=FALSE)##amine chemicals, but not amide
    Match.loss2= smartsSearchOB(sdf,"[#6][CX3](=O)[OX2H0][#6]",uniqueMatches=FALSE)##ester and aldehyde
    Match.loss3= smartsSearchOB(sdf,"[#6][CX3](=O)[#6]",uniqueMatches=FALSE)##ketone
    Match.loss4= smartsSearchOB(sdf,"[#6][OX2H]",uniqueMatches=FALSE)##hydroxyl in alcohol
    Match.loss<-c(Match.loss1,Match.loss2,Match.loss3,Match.loss4)
    if (max(Match.loss)==0){return(-1)}##matched
  }
  return (0)
}

#--------------------------------------------
#function to fill peaks
#--------------------------------------------
fillpeak<-function(xset.input,ppm,btw){
  ppm<-ppm/10^6
  idsave<-matrix(rep(0,length(xset.input@groupidx)*length(msfiles)*4),nrow=length(xset.input@groupidx)*length(msfiles),ncol=4)
  k<-1
  for (i in 1:length(xset.input@groupidx)){
    index<-unlist(xset.input@groupidx[i])
    sampleid<-xset.input@peaks[index,11]
    mz.value<-mean(xset.input@peaks[index,1])
    rt.value<-mean(xset.input@peaks[index,4])
    for (j in 1:length(unlist(phenoData(xset.input)))){
      index<-which(sampleid==j)
      if (length(index)<1){
        idsave[k,]<-c(i,j,mz.value,rt.value)##save groupidx,sampleid
        k<-k+1
      }}}
  index<-which(idsave[,1]==0)
  idsave<-idsave[-index,]
  newpeak<-matrix(rep(0,nrow(idsave)*5),nrow=nrow(idsave),ncol=5)
  kk<-1
  msfile<-filepaths(xset.input)
  minmz<-min(xset.input@peaks[,1])
  maxmz<-max(xset.input@peaks[,1])
  minrt<-min(xset.input@peaks[,4])
  maxrt<-max(xset.input@peaks[,4])
  if (length(idsave)>0){
    for (k in 1:length(unlist(phenoData(xset.input)))){
      print(c('filling peakID...',k,'of...',length(unlist(phenoData(xset.input)))))
      index<-which(idsave[,2]==k)
      if (length(index)==0){next}
      xraw<-xcmsRaw(msfiles[k])
      for (n in 1:length(index)){
        mz.value<-idsave[index[n],3]
        mzmin<-max(minmz,mz.value-mz.value*ppm)
        mzmax<-min(maxmz,mz.value+mz.value*ppm)
        rt.value<-idsave[index[n],4]
        rtmin<-max(minrt,rt.value-btw)
        rtmax<-min(maxrt,rt.value+btw)
        peak<-rawEIC(xraw,mzrange=cbind(mzmin,mzmax),rtrange=cbind(rtmin,rtmax))
        intensity<-max(peak$intensity)
        newpeak[kk,]<-c(mz.value,rt.value,intensity,idsave[index[n],1],k)##mz, rt, intensity,groupidx,sampleid
        kk<-kk+1
      }}}
  
  
  peakid<-length(xset.input@peaks[,1])
  for (i in 1:nrow(newpeak)){
    groupid<-newpeak[i,4]
    tempvalue<-unlist(xset.input@groupidx[groupid])
    peakid<-peakid+1
    xset.input@groupidx[groupid]<-list(c(tempvalue,peakid))}
  peak.combine<-matrix(0,ncol=11,nrow=nrow(newpeak))
  peak.combine[,1]<-newpeak[,1]
  peak.combine[,4]<-newpeak[,2]
  peak.combine[,9]<-newpeak[,3]
  peak.combine[,11]<-newpeak[,5]
  xset.input@peaks<-rbind(xset.input@peaks,peak.combine)
  return(xset.input)
}

#------------------------------------------------------------------
#score the adducts
#------------------------------------------------------------------
Score.Adducts<-function(mylib,Adducts.Find){
  mylib$adductscore<-rep(0,nrow(mylib))
  for (i in 1:nrow(mylib)){
    index<-which(Adducts.Find$ID==i)
    if(length(index)==0){next}##no adduct
    if (mylib$Adducts[i]==0){next}
  Adducts.temp<-strsplit(mylib$Adducts[i],';')
  adduct<-Adducts.temp[[1]]
  score.save<-rep(0,length(adduct))
  for (k in 1:length(adduct)){
    if(adduct[k]==Adducts.Find$Adducts[index[1]]){
    score.save[k]<-1}
    if (adduct[k]=='[M-Br+O]-'||adduct[k]=='[M]-'){#it is not possible to have M-H for M-Br+O, to avoid bias
      score.save[k]<-0.5
    }
  }
  mylib$adductscore[i]<-paste(score.save,collapse = ';')
  }
  return(mylib)
}

#-------------------------------------------------------------------
#Database matching
#------------------------------------------------------------------
DatabaseSearching<-function(cutoff,polarity,Database,mwoffset){#cutoff for intensity
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
    print(c('isotope...',i))
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
    }}
  
  Database$mz<-Database$mz+mwoffset
  IsotopeDB$mz<-IsotopeDB$mz+mwoffset
  
  setwd(path.out)
  msfiles<-list.files()
  mylib<-Library.new
  mylib$formula.pred<-rep(0,nrow(Library.new))
  mylib$isoscore<-rep(0,nrow(Library.new))
  mylib$mserror<-rep(0,nrow(Library.new))
  mylib$ms1score<-rep(0,nrow(Library.new))
  mylib$SMILES<-rep(0,nrow(Library.new))
  mylib$DbID<-rep(0,nrow(Library.new))
  mylib$Adducts<-rep(0,nrow(Library.new))
  for (i in 1:length(msfiles)){
    print(c('msfiles...',i))
    index<-which(Library.new$SampleID==i)
    if (length(index)==0){next}
    xrawdata<-xcmsRaw(msfiles[i])
    for (j in 1:length(index)){
      formula<-NULL
      formula<-formcal(number_input,Library.new,index[j],iso_list,xrawdata,ppm,ppm.ms2,IsotopeData,Database,Cal.Frag,adducts.input,IsotopeDB,mwoffset)
      if (length(formula)<2){next}
      max.form<-nrow(formula)
      mylib$formula.pred[index[j]]<-paste(formula[,1],collapse=';')
      mylib$isoscore[index[j]]<-paste(formula[,6],collapse=';')
      mylib$mserror[index[j]]<-paste(formula[,5],collapse=';')
      mylib$ms1score[index[j]]<-paste(formula[,8],collapse=';')
      mylib$SMILES[index[j]]<-paste(formula[,3],collapse=';')
      mylib$DbID[index[j]]<-paste(formula[,4],collapse=';')
      mylib$Adducts[index[j]]<-paste(formula[,2],collapse=';')
    }}
  
  ##-------------------------------------------------
  #neutral loss,Anal. Chem. 2014, 10724
  #-------------------------------------------------
  library(ChemmineOB)
  if (polarity==1){
    Neutral<-data.frame(frag=c('NH3','H2O','CH2O','CH3OH','H2S','HCl','NO2','HCOOH','H3PO4','CO2','CO'),
                        mz=c(17.0265,18.0106,30.0106,32.0262,33.9877,35.9767,45.9929,46.0055,97.97689,43.9898,27.9949))
    
  }
  if (polarity==-1){
    Neutral<-data.frame(frag=c('CH2O','H2O','NO2','HCOOH','CO2','CO','HBr','Br','HCl'),
                        mz=c(30.0106,18.0106,45.9929,46.0055,43.9898,27.9949,79.9262,78.9183,35.9767))
  }
  mylib.neutral<-Score.neutral(mylib,Cal.Frag,Neutral)
  
  ####---------------------------------------------
  #characteristic ion
  #------------------------------------------------
  if (polarity==1){
    Character.ion<-data.frame(frag=c('NH3'),
                              mz=c(17.0265))
    
  }
  if (polarity==-1){
    Character.ion<-data.frame(frag=c('Br','Br','SO3','PO3
                                     '),
                              mz=c(78.91885,80.91685,79.95735,78.95905))
  }
  mylib.charac<-Score.charac(mylib.neutral,Cal.Frag,Neutral)
  
  #------------------------------------------------
  #ion mode
  #------------------------------------------------
  mylib.ionmode<-mylib.charac
  mylib.ionmode$ionmodescore<-rep(0,nrow(mylib.ionmode))
  for (i in 1:nrow(mylib.charac)){
    print(c('ionmode',i))
    smiles<-mylib.ionmode$SMILES[i]##SMILES
    if (smiles==0){next}##no library compound
    smiles<-strsplit(smiles,';')
    smiles<-smiles[[1]]
    adduct<-mylib.ionmode$Adducts[i]##adducts
    if (adduct==0){next}##no library compound
    adduct<-strsplit(adduct,';')
    adduct<-adduct[[1]]
    score.save<-rep(0,length(smiles))
    for (j in 1:length(smiles)){
      score.save[j]<-pattern.ionmode(smiles[j],polarity,Ionsource,adduct[j])/2##score
    }
    mylib.ionmode$ionmodescore[i]<-paste(score.save,collapse=';')
  }
#-------------------------------
#score adducts
#------------------------------
myfinallib<-Score.Adducts(mylib.ionmode,Adducts.Find)##if the adduct is matched to the formula, plus 1
  
  return(myfinallib)
}

#--------------------------------
#Write the mass to ms files for Sirius searching
#-------------------------------
write.ms<-function (mydata, names, file.out, open = "w") 
{
  outfile <- file(description = file.out, open = open)
  for (i in 1:length(names)){
    temp<-paste(">",names[i],sep="")
    temp<-paste(c(temp,mydata[[i]]),collapse=" ")
    if (names[i]=='collision'||names[i]=='ms1peaks'){
      temp<-paste(">",names[i],sep="")
      writeLines(temp, outfile)
      temp.mz<-mydata[[i]]
      for (j in 1:nrow(temp.mz)){
        temp<-paste(temp.mz[j,],collapse=' ')
        writeLines(temp, outfile)
      }
    }else{
      writeLines(temp, outfile)}
  }
  close(outfile)
}

#----------------------------------
#write library to Sirius format
#----------------------------------
Sirius.build<-function(mylib,Cal.Frag,Isotope.Data){
  for (i in 1:nrow(mylib)){
    if (mylib$formula.pred[i]==0){next}
    fragment.id<-which(Cal.Frag$libraryid==i)
    if (length(fragment.id)==0){next}##no fragment
    mz<-mylib$mz[i]
    intensity<-10000
    Isotope.id<-which(Isotope.Data$ID==i)
    if (length(Isotope.id)>0){
      mz.isotope<-Isotope.Data$Isotope[Isotope.id]
      intensity.isotope<-Isotope.Data$intensity.iso[Isotope.id]
      intensity<-Isotope.Data$intensity[Isotope.id[1]]
    }
    ms1peak<-data.frame(ms=c(mz,mz.isotope),intensity=c(intensity,intensity.isotope))
    ms2peak<-ms1peak
    if (length(fragment.id)>0){
      mz.frag<-Cal.Frag$mz[fragment.id]
      intensity.frag<-Cal.Frag$intensity[fragment.id]
      ms2peak<-data.frame(ms=mz.frag,intensity=intensity.frag)
    }
    formula<-strsplit(mylib$formula.pred[i],';')
    formula<-unlist(formula)
    adduct<-strsplit(mylib$Adducts[i],';')
    adduct<-unlist(adduct)
    for (j in 1:length(formula)){
      mydata<-list(compound=paste(c(i,'_',j),collapse=''),formula=formula[j],parentmass=mz,
                   ionization=adduct[j],
                   collision=ms2peak,ms1peaks=ms1peak)
      write.ms(mydata, names(mydata), file.out = paste(c(i,'_',j,'.ms'),collapse=''))
    }}}

#---------------------------------------
#import Sirius data
#--------------------------------------
Import.MS2<-function(mylib,MS2files,Cal.Frag){
  mylib$MS2score<-rep(0,nrow(mylib))
  ms2scoresave<-NULL
  for (i in 1:nrow(mylib)){
    print(c(i,'of...',nrow(mylib)))
    index<-which(Cal.Frag$libraryid==i)
    if (length(index)==0){next}#no fragments
    formula<-mylib$formula.pred[i]
    if (formula=='0'){next}#no prediction
    formula<-strsplit(formula,';')
    formula<-unlist(formula)
    smiles<-mylib$SMILES[i]
    smiles<-strsplit(smiles,';')
    smiles<-unlist(smiles)
    treescore<-NULL
    for (j in 1:length(formula)){
    filename<-paste(c(i,'_',j,'.csv'),collapse = '')
    index2<-which(MS2files==filename)
    if (length(index2)==0){
      treescore<-c(treescore,0)#give a 0
      next
    }
    MS2file<-read.table(filename,header=TRUE,sep='\t',fill=TRUE)
    if (length(MS2file$smiles)==0){
      treescore<-c(treescore,0)#give a 0
    next}
    index.score<-NULL
      sdf2 <- smiles2sdf(as.character(smiles[j]))
      sdf2.ap<-sdf2ap(sdf2)
    for (k in 1:length(MS2file$smiles)){
      if (MS2file$smiles[k]==smiles[j]){
        index.score<-c(index.score,k)
        break
      }
      sdf1 <- try(smiles2sdf(as.character(MS2file$smiles[k])))
      if ('try-error' %in% class(sdf1)){next}
      sdf1.ap<-sdf2ap(sdf1)
      similarity<-cmp.similarity(sdf1.ap[1],sdf2.ap[1])
      if (similarity==1){
        index.score<-c(index.score,k)
        break
      }
    }
    if (length(index.score)==0){
      treescore<-c(treescore,0)#give a 0
      next}
    temp.diff<-min(0,3040+MS2file$score[index.score])##
    temp.tree<-exp(-0.5*temp.diff^2/56.8^2)
    treescore<-c(treescore,temp.tree)
    ms2scoresave<-c(ms2scoresave,temp.tree)
    }
    mylib$MS2score[i]<-paste(c(treescore),collapse = ';')}
  return(mylib)
}

#------------------------------------------------
#calculate final score
#------------------------------------------------
Finalscore<-function(mylib,weightK,precursor){
  mylib$allscore<-rep(0,nrow(mylib))
  for (i in 1:nrow(mylib)){
  formula<-mylib$formula.pred[i]
  if (formula[1]=='0'){next}
  formula<-unlist(strsplit(formula,';'))

  MS1score<-unlist(strsplit(mylib$ms1score[i],';'))
  
  if (mylib$adductscore[i]==0){#adduct score
    adductscore<-rep(0,length(formula))
  }else{
    adductscore<-unlist(strsplit(mylib$adductscore[i],';'))
  }
  
  if (mylib$chracaterscore[i]==0){#character score
    chascore<-rep(0,length(formula))
  }else{
    chascore<-unlist(strsplit(mylib$chracaterscore[i],';'))
  }
  
  if (mylib$neutralscore[i]==0){#character score
    neutralscore<-rep(0,length(formula))
  }else{
    neutralscore<-unlist(strsplit(mylib$neutralscore[i],';'))
    neutralscore<-neutralscore[-1]
  }
    
  if (mylib$ionmodescore[i]==0){#character score
    ionmodescore<-rep(0,length(formula))
  }else{
    ionmodescore<-unlist(strsplit(mylib$ionmodescore[i],';'))
  }
  
  if (mylib$MS2score[i]==0){#character score
    MS2score<-rep(0,length(formula))
  }else{
    MS2score<-unlist(strsplit(mylib$MS2score[i],';'))
  }
  
  if (mylib$isoscore[i]==0){#isotope score
    isoscore<-rep(0,length(formula))
  }else{
    isoscore<-unlist(strsplit(mylib$isoscore[i],';'))
  }
  
  MS1score<-as.numeric(MS1score)
  MS2score<-as.numeric(MS2score)
  ionmodescore<-as.numeric(ionmodescore)
  neutralscore<-as.numeric(neutralscore)
  chascore<-as.numeric(chascore)
  adductscore<-as.numeric(adductscore)
  
  temp.score<-MS1score*weightK[1]+MS2score*weightK[2]+ionmodescore*weightK[3]+
    neutralscore*weightK[4]+chascore*weightK[5]+adductscore*weightK[6]
  
  ####ion mode score is -1, put all score to 0
  index.ion<-which(ionmodescore<0)
  if (length(index.ion)>0){
    temp.score[index.ion]<-0
  }
  
  ####isotope score is <0.8, put all score to 0
  isoscore<-as.numeric(isoscore)
  index.iso<-which(isoscore<(log(0.5)))
  if (length(index.iso)>0){
    temp.score[index.iso]<-0
  }
  
  mylib$allscore[i]<-paste(c(temp.score),collapse=';')
  }
  
  index<-which(mylib$mz>min(precursor)-2.5)
  index2<-which(mylib$mz[index]<max(precursor)+2.5)
  mylib<-mylib[index[index2],]##only output the results with MS2 fragments
  return(mylib)
}

#-----------------------
#function to find cutoff
#-----------------------
Find.cut<-function(Target,Decoy){
  Targetscore<-NULL
  for (i in 1:nrow(Target)){
    temp<-strsplit(as.character(Target$allscore[i]),';')
    Targetscore<-c(Targetscore,unlist(temp))
  }
  Decoyscore<-NULL
  for (i in 1:nrow(Decoy)){
    temp<-strsplit(as.character(Decoy$allscore[i]),';')
    Decoyscore<-c(Decoyscore,unlist(temp))
  }
  Targetscore<-as.numeric(Targetscore)
  Decoyscore<-as.numeric(Decoyscore)
  ####select cutoff at 95% FDR
  cutoff<-0
  for (score in 1:1000){
    temp<--1+6*score/1000
    index.T<-length(which(Targetscore>temp))
    index.D<-length(which(Decoyscore>temp))
    if (index.D<0.05*(index.T+index.D)){#95% FDR
      cutoff<-temp
      break
    }
  }
  return(cutoff)
}

#-----------------------------------
#output chemicals with scores greater than cutoff
#-----------------------------------
Output<-function(Target,cutoff){
  index.del<-NULL
  Target$Final<-rep(0,nrow(Target))
  for (i in 1:nrow(Target)){
    score<-Target$allscore[i]
    temp<-strsplit(as.character(score),';')
    temp<-as.numeric(unlist(temp))
    index<-which(temp>cutoff)
    if (length(index)==0){
      index.del<-c(index.del,i)
      next
    }
    index<-which.max(temp)#only output the maximal formula
    formula<-strsplit(as.character(Target$formula.pred[i]),';')
    formula<-unlist(formula)
    Target$Final[i]<-paste(c(formula[index]),collapse = ';')
  }
  output<-Target[-index.del,]
return(output)
}


#---------------------------------------
#predict in silico retention time
#---------------------------------------
Predict.RT<-function(Target,Database,cutoff){
  marker<-NULL
  for (i in 1:nrow(Target)){
    score<-Target$allscore[i]
    temp<-strsplit(as.character(score),';')
    temp<-as.numeric(unlist(temp))
    index<-which(temp>cutoff)
    if (length(index)==0){next}
    ID<-strsplit(as.character(Target$DbID[i]),';')
    ID<-unlist(ID)
    ID<-as.numeric(ID[index])
    ID.max<-which.max(temp[index])
    ID<-ID[ID.max]
    LogP<-NULL
    Tpsa<-NULL
    for (j in 1:length(ID)){
      mol<-Database$smiles[ID[j]]
      mol<-parse.smiles(mol)[[1]]
      convert.implicit.to.explicit(mol)
      LogP<-c(LogP,get.xlogp(mol))
      Tpsa<-c(Tpsa,get.tpsa(mol))
    }
    RT<-rep(Target$rt[i],length(LogP))
    marker<-rbind(marker,cbind(RT,LogP,Tpsa))
  }
  marker<-data.frame(marker)
  
  #iterative RT
  R.coeff<-0
  kk<-0
  while(R.coeff<0.8&&kk<nrow(marker)*4/5){##could delete more than 2/3
    print(R.coeff)
    kk<-kk+1
    results<-lm(RT~LogP+Tpsa,data=marker)
    S.single<-summary(results)
    R.coeff<-S.single$r.squared
    res.single<-S.single$residuals
    if (R.coeff<0.8){
    index<-which.max(abs(res.single))##the maximal residual error
    marker<-marker[-index,]
    }}
  coeff<-S.single$coefficients
  marker$predictRT<-coeff[1,1]+marker$LogP*coeff[2,1]+marker$Tpsa*coeff[3,1]
  SD.rt<-sd(S.single$residuals)##SD of the regressions
 
  return(c(coeff[,1],SD.rt)) 
}

#-------------------------
#predict score for rt
#----------------------------
Score.RT<-function(Target,Database,RT.coeff){
  Target$allscore0<-rep(0,nrow(Target))
  SD.rt<-RT.coeff[4]
  for (i in 1:nrow(Target)){
    score<-Target$allscore[i]
    if (score==0){next}
    temp<-strsplit(as.character(score),';')
    temp<-as.numeric(unlist(temp))
    ID<-strsplit(as.character(Target$DbID[i]),';')
    ID<-unlist(ID)
    ID<-as.numeric(ID)
    LogP<-NULL
    Tpsa<-NULL
    for (j in 1:length(ID)){
      mol<-Database$smiles[ID[j]]
      mol<-parse.smiles(mol)[[1]]
      convert.implicit.to.explicit(mol)
      LogP<-c(LogP,get.xlogp(mol))
      Tpsa<-c(Tpsa,get.tpsa(mol))
      predictRT<-RT.coeff[1]+LogP*RT.coeff[2]+Tpsa*RT.coeff[3]
      Diff.rt<-predictRT-Target$rt[i]
      temp.score<-log(exp(-0.25*Diff.rt^2/SD.rt^2))
    }
    Target$rtscore[i]<-paste(temp.score,collapse=';')
    Target$allscore0[i]<-paste(temp+temp.score,collapse=';')
  }
  Target$allscore<-Target$allscore0
  return(Target)
}

#-------------------------------------------------
#chemical library
#--------------------------------------------------
BuildLib<-function(TSCA){
setwd(path)
library(webchem)
library(ChemmineR)
library(ChemmineOB)
TSCA<-read.table("TSCA2.csv",header=TRUE,sep=',')
TSCA$SMILES<-rep(0,nrow(TSCA))
CASID<-sapply(TSCA$casregno,CAS.read)
TSCA$CASRN<-CASID
for (i in 1:nrow(TSCA)){####query the SMILES
TSCA$SMILES[i]<-cir_query(TSCA$CASRN[i],'smiles')
}
smi.set<-as(TSCA$SMILES,'SMIset')


###-----------------------------
#store the Database to MySQL
#-------------------------------
setwd(path.db)
library(ChemmineR)
library(ChemmineOB)
library(RSQLite)
conn <- initDb("TSCADB.db")
idsmiles<-function(smileset){
  return(unlist(smileset@smilist))
}
setwd(path.db)
ids<-list.files()
mydb<-dbConnect(RSQLite::SQLite(),ids[1])

DBinfo<-DB.build(TSCA)##return ID, SMILES, LogP, formula, mz
dbWriteTable(mydb,'TSCA8',DBinfo)

#read data###
finaldb<- initDb("TSCADB_all.db")
dbListTables(mydb)
dbListFields(mydb,'TSCA1')
rawdata<-rbind(dbReadTable(mydb,"TSCA1"),
               dbReadTable(mydb,"TSCA2"),
               dbReadTable(mydb,"TSCA3"),
               dbReadTable(mydb,"TSCA4"),
               dbReadTable(mydb,"TSCA5"),
               dbReadTable(mydb,"TSCA6"),
               dbReadTable(mydb,"TSCA7"),
               dbReadTable(mydb,"TSCA8"))
dbWriteTable(finaldb,'all.TSCA',rawdata)}

#----------------------------------------
#produce unique ID
#----------------------------------------
UniqueID<-function(mylib){
  mylib$avg<-rowMeans(mylib[,3:(2+length(msfiles))])
  mylib<-mylib[order(mylib$avg,decreasing=TRUE),]##RANK FROM BIG TO LOW
  index.del<-NULL
  for (i in 2:nrow(mylib)){
    index.formula<-which(mylib$Final[1:(i-1)]==mylib$Final[i])
    if (length(index.formula)==0){next}
    index.rt<-which(abs(mylib$rt[index.formula]-mylib$rt[i])<0.2)
    if (length(index.rt)==0){next}
    index.del<-c(index.del,i)
  }
  if (length(index.del)>0){
    mylib<-mylib[-index.del,]   #This for loop eliminates any table rows with duplicate formulas
  } 
  
  mylib$SMILES<-as.character(mylib$SMILES)##defactor
  mylib$allscore0<-as.character(mylib$allscore0)
  for (i in 1:nrow(mylib)){
    score<-mylib$allscore0[i]
    score<-strsplit(as.character(score),';')
    score<-as.numeric(unlist(score))
    ID<-which.max(score)
    smiles<-strsplit(as.character(mylib$SMILES[i]),';')
    smiles<-unlist(smiles)
    mylib$allscore0[i]<-score[ID[1]]
    mylib$SMILES[i]<-smiles[ID[1]]
    print(i)
  }
  return(mylib)
}
