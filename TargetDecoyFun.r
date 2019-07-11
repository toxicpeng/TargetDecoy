

##########FUNCTION to plot lockmass shift########
plotlock<-function(xraw,Lockmass,ppm){
    ppmmatrix<-matrix(rep(0,length(xraw@scanindex)*length(Lockmass)),nrow=length(Lockmass),ncol=length(xraw@scanindex))
    #Lockmass<-Lockmass-0.0005485799*polarity
    for (i in 2:length(xraw@scanindex)){##mass calibration using the closest lockmass
         scanNum<-c(xraw@scanindex[i-1],xraw@scanindex[i])
         correctindex<-scanNum[1]:scanNum[2]
         for (j in 1:length(Lockmass)){
              mz.index<-which(abs(xraw@env$mz[correctindex]-Lockmass[j])<ppm*10^(-6)*Lockmass[j])
              if (length(mz.index)==0){
                 next}
         index2<-which.max(xraw@env$intensity[correctindex[mz.index]])##if multiple data exist, pick out the most abundant one, questionable?
         mz.index<-mz.index[index2]
         mz.index<-mz.index[1]
         mz.index<-scanNum[1]+mz.index-1
         ppmmatrix[j,i]<-10^6*(xraw@env$mz[mz.index]-Lockmass[j])/Lockmass[j]
         if (xraw@env$intensity[mz.index]<10000){ppmmatrix[j,i]<-0}
         }         
         }
    return(ppmmatrix)
}


##############function to identify lockmass, written by Steven 2017/12/20
findlock<-function(msfile,intthresh,stepmz){
  newpeak<-NULL
  p<-0
  for (k in 1:length(msfile)){
    xraw<-xcmsRaw(msfile[k],mslevel=1)
    for (mz.value in seq(from=min(xraw@mzrange),to=max(xraw@mzrange),by=stepmz)){
      mz.min<-mz.value-0.001
      mz.min<-max(mz.min, xraw@mzrange[1])
      mz.max<-mz.value+0.001
      mz.max<-min(mz.max,xraw@mzrange[2])
      rt.min<-(min(xraw@scantime)+180) #range of analysis is the entire chromatograph
      rt.max<-(max(xraw@scantime)-180) #the first and last 3 minutes are skipped to avoid sections with no analytes
      peaks<-rawEIC(xraw,mzrange=cbind(mz.min,mz.max),rtrange=cbind(rt.min, rt.max))
      peaklist<-unlist(peaks$intensity)
      newlist<-NULL
      for(i in 1:length(peaklist)){
        if (peaklist[i]>0){newlist<-c(newlist,peaklist[i])}
        else if (is.na(peaklist[i+1])==TRUE){next}
        else if (peaklist[i+1]>0){next} ##allows mz values through which do not have consecutive zero intensity values
        else {newlist<-c(newlist,peaklist[i])}}
      if (length(newlist)==0){next}
      minintensity<-min(newlist)
      if (minintensity>intthresh){
        newpeak<-rbind(newpeak,c(mz.value,minintensity,k))
        p<-p+1
        print(c("running...",k,"-",p))
      }
    }
  }
  return (newpeak)
}

form.read<-function(formulainput,charinput){
element.lib<-strsplit(as.character(formulainput),charinput)
element.lib<-unlist(element.lib)
if (length(element.lib)==1&&element.lib!=as.character(formulainput)){return(1)}##at the end, CH3Br
if (length(element.lib)==1){return(0)}##no such element
element.lib<-strsplit(element.lib[2],"")
element.lib<-as.numeric(unlist(element.lib))
value<-NULL
if (is.na(element.lib[1])){##put it to 1, when like CH4
value<-1
}else{
for (i in 1:length(element.lib)){###when1,2,need to read, C12H5,
if (is.na(element.lib[i])){break}
value<-c(value,element.lib[i])}
value<-sum(value*10^((length(value)-1):0))
}
return(value)}

form.parse<-function(rawformula){##formulae with two column, the second column is formulae
formulae.save<-NULL
for (i in 1:nrow(rawformula)){
if (rawformula[i,2]=='0'){
formulae.save<-rbind(formulae.save,rep(0,9))
next}
formula.input<-as.character(rawformula[i,2])
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
         mz.cal<-MzCal(form.ReadOne(Library[i,2]),0.1)
         mz.cal<-mz.cal[,1]
         for (j in 1:min(length(mz.cal),6)){
              Database[i,3+j]<-mz.cal[j]}
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
     Database[i,3]<-mz
     Database[i,4:ncol(Database)]<-Database[i,4:ncol(Database)]+mz
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

################funtion to calculate isotope distribution of proposed formulae########
IsoscoreCal<-function(xrawdata,rtinput,isoinput,mzinput){
rtinput<-rtinput*60
mz_tol<-5*10^(-6)
iso.pattern<-isoinput
iso.pattern[,2]<-isoinput[,2]/(sum(isoinput[,2]))
iso.intensity<-NULL
save.peak1<-getEIC(xrawdata,mzrange=cbind(mzinput-mz_tol*mzinput,mzinput+mz_tol*mzinput),rtrange=cbind(rtinput-20,rtinput+20),step=0.001)
save.peak1<-unlist(save.peak1@eic[[1]])
save.peak<-save.peak1[((length(save.peak1)/2)+1):length(save.peak1)]
index.first<-which.max(save.peak)###find the exact peak rt for the mass spec files
for (i in 1:nrow(iso.pattern)){
          mzvalue<-iso.pattern[i,1]
          iso.peak1<-getEIC(xrawdata,mzrange=cbind(mzvalue-mz_tol*mzvalue,mzvalue+mz_tol*mzvalue),rtrange=cbind(rtinput-20,rtinput+20),step=0.003)
          iso.peak1<-unlist(iso.peak1@eic[[1]])
          iso.peak<-iso.peak1[((length(iso.peak1)/2)+1):length(iso.peak1)]
          cor_score<-cor(save.peak,iso.peak)
          if ((is.na(cor_score))||cor_score<0.8){##otherwise, just noise, should not calculate
          iso.peak[]<-0
          }
          iso.intensity<-c(iso.intensity,iso.peak[index.first])##always using the exact time of the first peak
          }
if (max(iso.intensity)==0){
return(0)} ###return large value if no ions are detected
iso.intensity<-iso.intensity/max(sum(iso.intensity),0.1)###set the theory one as 1, if the maximal intensity is 0, means that it is not that possible, put a big values
iso.pattern<-iso.pattern[,2]
ratio<-iso.intensity/iso.pattern
index.noise<-which(ratio>3)##too high signal means noise
iso.intensity[index.noise]<-iso.pattern[index.noise]##for these noise, not included for data analysis
RMSE<-sum(abs(iso.pattern-iso.intensity))##calculate the similarity of isotopic peaks pattern
return(1-RMSE)
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
    
#########Merge the isotope peaks with mass error###############
isopattern_Merge<-function(isotope,ppm){
     if (nrow(isotope)<2){return(isotope)}##only one isotope
     isotope<-isotope[order(isotope[,2],decreasing = TRUE),]
     index.save<-NULL
     for (i in 2:nrow(isotope)){
         index<-which(abs(isotope[i,1]-isotope[1:(i-1),1])<ppm*isotope[i,1])###overlapping with previous peaks, delete
         if (length(index)>0){
         index.save<-c(index.save,i)
         isotope[index[1],2]<-isotope[index[1],2]+isotope[i,2]##addition of the relative abundance
         isotope[i,1]<-0}}
     if (length(index.save)==0){return(isotope)}    
     return(isotope[-index.save,])}

##############mass calibration using lockmass#########       
MassCal<-function(xraw,LockMass){
    Lockmass<-sort(LockMass)
    Lockmass<-Lockmass-0.0005485799*polarity
    for (i in 2:length(xraw@scanindex)){##mass calibration using the closest lockmass
         scanNum<-c(xraw@scanindex[i-1],xraw@scanindex[i])
         correctindex<-scanNum[1]:scanNum[2]
         massshift<-rep(0,length(Lockmass))
         for (k in 1:length(Lockmass)){##define the shift for each lock mass
         mz.index<-which(abs(xraw@env$mz[correctindex]-Lockmass[k])<10*10^(-6)*Lockmass[k])##Determine the index of any mz in the given scan range which fall within 10 ppm of the current Lockmass
         if (length(mz.index)==0){
         massshift[k]<-100 ##There will be a '100' for every Lockmass which has no mz within 10 ppm of it (for the given scanindex range)
         next}
         index2<-which.max(xraw@env$intensity[correctindex[mz.index]])##if multiple data exist, pick out the most abundant one, questionable?
         mz.index<-mz.index[index2]
         mz.index<-scanNum[1]+mz.index-1##change the index from just the 'correctindex' range to the entire 'scanindex' range (still the same mz)
         if (xraw@env$intensity[mz.index]<1e4){
         massshift[k]<-100
         next}###if the mass lock intensity too low, delete it
         massshift[k]<-(xraw@env$mz[mz.index]-Lockmass[k])/Lockmass[k]}##calculate the mass shift
         Lockmass1<-Lockmass                                     
         index1<-which(massshift==100)
         if (length(index1)>0){
         Lockmass1<-Lockmass1[-index1]
         massshift<-massshift[-index1]} ##Remove the massshift values of '100' from the list (and the corresponding Lockmasses)
         if (length(Lockmass1)==0){next}
         for (j in 1:length(Lockmass1)){
         if (length(Lockmass1)==1){
             xraw@env$mz[correctindex]<-xraw@env$mz[correctindex]*(1-massshift[1])
             next
         }
         if (j==1){
            segment<-c(0,mean(c(Lockmass1[j],Lockmass1[j+1])))
            index<-which(xraw@env$mz[correctindex]<=segment[2])##find the closest lockmass
            if (length(index)<1){next}
            xraw@env$mz[correctindex[index]]<-xraw@env$mz[correctindex[index]]*(1-massshift[j])##calibration
            next
         }
         else if (j==length(Lockmass1)){
            segment<-mean(c(Lockmass1[j-1],Lockmass1[j]))
            index<-which(xraw@env$mz[correctindex]>segment)##find the closest lockmass
            if (length(index)<1){next}
            xraw@env$mz[correctindex[index]]<-xraw@env$mz[correctindex[index]]*(1-massshift[j])##calibration
            next
         }
         else if(j<length(Lockmass1)&&j>1){
            segment<-c(mean(c(Lockmass1[j],Lockmass1[j-1])),mean(c(Lockmass1[j],Lockmass1[j+1])))
            index<-which(xraw@env$mz[correctindex]>segment[1])##find the closest lockmass
            if (length(index)<1){next}
            index1<-which(xraw@env$mz[correctindex[index]]<=segment[2])
            if (length(index1)<1){next}
            xraw@env$mz[correctindex[index[index1]]]<-xraw@env$mz[correctindex[index[index1]]]*(1-massshift[j])##calibration
         }}}
         return(xraw)}         
   
