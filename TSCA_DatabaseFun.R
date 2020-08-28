#####Obtaining secondary results for polymer/salt matches#####

path<-"C:/Users/Steven Desktop/Documents/MSc2/TargetDecoy"
path.db<-paste(path,"/SMILES_DATABASE", sep="")

setwd(path.db)
library(RSQLite)
library(tidyr)
library(stringr)
mydb<-dbConnect(RSQLite::SQLite(),'TSCADB_all.db')
Database<-dbReadTable(mydb,"TSCA")

setwd("C:/Users/Steven Desktop/Google Drive/MSc 3/Algorithm Paper")
oldresults<-read.table("PolymerSaltResults.csv",header=TRUE, sep=',', fill=TRUE)
names(oldresults)<-c("SMILES","DbID","allscore","CmpdID")

searchtable<-matrix(data=NA,nrow=nrow(oldresults),ncol=5)
colnames(searchtable)<-c("CmpdID","DbID","newAllscore0","newSMILES","newformula")
dotindex<-nrow(oldresults)
ranksearch<-1

for(i in 1:nrow(searchtable)){
  allscorelist<-as.numeric(unlist(strsplit(as.character(oldresults[i,3]),';')))
  allscore<-sort(allscorelist,decreasing=TRUE)[ranksearch]
  allscoreindex<-which(allscorelist==allscore)
    
  DbIDList<-as.numeric(unlist(strsplit(as.character(oldresults[i,2]),';')))
  DbID<-DbIDList[allscoreindex]
  
  newSMILES<-Database$smiles[DbID]
  if(length(newSMILES)>1){newSMILES<-paste(newSMILES,collapse=";")}
  j<-1
  while(str_detect(newSMILES,"\\.")==TRUE || is.na(newSMILES)==TRUE){
    ranksearch<-ranksearch+1
    if(ranksearch>length(allscorelist)){
      newSMILES<-NA
      break
      }
    allscore<-sort(allscorelist,decreasing=TRUE)[ranksearch]
    allscoreindex<-which(allscorelist==allscore)
    if(length(allscoreindex)>1){
      allscoreindex<-allscoreindex[j]
      j<-j+1
      }else{j<-1}
    DbIDList<-as.numeric(unlist(strsplit(as.character(oldresults[i,2]),';')))
    DbID<-DbIDList[allscoreindex]
    newSMILES<-Database$smiles[DbID]
  }
print(newSMILES)
searchtable[i,1]<-oldresults$CmpdID[i]
searchtable[i,2]<-DbID
searchtable[i,3]<-allscore
searchtable[i,4]<-newSMILES
searchtable[i,5]<-Database$formula[DbID]
    
ranksearch<-1
}

lowscoreindex<-which(searchtable[,3]<=0.27)
searchtable<-searchtable[-lowscoreindex,]
naindex<-which(is.na(searchtable[,4]))
searchtable<-searchtable[-naindex,]

setwd("C:/Users/Steven Desktop/Google Drive/MSc 3/Algorithm Paper")
write.table(searchtable,file="PolymerSalt_NewResults.csv",sep=',',row.names = FALSE)

#####Quick Score Difference calculator#####
library(stringr)
setwd("C:/Users/Steven Desktop/Google Drive/MSc 3/Algorithm Paper")
compoundtable<-read.table("test.csv",header=TRUE, sep=',', fill=FALSE)

scoretable<-matrix(data=NA,nrow=nrow(compoundtable),ncol=4)
colnames(scoretable)<-c("CmpdID","Score1","Score2","ScoreDiff")

for(i in 1:nrow(scoretable)){
  allscorelist<-strsplit(as.character(compoundtable[i,1]),';')
  allscorelist<-as.numeric(unlist(allscorelist))
  allscorelist<-sort(allscorelist,decreasing=TRUE)
  score1<-allscorelist[1]
  if(length(allscorelist<2)==TRUE){
    print(i)
    score2<-score1
    scorediff<-"Unique"
  }else{
  score2<-allscorelist[2]
  scorediff<-score1-score2
  }
  scoretable[i,1]<-compoundtable$CpdID[i]
  scoretable[i,2]<-score1
  scoretable[i,3]<-score2
  scoretable[i,4]<-scorediff
}
write.table(scoretable,file="ScoreDifferences.csv",sep=',',row.names = FALSE)


oldresults<-read.table("PolymerSaltResults.csv",header=TRUE, sep=',', fill=TRUE)
names(oldresults)<-c("SMILES","DbID","allscore","CmpdID")
scoretable2<-matrix(data=NA,nrow=nrow(oldresults),ncol=4)
colnames(scoretable2)<-c("CmpdID","Score1","Score2","ScoreDiff")

for(i in 1:nrow(scoretable2)){
  allscorelist<-strsplit(as.character(oldresults[i,3]),';')
  allscorelist<-as.numeric(unlist(allscorelist))
  allscorelist<-sort(allscorelist,decreasing=TRUE)
  score1<-allscorelist[1]
  if(length(allscorelist<2)==TRUE){
    print(i)
    score2<-score1
    scorediff<-"Unique"
  }else{
    score2<-allscorelist[2]
    scorediff<-score1-score2
  }
  scoretable2[i,1]<-oldresults$CmpdID[i]
  scoretable2[i,2]<-score1
  scoretable2[i,3]<-score2
  scoretable2[i,4]<-scorediff
}
write.table(scoretable2,file="ScoreDifferences2.csv",sep=',',row.names = FALSE)