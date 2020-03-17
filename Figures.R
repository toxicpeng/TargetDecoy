####ggplot histogram####
library(ggplot2)
library(plyr)
df<-data.frame(group=c(rep('target',length(Targetscore)),rep('decoy',length(Decoyscore))),cor=c(Targetscore,Decoyscore))
mu <- ddply(df, "group", summarise, grp.mean=mean(cor))
ggplot(df, aes(x=cor, color=group, fill=group)) +
  geom_histogram(aes(y=..density..), position="identity", alpha=0.5)+
  geom_density(alpha=0.6)+
  geom_vline(data=mu, aes(xintercept=grp.mean, color=group),
             linetype="dashed")+
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  labs(title="Weight histogram plot",x="score", y = "Density")+
  xlim(-5,5)+
  theme_classic()

###single###
ggplot(df, aes(x=cor)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666") 

###tune the weight#
index.save<-0
saveweight<-0
for (i1 in 1:3){
  print(i1)
  for (i2 in 1:3){
    for (i3 in 1:3){
      for (i4 in 1:3){
        for (i5 in 1:3){
          for (i6 in 1:3){
            weightK<-c(i1,i2,i3,i4,i5,i6)/3#weight for MS1, ms2, ionmode, neutral, characteristic, and adducts
            mylib.output.decoy<-Finalscore(mylib.output1,weightK)
            mylib.output.target<-Finalscore(mylib.output2,weightK)
            
            setwd(path)
            Target<-mylib.output.target
            Decoy<-mylib.output.decoy
            Targetscore<-NULL
            for (i in 5427:9029){
              temp<-strsplit(as.character(Target$allscore[i]),';')
              Targetscore<-c(Targetscore,unlist(temp))
            }
            Decoyscore<-NULL
            for (i in 5427:9029){
              temp<-strsplit(as.character(Decoy$allscore[i]),';')
              Decoyscore<-c(Decoyscore,unlist(temp))
            }
            Targetscore<-as.numeric(Targetscore)
            Decoyscore<-as.numeric(Decoyscore)
            ####select cutoff at 95% FDR
            cutoff<-0
            for (score in 1:1000){
              temp<-1+6*score/1000
              index.T<-length(which(Targetscore>temp))
              index.D<-length(which(Decoyscore>temp))
              if (index.D<0.05*(index.T+index.D)){#95% FDR
                cutoff<-temp
                break
              }
            }
            if (index.T>index.save){
              index.save<-index.T
              saveweight<-weightK
            }
          }
        }
      }
    }
  }
}

##
Targetscore<-NULL
for (i in 1:nrow(Target.rt)){
  if (Target.rt$allscore[i]==0){
    next
  }
  temp<-strsplit(as.character(Target.rt$allscore[i]),';')
  temp<-as.numeric(unlist(temp))
  Targetscore<-c(Targetscore,temp)
}
Decoyscore<-NULL
for (i in 1:nrow(Decoy.rt)){
  if (Decoy.rt$allscore[i]==0){
    next
  }
  temp<-strsplit(as.character(Decoy.rt$allscore[i]),';')
  temp<-as.numeric(unlist(temp))
  Decoyscore<-c(Decoyscore,temp)
}

#ROC curve#
library(ggplot2)
library(plotROC)
targetdata<-Targetscore
decoydata<-Decoyscore
df<-data.frame(ID=c(rep(1,length(targetdata)),rep(0,length(decoydata))),
               group=c(rep('target',length(targetdata)),rep('decoy',length(decoydata))),
               score=c(targetdata,decoydata))
basicplot <- ggplot(df, aes(d = ID, m = score)) + geom_roc()
styledplot <- basicplot + style_roc()
styledplot

score.silico<-NULL
score.expt<-NULL
score.sep<-NULL
for (i in 1:1000){
  score.temp<--10+i*15/1000
  score.sep<-c(score.sep,score.temp)
  index.decoy<-length(which(Decoyscore>score.temp))
  index.target<-length(which(Targetscore>score.temp))
  score.silico<-c(score.silico,index.decoy/index.target)
  index.1<-length(which(targetdata>score.temp))
  index.2<-length(which(decoydata>score.temp))
  score.expt<-c(score.expt,index.2/index.1)
}
Score<-data.frame(FDR=c(score.silico/1.2,score.expt),score=c(score.sep,score.sep),
                  group=factor(c(rep('Silico',length(score.sep)),rep('Expt',length(score.sep)))))
ggplot(Score, aes(x=score, y=FDR,group=group,colour=group)) + 
  geom_line(size=1)

ID<-60
smiles<-Uniqueid$SMILES[ID]
sdf<-smiles2sdf(smiles)
plot(sdf)


####Chlorinated Azo Dyes Heatmap with Dendrogram####

library(ggplot2)
library(tidyr)
library(ggdendro)
library(reshape2)
library(grid)
library(corrplot)
library(Hmisc)

chlor<-read.csv("E:/Steven/Target Decoy/analysis/dust analysis/DustCompoundsUnique.csv", header=TRUE)

# Scale data
chlor.scaled <- chlor[,c(1:27,43,44,47,48)]
chlor.scaled[, c(3:26)] <- scale(chlor.scaled[, 3:26])

# Run clustering
chlor.matrix <- as.matrix(chlor.scaled[, c(3:26)])
rownames(chlor.matrix) <- chlor.scaled$CpdID
chlor.dendro <- as.dendrogram(hclust(d = dist(x = chlor.matrix)))

# Create dendrogram
dendro.plot <- ggdendrogram(data = chlor.dendro, rotate = TRUE)+
            theme(axis.text.x = element_blank())

# Data reordering
chlor.long <- gather(data = chlor.scaled, key = Sample, value = Intensity, -c(1:2), -c(27:33))
# Extract the order of the tips in the dendrogram
chlor.order <- order.dendrogram(chlor.dendro)
# Order the levels according to their position in the cluster
chlor.long$CpdID <- factor(x = chlor.long$CpdID,
                               levels = chlor.scaled$CpdID[chlor.order], 
                               ordered = TRUE)

##Heatmap
heatmap.plot <- ggplot(data=chlor.long, mapping=aes(x=Sample, y=CpdID, fill=Intensity)) +
      geom_raster()+
      xlab(label="Sample")+
      scale_fill_gradient(name="Intensity", low="yellow", high="blue")+
      theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "top")

##All together
#grid.newpage()
#print(heatmap.plot, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
#print(dendro.plot, vp = viewport(x = 0.90, y = 0.455, width = 0.2, height = 0.91))

##Alternate (much faster) option for plotting correlation score heatmap
chlor_rows<-which(chlor.scaled[,30]==TRUE)
chlor.matrix <- as.matrix(chlor.scaled[chlor_rows, c(3:26)])
rownames(chlor.matrix) <- chlor.scaled[chlor_rows,31]
#melt_chlor<-melt(chlor.matrix)
#cast_chlor<-acast(melt_chlor, Var2~Var1)
#chlor.cormatrix<-rcorr(cast_chlor)
chlor.cormatrix<-rcorr(chlor.matrix)
colour1 <- colorRampPalette(c("red","white","blue"))
corrplot(chlor.cormatrix$r, method = "square", type = "full", order = "hclust" ,col = colour1(200), tl.col = "black", tl.srt = 45, p.mat = chlor.cormatrix$P, sig.level = 0.05, insig = "blank", addrect = 5)


####Faceted Bar Graph for All Chlorine-Containing Compounds####

library(ggplot2)
library(tidyr)
library(ggdendro)
library(reshape2)
library(ggthemes)

clcpds<-read.csv("C:/Users/skutarna/Google Drive/NTA/TargetDecoy/ChlorAll.csv", header=TRUE)
clcpds$avg<-log(clcpds$avg)

grid.plot <- ggplot(data=clcpds, mapping=aes(x=CpdID, y=avg))+
                  geom_bar(stat="identity")+
                  coord_cartesian(ylim=c(9,15))+
                  facet_wrap(~Class)+
                  xlab(label="Compound #")+
                  ylab(label="Log Average Intensity")+
                  scale_x_continuous(breaks = pretty(clcpds$CpdID, n = 25))+
                  theme(axis.text.x = element_text(size = 6))               
  
print(grid.plot)



####Chlorinated Azo Dyes by Colour and Dye Type####

library(ggplot2)
library(tidyr)
library(ggthemes)

cldyes<-read.csv("C:/Users/skutarna/Google Drive/NTA/TargetDecoy/KnownDyes.csv", header=TRUE)
cldyes$avg<-log(cldyes$avg)

color.plot <- ggplot(data=cldyes, mapping=aes(x=Colour, y=avg))+
  geom_bar(stat="identity")+
  coord_cartesian(ylim=c(9,15))+
  facet_grid(.~Type)+
  xlab(label="Dye Colour")+
  ylab(label="Log Average Intensity")
             
print(color.plot)


####Lockmass Calibration - Searching Space Reduction Figure####

library(ggplot2)
library(ggthemes)

data5ppm<-read.csv("C:/Users/skutarna/Documents/MSc 2/NTA/Algorithm Paper/formulae5ppm.csv", header=TRUE)
names(data5ppm)<-c("mz","results")

group10<-NULL
group10<-data.frame()
for (i in 1:length(data5ppm[,1])){
  tens<-signif(data5ppm[i,1],2)
  group10[i,1]<-tens
  group10[i,2]<-data5ppm[i,2]
}
names(group10)<-c("mz","results")
group10data<-aggregate(group10[,2],list(group10[,1]),max)
names(group10data)<-c("mz","results")

group100<-NULL
group100<-data.frame()
for (i in 1:length(data5ppm[,1])){
  hundreds<-signif(data5ppm[i,1],1)
  group100[i,1]<-hundreds
  group100[i,2]<-data5ppm[i,2]
}
names(group100)<-c("mz","results")
group100data<-aggregate(group100[,2],list(group100[,1]),max)
names(group100data)<-c("mz","results")

formula.plot <- ggplot(group10data, aes(mz,results))+
  geom_area(color="black")+
  stat_smooth()+
  xlab(label="m/z")+
  ylab(label="# of formulas")

print(formula.plot)



####Reference Compounds Figure - Before and After Calibration####

library(ggplot2)

setwd("E:/Steven/Target Decoy/data/20190426")
msfilesraw<-list.files()
msfilesraw<-msfilesraw[1:5]
xdataraw<-NULL
xdataraw<-list()
for (i in 1:length(msfilesraw)){
  rawlist<-xcmsRaw(msfiles[i])
  xdataraw[i]<-rawlist
}

setwd("E:/Steven/Target Decoy/caldata")
msfilescal<-list.files()
msfilescal<-msfilescal[1:5]
xdatacal<-NULL
xdatacal<-list()
for (i in 1:length(msfilescal)){
  callist<-xcmsRaw(msfilescal[i])
  xdatacal[i]<-callist
}

azoref<-read.csv("E:/Steven/Target Decoy/dust analysis/azoref.csv", header=TRUE)
azoref<-azoref$Ref
reflistraw<-NULL
ppm<-15/(10^6)

for (m in 1:length(msfilesraw)){
      xraw<-xdataraw[[m]]
      plist<-NULL
      plist<-matrix(data=NA,nrow=length(azoref),ncol=2)
    for (p in 1:length(azoref)){
      indexref<-which(abs((xraw@env$mz-azoref[p])/azoref[p])<ppm)
      indexref2<-which.max(xraw@env$intensity[indexref])
      indexref<-indexref[indexref2]
      indexref<-indexref[1]
      mzmatch<-xraw@env$mz[indexref]
      plist[p,1]<-mzmatch
      if (is.na(indexref)){
        plist[p,1]<-100
        next}
      shift<-(xraw@env$mz[indexref]-azoref[p])/azoref[p]##calculate the mass shift
      plist[p,2]<-shift*10^6
      print(p)
      }
  reflistraw<-rbind(reflistraw,plist)
  colnames(reflistraw)<-c("mz","shift")
  }
setwd("E:/Steven/Target Decoy/products analysis")
write.table(reflistraw, file="azoref in raw2.csv", sep = ',',row.names=FALSE,col.names=c("mz","ppm"))

 
reflistcal<-NULL
ppm<-15/(10^6)

for (m in 1:length(msfilescal)){
  xcal<-xdatacal[[m]]
  plist<-NULL
  plist<-matrix(data=NA,nrow=length(azoref),ncol=2)
  for (p in 1:length(azoref)){
    indexref<-which(abs((xcal@env$mz-azoref[p])/azoref[p])<ppm)
    indexref2<-which.max(xcal@env$intensity[indexref])
    indexref<-indexref[indexref2]
    indexref<-indexref[1]
    mzmatch<-xcal@env$mz[indexref]
    plist[p,1]<-mzmatch
    if (is.na(indexref)){
      plist[p,1]<-100
      next}
    shift<-(xcal@env$mz[indexref]-azoref[p])/azoref[p]##calculate the mass shift
    plist[p,2]<-shift*10^6
    print(p)
  }
  reflistcal<-rbind(reflistcal,plist)
  colnames(reflistcal)<-c("mz","shift")
}
setwd("E:/Steven/Target Decoy/products analysis")
write.table(reflistcal, file="azoref in caldata.csv", sep = ',',row.names=FALSE,col.names=c("mz","ppm"))




