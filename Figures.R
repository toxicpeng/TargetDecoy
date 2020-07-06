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

chlor<-read.csv("C:/Users/Steven Desktop/Documents/MSc2/TargetDecoy/analysis/dust analysis/ChlorAll.csv", header=TRUE)

# Scale data
chlor.scaled <- chlor[,c(1:27,43,46:48)]
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
chlor.scaled <- chlor[,c(1:27,43,46:48)]
chlor.scaled[, c(3:26)] <- scale(chlor.scaled[, 3:26])


chlor<-read.csv("C:/Users/Steven Desktop/Documents/MSc2/TargetDecoy/analysis/dust analysis/ChlorAllTransposed.csv", header=TRUE)

chlor.matrix <- as.matrix(chlor[,2:43])
rownames(chlor.matrix) <- chlor[,1]
#melt_chlor<-melt(chlor.matrix)
#cast_chlor<-acast(melt_chlor, Var2~Var1)
#chlor.cormatrix<-rcorr(cast_chlor)
chlor.cormatrix<-rcorr(chlor.matrix)
colour1 <- colorRampPalette(c("red","white","blue"))
par(xpd=TRUE)
corrplot(
  chlor.cormatrix$r, 
  method = "square", 
  type = "full", 
  order = "hclust",
  col = colour1(200), 
  tl.col = "black", 
  tl.srt = 90, 
  tl.cex = 0.5,
  tl.pos = "tl",
  cl.cex = 1.7, #Font size of color index labels
  cl.ratio = 0.25, #Width of color index + labels
  p.mat = chlor.cormatrix$P, 
  sig.level = 0.05, 
  insig = "blank", 
  addrect = 10, #Number of square groupings to outline 
  mar = c(1,1,1,1) #Outer graph margins
  )


#####Boxplot for All House Dust Results#####
library(scales)
library(ggthemes)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)


dustdata<-read.csv("C:/Users/Steven Desktop/Documents/MSc2/TargetDecoy/analysis/dust analysis/DustCompoundsUnique.csv", header=TRUE)

####Attempt at calculating my own boxplot values (Not currently used)####
dustcpds<-dustdata[,c(43:46,48,51)]
dustcpds$AvgIntNew<-NA
dustvalues<-dustdata[,3:26]
for(i in 1:nrow(dustvalues)){
  intset<-as.numeric(dustvalues[i,])
  zeroind<-which(intset==0)
  intset<-intset[-zeroind]
  upperbound<-quantile(intset)[4]+(1.5*IQR(intset))
  lowerbound<-quantile(intset)[2]-(1.5*IQR(intset))
  outlierupp<-which(intset>upperbound)
  outlierlow<-which(intset<lowerbound)
  outliers<-c(outlierupp,outlierlow)
  if(length(outliers)==0){
    dustcpds$AvgIntNew[i]<-dustcpds$AvgInt[i]
    next
    }
  intset<-intset[-outliers]
  averageint<-mean(intset)
  dustcpds$AvgIntNew[i]<-averageint
}
dustcpds$AvgIntNew<-log10(dustcpds$AvgIntNew)
dustcpds<-dustcpds[order(-dustcpds$AvgIntNew),]
CpdOrder<-seq(from=1, to=nrow(dustcpds), by=1)
####Actual Code####

# Is a grouping available?
# (Will return TRUE if an explicit group or a discrete variable with only one
# level existed when add_group() was called.)
NO_GROUP <- -1L
has_groups <- function(data) {
  # If no group aesthetic is specified, all values of the group column equal to
  # NO_GROUP. On the other hand, if a group aesthetic is specified, all values
  # are different from NO_GROUP (since they are a result of plyr::id()). NA is
  # returned for 0-row data frames.
  data$group[1L] != NO_GROUP
}


dustboxdata<-dustdata[,c(3:26,43:44,46:49,52)]
dustboxdata<-dustboxdata %>%
  gather(key="Sample",value="Intensity",-Formula,-MedInt,-Name,-Class, -ClAzoDye, -CpdID,-IUPAC)
zeroindex<-which(dustboxdata$Intensity==0)
dustboxdata<-dustboxdata[-zeroindex,]
dustboxdata$Formula<-as.character(dustboxdata$Formula)

##Boxplot of all dust compounds, ordered by median intensity
dustplotall<-ggplot(data = dustboxdata, mapping=aes(x=CpdID, y=Intensity, group=CpdID))+
  geom_boxplot(size=0.1, outlier.size=0.1)+
  scale_y_log10(breaks=trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
  scale_x_continuous(breaks=seq(from = 0, to = 600, by = 100),limits=c(0,600))+
  theme_classic()+
  xlab("Compound Rank by Median Intensity")+
  ylab("Signal Intensity")
print(dustplotall)

##Exact same graph as above, but in Tufte style
dustplotall<-ggplot(data = dustboxdata, mapping=aes(x=CpdID, y=Intensity, group=CpdID))+
  theme_classic()+
  geom_tufteboxplot()+
  scale_y_log10(breaks=trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
  scale_x_continuous(breaks=seq(from = 0, to = 600, by = 100),limits=c(0,50))+
  xlab("Compound Rank by Median Intensity")+
  ylab("Signal Intensity")
print(dustplotall)

##Just the top 50 compounds (by intensity)
top50<-which(dustboxdata$CpdID<=50)
dusttopdata<-dustboxdata[top50,]
dustplottop<-ggplot(data = dusttopdata, mapping=aes(x=CpdID, y=Intensity, group=CpdID))+
  geom_boxplot(size=0.1, outlier.size=1)+
  scale_y_log10(breaks=trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
  scale_x_continuous(breaks=seq(from = 0, to = 50, by = 10))+
  theme_classic()+
  xlab("Compound Rank by Median Intensity")+
  ylab("Signal Intensity")
print(dustplottop)
  

##Chlorine-containing compounds only
Cldustdata<-read.csv("C:/Users/Steven Desktop/Documents/MSc2/TargetDecoy/analysis/dust analysis/ChlorAllLog.csv", header=TRUE)

Clboxdata<-Cldustdata[,c(3:27,44,46:50,53)]
Clboxdata<-Clboxdata %>%
  gather(key="Sample",value="Intensity",-Formula,-MedInt,-Name,-Class, -ClAzoDye, -CpdID,-ClID, -IUPAC)
zeroindex<-which(Clboxdata$Intensity==0)
Clboxdata<-Clboxdata[-zeroindex,]
Clboxdata$Formula<-as.character(Clboxdata$Formula)

ClSorted<-Clboxdata[order(Clboxdata$ClID),]

dustplotCl<-ggplot(data = ClSorted, mapping=aes(x=ClID, y=Intensity, group=ClID))+
  geom_boxplot(outlier.size=1)+
  #scale_y_log10(labels = trans_format("log10", math_format(10^.x)))+ #***For some reason, the data is being log-transformed before the boxplot calculations are made.  This is why some of the boxes appear to be out of order. 
  scale_y_continuous(breaks = c(4,5,6,7,8), labels = c(expression(10^4),expression(10^5),expression(10^6),expression(10^7),expression(10^8)))+
  scale_x_continuous(breaks=seq(from = 0, to = 50, by = 10))+
  theme_Steven()+
  xlab("Compound Rank by Median Intensity")+
  ylab("Signal Intensity")
print(dustplotCl)



##Chlorinated Azo Dyes only
ClAzo<-which(dustboxdata$ClAzoDye==TRUE)
dustClAzodata<-dustboxdata[ClAzo,]
ClAzoSorted<-dustClAzodata[order(dustClAzodata$CpdID),]
newcount<-1
for(i in 1:nrow(ClAzoSorted)){
  if(ClAzoSorted$CpdID[i]<30){next}
  AzoID<-which(ClAzoSorted$CpdID==ClAzoSorted$CpdID[i])
  ClAzoSorted$CpdID[AzoID]<-newcount
  newcount<-newcount+1  
}

dustplotClAzo<-ggplot(data = ClAzoSorted, mapping=aes(x=CpdID, y=Intensity, group=CpdID))+
  geom_boxplot(outlier.size=1)+
  scale_y_log10(breaks=trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
  scale_x_continuous(breaks=seq(from = 0, to = 600, by = 100))+
  theme_classic()+
  xlab("Compound Rank by Median Intensity")+
  ylab("Signal Intensity")
print(dustplotClAzo)





####Figures comparing abundance of Br azo dyes and their chlorine analogues (boxplot and scatterplots)####

library(scales)
library(ggthemes)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)


dyedata<-read.csv("C:/Users/Steven Desktop/Google Drive/MSc 3/Algorithm Paper/Cl_vs_Br_Azo_Dyes_PeakAreas.csv", header=TRUE)

indexDB373<-which(dyedata$DB373<1000000)
indexCl373<-which(dyedata$Cl373<1000000)
dyedata373<-dyedata[-c(indexDB373,indexCl373),2:3]
dyedata373<-log10(dyedata373)

trendline373<-lm(dyedata373$Cl373~dyedata373$DB373, data=dyedata373) #Can use 'summary.lm()' after this line to obtain R^2 and p-value
confribbon<-predict(trendline373,interval = "confidence", level=0.95)
trendline373<-predict(trendline373, level=0.95)

dyedata373<-as.data.frame(10^(dyedata373))
trendline373<-10^(trendline373)
confribbon<-10^(confribbon)

dyes373plot<-ggplot(data = dyedata373, mapping=aes(x=DB373, y=Cl373))+
  geom_line(aes(y=trendline373, x=DB373), colour = "red", size=1)+
  geom_ribbon(aes(ymin=confribbon[,2],ymax=confribbon[,3]), fill="red", alpha=0.1)+
  geom_point(colour="blue", size=2)+
  scale_x_log10(breaks=c(1e7,1e8),labels = c(expression(10^7),expression(10^8)))+
  scale_y_log10(breaks=c(1e7,1e8,1e9),labels = c(expression(10^7),expression(10^8),expression(10^9)))+
  annotation_logticks(sides="lb")+
  xlab("Br-DB373")+
  ylab("Cl-DB373")+
  theme_Steven()

print(dyes373plot)


indexDV93<-which(dyedata$DV93<1000000)
indexCl93<-which(dyedata$Cl93<1000000)
dyedata93<-dyedata[-c(indexDV93,indexCl93),4:5]
dyedata93<-log10(dyedata93)

trendline93<-lm(dyedata93$Cl93~dyedata93$DV93, data=dyedata93) #Can use 'summary.lm()' after this line to obtain R^2 and p-value
confribbon<-predict(trendline93,interval = "confidence", level=0.95)
trendline93<-predict(trendline93, level=0.95)

dyedata93<-as.data.frame(10^(dyedata93))
trendline93<-10^(trendline93)
confribbon<-10^(confribbon)

dyes93plot<-ggplot(data = dyedata93, mapping=aes(x=DV93, y=Cl93))+
  geom_line(aes(y=trendline93, x=DV93), colour = "red", size=1)+
  geom_ribbon(aes(ymin=confribbon[,2],ymax=confribbon[,3]), fill="red", alpha=0.1)+
  geom_point(colour="blue", size=2)+
  scale_x_log10(breaks=c(1e7,1e8,1e9),labels = c(expression(10^7),expression(10^8),expression(10^9)))+
  scale_y_log10(breaks=c(1e8,1e9),labels = c(expression(10^8),expression(10^9)))+
  annotation_logticks(sides="lb")+
  xlab("Br-DV93")+
  ylab("Cl-DV93")+
  theme_Steven()

print(dyes93plot)


indexBDNA<-which(dyedata$BDNA<100000)
indexCDNA<-which(dyedata$CDNA<100000)
dyedataDNA<-dyedata[-c(indexBDNA,indexCDNA),6:7]
dyedataDNA<-log10(dyedataDNA)

trendlineDNA<-lm(dyedataDNA$CDNA~dyedataDNA$BDNA, data=dyedataDNA) #Can use 'summary.lm()' after this line to obtain R^2 and p-value
confribbon<-predict(trendlineDNA,interval = "confidence", level=0.95)
trendlineDNA<-predict(trendlineDNA, level=0.95)

dyedataDNA<-as.data.frame(10^(dyedataDNA))
trendlineDNA<-10^(trendlineDNA)
confribbon<-10^(confribbon)

dyesDNAplot<-ggplot(data = dyedataDNA, mapping=aes(x=BDNA, y=CDNA))+
  geom_line(aes(y=trendlineDNA, x=BDNA), colour = "red", size=1)+
  geom_ribbon(aes(ymin=confribbon[,2],ymax=confribbon[,3]), fill="red", alpha=0.1)+
  geom_point(colour="blue", size=2)+
  scale_x_log10(breaks=c(1e6),labels = c(expression(10^6)))+
  scale_y_log10(breaks=c(1e6,1e7),labels = c(expression(10^6),expression(10^7)))+
  annotation_logticks(sides="lb")+
  xlab("Br-DNA")+
  ylab("Cl-DNA")+
  theme_Steven()

print(dyesDNAplot)


BrCldata<-dyedata
BrCldata<-BrCldata[-c(1:2),]
BrClratios<-data.frame(DB373=(BrCldata[,3]/BrCldata[,2]), DV93=(BrCldata[,5]/BrCldata[,4]), DNA=(BrCldata[,7]/BrCldata[,6]))

BrClratios.long<-gather(BrClratios, Compound, Ratios)
BrClratios.long$Compound<-factor(BrClratios.long$Compound, levels = c("DB373","DV93","DNA"))



boxplotBrCl<-ggplot(data = BrClratios.long, mapping=aes(x=Compound, y=Ratios, group=Compound))+
  geom_boxplot(outlier.size=2,lwd=1,width=0.5)+
  scale_y_continuous(name="Cl:Br Ratios", trans="log10", breaks=c(1,10,100))+
  #coord_trans(y="log10")+
  theme_Steven()+
  annotation_logticks(sides="l", scaled=TRUE, size=1)
  #theme(axis.ticks=element_line())
  #ylab("Cl:Br Ratios")
print(boxplotBrCl)

####Multi-Part Bar Graph for All Chlorine-Containing Compounds####

library(ggplot2)
library(tidyr)
library(ggdendro)
library(reshape2)
library(ggthemes)

clcpds<-read.csv("C:/Users/Steven Desktop/Documents/MSc2/TargetDecoy/analysis/dust analysis/ChlorAll.csv", header=TRUE)
clcpds$avg<-log10(clcpds$avg)

grid.plot <- ggplot(data=clcpds, mapping=aes(x=CpdID, y=avg))+
                  geom_bar(stat="identity")+
                  coord_cartesian(ylim=c(3,7))+
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





####Custom ggplot2 themes####
theme_Steven <- function () { 
  theme_classic(base_size=20, base_family="sans") %+replace% 
    theme(
      axis.text=element_text(color="black")
    )
}
