#-----------------------------
#ggplot histgram
#-----------------------------
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

######single############
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

