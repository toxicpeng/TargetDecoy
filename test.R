data(sdfsample)
sdf1 <- smiles2sdf(as.character(MS2file$smiles[4]))
sdf2 <- smiles2sdf(as.character(smiles[4]))
sdfsample[[1]]<-sdf1
sdfsample[[2]]<-sdf2
K <- sd2gram(sdfsample)
heatmap(K,Rowv=NA,Colv=NA,scale="none")