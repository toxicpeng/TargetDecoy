library(rcdk)


mysmiles<-'CCCCOC(=O)CC(CC(=O)OCCCC)(C(=O)OCCCC)O'
cpd<-parse.smiles(mysmiles)
get.xlogp(cpd[[1]])