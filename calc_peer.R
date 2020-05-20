args <- commandArgs(TRUE)
tissue <- as.numeric(args[1])
pc<-as.numeric(args[2])
library('data.table')
tissue_info <- fread("######") ###load a file contains the names of 48 tissues


tissue=tissue_info[tissue]

library(data.table)
library(peer)
dt <- data.frame(fread(paste(tissue, "_expression.txt", sep=""), header = T))
dt <- dt[,-c(1:4)]; dt <- t(as.matrix(dt))
#for (pc in c(80,90,100)) {
#for (pc in 1 2 5 10 15 20 25 30 35 40 50 60 70 80 90 100) {
  model = PEER()
  PEER_setPhenoMean(model,as.matrix(dt))
  dim(PEER_getPhenoMean(model))

  PEER_setAdd_mean(model, TRUE)
  PEER_setNk(model,pc)  ## pc hidden confounders
  PEER_getNk(model)

  PEER_update(model)
  #PEER_setNMax_iterations(model, 10000)

  factors = PEER_getX(model)
  factors=factors[,-1]
  #dim(factors)
  residuals = PEER_getResiduals(model)
  #dim(residuals)
  #write.table(factors, paste("peer_factor_",tissue,'_', pc, ".txt", sep=""), quote=F, col.names=F, row.names=F)
write.table(residuals, paste(tissue,'_peer', pc, ".txt", sep=""), quote=F, col.names=F, row.names=F)
  print (pc)
#}
