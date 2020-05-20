args <- commandArgs(TRUE)
chr_no <- as.integer(args[1])

#library(data.table)

bedfile<-data.frame(read.table(paste(".../chr", chr_no, ".bed", sep=""))) ### read the genotype data

MAF1=bedfile[,9]
MAF2=bedfile[,10]

idx1=substr(MAF1,1,3)
idxx1=which(idx1=='MAF')
idx2=substr(MAF2,1,3)
idxx2=which(idx2=='MAF')
#MAF1=MAF1[idxx1]
#MAF2=MAF2[idxx2]
MAF1=as.character(MAF1)
MAF2=as.character(MAF2)
maf1=substr(MAF1,5,nchar(MAF1))
maf2=substr(MAF2,5,nchar(MAF2))
maf1=as.numeric(maf1)
maf2=as.numeric(maf2)

maf1[idxx2]=0
maf2[idxx1]=0

idx1=which(maf1>0.05)
idx2=which(maf2>0.05)
idx=union(idx1,idx2)

maf2=MAF2[idxx1]
maf2=substr(maf2,4,nchar(maf2))
maf2=as.numeric(maf2)
idxx=which(maf2<0.99)

idx=setdiff(idx,idxx1[idxx])

idx=sort(idx)

write.table(bedfile[idx,], paste("./chr", chr_no, ".bed", sep=""), quote=F, col.names=F, row.names=F)

write.table(bedfile[idx,3], paste("./snp_list", chr_no, ".txt", sep=""), quote=F, col.names=F, row.names=F)




