args <- commandArgs(TRUE)
args <- as.numeric(args)
chr=args[1]
d=read.table('gene.bed')
library(data.table)
library(doParallel)
numCore = 10

idx=which(d[,1]==chr)
d=d[idx,]

registerDoParallel(cores=numCore)
N=nrow(d)
resBMM <- foreach(i=1:N, .combine=rbind, .errorhandling = 'remove')%dopar%
{
  res <- data.frame()
name=d[i,4]  
name=as.character(name)
file=paste0('/net/mulan/disk2/yuef/data/GTEX/GTEx_v7/qc/genotype/',name)
if(file.exists(file))
{
  snp_raw=fread(file)
  snp_raw=snp_raw[,7:641]
  sfile=paste0('/net/mulan/disk2/yuef/data/GTEX/GTEx_v7/qc/genotype1/',name,'.RData')
  save(snp_raw,file=sfile)
}
print(i)
res <- data.frame(iter=i)
return(res)
}


