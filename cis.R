args <- commandArgs(TRUE)
tissue <- as.character(args[1])

library(data.table)
gene_info <- data.frame(fread(paste(tissue, "_gene_info.txt", sep=""), header=T))
for (i in 1:nrow(gene_info)) {
  gene_info[i,2] <- max(0, gene_info[i,2]-10^6)
  gene_info[i,3] <- gene_info[i,3]+10^6
}
write.table(gene_info, paste(tissue, "_gene_info1.bed", sep=""), quote=F, col.names=F, row.names=F)
