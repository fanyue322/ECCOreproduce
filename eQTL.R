
args <- commandArgs(TRUE)
peer <- as.integer(args[1])
###################
### Read files ####
###################
  for(iter in 1:10000)
  {
    M <- M_matrix[,iter]
    res1=M
    
    load(‘#########’) ###load the genotype data of the cis-SNP of iterth gene
    snp_raw <- data.frame(snp_raw)
    A <- t(snp_raw[,geno_id %in% common])
    
    snps = SlicedData$new();
    A=as.matrix(A)
    A=t(A)
    snps$CreateFromMatrix(A)
    
    gene=SlicedData$new();
    res1=as.matrix(res1)
    #M=sample(M)
    res1=t(res1)
    gene$CreateFromMatrix(res1)
    
    useModel = modelLINEAR; 
    pvOutputThreshold = 5e-1;
    errorCovariance = numeric();
    output_file_name = tempfile();
    me = Matrix_eQTL_engine(
      snps = snps,
      gene = gene,
      #  cvrt = cvrt,
      output_file_name = output_file_name,
      pvOutputThreshold = pvOutputThreshold,
      useModel = useModel,
      errorCovariance = errorCovariance,
      verbose = TRUE,
      pvalue.hist = TRUE,
      min.pv.by.genesnp = FALSE,
      noFDRsaveMemory = FALSE);
    
    me=me$all$eqtls[1,]
    time=1
    summary <- rbind(summary, c(iter, pc, me[6], me[4],time))
  }
  
  
  write.table(summary, paste("eQTL_", peer, ".txt", sep=""), quote=F, col.names=F, row.names=F)
