###################
### Read files ####
###################
iv_snp=c()
peer=0
  for(iter in iteration)
  {
    
    M <- M_matrix[,iter]
    

    load(‘#########’) ###load the genotype data of the cis-SNP of iterth gene
    snp_raw <- data.frame(snp_raw)
    A <- t(snp_raw[,geno_id %in% common])
    
   
    snps = SlicedData$new();
    A=as.matrix(A)
    A=t(A)
    snps$CreateFromMatrix(A)
    
    gene=SlicedData$new();
    res1=as.matrix(M)
    res1=t(res1)
    gene$CreateFromMatrix(res1)
    
    useModel = modelLINEAR; 
    pvOutputThreshold = 5e-1;
    errorCovariance = numeric();
    output_file_name = tempfile();
    me = Matrix_eQTL_engine(
      snps = snps,
      gene = gene,
      output_file_name = output_file_name,
      pvOutputThreshold = pvOutputThreshold,
      useModel = useModel,
      errorCovariance = errorCovariance,
      verbose = TRUE,
      pvalue.hist = TRUE,
      min.pv.by.genesnp = FALSE,
      noFDRsaveMemory = FALSE);
    time=as.numeric(me$time.in.sec)
    me=me$all$eqtls[1,]
    tmp=as.character(me$snps)
    onesnp=as.numeric(substr(tmp,4,nchar(tmp)))
    
    iv_snp=rbind(iv_snp,A[onesnp,])
    num_snp=rbind(num_snp,c(iter,num1,num2))
    
    
    l1 <- lm(M ~ A[onesnp,])
    l2 <- lm(Y ~ M)
    l3 <- lm(Y ~  A[onesnp,])


    summary <- rbind(summary, c(iter,summary(l1)$coeff[2,4], summary(l1)$coeff[2,1], summary(l3)$coeff[2,1],summary(l3)$coeff[2,4], summary(l3)$coeff[2,1]/summary(l1)$coeff[2,1], summary(l2)$coeff[2,1],summary(l2)$coeff[2,4],time))
    
    print(iter)
  }
  
  write.table(summary, paste("mr2_sim_", pc, ".txt", sep=""), quote=F, col.names=F, row.names=F)
  file=paste("mr_sim.RData", sep="")
  save(iv_snp,file=file)
