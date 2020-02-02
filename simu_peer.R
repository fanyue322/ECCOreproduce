args <- commandArgs(TRUE)
peer <- as.integer(args[1])
###################
### Read files ####
###################
  for(iter in 1:N)
  {
    ind=num_snp[iter,1]
    res1=gene_peer[,ind]
    M <- M_matrix[,ind]
    A=iv_snp[iter,]
    
    
    t1=Sys.time()
    l1 <- lm(M ~ A)
    
    l2 <- lm(Y ~ res1)
    
    l3 <- lm(Y ~ A)
    liv=ivreg(Y~M,~A)
    t2=Sys.time()
    time=as.numeric(t2-t1)
    
    summary <- rbind(summary, c(iter,summary(l1)$coeff[2,4], summary(l1)$coeff[2,1], summary(l3)$coeff[2,1],summary(l3)$coeff[2,4], summary(l3)$coeff[2,1]/summary(l1)$coeff[2,1],summary(liv)$coeff[2,4], summary(l2)$coeff[2,1],summary(l2)$coeff[2,4],time))
    
    
  }
  
  write.table(summary, paste("mr2_sim_peer", peer, ".txt", sep=""), quote=F, col.names=F, row.names=F)
  
