n = 491
n_genes = 10000
n_conf = 10
pve1 = 0.03; pve2 = 0.5; pve3 = 0.25; 
C <- matrix(rnorm(n*n_conf, 0, sd = 1), nrow=n)
M_matrix <- c(); Y_matrix <- c()
g_matrix <- c(); beta_M1_matrix <- c()
IV_mark<-c()



gene <- fread(paste("/net/mulan/disk2/yuef/data/GTEX/GTEx_v7/expression/", tissue, "_expression.txt", sep=""), header=T)
gene_id <- colnames(gene)[-c(1:4)]
gene_names <- gene$gene_id
common <- gene_id
num=length(gene_names)
idx=sample(1:num)
num=1:n_genes


gene_names=gene_names[num]

for(i in 1:n_genes)
{
  #idx=sample(1:10)[1:5] #heterogeneous confounding scenario
  file=paste("/net/mulan/disk2/yuef/data/GTEX/GTEx_v7/qc/", tissue, "/", gene_names[i], sep="")
  snp_raw <- fread(file)
  snp_raw <- data.frame(snp_raw[,7:641])
  A <- t(snp_raw[,geno_id %in% common])
  num_snp=ncol(A)
  
 
  j=sample(1:num_snp,1)
  g1 <- A[,j]

  beta_g1 <- rnorm(1, 0, sd = 1)
  beta_g1 <- beta_g1/(sd(g1 * beta_g1)) * sqrt(pve1)

  beta_C <- rnorm(n_conf, 0, sd = 1)
  beta_C <- beta_C/(sd(C %*% beta_C)) * sqrt(pve2)
  #  beta_C[idx] <- beta_C[idx]/(sd(C[,idx] %*% beta_C[idx])) * sqrt(pve2)# heterogeneous confounding scenario
  e1 <- rnorm(n, 0, sd = sqrt(1-pve1-pve2))
  
  g_matrix <- cbind(g_matrix, g1 * beta_g1 + e1)
  M_matrix <- cbind(M_matrix, g1 * beta_g1 + C %*% beta_C + e1)
 # M_matrix <- cbind(M_matrix, g1 * beta_g1 + C[,idx] %*% beta_C[idx] + e1)# heterogeneous confounding scenario
  IV_mark<-c(IV_mark,j)
}

gamma <- array(0, n_genes)
gamma[1:100] <- rnorm(100, 0, 1)
gamma <- gamma/(sd(g_matrix %*% gamma)) * sqrt(pve_M)
e2 <- rnorm(n, 0, sd = sqrt(1-pve_M))
Y <- g_matrix %*% gamma + e2

kmatrix <- M_matrix %*% t(M_matrix)
kmatrix <- as.matrix(kmatrix)
k.pca <- prcomp(kmatrix,
                center = TRUE,
                scale = TRUE)
pc_matrix <- predict(k.pca,
                     newdata=kmatrix)
pc_matrix <- scale(pc_matrix, center = TRUE, scale = TRUE)

save(pc_matrix,Y,A,M_matrix,file='sim_pc10_sparse.RData')

gamma <- array(0, n_genes)
gamma[1:10000] <- rnorm(10000, 0, 1)
gamma <- gamma/(sd(g_matrix %*% gamma)) * sqrt(pve_M)
e2 <- rnorm(n, 0, sd = sqrt(1-pve_M))
Y <- g_matrix %*% gamma + e2

kmatrix <- M_matrix %*% t(M_matrix)
kmatrix <- as.matrix(kmatrix)
k.pca <- prcomp(kmatrix,
                center = TRUE,
                scale = TRUE)
pc_matrix <- predict(k.pca,
                     newdata=kmatrix)
pc_matrix <- scale(pc_matrix, center = TRUE, scale = TRUE)

save(pc_matrix,Y,A,M_matrix,file='sim_pc10_poly.RData')
