load("Sugarcane_data.RData")

## Preparing Gmatrix - allele dosage

#Get allele frequencies
allele_freq <- as.data.frame(rowSums(SNP_matrix, na.rm = T)/rowSums(Ploidy_matrix, na.rm = T))[,1]

#Scale SNP genotypes between 0 and 2
Wtmp = 2*SNP_matrix/(Ploidy_matrix) 

#Center SNP matrix
W = Wtmp - (2*allele_freq)

#Get G matrix
WWt = t(W) %*% W
denom = 2*allele_freq * (1 - allele_freq)
denom = sum(denom)

Gmatrix = WWt / denom

#If the matrix is not invertible - add very small constant (10e-3) to diagonals
Gmatrix = Gmatrix + 0.001*diag(nrow(Gmatrix))

## Preparing digenic dominance matrix 

#Get matrix with number of possible pair combinations (C):
C_matrix = apply(Ploidy_matrix,1:2, function(x) length(combn(x,2))/2)

#Get digenic dominance matrix (Q):

Q = (allele_freq^2 * C_matrix) - 
  (Ploidy_matrix - 1) * allele_freq * SNP_matrix + 
  0.5 * (SNP_matrix) * (SNP_matrix-1) 

denomDom = sum(C_matrix[,1]*allele_freq^2*(1- allele_freq)^2) 
D = t(Q)%*% Q
Dmatrix = D/denomDom

#If the matrix is not invertible - add very small constant (10e-3) to diagonals
Dmatrix = Dmatrix + 0.001*diag(nrow(Dmatrix))

