load("SweetPotato_data.RData")

#Additive matrix (hexa)
allele_freq = colSums(SNP_matrix)/(6*nrow(SNP_matrix))
W = t(SNP_matrix) - 6*allele_freq

WWt = crossprod(W)
denom = sum(6*allele_freq * (1 - allele_freq))
Gmatrix = WWt/denom

#Digenic dominance matrix (hexa)
C_matrix = matrix(length(combn(6,2))/2, nrow = nrow(SNP_matrix), ncol = ncol(SNP_matrix))
Ploidy_matrix = matrix(6, nrow = nrow(SNP_matrix), ncol = ncol(SNP_matrix))

Q = (allele_freq^2 * C_matrix) - 
  (Ploidy_matrix - 1) * allele_freq * SNP_matrix + 
  0.5 * (SNP_matrix) * (SNP_matrix-1) 

D = t(Q)%*% Q
denomDom = sum(C_matrix[,1]*allele_freq^2*(1- allele_freq)^2) 
Dmatrix = D/denomDom