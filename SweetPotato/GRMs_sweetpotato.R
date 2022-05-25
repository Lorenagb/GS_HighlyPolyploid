load("SweetPotato_data.RData")

#Additive matrix (hexa)
allele_freq = colSums(SNP_matrix)/(6*nrow(SNP_matrix))
W = t(SNP_matrix) - 6*allele_freq

WWt = crossprod(W)
denom = sum(6*allele_freq * (1 - allele_freq))
Gmatrix = WWt/denom

#If the matrix is not invertible - add very small constant (10e-3) to diagonals
Gmatrix = Gmatrix + 0.001*diag(nrow(Gmatrix))

#Digenic dominance matrix (hexa)
C_matrix = matrix(length(combn(6,2))/2, nrow = nrow(t(SNP_matrix)), ncol = ncol(t(SNP_matrix)))
Ploidy_matrix = matrix(6, nrow = nrow(t(SNP_matrix)), ncol = ncol(t(SNP_matrix)))

Q = (allele_freq^2 * C_matrix) - 
  (Ploidy_matrix - 1) * allele_freq * t(SNP_matrix) + 
  0.5 * t(SNP_matrix) * (t(SNP_matrix)-1) 

D = crossprod(Q)
denomDom = sum(C_matrix[,1]*allele_freq^2*(1- allele_freq)^2) 
Dmatrix = D/denomDom

#If the matrix is not invertible - add very small constant (10e-3) to diagonals
Dmatrix = Dmatrix + 0.001*diag(nrow(Dmatrix))
