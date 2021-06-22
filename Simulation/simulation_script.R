library(AlphaSimR)
library(dplyr)

###################################### Simulate Population 1 #############################################

founderPop = runMacs2(nInd = 50, nChr = 15, segSites = 1000, ploidy = 6L, 
                      Ne = 50, bp = 2e+07, genLen = 1.43, mutRate = 2e-9)
#High Dominance
SP = SimParam$
  new(founderPop)$
  addTraitAD(250,meanDD=1,varDD=0.2)$
  addSnpChip(750)

SP$quadProb = 0.15

pop = newPop(founderPop)
F1 = randCross(pop,300)  

F1 = setPheno(F1, H2 = 0.5)

Pop1_SNP_matrix_HD = pullSnpGeno(F1)
Pop1_pheno_HD = F1@pheno

#Low Dominance
SP = SimParam$
  new(founderPop)$
  addTraitAD(250,meanDD=0.3,varDD=0.2)$
  addSnpChip(750)

SP$quadProb = 0.15

pop = newPop(founderPop)
F1 = randCross(pop,300)  

F1 = setPheno(F1, H2 = 0.5)

Pop1_SNP_matrix_LD = pullSnpGeno(F1)
Pop1_pheno_LD = F1@pheno

###################################### Simulate Population 2 #############################################

#Function to select simplex and nuliplex markers

getSimplexSNP = function(pop,sp,nSelSNPs,
                         simplex1_importance,nuliplex1_importance,
                         simplex0_importance,nuliplex0_importance,
                         alt_importance){
  
  tmp2 = (tmp == 1) #simplex alt
  tmp3 = (tmp == 0) #nuliplex alt
  tmp4 = (tmp == 5) #simplex ref
  tmp5 = (tmp == 6) #nuliplex ref
  
  #Number of simplex genotypes (alternative allele) per marker
  simplex = colSums(tmp2)
  simplex = as.data.frame(simplex)
  #Number of nuliplex genotypes (alternative allele) per marker
  simplex$nulli = colSums(tmp3)
  #Apply selection criteria 
  simplex$ind = simplex1_importance*simplex$simplex + nuliplex1_importance*simplex$nulli
  simplex = simplex %>% arrange(-ind)
  
  #Number of simplex genotypes (reference allele) per marker
  simplex2 = colSums(tmp4)
  simplex2 = as.data.frame(simplex2)
  #Number of nuliplex genotypes (reference allele) per marker
  simplex2$nulli = colSums(tmp5)
  #Apply selection criteria 
  simplex2$ind = simplex0_importance*simplex2$simplex + nuliplex0_importance*simplex2$nulli
  simplex2 = simplex2 %>% arrange(-ind)
  
  #Combine 
  nAlt = nSelSNPs*alt_importance
  nRef = nSelSNPs*(alt_importance-1)
  
  tmp = tmp[,c(rownames(simplex)[1:nAlt],rownames(simplex2)[1:nRef])]
}


#Simulation

founderPop = runMacs2(nInd = 50, nChr = 15, segSites = 5000, ploidy = 6L, 
                      Ne = 50, bp = 2e+07, genLen = 1.43, mutRate = 2e-7)
#High Dominance
SP = SimParam$
  new(founderPop)$
  addTraitAD(250,meanDD=1,varDD=0.2)$
  addSnpChip(4750)

SP$quadProb = 0.15

pop = newPop(founderPop)
F1 = randCross(pop,300)  

F1 = setPheno(F1, H2 = 0.5)

#Get the SNPs of interest
#The importance values were chosen ad-hoc, in order to mimic the genotype frequencies found in the Sweet Potato dataset
Pop2_SNP_matrix_HD = getSimplexSNP(pop = F1, sp = SP, nSelSNPs = 11250,
                                   simplex1_importance = 2, nuliplex1_importance = 1,
                                   simplex0_importance = 1.5, nuliplex0_importance = 1,
                                   alt_importance = 0.8)
Pop2_pheno_HD = F1@pheno

#Low Dominance
SP = SimParam$
  new(founderPop)$
  addTraitAD(250,meanDD=0.3,varDD=0.2)$
  addSnpChip(4750)

SP$quadProb = 0.15

F1 = setPheno(F1, H2 = 0.5)

#Get the SNPs of interest
#The importance values were chosen ad-hoc, in order to mimic the genotype frequencies found in the Sweet Potato dataset
Pop2_SNP_matrix_LD = getSimplexSNP(pop = F1, sp = SP, nSelSNPs = 11250,
                                   simplex1_importance = 2, nuliplex1_importance = 1,
                                   simplex0_importance = 1.5, nuliplex0_importance = 1,
                                   alt_importance = 0.8)
Pop2_pheno_LD = F1@pheno



