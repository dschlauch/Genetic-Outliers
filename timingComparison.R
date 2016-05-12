library(data.table)

speedTest <- function(genotypesSubpop){
    
    numSamples <- ncol(genotypesSubpop)
    numVariants <- nrow(genotypesSubpop)
    sumVariants <- rowSums(genotypesSubpop)
  
    # remove < n variants
    sumVariants <- rowSums(genotypesSubpop)
    genotypesSubpop <- as.matrix(genotypesSubpop)
    
    print("Number of used variants")
    print(nrow(genotypesSubpop))
    numFilteredVariants <- nrow(genotypesSubpop)
    sumFilteredVariants <- rowSums(genotypesSubpop)

    totalPossiblePairs <- choose(numSamples,2)
    totalPairs <- choose(sumFilteredVariants,2)
    weights <- totalPossiblePairs/totalPairs
    p <- 1/weights
    
    var_s_hap <- sum(weights-1)/(numFilteredVariants^2)
    
    # Calculate expected values conditional on kinship
    pkweightsMean <- mean(((sumFilteredVariants-2)/numSamples)*weights)
    
    s_matrix_numerator <- t(genotypesSubpop*weights)%*%genotypesSubpop
    s_matrix_denominator <- numFilteredVariants
    s_matrix_hap <- s_matrix_numerator/s_matrix_denominator
}

numSamples <- 2000
numVariants <- 100000

genotypes <- data.table(matrix(rbinom(numSamples*numVariants,1,.2), ncol=numSamples))

print(system.time(speedTest(genotypes)))
# print(system.time(princomp(genotypes)))
# print(system.time(eigen(cor(genotypes))))
# print(system.time(svd(cor(genotypes))))
print(system.time(cor(genotypes)))
