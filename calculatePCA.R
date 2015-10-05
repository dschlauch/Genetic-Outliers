macDir <- "counts_20_50"
outputDir <- "output_20_50"
aveAlleles <- 30

# get the total allele counts for each sample (haplotype)
counts <- Reduce('+',sapply(file.path('~/1000GP',macDir, list.files(path=file.path('~/1000GP',macDir), pattern="*.csv")), read.csv, header=F))

allU <- read.csv(file.path('~/1000GP',outputDir,'all_U.csv'), row.names=1)
allV <- read.csv(file.path('~/1000GP',outputDir,'all_V.csv'), row.names=1)

# estimated number of loci
# average alleles per loci=denominator
numLoci <- sum(diag(as.matrix(allV)))/aveAlleles

numX1Y0 <- counts-allV
numX0Y1 <- t(numX1Y0)
numX1Y1 <- allV
numX0Y0 <- numLoci - numX1Y0 - numX0Y1 + numX1Y1

alleleProportions <- counts/numLoci

valueX1Y0 <- (1-alleleProportions)%*%t(0-alleleProportions)
valueX0Y1 <- (0-alleleProportions)%*%t(1-alleleProportions)
valueX1Y1 <- (1-alleleProportions)%*%t(1-alleleProportions)
valueX0Y0 <- (0-alleleProportions)%*%t(0-alleleProportions)

varcovMatrixHap <- (1/(numLoci-1))*(numX0Y0*valueX0Y0 + numX0Y1*valueX0Y1 + numX1Y0*valueX1Y0 + numX1Y1*valueX1Y1)
varcovMatrix <- (varcovMatrixHap[c(T,F),c(T,F)] + varcovMatrixHap[c(T,F),c(F,T)] + varcovMatrixHap[c(F,T),c(T,F)] + varcovMatrixHap[c(F,T),c(F,T)])/4
