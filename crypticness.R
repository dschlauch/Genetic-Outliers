library(readr)
library(ggplot2)
dat <- read_table()
read_lines("./data/1000GP_Phase3_chr22.hap.gz", skip = 0, n_max = 100, locale = default_locale(),           progress = interactive())


#### Simulated data
simn <- 3000
numVariants <- 100000
sim.w.numerator <- 2*(2*simn-1)
simMAF <- runif(numVariants,0,.5)
genotypes <- sapply(1:simn,function(x){
    rbinom(numVariants,1,simMAF)
})

# Real data
sample <- read.table("~/1000GP/data/1000GP_Phase3.sample", sep=" ", header=T)
sampleIDs <- as.character(sample[,1])
pop <- as.character(sample[,2])
group <- as.character(sample[,3])
sex <- as.character(sample[,4])
hap.pop <- rep(pop,each=2)

calculateSMatrix <- function(subpop="CEU", numberOfLines=5000, minVariants=5){
    filename <- "./data/combinedFiltered1000.gz"
    con <- file(filename, "rt")
    system.time(genotypes <- apply(do.call(cbind, strsplit(readLines(con, numberOfLines)," ")), 1,as.numeric)[,hap.pop=="STU"])
    system.time(genotypes <- apply(do.call(cbind, strsplit(readLines(con, numberOfLines)," ")), 1,as.numeric)[,hap.pop=="STU"])
    
    numSamples <- ncol(genotypes)
    sumVariants <- rowSums(genotypes)
    # reverse so that MAF<.5
    genotypes[sumVariants>(numSamples/2),] <- 1-genotypes[sumVariants>(numSamples/2),]
    sumVariants <- rowSums(genotypes)
    # remove < n variants
    genotypes <- genotypes[sumVariants>minVariants,]
    numFilteredVariants <- nrow(genotypes)
    sumFilteredVariants <- rowSums(genotype)s
    
    totalPossiblePairs <- choose(numSamples,2)
    totalPairs <- choose(sumFilteredVariants,2)
    weights <- totalPossiblePairs/totalPairs
    s.i.j.numerator <- t(genotypes*weights)%*%genotypes
    s.i.j.denominator <- numFilteredVariants
    s.i.j <- s.i.j.numerator/s.i.j.denominator
    
    print(mean(s.i.j[row(s.i.j)!=col(s.i.j)]))
    print(median(s.i.j[row(s.i.j)!=col(s.i.j)]))
    topValues <- sort(s.i.j[row(s.i.j)>col(s.i.j)], decreasing=T)[1:40]
    print(topValues)
    print(which(s.i.j==topValues[1], arr.ind=T))
    print(which(s.i.j==topValues[2], arr.ind=T))
    print(which(s.i.j==topValues[3], arr.ind=T))
    qplot(s.i.j[row(s.i.j)!=col(s.i.j)], binwidth=.02) + ggtitle(paste0("Distribution of s, population: ",subpop)) + xlab("s")
}

calculateSMatrix("CEU", numberOfLines=30000, minVariants=30)
