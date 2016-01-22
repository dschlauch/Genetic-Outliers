library(readr)
library(ggplot2)


#### Simulated data
# simn <- 3000
# numVariants <- 100000
# sim.w.numerator <- 2*(2*simn-1)
# simMAF <- runif(numVariants,0,.5)
# genotypes <- sapply(1:simn,function(x){
#     rbinom(numVariants,1,simMAF)
# })

# Real data
sample <- read.table("~/1000GP/data/1000GP_Phase3.sample", sep=" ", header=T)
sampleIDs <- as.character(sample[,1])
pop <- as.character(sample[,2])
group <- as.character(sample[,3])
sex <- as.character(sample[,4])
hap.pop <- rep(pop,each=2)
hap.sampleIDs <- rep(as.character(sample[,1]),each=2)

source('~/1000GP/s_matrix_functions.R')
# calculateSMatrix(c("CEU"), numberOfLines=10695, minVariants=20)
# calculateSMatrix(c("CEU","CHB"), numberOfLines=10695, minVariants=20)
# sapply(unique(pop),calculateSMatrix, numberOfLines=10695, minVariants=20)


library(foreach)
library(doParallel)

num_cores <- 4 #detectCores() - 4

# Initiate cluster
if(!is.na(num_cores)){
    cl <- makeCluster(num_cores)
    registerDoParallel(cl)
}

res <- foreach(i=unique(pop),.packages=c("ggplot2")) %dopar% {
    calculateSMatrix(i, numberOfLines=80695, minVariants=10)
}

print(Sys.time()-strt)
if(!is.na(num_cores)){
    stopCluster(cl)
}
