library(ggplot2)


genotypeFile <- "./data/1000GP_Phase3_chr10.hap.gz"
minVariants <- 10
args<-commandArgs(TRUE)
if(length(args)!=0){
    genotypeFile <- args[1]
    numberOfLines <- as.numeric(args[2])
    minVariants <- as.numeric(args[3])
}

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
# source('~/1000GP/all_samples_HCL.R')
# calculateSMatrix(c("CEU"), numberOfLines=695, minVariants=20)
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

res <- foreach(pop_i=unique(pop),.packages=c("ggplot2")) %dopar% {
    calculateSMatrix(pop_i, filename=genotypeFile, numberOfLines=numberOfLines, minVariants=minVariants)
}

res <- calculateSMatrix(unique(pop), filename=genotypeFile, numberOfLines=numberOfLines, minVariants=minVariants)
res <- calculateSMatrix("CEU", filename=genotypeFile, numberOfLines=numberOfLines, minVariants=minVariants)

res <- readRDS('~/1000GP/plots/s_distributions/plotdata/allSamples_sij_80695.rds')
allSamplesGSM <- res[['s_i_j']]
varcovMat <- res[['varcovMat']]

pdf('~/1000GP/plots/s_matrix.pdf', width=9, height=9)
plotHeatmap(allSamplesGSM, title="s GSM")
dev.off()

pdf('~/1000GP/plots/varcov_matrix.pdf', width=9, height=9)
plotHeatmap(varcovMat,title="varcov GSM")
dev.off()

print(Sys.time()-strt)
if(!is.na(num_cores)){
    stopCluster(cl)
}