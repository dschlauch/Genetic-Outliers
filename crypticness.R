library(ggplot2)

genotypeFile <- "./data/1000GP_Phase3_chr10.hap.gz"
genotypeFile <- "./data/combinedFiltered1000.gz"
numberOfLines <- 1000
minVariants <- 10
numCores <- 4
args<-commandArgs(TRUE)
if(length(args)!=0){
    genotypeFile <- args[1]
    numberOfLines <- as.numeric(args[2])
    minVariants <- as.numeric(args[3])
    numCores <- as.numeric(args[4])
    outputDir <- args[5]
    dir.create(paste0("~/1000GP/plots/s_distributions/",outputDir))
    dir.create(paste0("~/1000GP/plots/s_distributions/",outputDir,"/plotdata"))
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


library(foreach)
library(doParallel)

# Initiate cluster
if(!is.na(numCores)){
    cl <- makeCluster(numCores)
    registerDoParallel(cl)
}

strt <- Sys.time()
# Calculate all the s matrices and save
res <- foreach(pop_i=unique(pop),.packages=c("ggplot2")) %dopar% {
    result <- calculateSMatrix(pop_i, filename=genotypeFile, numberOfLines=numberOfLines, minVariants=minVariants)
    saveRDS(result, paste0("./plots/s_distributions/",outputDir,"/plotdata/",pop_i, "_data.rds"))
}

print(Sys.time()-strt)
if(!is.na(numCores)){
    stopCluster(cl)
}
# Read in all the results. plot histograms.  Calculated structure p.value
sapply(unique(pop),function(pop_i){
    result <- readRDS(paste0("./plots/s_distributions/",outputDir,"/plotdata/",pop_i, "_data.rds"))
    plotFromGSM(pop_i, result$s_matrix_dip, result$var_s_dip, result$weightsMean, sampleIDs[pop%in%pop_i], "diploid", outputDir=outputDir)
    s_vector <- result$s_matrix_dip[row(result$s_matrix_dip)>col(result$s_matrix_dip)]
    topValuesKinship <- (sort(s_vector, decreasing=T)-1)/(result$weightsMean-1)
    btest <- binom.test(sum(s_vector>mean(s_vector)), length(s_vector), alternative="less")
    btest$p.value
})

# # Same as function above, merge when ready
# p.values <- sapply(unique(pop),function(pop_i){
#     result <- readRDS(paste0("./plots/s_distributions/",outputDir,"/plotdata/",pop_i, "_sij.rds"))
#     s_vector <- result$s_matrix_dip[row(result$s_matrix_dip)>col(result$s_matrix_dip)]
#     
# })
# res <- calculateSMatrix(unique(pop), filename=genotypeFile, numberOfLines=numberOfLines, minVariants=minVariants)
# res <- calculateSMatrix("CEU", filename=genotypeFile, numberOfLines=numberOfLines, minVariants=minVariants)

# res <- readRDS('~/1000GP/plots/s_distributions/",outputDir,"plotdata/allSamples_sij_80695.rds')
# allSamplesGSM <- res[['s_i_j']]
# varcovMat <- res[['varcovMat']]
# 
# pdf('~/1000GP/plots/s_matrix.pdf', width=9, height=9)
# plotHeatmap(allSamplesGSM, title="s GSM")
# dev.off()
# 
# pdf('~/1000GP/plots/varcov_matrix.pdf', width=9, height=9)
# plotHeatmap(varcovMat,title="varcov GSM")
# dev.off()
