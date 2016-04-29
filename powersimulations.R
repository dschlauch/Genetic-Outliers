# 4/29/16
# For running power simulations
#

# Read in command line args -----------------------------------------------

args<-commandArgs(TRUE)
if(length(args)!=0){
    numSimulatedSamples <- as.numeric(args[1])
    nVariants <- as.numeric(args[2])
    numSimulations <- as.numeric(args[3])
    minVariants <- as.numeric(args[4])
    outputDir <- args[5]
    numCores <- as.numeric(args[6])
    dir.create(paste0("~/1000GP/plots/s_distributions/",outputDir))
    dir.create(paste0("~/1000GP/plots/s_distributions/",outputDir,"/plotdata"))
}


library(ggplot2)
library(reshape2)
library(data.table)
library(gtools)
library(grid)
source('~/1000GP/s_matrix_functions.R')

kinshipCoefs <- seq(0,.0625,.0025)
system.time(sapply(kinshipCoefs, function(cok){
    system.time(results <- homogeneousSimulations(numSimulatedSamples=20, nVariants=10000, cok=cok, numSimulations=100, minVariants=5, numCores, outputDir=outputDir))
    saveRDS(results, paste0("./plots/s_distributions/",outputDir,"/plotdata/Simulated_",cok,"data.rds"))
    print(cok)
}))

combinedSimTable <- data.table(do.call(rbind,lapply(kinshipCoefs, function(cok){
    simResults <- readRDS(paste0("./plots/s_distributions/",outputDir,"/plotdata/Simulated_",cok,"data.rds"))
    simTable <- getPopResults(simResults)
    simTable$cok <- cok
    simTable
})))

powerCurve <- combinedSimTable[, mean(crypticSig=="NO"), by=cok ]
plot(powerCurve)
lines(powerCurve)
# 
# pValues <- sapply(simResults, function(res){
#     #     s_vector <- sort(res$s_matrix_hap[row(res$s_matrix_hap)>col(res$s_matrix_hap)], decreasing=T)
#     s_vector <- res$s_matrix_hap[row(res$s_matrix_hap)>col(res$s_matrix_hap)]
#     crypticPValue <- 1-pnorm((s_vector-1)/sqrt(res$var_s_hap))
#     crypticPValue
# })