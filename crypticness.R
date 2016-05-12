library(ggplot2)
library(reshape2)
library(data.table)
library(gtools)
library(grid)


# Analysis parameters -----------------------------------------------------
genotypeFile <- "./data/combinedFiltered1000.gz"
numberOfLines <- 30000
minVariants <- 10
numCores <- 4
args<-commandArgs(TRUE)
outputDir <- '.'
# outputDir <- 'filtered40_LD10'
# outputDir <- 'filtered40_TGP2261_LD10'

# gazalFilter <-"TGP2261"
gazalFilter <- "NA"
ldPrune <- 10
if(length(args)!=0){
    genotypeFile <- args[1]
    numberOfLines <- as.numeric(args[2])
    minVariants <- as.numeric(args[3])
    numCores <- as.numeric(args[4])
    outputDir <- args[5]
    gazalFilter <- args[6]
    ldPrune <- as.numeric(args[7])
    dir.create(paste0("~/1000GP/plots/s_distributions/",outputDir))
    dir.create(paste0("~/1000GP/plots/s_distributions/",outputDir,"/plotdata"))
}

source('~/1000GP/read1000GPsupportFiles.R')
source('~/1000GP/s_matrix_functions.R')


# Run all for HCL ---------------------------------------------------------
system.time(results <- calculateSMatrix("All", filename=genotypeFile, numberOfLines=numberOfLines, minVariants=minVariants, qcFilter=qcFilter, ldPrune, outputDir))

# Run some for HCL ---------------------------------------------------------
system.time(results <- calculateSMatrix(c("STU","ITU"), filename=genotypeFile, numberOfLines=numberOfLines, minVariants=minVariants, qcFilter=qcFilter, ldPrune, outputDir))
system.time(results <- calculateSMatrix(c("IBS","TSI"), filename=genotypeFile, numberOfLines=numberOfLines, minVariants=minVariants, qcFilter=qcFilter, ldPrune, outputDir))


# Calculate all the s matrices and save
system.time(results <- calculateSMatrix("each", filename=genotypeFile, numberOfLines=numberOfLines, minVariants=minVariants, qcFilter=qcFilter, ldPrune, outputDir))
saveRDS(results, paste0("./plots/s_distributions/",outputDir,"/plotdata/all_data.rds"))

# Simulated results
kinshipCoefs <- seq(0,.0625,.0025)
sapply(kinshipCoefs, function(cok){
    system.time(simResults <- homogeneousSimulations(numSimulatedSamples=100, nVariants=100000, cok=NA, numSimulations=30, minVariants=5, outputDir=outputDir))
    saveRDS(results, paste0("./plots/s_distributions/",outputDir,"/plotdata/Simulated_",cok,"data.rds"))
    print(cok)
})
saveRDS(results, paste0("./plots/s_distributions/",outputDir,"/plotdata/Simulated_data.rds"))
saveRDS(results, paste0("./plots/s_distributions/",outputDir,"/plotdata/Simulated_data_cok0625.rds"))

source('~/1000GP/process1000GPResults.R')
