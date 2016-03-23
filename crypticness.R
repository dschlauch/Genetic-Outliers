library(ggplot2)
library(reshape2)
library(data.table)
library(gtools)

source('~/1000GP/read1000GPsupportFiles.R')
source('~/1000GP/s_matrix_functions.R')

# Calculate all the s matrices and save
system.time(results <- calculateSMatrix(pop_i, filename=genotypeFile, numberOfLines=numberOfLines, minVariants=minVariants, qcFilter=qcFilter, ldPrune))
names(results) <- unique(pop)
saveRDS(results, paste0("./plots/s_distributions/",outputDir,"/plotdata/all_data.rds"))

source('~/1000GP/process1000GPResults.R')
