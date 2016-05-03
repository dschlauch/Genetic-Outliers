# 4/29/16
# For running power simulations
#
library(ggplot2)
library(reshape2)
library(data.table)
library(gtools)
library(grid)
source('~/1000GP/s_matrix_functions.R')

numSimulatedSamples <- 20
nVariants <- 10000
numSimulations <- 10
minVariants <- 5
outputDir <- '.'
numCores <- 4
kinshipStart <- 0
kinshipEnd <- .034
kinshipSeq <- .001
kinshipCoefs <- seq(kinshipStart,kinshipEnd,kinshipSeq)


# Read in command line args -----------------------------------------------

args<-commandArgs(TRUE)
if(length(args)!=0){
    numSimulatedSamples <- as.numeric(args[1])
    nVariants <- as.numeric(args[2])
    numSimulations <- as.numeric(args[3])
    minVariants <- as.numeric(args[4])
    outputDir <- args[5]
    numCores <- as.numeric(args[6])
    kinshipStart <- as.numeric(args[7])
    kinshipEnd <- as.numeric(args[8])
    kinshipSeq <- as.numeric(args[9])
    kinshipCoefs <- seq(kinshipStart,kinshipEnd,kinshipSeq)
    dir.create(paste0("~/1000GP/plots/s_distributions/",outputDir))
    dir.create(paste0("~/1000GP/plots/s_distributions/",outputDir,"/plotdata"))
}




# system.time(sapply(kinshipCoefs, function(cok){
#     system.time(results <- homogeneousSimulations(numSimulatedSamples=numSimulatedSamples, nVariants=nVariants, cok=cok, numSimulations=numSimulations, minVariants=minVariants, numCores=numCores, outputDir=outputDir))
#     saveRDS(results, paste0("./plots/s_distributions/",outputDir,"/plotdata/Simulated_",cok,"data.rds"))
#     print(cok)
# }))


powerSimulationResult <- lapply(kinshipCoefs, function(cok){
    system.time(simResults <- homogeneousSimulations(numSimulatedSamples=numSimulatedSamples, nVariants=nVariants, cok=cok, numSimulations=numSimulations, minVariants=minVariants, numCores=numCores, outputDir=outputDir))
#     saveRDS(results, paste0("./plots/s_distributions/",outputDir,"/plotdata/Simulated_",cok,"data.rds"))
    print(cok)
#     simResults <- readRDS(paste0("./plots/s_distributions/",outputDir,"/plotdata/Simulated_",cok,"data.rds"))
    nSamples <- ncol(simResults[[1]]$s_matrix_hap)
    varS <- simResults[[1]]$var_s_hap
    expectedS <- (1-cok) + cok*(.1*choose(nSamples,2)/(choose(.1*(nSamples-2)+2,2))) 
    simTable <- getPopResults(simResults)
    simTable$cok <- cok
    resList <- list(combinedSimTable=simTable, expectedS=expectedS, nSamples=nSamples, varS=varS)
    saveRDS(resList, paste0("./plots/s_distributions/",outputDir,"/plotdata/Simulated_",cok,"_results.rds"))
    resList
})
saveRDS(powerSimulationResult, paste0("./plots/s_distributions/",outputDir,"/plotdata/SimulationResults.rds"))
powerSimulationResult <- readRDS(paste0("./plots/s_distributions/",outputDir,"/plotdata/SimulationResults.rds"))

combinedSimTable <- do.call(rbind, lapply(powerSimulationResult, '[[',1))
expectedSvector <- unlist(lapply(powerSimulationResult, '[[',2))
nSamples <- unlist(lapply(powerSimulationResult, '[[',3))[1]
varSvector <- unlist(lapply(powerSimulationResult, '[[',4))
fwalpha <- .05
bonferroniAlpha <- fwalpha/choose(nSamples,2)
cutoff <- qnorm(1-bonferroniAlpha, mean=1,sd=sqrt(varSvector))
simPower <- fwalpha + (1-fwalpha)*pnorm(expectedSvector-cutoff, sd=sqrt(varSvector))

powerCurve <- combinedSimTable[, mean(crypticSig=="NO"), by=cok ]
powerCurve$typeIIError <- 1-simPower
pdf("./plots/powerCurve.pdf", width=10)
ggplot(powerCurve) + 
    geom_line(aes(x=cok,y=typeIIError, col="Expected")) +
    geom_point(aes(x=cok,y=V1, col="Simulated"), size=4) +
    scale_colour_manual(values=c("blue","red")) +
    geom_hline(yintercept=.95, linetype="dotted") + 
    geom_hline(yintercept=0) +
    geom_vline(xintercept=0) +  
    guides(colour = guide_legend(override.aes = list(shape=c(NA,16),linetype=c(1,0)))) + 
    ggtitle(expression(paste("Type II error vs ", phi, " ,",alpha[fw]==.05))) + xlab(expression(phi)) + ylab("Type II Error") +
    theme_bw() +
    theme(plot.title = element_text(size=40), axis.title.x = element_text(size = 30), axis.title.y = element_text(size = 30), 
          axis.text.x=element_text(size=20), axis.text.y=element_text(size=20),
          legend.title=element_blank(),legend.position=c(0.04, .06),legend.justification=c(0,0), 
          legend.key.size = unit(1.5, "cm"), legend.text=element_text(size=20),
          legend.background = element_rect(fill=alpha('white', 0.5)),legend.key = element_rect(colour = NA))
dev.off()
# 
# pValues <- sapply(simResults, function(res){
#     #     s_vector <- sort(res$s_matrix_hap[row(res$s_matrix_hap)>col(res$s_matrix_hap)], decreasing=T)
#     s_vector <- res$s_matrix_hap[row(res$s_matrix_hap)>col(res$s_matrix_hap)]
#     crypticPValue <- 1-pnorm((s_vector-1)/sqrt(res$var_s_hap))
#     crypticPValue
# })