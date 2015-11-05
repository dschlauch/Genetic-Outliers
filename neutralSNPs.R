library(ggplot2)
library(foreach)
library(doParallel)


runAlleleFreqSim <- function(populationSize, numGen, numSim, alleleFreq){
  # Initiate cluster
  cores <- 4
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  #start time
  strt  <- Sys.time()
  res <- foreach(i=1:cores,.packages=c()) %dopar% {
    generation <- function(alleles){
      fatherHap <- rbinom(length(alleles), 1, sample(alleles[1:(populationSize/2)], populationSize, replace=T)/2)
      motherHap <- rbinom(length(alleles), 1, sample(alleles[(populationSize/2+1):populationSize], populationSize, replace=T)/2)
      fatherHap + motherHap
    }
    res <- replicate(numSim/cores,{
      alleles <- rbinom(populationSize, 2, alleleFreq)
      for(i in 1:numGen){
        alleles <- generation(alleles)
      }
      mean(alleles)/2
    })
  }
  print(Sys.time()-strt)
  stopCluster(cl)
  unlist(res)
}

populationSize <- 10000
numGen <- 1000
numSim <- 1000
alleleFreq = .1


res <- runAlleleFreqSim(populationSize, numGen, numSim, alleleFreq)

ggplot(data.frame(res), aes(res))+ geom_histogram(binwidth = .01, color="red", fill="blue") +
  ggtitle(paste0("Allele frequency for population of ", populationSize, " after ", numGen, " generations (", numSim," simulations)")) + 
  geom_vline(xintercept = alleleFreq) +
  annotate("text", x=alleleFreq-.001, y=20, label="Starting Allele Frequency", color="red",angle = 90) +
  xlab("Final allele frequency")

