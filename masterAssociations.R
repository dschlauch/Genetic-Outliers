

chr <- "filtered100"
numRows <- 100000
varianceEstimate <- "adjVar" 
subpops <- c("all","all")
#subpops <- c("TSI","IBS")
args<-commandArgs(TRUE)
correctionMethods <- c("subpop","varcov","jaccard", "uncorrected")
phenotypeVariable <- "continuous" #"binary"
if(length(args)!=0){
    chr <- as.numeric(args[1])
    subpops <- c(args[2], args[3])
    numRows <- as.numeric(args[4])
    varianceEstimate <- args[5] # "adjVar" or "ATT"
    phenotypeVariable <- args[6]
    correctionMethods <- args[-1:-6] # "uncorrected","varcov","jaccard", and "subpop"
}

source('~/1000GP/calculatePCA.R')
source('~/1000GP/associations-copy.R')

library(foreach)
library(doParallel)

cl <- makeCluster(4)
registerDoParallel(cl)
strt  <- Sys.time()
foreach(i=1:length(correctionMethods),.packages=c("ggplot2","MASS")) %dopar% {
    print(args)
    print(correctionMethods[i])
    runAndPlot(chr=chr,correctMethod=correctionMethods[i], ATT=varianceEstimate, subpops=subpops, numEigenVectors=2, numRows=numRows, phenotypeVariable)
    gc()
}
print(Sys.time()-strt)
stopCluster(cl)

# 
# sapply(correctionMethods, function(method){
#     print(args)
#     runAndPlot(chr=chr,correctMethod=method, ATT=varianceEstimate, subpops=subpops, numEigenVectors=2, numRows=numRows, phenotypeVariable)
#     gc()
# })
# 
# runAndPlot(chr=chr,correctMethod="uncorrected", ATT="adjVar", subpops=subpops, numRows=numRows)
# gc()
# runAndPlot(chr=chr,correctMethod="uncorrected", ATT="ATT", subpops=subpops, numRows=numRows)
# gc()
# runAndPlot(chr=chr,correctMethod="varcov", ATT="adjVar", subpops=subpops, numRows=numRows)
# gc()
# runAndPlot(chr=chr,correctMethod="varcov", ATT="ATT", subpops=subpops, numRows=numRows)
# gc()
# runAndPlot(chr=chr,correctMethod="jaccard", ATT="adjVar", subpops=subpops, numEigenVectors=1, numRows=numRows)
# gc()
# runAndPlot(chr=chr,correctMethod="jaccard", ATT="ATT", subpops=subpops, numRows=numRows)
# gc()
# runAndPlot(chr=chr,correctMethod="subpop", ATT="adjVar", subpops=subpops, numRows=numRows)
# gc()
# runAndPlot(chr=chr,correctMethod="subpop", ATT="ATT", subpops=subpops, numRows=numRows)
# gc()
# 
# 
# runAndPlot(chr=chr,correctMethod="subpop", ATT="ATT", subpops="all", numRows=numRows)
# gc()
# runAndPlot(chr=chr,correctMethod="subpop", ATT="adjVar", subpops="all", numRows=numRows)
# gc()
# runAndPlot(chr=chr,correctMethod="jaccard", ATT="adjVar", subpops="all", numRows=numRows)
# gc()
# runAndPlot(chr=chr,correctMethod="varcov", ATT="adjVar", subpops="all", numRows=numRows)
# gc() 

