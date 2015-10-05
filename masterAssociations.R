source('~/1000GP/calculatePCA.R')
source('~/1000GP/associations-copy.R')

chr <- 7
#subpops <- c("TSI","IBS")
numRows <- 1000
varianceEstimate <- "adjVar" 
subpops <- c("all","all")
args<-commandArgs(TRUE)
correctionMethods <- c("jaccard")
if(length(args)!=0){
    chr <- as.numeric(args[1])
    subpops <- c(args[2], args[3])
    numRows <- as.numeric(args[4])
    varianceEstimate <- args[5] # "adjVar" or "ATT"
    correctionMethods <- args[-1:-5] # "uncorrected","varcov","jaccard", and "subpop"
}
sapply(correctionMethods, function(method){
    print(args)
    runAndPlot(chr=chr,correctMethod=method, ATT=varianceEstimate, subpops=subpops, numRows=numRows)
    gc()
})
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

