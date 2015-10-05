setwd('~/1000GP/')
library(data.table)

minAC<-20
maxAC<-500
outputDir <- paste('./output_',minAC,'_',maxAC, sep="")

V <- Reduce('+', lapply(1:22, function(chr){
    as.data.frame(fread(paste('gunzip -c ',outputDir,'/chr',chr,'_V.csv.gz', sep=""), sep=","))
}))
write.csv(V, file=paste(outputDir,"/all_V.csv",sep=""))
rm(list = ls())

minAC<-20
maxAC<-500
outputDir <- paste('./output_',minAC,'_',maxAC, sep="")
U <- Reduce('+', lapply(1:22, function(chr){
    as.data.frame(fread(paste('gunzip -c ',outputDir,'/chr',chr,'_U.csv.gz', sep=""), sep=","))
}))
write.csv(U, file=paste(outputDir,"/all_U.csv",sep=""))
rm(list = ls())