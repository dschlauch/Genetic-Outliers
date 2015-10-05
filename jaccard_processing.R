library(data.table)
library(gplots)

args<-commandArgs(TRUE)
minAC <- args[1]
maxAC <- args[2]

outputDir <- paste('./output_',minAC,'_',maxAC, sep="")

sample <- read.table("./data/1000GP_Phase3.sample", sep=" ", header=T)

V <- Reduce('+', lapply(1:22, function(chr){
    as.data.frame(fread(paste('gunzip -c ',outputDir,'/chr',chr,'_V.csv.gz', sep=""), sep=","))
}))
U <- Reduce('+', lapply(1:22, function(chr){
    as.data.frame(fread(paste('gunzip -c ',outputDir,'/chr',chr,'_U.csv.gz', sep=""), sep=","))
}))

jaccardMat <- V/U

sampleIDs <- unlist(lapply(as.character(sample[,1]),function(x){c(x,x)}))
pop <- unlist(lapply(as.character(sample[,2]),function(x){c(x,x)}))
group <- unlist(lapply(as.character(sample[,3]),function(x){c(x,x)}))
sex <- unlist(lapply(as.character(sample[,4]),function(x){c(x,x)}))

numSamps <- nrow(sample)

combinedJaccard <- (jaccardMat[seq(1,nrow(jaccardMat),2),seq(1,nrow(jaccardMat),2)] + 
                    jaccardMat[seq(2,nrow(jaccardMat),2),seq(1,nrow(jaccardMat),2)] +
                    jaccardMat[seq(1,nrow(jaccardMat),2),seq(2,nrow(jaccardMat),2)] + 
                    jaccardMat[seq(2,nrow(jaccardMat),2),seq(2,nrow(jaccardMat),2)])/4

rownames(combinedJaccard) <- sample[,1] # Sample ID

#subset <- seq(1,nrow(combinedJaccard),12)
subset <- 1:nrow(combinedJaccard)

colnames(combinedJaccard) <- sample[,2] # Pop name
png(paste(outputDir,"/dendro_pop.png",sep=""), width=32000)
plot(hclust(dist(combinedJaccard[subset,subset]),"ave"),labels = colnames(combinedJaccard)[subset], xlab="", cex=.7, main="1000GP : Population Dendrogram")
dev.off()

colnames(combinedJaccard) <- sample[,3] # Super pop name
png(paste(outputDir,"/dendro_super.png",sep=""), width=32000)
plot(hclust(dist(combinedJaccard[subset,subset]),"ave"),labels = colnames(combinedJaccard)[subset], xlab="", cex=.7, main="1000GP : Population Dendrogram")
dev.off()

diag(combinedJaccard) <- 0
png(paste(outputDir,"/heatmap_group.png",sep=""), width=12000, height=12000)
heatmap.2(10*as.matrix(combinedJaccard),dendrogram="none", trace="none", density.info="none", key=F,col=bluered)
dev.off()

colnames(combinedJaccard) <- sample[,2] # Pop name
write.csv(combinedJaccard, file=paste(outputDir,"/combined_jaccard.csv",sep=""))
write.csv(V, file=paste(outputDir,"/all_V.csv",sep=""))
write.csv(U, file=paste(outputDir,"/all_U.csv",sep=""))
