library(gplots)
colorCodes <- c(AFR="red", AMR="green", EUR="blue", SAS="yellow", EAS="black")
colorCodesPop <- c("dodgerblue2","#E31A1C", # red
                          "green4",
                          "#6A3D9A", # purple
                          "#FF7F00", # orange
                          "black","gold1",
                          "skyblue2","#FB9A99", # lt pink
                          "palegreen2",
                          "#CAB2D6", # lt purple
                          "#FDBF6F", # lt orange
                          "gray70", "khaki2",
                          "maroon","orchid1","deeppink1","blue1","steelblue4",
                          "darkturquoise","green1","yellow4","yellow3",
                          "darkorange4","brown","black")
names(colorCodesPop)<-unique(pop)
labelCol <- function(x) {
  if (is.leaf(x)) {
    ## fetch label
    label <- attr(x, "label")
    #     code <- substr(label, 1, 1)
    ## use the following line to reset the label to one letter code
    # attr(x, "label") <- code
    attr(x, "nodePar") <- list(lab.col=colorCodes[label])
  }
  return(x)
}

# Load in data
# jaccardMatrix <- read.csv('~/1000GP/output_0_20/combined_jaccard.csv', row.names=1)
# 
# ### hierarchical clustering
# 
subset <- sample(2504,100)
subset <- -related
subset <- rep(T,2504)
jm <- as.matrix(jaccardMatrix)[subset,subset]
diag(jm)<-0
rownames(jm) <- group[subset]
colnames(jm) <- pop[subset]
hc <- hclust(as.dist(max(jm)-jm),method="average")
d <- dendrapply(as.dendrogram(hc, hang=max(jm)*.1), labelCol)
plot(d, main="Hierarchical Clustering across superpopulations")
# 
# jm[jm>.005]<-.005
# dev.off()
# heatmap.2(jm, Rowv=d, Colv=d, dendrogram="none", trace="none",
#           labRow="",labCol="",key=FALSE, RowSideColors=colorCodes[rownames(jm)], ColSideColors=colorCodesPop[colnames(jm)],
# #           lwid = c(1,10),lhei = c(.01,5), margins = c(5,10))
# )
# legend(.1,.5, names(colorCodes), inset=c(-0.01,0),  lty=1, lwd=10, col=colorCodes, cex = 0.75)
# legend("topright", names(colorCodesPop[1:13]),  lty=1, lwd=10, col=colorCodesPop[1:13], cex = 0.5, horiz = T,  inset=c(0,.05))
# legend("topright", names(colorCodesPop[14:26]),  lty=1, lwd=10, col=colorCodesPop[14:26], cex = 0.5, horiz = T,  inset=c(0,.10))

plotHeatmap <-  function(x, subset=NA, title="GSM"){
    if(is.na(subset)){
        subset <- rep(T,nrow(x))
    }
    simMat <- as.matrix(x)[subset,subset]
    diag(simMat)<-0
    rownames(simMat) <- group[subset]
    colnames(simMat) <- pop[subset]
    hc <- hclust(as.dist(1-simMat),method="average")
    d <- dendrapply(as.dendrogram(hc, hang=max(simMat)*.1), labelCol)
#     plot(d, main="Hierarchical Clustering across superpopulations", ylim=c(1-max(simMat)*1.25,1))
    
#     simMat[simMat>.005]<-.005
    simMat[simMat>quantile(simMat,.99)]<-quantile(simMat,.99)
    heatmap.2(simMat, Rowv=d, Colv=d, dendrogram="none", trace="none", main=title, col="bluered",
              labRow="",labCol="",key=FALSE, RowSideColors=colorCodes[rownames(simMat)], ColSideColors=colorCodesPop[colnames(simMat)],
              #           lwid = c(1,10),lhei = c(.01,5), margins = c(5,10))
    )
    legend("left", names(colorCodes), inset=c(0.1,0),  lty=1, lwd=10, col=colorCodes, cex = 0.75)
    legend("topright", names(colorCodesPop[1:13]),  lty=1, lwd=10, col=colorCodesPop[1:13], cex = 0.5, horiz = T,  inset=c(0,.05))
    legend("topright", names(colorCodesPop[14:26]),  lty=1, lwd=10, col=colorCodesPop[14:26], cex = 0.5, horiz = T,  inset=c(0,.10))
}


## Usage
## Commented so can source from crypticness.R
# plotHeatmap(jaccardMatrix,title="Rare jaccard GSM")
# plotHeatmap(varcovMatrix,title="Varcov GSM")
# # for removing possible related individuals
# related <- which(jm>.1,arr.ind=T)[,1]
