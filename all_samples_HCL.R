colorCodes <- c(AFR="red", AMR="green", EUR="blue", SAS="yellow", EAS="black")
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


### hierarchical clustering

subset <- sample(2504,500)
subset <- -related
subset <- rep(T,2504)
jm <- as.matrix(jaccardMatrix)[subset,subset]
diag(jm)<-0
rownames(jm) <- group[subset]
colnames(jm) <- pop[subset]
hc <- hclust(as.dist(1-jm),method="average")
d <- dendrapply(as.dendrogram(hc, hang=max(jm)*.1), labelCol)
plot(d, main="Hierarchical Clustering across superpopulations", ylim=c(1-max(jm)*1.25,1))

jm[jm>.005]<-.005
heatmap.2(1-jm, dendrogram="none", trace="none",col=bluered)
# for removing possible related individuals
related <- which(jm>.1,arr.ind=T)[,1]
