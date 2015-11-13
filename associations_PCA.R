# after sourcing files from masterAssociations.R
library(ggplot2)
library(reshape2)
require(gridExtra)

getSampleSubset <- function(subpops){
  if("all"%in%subpops){
    sampleSubset <- rep(T,length(pop))
  } else {
    sampleSubset <- pop %in% subpops
  }
}
generateggplotDF <- function(jaccard.correction, varcov.correction, sampleSubset, components=c(1,2)){
  components.df <- data.frame(row.names=sampleIDs[sampleSubset])
  components.df$JaccardVector1 <- jaccard.correction[,components[1]]
  components.df$JaccardVector2 <- jaccard.correction[,components[2]]
  components.df$VarcovVector1 <- varcov.correction[,components[1]]
  components.df$VarcovVector2 <- varcov.correction[,components[2]]
  components.df$pop <- pop[sampleSubset]
  components.df$group <- group[sampleSubset]
  components.df
}
plotPCAsSideBySide <-  function(components.df){
  jaccardPlot <- ggplot(components.df, aes(x = JaccardVector1, y = JaccardVector2)) + ggtitle("Weighted Jaccard") + geom_point(aes(color=factor(pop))) + ylab("Component 2")+ xlab("Component 1") + guides(colour=FALSE)
  varcovPlot <-  ggplot(components.df, aes(x = VarcovVector1, y = VarcovVector2))   + ggtitle("Variance") + geom_point(aes(color=factor(pop))) + ylab("Component 2")+ xlab("Component 1")
  grid.arrange(jaccardPlot, varcovPlot, ncol=2,top = "Principal Component Plots")
}
withinVsBetween <-  function(components.df){
  ind <- t(combn(nrow(components.df),2))
  comparison <- t(combn(pop[sampleSubset],2))
  inGroup <- comparison[,1]==comparison[,2]
  jaccardDist <- apply(ind, 1, function(x) sqrt((components.df[x[1], 1]-components.df[x[2], 1])^2 + (components.df[x[1], 2]-components.df[x[2], 2])^2))
  varianceDist <- apply(ind, 1, function(x) sqrt((components.df[x[1], 3]-components.df[x[2], 3])^2 + (components.df[x[1], 4]-components.df[x[2], 4])^2))
  jaccardRatio <- mean(jaccardDist[inGroup])/mean(jaccardDist[!inGroup])
  varianceRatio <-mean(varianceDist[inGroup])/mean(varianceDist[!inGroup])
  c("jaccard"=jaccardRatio, "variance"=varianceRatio)  
}
####################################
####################################
####################################


sampleSubset <- getSampleSubset("all")
sampleSubset <- getSampleSubset(c("TSI","IBS"))
sampleSubset <- getSampleSubset(c("STU","ITU"))
sampleSubset <- getSampleSubset(c("CEU","YRI"))
sampleSubset <- getSampleSubset(c("CHB","CHS"))
sampleSubset <- getSampleSubset(c("GBR","FIN"))
sampleSubset <- getSampleSubset(c("GIH","ITU"))
sampleSubset <- getSampleSubset(c("ESN","LWK"))

jaccard.correction <- eigen(jaccardMatrix[sampleSubset,sampleSubset])$vectors
varcov.correction <- eigen(varcovMatrix[sampleSubset,sampleSubset])$vectors
components.df <- generateggplotDF(jaccard.correction, varcov.correction, sampleSubset, c(3,2))
plotPCAsSideBySide(components.df)
withinVsBetween(components.df)

