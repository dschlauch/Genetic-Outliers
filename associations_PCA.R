# after sourcing files from masterAssociations.R
library(ggplot2)
library(reshape2)
library(data.table)
require(gridExtra)
library(ggrepel)

getSampleSubset <- function(subpops){
  if("all"%in%subpops){
    sampleSubset <- rep(T,length(pop))
  } else {
    sampleSubset <- pop %in% subpops
  }
}
generateggplotDT <- function(jaccard.correction, varcov.correction, sampleSubset, components=c(1,2)){
  components.df <- data.frame(samples=rep(sampleIDs[sampleSubset],2))
  components.df$Vector1 <- c(jaccard.correction[,components[1]],varcov.correction[,components[1]])
  components.df$Vector2 <- c(jaccard.correction[,components[2]],varcov.correction[,components[2]])
#   components.df$VarcovVector1 <- 
#   components.df$VarcovVector2 <- varcov.correction[,components[2]]
  components.df$pop <- rep(pop[sampleSubset],2)
  components.df$group <- rep(group[sampleSubset],2)
  components.df$method <- c(rep("Ours",sum(sampleSubset)),rep("varcov",sum(sampleSubset)))
  data.table(components.df)
}
plotPCAsSideBySide <-  function(components.dt){
    plot_labeller <- function(variable,value){
        list(Ours="Our Method",varcov="Variance-Covariance")[value]
    }
  ggplot(components.dt, aes(x = Vector1, y = Vector2)) + 
      ggtitle("Principal Component Plots") + ylab("Component 2")+ xlab("Component 1") + guides(colour=FALSE) + theme_bw() + facet_grid(.~method, labeller=plot_labeller) +
      geom_point(aes(color=factor(pop)), size=5, alpha=.8) + geom_text_repel(data=subset(components.dt, Vector1>.1 | abs(Vector2)>.1 ),aes( label=samples)) 
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

sampleSubset <- T
sampleSubset <- getSampleSubset("all")
sampleSubset <- getSampleSubset(c("TSI","IBS"))
sampleSubset <- getSampleSubset(c("ITU","STU"))
sampleSubset <- getSampleSubset(c("CEU","YRI"))
sampleSubset <- getSampleSubset(c("CHB","CHS"))
sampleSubset <- getSampleSubset(c("GBR","FIN"))
sampleSubset <- getSampleSubset(c("GIH","ITU"))
sampleSubset <- getSampleSubset(c("ESN","LWK"))
sampleSubset <- getSampleSubset(c("PJL","BEB"))
sampleSubset <- getSampleSubset(c("CDX","CHB"))
sampleSubset <- getSampleSubset(c("CDX","CHS"))
sampleSubset <- getSampleSubset(c("PUR","ACB"))
sampleSubset <- getSampleSubset(c("KHV","CHS"))
sampleSubset <- getSampleSubset(c("ITU","BEB"))

jaccard.correction <- eigen(jaccardMatrix[sampleSubset,sampleSubset])$vectors
varcov.correction <- eigen(varcovMatrix[sampleSubset,sampleSubset])$vectors

jaccard.correction <- eigen(jaccardMatrix)$vectors
varcov.correction <- eigen(varcovMatrix)$vectors
components.dt <- generateggplotDT(jaccard.correction, varcov.correction, sampleSubset, c(2,3))
plotPCAsSideBySide(components.dt)
withinVsBetween(components.df)

