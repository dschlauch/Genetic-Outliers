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
generateggplotDT <- function(jaccard.correction, varcov.correction, sampleNames, components=c(1,2)){
  components.df <- data.frame(samples=rep(sampleNames,2))
  components.df$Vector1 <- c(jaccard.correction[,components[1]],varcov.correction[,components[1]])
  components.df$Vector2 <- c(jaccard.correction[,components[2]],varcov.correction[,components[2]])
#   components.df$VarcovVector1 <- 
#   components.df$VarcovVector2 <- varcov.correction[,components[2]]
  components.df$pop <- rep(pop[sampleIDs%in%sampleNames],2)
  components.df$group <- rep(group[sampleIDs%in%sampleNames],2)
  components.df$method <- c(rep("Ours",length(sampleNames)),rep("varcov",length(sampleNames)))
  data.table(components.df)
}

plotPCAsSideBySide <-  function(components.dt, super=T, labelOutliers=5){
    if (super){
        colorVar <- "group"
        legendTitle <- "Super\nPopulation"
    } else {
        colorVar <- "pop"
        legendTitle <- "Population"
    }
    plot_labeller <- function(variable,value){
        list(Ours="Our Method",varcov="Variance-Covariance")[value]
    }
    components.dt$rank <-c(rank(-(components.dt[method=="Ours"]$Vector1)^2 - (components.dt[method=="Ours"]$Vector2)^2),
                           rank(-(components.dt[method=="varcov"]$Vector1)^2 - (components.dt[method=="varcov"]$Vector2)^2))
    ggplot(components.dt, aes(x = Vector1, y = Vector2)) + 
        ggtitle("Principal Component Plots") + ylab("Component 2")+ xlab("Component 1") + theme_bw() + facet_grid(.~method, labeller=plot_labeller) +
        geom_point(aes_string(color=colorVar), size=3, alpha=.7) + 
        guides(colour=guide_legend(title=legendTitle)) + 
        geom_text_repel(data=subset(components.dt, rank<=labelOutliers ),aes( label=samples)) +
        scale_colour_brewer(palette="Set1")
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


results <- readRDS("~/1000GP/plots/s_distributions/plotdata/STU_ITU_data.rds")
results <- readRDS("~/1000GP/plots/s_distributions/plotdata/CEU_YRI_data.rds")
results <- readRDS("~/1000GP/plots/s_distributions/plotdata/IBS_TSI_data.rds")

jaccardMatrix <- results$s_matrix_dip
varcovMatrix <- results$varcovMat

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

jaccard.correction <- eigen(apply(jaccardMatrix,2,scale),symmetric=T)$vectors
varcov.correction <- eigen(apply(varcovMatrix,2,scale),symmetric=T)$vectors
varcov.correction[,3] <- -varcov.correction[,3] 

components.dt <- generateggplotDT(jaccard.correction, varcov.correction, rownames(jaccardMatrix), c(1,2))

pdf("./plots/PCA_all.pdf", height=5, width=10)
plotPCAsSideBySide(components.dt, super=F, labelOutliers=0)
dev.off()

withinVsBetween(components.df)

