# after sourcing files from masterAssociations.R
library(ggplot2)
library(reshape2)
library(data.table)
library(gridExtra)
library(ggrepel)
library(rms)

getSampleSubset <- function(subpops){
  if("all"%in%subpops){
    sampleSubset <- rep(T,length(pop))
  } else {
    sampleSubset <- pop %in% subpops
  }
}
generateggplotDT <- function(jaccard.correction, varcov.correction, sampleNames, components=c(1,2,3)){
    components.df <- data.frame(samples=rep(sampleNames,2))
    components.df$Vector1 <- c(jaccard.correction[,components[1]],varcov.correction[,components[1]])
    components.df$Vector2 <- c(jaccard.correction[,components[2]],varcov.correction[,components[2]])
    components.df$Vector3 <- c(jaccard.correction[,components[3]],varcov.correction[,components[3]])
    #   components.df$VarcovVector1 <- 
    #   components.df$VarcovVector2 <- varcov.correction[,components[2]]
    components.df$pop <- rep(pop[sampleIDs%in%sampleNames],2)
    components.df$group <- rep(group[sampleIDs%in%sampleNames],2)
    components.df$method <- c(rep("STEGO",length(sampleNames)),rep("varcov",length(sampleNames)))
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
        list(STEGO="STEGO",varcov="Variance-Covariance")[value]
    }
    ranks <- c(rank(-(components.dt[method=="STEGO"]$Vector1)^2 - (components.dt[method=="STEGO"]$Vector2)^2),
                  rank(-(components.dt[method=="varcov"]$Vector1)^2 - (components.dt[method=="varcov"]$Vector2)^2))
    components.dt$rank <- sapply(1:(length(ranks)/2), function(i){
        min(ranks[i],ranks[i+length(ranks)/2])
    })
    
#     ranks <- c(rep(,2))
#     components.dt$rank <- c(rank(-(components.dt[method=="STEGO"]$Vector1)^2 - (components.dt[method=="STEGO"]$Vector2)^2),
#                            rank(-(components.dt[method=="varcov"]$Vector1)^2 - (components.dt[method=="varcov"]$Vector2)^2))
    
    ggplot(components.dt, aes(x = Vector1, y = Vector2)) + 
        ggtitle("Principal Component Plots") + ylab("Component 2")+ xlab("Component 1") + theme_bw() + facet_grid(.~method, labeller=plot_labeller) +
        geom_point(aes_string(color=colorVar), size=3, alpha=.7) + 
        guides(colour=guide_legend(title=legendTitle)) + 
        geom_text_repel(data=subset(components.dt, rank<=labelOutliers ),aes( label=samples)) +
        scale_colour_brewer(palette="Set1")
}


withinVsBetween <-  function(components.dt){
    totalmseSTEGO <- mean(components.dt[method=="STEGO"]$Vector1^2+components.dt[method=="STEGO"]$Vector2^2+components.dt[method=="STEGO"]$Vector3^2)
    Vector1DifStego <- components.dt[method=="STEGO",Vector1-mean(Vector1),by=pop]$V1
    Vector2DifStego <- components.dt[method=="STEGO",Vector2-mean(Vector2),by=pop]$V1
    Vector3DifStego <- components.dt[method=="STEGO",Vector3-mean(Vector3),by=pop]$V1
    stegoMSE <- mean(Vector1DifStego^2+Vector2DifStego^2+Vector3DifStego^2)
    
    totalmseVARCOV <- mean(components.dt[method=="varcov"]$Vector1^2+components.dt[method=="varcov"]$Vector2^2+components.dt[method=="varcov"]$Vector3^2)
    Vector1Difvarcov <- components.dt[method=="varcov",Vector1-mean(Vector1),by=pop]$V1
    Vector2Difvarcov <- components.dt[method=="varcov",Vector2-mean(Vector2),by=pop]$V1
    Vector3Difvarcov <- components.dt[method=="varcov",Vector3-mean(Vector3),by=pop]$V1
    varcovMSE <- mean(Vector1Difvarcov^2+Vector2Difvarcov^2+Vector3Difvarcov^2)
    c("stego"=stegoMSE/totalmseSTEGO, "varcov"=varcovMSE/totalmseVARCOV)  
}
withinVsBetween2 <-  function(components.dt){
    pair <- unique(components.dt$pop)
    totalmseSTEGO <- mean(components.dt[method=="STEGO"]$Vector1^2+components.dt[method=="STEGO"]$Vector2^2)
    Vector1DifStego <- components.dt[method=="STEGO",Vector1-mean(Vector1),by=pop]$V1
    Vector2DifStego <- components.dt[method=="STEGO",Vector2-mean(Vector2),by=pop]$V1
    stegoMSE <- mean(Vector1DifStego^2+Vector2DifStego^2)
    
    totalmseVARCOV <- mean(components.dt[method=="varcov"]$Vector1^2+components.dt[method=="varcov"]$Vector2^2)
    Vector1Difvarcov <- components.dt[method=="varcov",Vector1-mean(Vector1),by=pop]$V1
    Vector2Difvarcov <- components.dt[method=="varcov",Vector2-mean(Vector2),by=pop]$V1
    varcovMSE <- mean(Vector1Difvarcov^2+Vector2Difvarcov^2)
    c("stego"=stegoMSE/totalmseSTEGO, "varcov"=varcovMSE/totalmseVARCOV)  
}



firstTwoCompRSq <- function(components.dt){
    c("STEGO"=summary(glm((as.numeric(as.factor(components.dt[method=="STEGO",pop]))-1)~components.dt[method=="STEGO",Vector1]+components.dt[method=="STEGO",Vector2],family=binomial(link='logit')))$r.squared,
    "varcov"=summary(glm((as.numeric(as.factor(components.dt[method=="varcov",pop]))-1)~components.dt[method=="varcov",Vector1]+components.dt[method=="varcov",Vector2],family=binomial(link='logit')))$r.squared)
}
firstThreeCompRSq <- function(components.dt){
    c("STEGO"=summary(lm((as.numeric(as.factor(components.dt[method=="STEGO",pop]))-1)~components.dt[method=="STEGO",Vector1]+components.dt[method=="STEGO",Vector2]+components.dt[method=="STEGO",Vector3]))$r.squared,
      "varcov"=summary(lm((as.numeric(as.factor(components.dt[method=="varcov",pop]))-1)~components.dt[method=="varcov",Vector1]+components.dt[method=="varcov",Vector2]+components.dt[method=="varcov",Vector3]))$r.squared)
}

firstThreeCompRSq <- function(components.dt){
    xx <- as.matrix(data.frame(aa=components.dt[method=="STEGO",Vector1],bb=components.dt[method=="STEGO",Vector2],cc=components.dt[method=="STEGO",Vector3]))
    zz <- as.matrix(data.frame(aa=components.dt[method=="varcov",Vector1],bb=components.dt[method=="varcov",Vector2],cc=components.dt[method=="varcov",Vector3]))
    c("STEGO"=lrm(as.factor(components.dt[method=="STEGO",pop])~xx)$stats['R2'],
      "varcov"=lrm(as.factor(components.dt[method=="STEGO",pop])~zz)$stats['R2'])
}
firstCompRSq <- function(components.dt){
    c("STEGO"=summary(lm(components.dt[method=="STEGO",Vector1]~components.dt[method=="STEGO",pop]))$r.squared,
      "varcov"=summary(lm(components.dt[method=="varcov",Vector1]~components.dt[method=="varcov",pop]))$r.squared,
      "STEGO2"=summary(lm(components.dt[method=="STEGO",Vector2]~components.dt[method=="STEGO",pop]))$r.squared,
      "varcov2"=summary(lm(components.dt[method=="varcov",Vector2]~components.dt[method=="varcov",pop]))$r.squared)
}
####################################
####################################
####################################


# results <- readRDS("~/1000GP/plots/s_distributions/plotdata/STU_ITU_data.rds")
# results <- readRDS("~/1000GP/plots/s_distributions/plotdata/CEU_YRI_data.rds")
# results <- readRDS("~/1000GP/plots/s_distributions/plotdata/IBS_TSI_data.rds")
results <- readRDS("~/1000GP/allVariants_CHR14/GBR-FIN.rds")



sampleSubset <- T
# sampleSubset <- getSampleSubset("all")
# sampleSubset <- getSampleSubset(c("TSI","IBS"))
# sampleSubset <- getSampleSubset(c("ITU","STU"))
# sampleSubset <- getSampleSubset(c("CEU","YRI"))
# sampleSubset <- getSampleSubset(c("CHB","CHS"))
# sampleSubset <- getSampleSubset(c("GBR","CEU"))
# sampleSubset <- getSampleSubset(c("GBR","FIN"))
# sampleSubset <- getSampleSubset(c("GIH","ITU"))
# sampleSubset <- getSampleSubset(c("ESN","LWK"))
# sampleSubset <- getSampleSubset(c("PJL","BEB"))
# sampleSubset <- getSampleSubset(c("CDX","CHB"))
# sampleSubset <- getSampleSubset(c("CDX","CHS"))
# sampleSubset <- getSampleSubset(c("PUR","ACB"))
# sampleSubset <- getSampleSubset(c("KHV","CHS"))
# sampleSubset <- getSampleSubset(c("ITU","BEB"))
# 
# sampleSubset <- qcFilter&sampleSubset

pairwiseResults <- readRDS("./data/pairwiseResults1000.rds")[1:57]


directory <- paste0("./block19_",0:4)
# directory <- "./block19_cor_minVariants2/"
# directory <- "./block19_minvariants5/"
directory <- "./allVariants_CHR14/"

pairwiseResults <- lapply(directory, function(dir) lapply(list.files(dir, full.names = T), readRDS))

combineResults <- mapply(function(resA, resB, resC, resD, resE){
    list(simMat=(resA$simMat+resB$simMat+resC$simMat+resD$simMat+resE$simMat)/5, 
         s_matrix_dip=(resA$s_matrix_dip+resB$s_matrix_dip+resC$s_matrix_dip+resD$s_matrix_dip+resE$s_matrix_dip)/5)
}, pairwiseResults[[1]], pairwiseResults[[2]], pairwiseResults[[3]], pairwiseResults[[4]], pairwiseResults[[1]], SIMPLIFY = FALSE)

pairwiseResults <- combineResults
names(pairwiseResults) <- sub("-","_",sub(".rds","",list.files(directory[1])))

results <- pairwiseResults[["STU_ITU"]]
results <- pairwiseResults[["BEB_STU"]]
results <- pairwiseResults[["IBS_TSI"]]
results <- pairwiseResults[["CHS_CHB"]]
results <- pairwiseResults[["GBR_CEU"]]
results <- pairwiseResults[["ESN_YRI"]]
results <- pairwiseResults[["CDX_KHV"]]
results <- pairwiseResults[["CLM_MXL"]]
results <- pairwiseResults[["PJL_STU"]]
results <- pairwiseResults[["ACB_ESN"]]
results <- pairwiseResults[["ACB_ASW"]]
results <- pairwiseResults[["ACB_LWK"]]
results <- pairwiseResults[["ACB_YRI"]]
results <- pairwiseResults[["PEL_MXL"]]
results <- pairwiseResults[["CLM_PEL"]]
results <- pairwiseResults[["PUR_CLM"]]
results <- pairwiseResults[["GBR_FIN"]]
results <- pairwiseResults[["CEU_TSI"]]

jaccardMatrix <- results$s_matrix_dip
simMat <- results$simMat
jaccard.correction <- eigen(apply(jaccardMatrix[sampleSubset,sampleSubset],2,scale),symmetric=T)$vectors
varcov.correction <- eigen(apply(simMat[sampleSubset,sampleSubset],2,scale),symmetric=T)$vectors
# varcov.correction[,3] <- -varcov.correction[,3] 

components.dt <- generateggplotDT(jaccard.correction, varcov.correction, rownames(jaccardMatrix[sampleSubset,sampleSubset]), c(1,2,3))


components.dt$rank <- c(rank(-(components.dt[method=="STEGO"]$Vector1)^2 - (components.dt[method=="STEGO"]$Vector2)^2))
withinVsBetween2(components.dt)
# names(components.dt)[2:4] <- paste0("Vector",c(2,3,1))
# pdf(paste0("./plots/s_distributions/",outputDir,"/PCA_.pdf", height=5, width=10)
plotPCAsSideBySide(components.dt[rank>11], super=F, labelOutliers=4)
# dev.off()

# withinVsBetween(components.dt)


# Get the FSTs
# pairwiseFSTs <- unlist(lapply(pairwiseResults, function(pairresult){
#     pairresult$FST
# }))
# sort(pairwiseFSTs)



pairwiseWithinVsBetween <- lapply(pairwiseResults, function(pairresult){
    jaccardMatrix <- pairresult$s_matrix_dip
    simMat <- pairresult$simMat
    
    jaccard.correction <- eigen(apply(jaccardMatrix[sampleSubset,sampleSubset],2,scale),symmetric=T)$vectors
    sim.correction <- eigen(apply(simMat[sampleSubset,sampleSubset],2,scale),symmetric=T)$vectors
    # varcov.correction[,3] <- -varcov.correction[,3] 
    
    components.dt <- generateggplotDT(jaccard.correction, sim.correction, rownames(jaccardMatrix[sampleSubset,sampleSubset]), c(1,2,3))
    asd <- NULL
    tryCatch(return(withinVsBetween(components.dt)), 
             error = function(err){return(c(STEGO.R2=1,varcov.R2=1))} ,
             finally = cat(".")
    )
    return(asd)

})

pairwiseWithinVsBetween <- lapply(pairwiseWithinVsBetween, function(x){
    if(is.null(x)){
        return(c(STEGO.R2=1,varcov.R2=1))
    } else {return(x)}
    }
)
pairwiseWithinVsBetweenDT <- data.table(do.call(rbind,pairwiseWithinVsBetween),keep.rownames="pair")
# pairwiseWithinVsBetweenDT$FST <-pairwiseFSTs[pairwiseWithinVsBetweenDT$pair]
pairwiseWithinVsBetweenDT$pair <- factor(pairwiseWithinVsBetweenDT$pair, levels=pairwiseWithinVsBetweenDT$pair[order(pairwiseWithinVsBetweenDT$varcov)])
names(pairwiseWithinVsBetweenDT) <- c("pair","STEGO","varcov")
pairwiseWithinVsBetweenDT <- pairwiseWithinVsBetweenDT[order(varcov)]
pairwisePlot <- ggplot(melt(pairwiseWithinVsBetweenDT, id='pair', variable.name="Method")) + 
    theme_bw() +
    theme(plot.title = element_text(size=22)) + 
    geom_point(aes(x=pair,y=value, color=Method, shape=Method),size=3) + 
    scale_color_manual(values = c("STEGO" = 'red','varcov' = 'blue'), labels=c("STEGO", "PCA")) +
    scale_shape_manual(values = c('STEGO' = 17, 'varcov' = 16), labels=c("STEGO", "PCA")) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5)) + 
    ggtitle("Comparison of recent ancestry separation: STEGO vs PCA") + xlab("Population Pair") + ylab("Ratio of within-population variance to total variance")
# png("~/1000GP/plots/pairwise_variances.png", width=1000)
pairwisePlot
# dev.off()

# ggplot(pairwiseWithinVsBetweenDT) + geom_point(aes(x=pair,y=STEGO),color="red",size=3) + geom_point(aes(x=pair,y=varcov),size=3)+ geom_point(aes(x=pair,y=STEGO2),color="blue",size=3) + geom_point(aes(x=pair,y=varcov2),color="grey",size=3) +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1))
t.test(pairwiseWithinVsBetweenDT$STEGO-pairwiseWithinVsBetweenDT$varcov)
binom.test(sum(pairwiseWithinVsBetweenDT$STEGO-pairwiseWithinVsBetweenDT$varcov<0),nrow(pairwiseWithinVsBetweenDT),p=.5)
pairwiseWithinVsBetweenDT <- pairwiseWithinVsBetweenDT[order(varcov)]
colMeans(pairwiseWithinVsBetweenDT[,-1,with=F])

sum(pairwiseWithinVsBetweenDT$STEGO<pairwiseWithinVsBetweenDT$varcov)



pairwiseWithinVsBetweenDT0 <- pairwiseWithinVsBetweenDT
# Merging multiple sources
pairwiseWithinVsBetweenDT_all <- merge(pairwiseWithinVsBetweenDT0,pairwiseWithinVsBetweenDT1,by="pair")
ggplot(pairwiseWithinVsBetweenDT_all) + 
    geom_point(aes(x=pair,y=STEGO.x),color="red",size=3) + geom_point(aes(x=pair,y=varcov.x),color="grey",size=3)+ 
    geom_point(aes(x=pair,y=STEGO.y),color="blue",size=3) + geom_point(aes(x=pair,y=varcov.y),size=3)+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

pairwiseWithinVsBetweenDT_all$STEGO <- (pairwiseWithinVsBetweenDT_all$STEGO.x+pairwiseWithinVsBetweenDT_all$STEGO.y)/2
pairwiseWithinVsBetweenDT_all$varcov <- (pairwiseWithinVsBetweenDT_all$varcov.x +pairwiseWithinVsBetweenDT_all$varcov.y)/2
sum(pairwiseWithinVsBetweenDT_all$STEGO<pairwiseWithinVsBetweenDT_all$varcov)
binom.test(sum(pairwiseWithinVsBetweenDT_all$STEGO-pairwiseWithinVsBetweenDT_all$varcov<0),nrow(pairwiseWithinVsBetweenDT_all),p=.5)
