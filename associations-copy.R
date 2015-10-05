library(MASS)
library(ggplot2)

chr <- 1
numRows <- 1000

# args<-commandArgs(TRUE)
# chr <- as.numeric(args[1])
# numRows <- as.numeric(args[2])

sample <- read.table("~/1000GP/data/1000GP_Phase3.sample", sep=" ", header=T)
sampleIDs <- as.character(sample[,1])
pop <- as.character(sample[,2])
group <- as.character(sample[,3])
sex <- as.character(sample[,4])



getResultDF <- function(filename, phenotype, numRows, linesAtATime=1000, correctionVec=NULL, varFactor=T, sampleSubset=NULL){
    
    if(is.null(correctionVec)){
        correction <- function(geno){
            fitted <- mean(geno)
            list("corrected"=(geno-fitted),"fitted"=rep(fitted, length(geno)))
        }
    } else {
        x <- cbind(1, correctionVec)
        correctionHat  <- x %*% ginv(t(x)%*%x) %*% t(x)
        correction <- function(geno){
            fitted <- correctionHat %*% geno
            # Dan's correction
            fitted[fitted<.001]<-.001
            fitted[fitted>.999]<-.999
            list("corrected"=(geno - fitted),"fitted"=fitted)
        }
    }
    # Unless specified, use all the samples
    if(is.null(sampleSubset)){
        sampleSubset <- rep(T, numSamples)
    } else {
        numSamples = sum(sampleSubset)
    }
    phenotypeMatrix <- do.call(cbind, phenotype)
    phenotypeMatrixCorrected <- correction(phenotypeMatrix)
#     varY <- sqrt(phenotypeMatrixCorrected$fitted*(1-phenotypeMatrixCorrected$fitted))
    con <- file(filename, "rt")
    # fix this to stop at end of file
    resultMatrix <- t(do.call(cbind, lapply(1:(numRows/linesAtATime), function(i){
        genotypes <- apply(do.call(rbind, strsplit(readLines(con, linesAtATime)," ")), 1,as.numeric)
        alleles <- genotypes[c(T,F),] + genotypes[c(F,T),]
        
        # Subset of samples
        alleles <- alleles[sampleSubset,]
        
        numAlleles <- colSums(alleles)
        # apply correction here
        allelesCorrected <- correction(alleles)
        alleles <- allelesCorrected$corrected
        phenotypeMatrix <- phenotypeMatrixCorrected$corrected
        if(is.null(correctionVec)|!varFactor){
            varianceFactor <- matrix(numSamples, nrow=ncol(phenotypeMatrix), ncol=ncol(alleles))
        } else {
            # New variance factor
            xProb <- allelesCorrected$fitted/2
            xVariances <- 2*xProb*(1-xProb) #variance of Binomial(2,xProb)
            yVariances <- phenotypeMatrixCorrected$fitted*(1-phenotypeMatrixCorrected$fitted)
            
            zvar <- t(yVariances)%*%xVariances
            numerator <- zvar
            denominator <- colSums(yVariances)%*%t(colSums(xVariances))
            varRsq <- numerator/denominator
            varianceFactor <- 1/varRsq
            
        }
        print(dim(varianceFactor))
        # HW sandwich estimator.
        # speed up with matrix operations?
#         xTx.inv_all <- 1/colSums(alleles^2) # sapply(colSums(alleles^2),ginv)
#         x.hat <- t(sweep(alleles, MARGIN=2, xTx.inv_all, '*'))
#         betas <- x.hat%*%phenotypeMatrix
#         sigma <- (phenotypeMatrix-alleles%*%betas)^2  # The diagonal of the Sigma matrix
#         seHWs <- sqrt(xTx.inv_all^2*t(alleles^2)%*%sigma)
#         chisqHW2 <- t((betas/seHWs)^2)
        
#         chisqHW <- apply(phenotypeMatrix, 2, function(phenoCol){
#             apply(alleles, 2, function(alleleCol){
#                 xTx.inv <- ginv(t(alleleCol)%*%alleleCol)
#                 betaOLS <- xTx.inv%*%t(alleleCol)%*%phenoCol
#                 sigma <- diag(c(phenoCol - alleleCol%*%betaOLS)^2)
#                 seHW <- sqrt(xTx.inv%*%t(alleleCol)%*%sigma%*%alleleCol%*%xTx.inv)
#                 chisq <- (betaOLS/seHW)^2
#                 chisq
#             })
#         })
#         mod <- lm(y~x)
#         p_lm <- summary(mod)$coef[2,4]
#         sand_SE <- sqrt(vcovHC(mod, type="HC0")["x","x"])
#         z_sand <- coef(mod)["x"]/sand_SE
#         p_sand <- pchisq(z_sand^2, 1, lower.tail=FALSE)
#         kendall <- cor.test(x, y, method = "kendall")$p.value

        # Change to t-distribution 8/26/15
#         rbind(numAlleles, varianceFactor*(cor(phenotypeMatrix, alleles)^2))
        r.values <- cor(phenotypeMatrix, alleles)
        negLogPValue <- -log(1-pt(r.values*sqrt((varianceFactor-2)/(1-r.values^2)), varianceFactor))
        rbind(numAlleles, negLogPValue, varianceFactor)

#         rbind(numAlleles, chisqHW2)
        
    })))
    close(con)
    colnames(resultMatrix)[1] <- "NumAlleles"
    resultMatrix <- data.frame(resultMatrix)
    #Filter out 0,1 allele SNPs
    filter <- resultMatrix[,1]>(numSamples*2/50) & (numSamples*2-resultMatrix[,1])>(numSamples*2/50)
    resultMatrix <- resultMatrix[filter,]
    
    obsPValues <- rep(NA, nrow(resultMatrix))
    rareIndices <- resultMatrix$NumAlleles>2 & resultMatrix$NumAlleles<=(numSamples*2*.04)
    lowIndices <- resultMatrix$NumAlleles>(numSamples*2*.04) & resultMatrix$NumAlleles<=(numSamples*2*.1)
    commonIndices <- resultMatrix$NumAlleles>(numSamples*2*.1)
    resultMatrix <- resultMatrix[,-1]

    # 1st half of matrix is test statistics
    # 2nd half of matrix is effective samples
    degFreedom <- resultMatrix[,((ncol(resultMatrix))/2+1):(ncol(resultMatrix))]
    resultMatrix <- resultMatrix[1:(ncol(resultMatrix)/2)]
    nIt <- ncol(resultMatrix)/3 # number of phenotype risks (uniform, gradual, sharp)
    combinedResultList <- lapply(list(list(rm=resultMatrix[,1:nIt],df=degFreedom[,1:nIt]), 
                                      list(rm=resultMatrix[,(nIt+1):(2*nIt)],df=degFreedom[,(nIt+1):(2*nIt)]), 
                                      list(rm=resultMatrix[,(2*nIt+1):(3*nIt)],df=degFreedom[,(2*nIt+1):(3*nIt)])), 
                                 function(resultsList){
        results <- resultsList$rm
        df      <- resultsList$df
        rareobsPValues <- apply(apply(results[rareIndices,], 2, sort, decreasing=TRUE),1,function(x){
            mean(x)
#             -log(1-pchisq(mean(x), 1))
#             mean(-log(2*pt(abs(x))))
        })
        lowobsPValues <- apply(apply(results[lowIndices,], 2, sort, decreasing=TRUE),1,function(x){
#             -log(1-pchisq(mean(x), 1))
            mean(x)
        })
        commonobsPValues <- apply(apply(results[commonIndices,], 2, sort, decreasing=TRUE),1,function(x){
#             -log(1-pchisq(mean(x), 1))
            mean(x)
        })
        
        rareobsPValues <- data.frame(cbind(rareobsPValues, -log(1:length(rareobsPValues)/length(rareobsPValues))))
        lowobsPValues <- data.frame(cbind(lowobsPValues, -log(1:length(lowobsPValues)/length(lowobsPValues))))
        commonobsPValues <- data.frame(cbind(commonobsPValues, -log(1:length(commonobsPValues)/length(commonobsPValues))))
        rareobsPValues <- cbind(rareobsPValues, "rare")
        lowobsPValues <- cbind(lowobsPValues, "low")
        commonobsPValues <- cbind(commonobsPValues, "common")
        colnames(rareobsPValues) <- c("Observed", "Expected", "MAF")
        colnames(lowobsPValues) <- c("Observed", "Expected", "MAF")
        colnames(commonobsPValues) <- c("Observed", "Expected", "MAF")
        
        do.call(rbind, list(rareobsPValues,lowobsPValues,commonobsPValues))
    })
    uniformResult <- cbind(combinedResultList[[1]], phenotype="Uniform")
    gradualResult <- cbind(combinedResultList[[2]], phenotype="Gradual")
    sharpResult <- cbind(combinedResultList[[3]], phenotype="Sharp")
    do.call(rbind, list(uniformResult, gradualResult, sharpResult))
}

# Make non-genetic risk
numSamples <- length(sampleIDs)
superPops <- levels(factor(group))
superRisks <- runif(length(superPops),0,.5)
names(superRisks) <- superPops

pops <- levels(factor(pop))
popRisks <- runif(length(pops),0,.2)
names(popRisks) <- pops

sharpRisk <-as.numeric(pop%in%sample(pops,2))

gradualRisk <- replicate(100, rbinom(numSamples, size=1, prob=(superRisks[group] + popRisks[pop])))
sharpRisk   <- replicate(100, rbinom(numSamples, size=1, prob=as.numeric(pop%in%sample(pops,1))/5))
randomRisk  <- replicate(100, rbinom(numSamples, size=1, prob=.25))

superpop.correction <- cbind(as.numeric(group%in%"AFR"),as.numeric(group%in%"AMR"),as.numeric(group%in%"EAS"),as.numeric(group%in%"EUR"))
subpop.correction <- sapply(pops, function(x){as.numeric(pop %in% x)})[,-1]
## End of day on 7/8/15
## TODO: implement correction for sharp phenotype (from jaccard or something)
jaccardMatrix <- read.csv('~/1000GP/output_0_20/combined_jaccard.csv', row.names=1)
diag(jaccardMatrix) <- 1 # This is to undo the 'set diagonal to zero' that was done for visualization in jaccard processing.
jaccard.correction <- eigen(jaccardMatrix)$vectors[,1:10]

# Comparison to PCA
varcov.correction <- eigen(varcovMatrix)$vectors[,1:10]

########### BLOCK COMMENT below replaced by better code 8/6/15

# # make the NULL phenotype plots
# # run it n times, sort each by pvalue and average over all runs
# results.uncorrected  <- getResultDF(paste('~/1000GP/data/1000GP_Phase3_chr', chr, '.hap.gz',sep=''), phenotype= list(randomRisk, gradualRisk, sharpRisk), numRows=50000, linesAtATime=1000, correctionVec = NULL, varFactor=T)
# results.jaccard.cor  <- getResultDF(paste('~/1000GP/data/1000GP_Phase3_chr', chr, '.hap.gz',sep=''), phenotype= list(randomRisk, gradualRisk, sharpRisk), numRows=50000, linesAtATime=1000, correctionVec = jaccard.correction, varFactor=T)
# results.varcov.cor  <- getResultDF(paste('~/1000GP/data/1000GP_Phase3_chr', chr, '.hap.gz',sep=''), phenotype= list(randomRisk, gradualRisk, sharpRisk), numRows=50000, linesAtATime=1000, correctionVec = varcov.correction, varFactor=T)
# results.superpop.cor <- getResultDF(paste('~/1000GP/data/1000GP_Phase3_chr', chr, '.hap.gz',sep=''), phenotype= list(randomRisk, gradualRisk, sharpRisk), numRows=50000, linesAtATime=1000, correctionVec = superpop.correction, varFactor=T)
# results.subpop.cor   <- getResultDF(paste('~/1000GP/data/1000GP_Phase3_chr', chr, '.hap.gz',sep=''), phenotype= list(randomRisk, gradualRisk, sharpRisk), numRows=50000, linesAtATime=1000, correctionVec = subpop.correction, varFactor=T)
# results.uncorrectedATT  <- getResultDF(paste('~/1000GP/data/1000GP_Phase3_chr', chr, '.hap.gz',sep=''), phenotype= list(randomRisk, gradualRisk, sharpRisk), numRows=50000, linesAtATime=1000, correctionVec = NULL, varFactor=F)
# results.jaccard.corATT  <- getResultDF(paste('~/1000GP/data/1000GP_Phase3_chr', chr, '.hap.gz',sep=''), phenotype= list(randomRisk, gradualRisk, sharpRisk), numRows=50000, linesAtATime=1000, correctionVec = jaccard.correction, varFactor=F)
# results.varcov.corATT  <- getResultDF(paste('~/1000GP/data/1000GP_Phase3_chr', chr, '.hap.gz',sep=''), phenotype= list(randomRisk, gradualRisk, sharpRisk), numRows=50000, linesAtATime=1000, correctionVec = varcov.correction, varFactor=F)
# results.superpop.corATT <- getResultDF(paste('~/1000GP/data/1000GP_Phase3_chr', chr, '.hap.gz',sep=''), phenotype= list(randomRisk, gradualRisk, sharpRisk), numRows=50000, linesAtATime=1000, correctionVec = superpop.correction, varFactor=F)
# results.subpop.corATT   <- getResultDF(paste('~/1000GP/data/1000GP_Phase3_chr', chr, '.hap.gz',sep=''), phenotype= list(randomRisk, gradualRisk, sharpRisk), numRows=50000, linesAtATime=1000, correctionVec = subpop.correction, varFactor=F)
# 
# # pairwise related populations
# sampleSubset <- pop %in% c("CHB","CHS")
# gradualRisk <- replicate(100, rbinom(sum(sampleSubset), size=1, prob=(.2+as.numeric(pop%in%"CHB")[sampleSubset]/10)))
# sharpRisk   <- replicate(100, rbinom(sum(sampleSubset), size=1, prob=as.numeric(pop%in%"CHB")[sampleSubset]/5))
# randomRisk  <- replicate(100, rbinom(sum(sampleSubset), size=1, prob=.25))
# jaccard.correction <- eigen(jaccardMatrix[sampleSubset,sampleSubset])$vectors[,1]
# results.jaccard.cor.CHBCHS  <- getResultDF(paste('~/1000GP/data/1000GP_Phase3_chr', chr, '.hap.gz',sep=''), phenotype= list(randomRisk, gradualRisk, sharpRisk), numRows=50000, linesAtATime=1000, correctionVec = jaccard.correction, varFactor=T, sampleSubset=sampleSubset)
# results.uncorrected.CHBCHS  <- getResultDF(paste('~/1000GP/data/1000GP_Phase3_chr', chr, '.hap.gz',sep=''), phenotype= list(randomRisk, gradualRisk, sharpRisk), numRows=50000, linesAtATime=1000, correctionVec = NULL, varFactor=T, sampleSubset=sampleSubset)
# results.varcov.cor.CHBCHS  <- getResultDF(paste('~/1000GP/data/1000GP_Phase3_chr', chr, '.hap.gz',sep=''), phenotype= list(randomRisk, gradualRisk, sharpRisk), numRows=50000, linesAtATime=1000, correctionVec = varcov.correction, varFactor=T, sampleSubset=sampleSubset)
# 
# #write.csv(results, file=paste("~/1000GP/output_associations/chr", chr, ".csv", sep=""))
# 
# png(paste("~/1000GP/plots/chr", chr, "_uncorrected.png",sep=""), width=800)
# ggplot(results.uncorrected, aes(Expected, Observed, color=MAF))+ geom_abline(intercept = 0) + ylim(0, 30) + ggtitle("Nongenetic Phenotype Risk (uncorrected)") + geom_point() + facet_wrap(~ phenotype)
# dev.off()
# png(paste("~/1000GP/plots/chr", chr, "_jaccard_corrected.png",sep=""), width=800)
# ggplot(results.jaccard.cor, aes(Expected, Observed, color=MAF))+ geom_abline(intercept = 0) + ylim(0, 15) + ggtitle("Nongenetic Phenotype Risk (corrected by jaccard)") + geom_point() + facet_wrap(~ phenotype)
# dev.off()
# png(paste("~/1000GP/plots/chr", chr, "_varcov_corrected.png",sep=""), width=800)
# ggplot(results.varcov.cor, aes(Expected, Observed, color=MAF))+ geom_abline(intercept = 0) + ylim(0, 15) + ggtitle("Nongenetic Phenotype Risk (corrected by variance-covariance matrix)") + geom_point() + facet_wrap(~ phenotype)
# dev.off()
# png(paste("~/1000GP/plots/chr", chr, "_superpop_corrected.png",sep=""), width=800)
# ggplot(results.superpop.cor, aes(Expected, Observed, color=MAF))+ geom_abline(intercept = 0) + ylim(0, 15) + ggtitle("Nongenetic Phenotype Risk (corrected by superpop)") + geom_point() + facet_wrap(~ phenotype)
# dev.off()
# png(paste("~/1000GP/plots/chr", chr, "_subpop_corrected.png",sep=""), width=800)
# ggplot(results.subpop.cor, aes(Expected, Observed, color=MAF))+ geom_abline(intercept = 0) + ylim(0, 15) + ggtitle("Nongenetic Phenotype Risk (corrected by subpop)") + geom_point() + facet_wrap(~ phenotype)
# dev.off()
# png(paste("~/1000GP/plots/chr", chr, "_uncorrectedATT.png",sep=""), width=800)
# ggplot(results.uncorrectedATT, aes(Expected, Observed, color=MAF))+ geom_abline(intercept = 0) + ylim(0, 30) + ggtitle("Nongenetic Phenotype Risk (uncorrected)") + geom_point() + facet_wrap(~ phenotype)
# dev.off()
# png(paste("~/1000GP/plots/chr", chr, "_jaccard_correctedATT.png",sep=""), width=800)
# ggplot(results.jaccard.corATT, aes(Expected, Observed, color=MAF))+ geom_abline(intercept = 0) + ylim(0, 15) + ggtitle("Nongenetic Phenotype Risk (corrected by jaccard)") + geom_point() + facet_wrap(~ phenotype)
# dev.off()
# png(paste("~/1000GP/plots/chr", chr, "_varcov_correctedATT.png",sep=""), width=800)
# ggplot(results.varcov.corATT, aes(Expected, Observed, color=MAF))+ geom_abline(intercept = 0) + ylim(0, 15) + ggtitle("Nongenetic Phenotype Risk (corrected by variance-covariance matrix)") + geom_point() + facet_wrap(~ phenotype)
# dev.off()
# png(paste("~/1000GP/plots/chr", chr, "_superpop_correctedATT.png",sep=""), width=800)
# ggplot(results.superpop.corATT, aes(Expected, Observed, color=MAF))+ geom_abline(intercept = 0) + ylim(0, 15) + ggtitle("Nongenetic Phenotype Risk (corrected by superpop)") + geom_point() + facet_wrap(~ phenotype)
# dev.off()
# png(paste("~/1000GP/plots/chr", chr, "_subpop_correctedATT.png",sep=""), width=800)
# ggplot(results.subpop.corATT, aes(Expected, Observed, color=MAF))+ geom_abline(intercept = 0) + ylim(0, 15) + ggtitle("Nongenetic Phenotype Risk (corrected by subpop)") + geom_point() + facet_wrap(~ phenotype)
# dev.off()
# png(paste("~/1000GP/plots/chr", chr, "_jaccard_correctedCHBCHS.png",sep=""), width=800)
# ggplot(results.jaccard.cor.CHBCHS, aes(Expected, Observed, color=MAF))+ geom_abline(intercept = 0) + ylim(0, 15) + ggtitle("Nongenetic Phenotype Risk (corrected by jaccard)") + geom_point() + facet_wrap(~ phenotype)
# dev.off()
# png(paste("~/1000GP/plots/chr", chr, "_varcov_correctedCHBCHS.png",sep=""), width=800)
# ggplot(results.varcov.cor.CHBCHS, aes(Expected, Observed, color=MAF))+ geom_abline(intercept = 0) + ylim(0, 15) + ggtitle("Nongenetic Phenotype Risk (corrected by variance-covariance matrix)") + geom_point() + facet_wrap(~ phenotype)
# dev.off()
# png(paste("~/1000GP/plots/chr", chr, "_uncorrectedCHBCHS.png",sep=""), width=800)
# ggplot(results.uncorrected.CHBCHS, aes(Expected, Observed, color=MAF))+ geom_abline(intercept = 0) + ylim(0, 15) + ggtitle("Nongenetic Phenotype Risk (uncorrected)") + geom_point() + facet_wrap(~ phenotype)
# dev.off()

runAndPlot <- function(chr=1, correctMethod="uncorrected", ATT="ATT", subpops="all", numEigenVectors=2, numRows=5000){
    print(correctMethod)
    print(ATT)
    print(subpops)
    if("all"%in%subpops){
        sampleSubset <- rep(T,length(pop))
        gradualRisk <- replicate(100, rbinom(numSamples, size=1, prob=(superRisks[group] + popRisks[pop])))
        sharpRisk   <- replicate(100, rbinom(numSamples, size=1, prob=as.numeric(pop%in%sample(pops,1))/5))
        randomRisk  <- replicate(100, rbinom(numSamples, size=1, prob=.25))
    } else {
        sampleSubset <- pop %in% subpops
        gradualRisk <- replicate(100, rbinom(sum(sampleSubset), size=1, prob=(.2+as.numeric(pop%in%subpops[1])[sampleSubset]/10)))
        sharpRisk   <- replicate(100, rbinom(sum(sampleSubset), size=1, prob=as.numeric(pop%in%subpops[1])[sampleSubset]/5))
        randomRisk  <- replicate(100, rbinom(sum(sampleSubset), size=1, prob=.25))
    }
    if(ATT=="ATT"){
        varFactor<-F
    } else {
        varFactor<-T
    }
    
    
    if(correctMethod=="uncorrected"){
        correction<-NULL
        varFactor <- F
    } else if(correctMethod=="jaccard"){
        correction <- eigen(jaccardMatrix[sampleSubset,sampleSubset])$vectors[,1:numEigenVectors,drop=F]
    } else if(correctMethod=="varcov"){
        correction <- eigen(varcovMatrix[sampleSubset,sampleSubset])$vectors[,1:numEigenVectors,drop=F]
    } else if(correctMethod=="superpop"){
        correction <- cbind(as.numeric(group%in%"AFR"),as.numeric(group%in%"AMR"),as.numeric(group%in%"EAS"),as.numeric(group%in%"EUR"))[sampleSubset,]
    } else if(correctMethod=="subpop"){
        correction <- sapply(pops, function(x){as.numeric(pop %in% x)})[sampleSubset,-1]
    }
    
    results  <- getResultDF(filename=paste('~/1000GP/data/1000GP_Phase3_chr', chr, '.hap.gz',sep=''), phenotype= list(randomRisk, gradualRisk, sharpRisk), numRows=numRows, linesAtATime=1000, correctionVec = correction, varFactor=varFactor, sampleSubset=sampleSubset)
    png(paste("~/1000GP/plots/chr", chr,"_",correctMethod,"_", ATT, "_",paste(subpops, collapse=""),".png",sep=""), width=1200)
    print(ggplot(results, aes(Expected, Observed, color=MAF))+ geom_abline(intercept = 0) + ylim(0, 20) + ggtitle(paste("Nongenetic Phenotype Risk (",correctMethod,") - ",paste(subpops, collapse="-"),sep="")) + geom_point() + facet_wrap(~ phenotype))
    write.csv(results, file=paste("~/1000GP/output_associations/", chr,"_",correctMethod,"_", ATT, "_",paste(subpops, collapse=""),".csv", sep=""))
    dev.off()
    print(paste0("plot created: ","~/1000GP/plots/chr", chr,"_",correctMethod,"_", ATT, "_",paste(subpops, collapse=""),".png"))
}
