library(foreach)
library(doParallel)
library(tools)
library(hierfstat)

homogeneousSimulations <- function(numSimulatedSamples=200, nVariants=50000, numSimulations=10, cok=.0625, filename="./data/combinedFiltered1000.gz", 
                                   numberOfLines=-1, minVariants=5, qcFilter=NULL, ldPrune=1, numCores=4, outputDir='.'){
#     print("Starting read file")
#     system.time(genotypes <- fread(paste('zcat',filename), sep=" ", nrows=numberOfLines, header=F))
#     print("Finished read file")
#     nVariants <- 10000
    registerDoParallel(makeCluster(4))
    results <- list()
#     genotypes <- genotypes[, 1:numSimulatedSamples, with=F]
#     genotypes[,Related := NA ] 
    results <- foreach(i=1:numSimulations,.export=c("generateSResultsFromGenotypes"), .packages=c("data.table","gtools")) %dopar% {
        print(i)
        
        genotypes <- data.table(matrix(rbinom(numSimulatedSamples*nVariants,1,.1), ncol=numSimulatedSamples))
        if(is.na(cok)){
            
        } else {
            ibdVariants <- rbinom(nVariants, 1, cok)==1
            genotypes[[numSimulatedSamples]][ibdVariants] <- genotypes[[1]][ibdVariants]
    #         genotypes[, Related := relatedSample]
        }
        names(genotypes)[1:numSimulatedSamples] <- paste0("Sample",1:numSimulatedSamples)
        genotypes <- pruneGenotypes(genotypes, ldPrune)
        generateSResultsFromGenotypes(paste0("Simulated",i), genotypes, minVariants, scaleBySampleAF=T, outputDir, saveResult=F)
    }
    names(results) <- paste0("Simulated",1:numSimulations)
    results
}
getPopResults <- function(results){
    as.data.table(t(sapply(names(results), function(pop_i){
        s_vector <- sort(results[[pop_i]]$s_matrix_hap[row(results[[pop_i]]$s_matrix_hap)>col(results[[pop_i]]$s_matrix_hap)], decreasing=T)
        topKinship <- (s_vector[1]-1)/(results[[pop_i]]$pkweightsMean-1)
        btest <- binom.test(sum(s_vector>mean(s_vector)), length(s_vector), alternative="less")
        structureKSTest <- ks.test((s_vector-1)/sqrt(results[[pop_i]]$var_s_hap), "pnorm", alternative = c("less"))$p.value
        crypticSig <- ifelse((s_vector[1]-1)/sqrt(results[[pop_i]]$var_s_hap) > qnorm(1-.005/length(s_vector)), "YES+",
                             ifelse((s_vector[1]-1)/sqrt(results[[pop_i]]$var_s_hap) > qnorm(1-.025/length(s_vector)),"YES","NO"))
        structureSig <- ifelse(structureKSTest<.01, "YES+",ifelse(structureKSTest<.05,"YES","NO"))
        c(structurePValue=btest$p.value, var_s=results[[pop_i]]$var_s_hap, sampleVariance=var(s_vector),
          structureKSTest=structureKSTest, closestRelatives=topKinship, crypticSig=crypticSig, structureSig=structureSig)
    })), keep.rownames=T)
}

calculateSMatrix <- function(subpop="each", 
                             filename="./data/combinedFiltered1000.gz", 
                             numberOfLines=-1, 
                             minVariants=5, 
                             qcFilter=NULL, 
                             ldPrune=1, 
                             computeFST=T,
                             outputDir='.'){
    
    print("Starting read file")
    freadString <- ifelse(file_ext(filename)=="gz", paste('zcat',filename), filename)
    system.time(genotypes <- fread(freadString, sep=" ", nrows=numberOfLines, header=F))
    print("Finished read file")

    names(genotypes) <- hap.sampleIDs
    if (subpop=="All"){
        genotypes <- pruneGenotypes(genotypes, ldPrune)
        generateSResultsFromGenotypes("Allsamples", genotypes, minVariants, ldPrune)
        stop("Finished running all")
    }
    if(length(subpop)>1){
        print(subpop)
        if(is.null(qcFilter)){
            filterhap <- hap.pop%in%subpop
            filterdip <- pop%in%subpop 
        } else {
            filterhap <- hap.pop%in%subpop & rep(qcFilter,each=2)
            filterdip <- pop%in%subpop & qcFilter
        }
        hapsampleNames <- hap.sampleIDs[filterhap]
        dipsampleNames <- sampleIDs[filterdip]
        genotypes <- pruneGenotypes(genotypes[,filterhap, with=F], ldPrune)
        return(generateSResultsFromGenotypes(subpop, genotypes, minVariants=minVariants, outputDir=outputDir))
#         stop("Finished running combined")
    }
    
    results <- list()
    if (subpop=="each"){
        for(subpop in unique(pop)){
            print(gc())
            print(subpop)
            if(is.null(qcFilter)){
                filterhap <- hap.pop%in%subpop
                filterdip <- pop%in%subpop 
            } else {
                filterhap <- hap.pop%in%subpop & rep(qcFilter,each=2)
                filterdip <- pop%in%subpop & qcFilter
            }
            hapsampleNames <- hap.sampleIDs[filterhap]
            dipsampleNames <- sampleIDs[filterdip]
            genotypes <- pruneGenotypes(genotypes[,filterhap, with=F], ldPrune)
            results[[subpop]] <- generateSResultsFromGenotypes(subpop, genotypes, qcFilter, minVariants, outputDir)
        }
        names(results) <- unique(pop)
    }
    if (subpop=="pairwise"){
        for(continent in unique(group)){
            pairs <- combn(unique(pop[group==continent]),2)
            for(i in seq_len(ncol(pairs))){
                pair <- pairs[,i]
                print(gc())
                print(pair)
                if(is.null(qcFilter)){
                    filterhap <- hap.pop%in%pair
                    filterdip <- pop%in%pair 
                } else {
                    filterhap <- hap.pop%in%pair & rep(qcFilter,each=2)
                    filterdip <- pop%in%pair & qcFilter
                }
                hapsampleNames <- hap.sampleIDs[filterhap]
                dipsampleNames <- sampleIDs[filterdip]
                genotypesSubpop <- pruneGenotypes(genotypes[,filterhap, with=F], ldPrune)
                
                
                results[[paste(pair,collapse="_")]] <- generateSResultsFromGenotypes(subpop=pair, genotypesSubpop=genotypesSubpop, minVariants=minVariants, outputDir=outputDir)

                if(computeFST){
                    fstData <- as.data.frame(t(as.matrix(as.data.frame(genotypesSubpop)[,c(T,F)] + as.data.frame(genotypesSubpop)[,c(F,T)])))
                    fstData <- cbind(pop = pop[filterdip], fstData)
                    results[[paste(pair,collapse="_")]]$FST <- genet.dist(fstData,diploid=F,method="Fst")
                } 
            }
        }
    }
    
    results    
}

pruneGenotypes <-  function(genotypesSubpop, ldPrune=1){
    
    names(genotypesSubpop) <- make.unique(names(genotypesSubpop))
    numSamples <- ncol(genotypesSubpop)
    numVariants <- nrow(genotypesSubpop)
    sumVariants <- rowSums(genotypesSubpop)
    
    
    # reverse so that MAF<.5
    genotypesSubpop[sumVariants>(numSamples/2),] <- 1-genotypesSubpop[sumVariants>(numSamples/2),]
    sumVariants <- rowSums(genotypesSubpop)
    
    # Intelligently LD prune
    numblocks <- numVariants/ldPrune +1
    blocks <- rep(1:numblocks, each=ldPrune)[1:numVariants]
    
    system.time(runningWhichMax <- running(sumVariants,width=ldPrune,fun=which.max, by=ldPrune))
    prunedIndices <- runningWhichMax + seq(0,ldPrune*(length(runningWhichMax)-1),ldPrune)
    system.time(genotypesSubpop <- genotypesSubpop[prunedIndices])
    genotypesSubpop
}

generateSResultsFromGenotypes <- function(subpop, genotypesSubpop, minVariants=5, scaleBySampleAF=F, outputDir=".", saveResult=T, varcov=T){
    
#     names(genotypesSubpop) <- make.unique(names(genotypesSubpop))
    numSamples <- ncol(genotypesSubpop)
    numVariants <- nrow(genotypesSubpop)
    sumVariants <- rowSums(genotypesSubpop)
#     
#     
#     # reverse so that MAF<.5
    invertMinorAllele <- sumVariants>(numSamples/2)
    genotypesSubpop[invertMinorAllele] <- 1-genotypesSubpop[invertMinorAllele]
    sumVariants <- rowSums(genotypesSubpop)
#     
#     # Intelligently LD prune
#     numblocks <- numVariants/ldPrune +1
#     blocks <- rep(1:numblocks, each=ldPrune)[1:numVariants]
#     
#     system.time(runningWhichMax <- running(sumVariants,width=ldPrune,fun=which.max, by=ldPrune))
#     prunedIndices <- runningWhichMax + seq(0,ldPrune*(length(runningWhichMax)-1),ldPrune)
#     system.time(genotypesSubpop <- genotypesSubpop[prunedIndices])
    
    # remove < n variants
    sumVariants <- rowSums(genotypesSubpop)
    genotypesSubpop <- genotypesSubpop[sumVariants>minVariants]
    genotypesSubpop <- as.matrix(genotypesSubpop)

    print("Number of used variants")
    print(nrow(genotypesSubpop))
    numFilteredVariants <- nrow(genotypesSubpop)
    sumFilteredVariants <- rowSums(genotypesSubpop)
    varcovMat <- NULL
    if(varcov){
#         varcovMat <- cov(t(scale(t(genotypesSubpop[,c(T,F)] + genotypesSubpop[,c(F,T)]))),use="pairwise.complete.obs")
#         varcovMat <- cov(genotypesSubpop[,c(T,F)] + genotypesSubpop[,c(F,T)],use="pairwise.complete.obs")
        varcovMat <- cor(genotypesSubpop[,c(T,F)] + genotypesSubpop[,c(F,T)],use="pairwise.complete.obs")
    }
    totalPossiblePairs <- choose(numSamples,2)
    totalPairs <- choose(sumFilteredVariants,2)
    weights <- totalPossiblePairs/totalPairs
    p <- 1/weights

    # var_s_hap <- sum((1-p)/p)/(numFilteredVariants^2)
    # recalculated variance 3/23/16
    var_s_hap <- sum(weights-1)/(numFilteredVariants^2)
    
    print("variance of s (haploid)")
    print(var_s_hap)
    
    # Calculate expected values conditional on kinship
    pkweightsMean <- mean(((sumFilteredVariants-2)/numSamples)*weights)
    kinships <- seq(0,.25,.001)
    kinshipExpectation <- 1+kinships*(pkweightsMean-1)
        
    s_matrix_numerator <- t(genotypesSubpop*weights)%*%genotypesSubpop
    s_matrix_denominator <- numFilteredVariants
    s_matrix_hap <- s_matrix_numerator/s_matrix_denominator
    if (scaleBySampleAF){
        numAllelesPerSample <- colSums(genotypesSubpop)
        relativeAF <- sqrt(mean(numAllelesPerSample)/numAllelesPerSample)
#         s_matrix_hap <- s_matrix_hap*tcrossprod(relativeAF)
    }
    colnames(s_matrix_hap) <- colnames(genotypesSubpop)
    rownames(s_matrix_hap) <- colnames(genotypesSubpop)
    
    print(mean(s_matrix_hap[row(s_matrix_hap)!=col(s_matrix_hap)]))
    print(median(s_matrix_hap[row(s_matrix_hap)!=col(s_matrix_hap)]))
    
    estimatedKinship <- (s_matrix_hap-1)/(pkweightsMean-1)
    
    # Collapse to diploid
    s_matrix_dip <- (s_matrix_hap[c(T,F),c(T,F)] + s_matrix_hap[c(F,T),c(T,F)] +s_matrix_hap[c(T,F),c(F,T)] + s_matrix_hap[c(F,T),c(F,T)])/4
    colnames(s_matrix_dip) <- colnames(genotypesSubpop)[c(T,F)]
    rownames(s_matrix_dip) <- colnames(genotypesSubpop)[c(T,F)]
    # very lazy variance estimate...
    var_s_dip <- var_s_hap/4
    

    popResult <- list(s_matrix_dip=s_matrix_dip, s_matrix_hap=s_matrix_hap, pkweightsMean=pkweightsMean, var_s_dip=var_s_dip, var_s_hap=var_s_hap, varcovMat=varcovMat)
    if(saveResult){
        saveRDS(popResult, paste0("./plots/s_distributions/",outputDir,"/plotdata/",paste(subpop, collapse = "_"),"_data.rds"))
    }
    popResult

}

plotFromGSM <- function(subpop, gsm, var_s, pkweightsMean, plotname="", outputDir=".",alphaCutoff=.01){
    print(mean(gsm[row(gsm)!=col(gsm)]))
    print(median(gsm[row(gsm)!=col(gsm)]))
    num_comparisons_dip <- choose(ncol(gsm),2)
    sample_IDs <- rownames(gsm)
    var_s <- var(gsm[row(gsm)!=col(gsm)])
    bonferroni_cutoff_dip <- qnorm((1-alphaCutoff/2)^(1/num_comparisons_dip), sd=sqrt(var_s)) + 1
    
    topValuesDip <- sort(gsm[row(gsm)>col(gsm)], decreasing=T)
    topValuesKinship <- (topValuesDip-1)/(pkweightsMean-1)
    
    ks.pvalue <- ks.test((topValuesDip-1)/sd(topValuesDip), "pnorm", alternative = c("less"))$p.value
    ksString <- ifelse(ks.pvalue<.001,"p<.001",paste0("p=",round(ks.pvalue,3)))
    # Display only those that are above the cutoff and among the top 5
    label_cutoff <- max(bonferroni_cutoff_dip, topValuesDip[1])
 
    pairs <- outer(sample_IDs, sample_IDs, paste)
    plotData <- data.frame(values=gsm[row(gsm)>col(gsm)], pairs=paste0("  ",pairs[row(pairs)>col(pairs)]))
    minDip <- min(plotData$values)
    maxDip <- max(plotData$values)#ifelse(max(plotData$values)>bonferroni_cutoff_dip,max(plotData$values),NA)
    xmin <- min(minDip, 1-(maxDip-1)*.5)
    xmax <- maxDip + (maxDip-minDip)*.4
    dipPlot <- ggplot(plotData, aes(values)) + 
        geom_histogram(color="blue",binwidth=.008,fill=I("blue")) + 
        ggtitle(paste0(subpop,collapse="_"))  + xlab("Similarity score") + 
#         scale_x_continuous(expand=c(.4,0))+#
        xlim(xmin, xmax) +
        theme_classic() +
        theme(plot.title = element_text(size=40), axis.title.x = element_text(size = 10), axis.title.y = element_blank(), axis.text.y=element_blank()) + 
        
        geom_vline(xintercept = bonferroni_cutoff_dip, color="red", linetype="longdash") + 
        geom_vline(data=subset(plotData, (values == maxDip & values>bonferroni_cutoff_dip)),aes(xintercept = values), color="blue", linetype="dotted") + 
        geom_vline(xintercept = 1, color="darkgrey", linetype=2) + 
        
        geom_text(data=subset(plotData, values >= label_cutoff), aes(values,label=pairs), y=0, angle = 80, hjust=0, size=5) +
        geom_text(data=subset(plotData, (values == maxDip & values>bonferroni_cutoff_dip)), x=maxDip, y=Inf, label=paste0("hat(phi)==", round(topValuesKinship[1],3),"  "),parse = TRUE, color="blue", angle = 0, size = 6, vjust = 2, hjust = 0) +
        
#         annotate("text", x=bonferroni_cutoff_dip, y=Inf, label=paste0("alpha==",format(alphaCutoff/num_comparisons_dip, digits=1)),parse = TRUE, color="red", angle = 0, size = 6, vjust = 2, hjust = 1) +
#         annotate("text", x=bonferroni_cutoff_dip, y=Inf, label=paste0("alpha "),parse = TRUE, color="red", angle = 0, size = 6, vjust = 2, hjust = 1.5) +
        annotate("text", x=xmin, y=Inf, label=ksString, color="black", angle = 0, size=6, vjust=1.5, hjust = 0) 
    
    
    
    pdf(paste0("./plots/s_distributions/",outputDir,"/",paste0(subpop,collapse="_"), plotname, ".pdf"), width=4, height=4)
    print(dipPlot)    
    dev.off()
    
}

getPopResults <- function(results, var_s="var_s_hap", s_matrix="s_matrix_hap", varianceMethod="scale"){
    dt <- as.data.table(t(sapply(names(results), function(pop_i){
        s_vector <- sort(results[[pop_i]][[s_matrix]][row(results[[pop_i]][[s_matrix]])>col(results[[pop_i]][[s_matrix]])], decreasing=T)
        topKinship <- (s_vector[1]-1)/(results[[pop_i]]$pkweightsMean-1)
        btest <- binom.test(sum(s_vector>mean(s_vector)), length(s_vector), alternative="less") 
        if(varianceMethod=="scale"){
            sdCalc <-  sd(s_vector)
        } else {
            sdCalc <-  sqrt(results[[pop_i]][[var_s]])         
        }
        structureKSTest <- ks.test((s_vector-1)/sdCalc, "pnorm", alternative = c("less"))$p.value
        crypticSig <- ifelse((s_vector[1]-1)/sdCalc > qnorm(1-.005/length(s_vector)), "YES+",
                             ifelse((s_vector[1]-1)/sdCalc > qnorm(1-.025/length(s_vector)),"YES","NO"))
        structureSig <- ifelse(structureKSTest<.01, "YES+",ifelse(structureKSTest<.05,"YES","NO"))
        c(structurePValue=btest$p.value, var_s=results[[pop_i]][[var_s]], sampleVariance=var(s_vector),
          structureKSTest=structureKSTest, closestRelatives=topKinship, crypticSig=crypticSig, structureSig=structureSig)
    })), keep.rownames=T)
    dt
}

getScoresFromSimResults <- function(simResults, diploid=F, relatedPair=NA, method="z-score"){
    s_matrix <- ifelse(diploid, "s_matrix_dip", "s_matrix_hap")
    var_s <-    ifelse(diploid, "var_s_dip", "var_s_hap")
    if (method=="z-score"){
        score <- identity
    } else if (method=="p-value"){
        score <- function(x){ 1-pnorm(x)}
    }else if (method=="neglog-p-value"){
        score <- function(x){ -log(1-pnorm(x))}
    }
    scores <- sapply(simResults, function(res){
        s_vector <- res[[s_matrix]][row(res[[s_matrix]])>col(res[[s_matrix]])]
        crypticPValue <- (s_vector-1)/sqrt(res[[var_s]])
        score(crypticPValue)
    })
    if (is.na(relatedPair)){
        scores <- apply(scores, 2, sort)       
    } else { 
        related <-scores[relatedPair,]
        scores <- rbind(related, apply(scores[-relatedPair,], 2, sort))
    }
    scores
}
