calculateSMatrix <- function(subpop="CEU", filename="./data/combinedFiltered1000.gz", numberOfLines=40695, minVariants=5, qcFilter=NULL, ldPrune=1){
    
    print("Starting read file")
    system.time(genotypes <- fread(paste('zcat',filename), sep=" ", nrows=numberOfLines, header=F))
    print("Finished read file")
    
    names(genotypes) <- hap.sampleIDs
    
    # lapply runs into memory issues for large datasets, use for loop. *cringe*
    results <- list()
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
        results[[subpop]] <- generateSResultsFromGenotypes(subpop, genotypes[,filterhap, with=F], qcFilter, minVariants, ldPrune)
    }
    names(results) <- unique(pop)
    results    
}

generateSResultsFromGenotypes <- function(subpop, genotypesSubpop, qcFilter, minVariants, ldPrune=1){
    
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

    # remove < n variants
    sumVariants <- rowSums(genotypesSubpop)
    genotypesSubpop <- genotypesSubpop[sumVariants>minVariants,]
    genotypesSubpop <- as.matrix(genotypesSubpop)

    print("Number of used variants")
    print(nrow(genotypesSubpop))
    numFilteredVariants <- nrow(genotypesSubpop)
    sumFilteredVariants <- rowSums(genotypesSubpop)
    varcovMat <- cov(genotypesSubpop[,c(T,F)] + genotypesSubpop[,c(F,T)])
    
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
    saveRDS(popResult, paste0("./plots/s_distributions/",outputDir,"/plotdata/",subpop,"_data.rds"))
    popResult

}
plotFromGSM <- function(subpop, gsm, var_s, pkweightsMean, plotname="", outputDir=".",alphaCutoff=.01){
    print(mean(gsm[row(gsm)!=col(gsm)]))
    print(median(gsm[row(gsm)!=col(gsm)]))
    num_comparisons_dip <- choose(ncol(gsm),2)
    sample_IDs <- rownames(gsm)
    bonferroni_cutoff_dip <- qnorm((1-alphaCutoff)^(1/num_comparisons_dip), sd=sqrt(var_s)) + 1
    
    topValuesDip <- sort(gsm[row(gsm)>col(gsm)], decreasing=T)
    topValuesKinship <- (topValuesDip-1)/(pkweightsMean-1)
    
    # Display only those that are above the cutoff and among the top 5
    label_cutoff <- max(bonferroni_cutoff_dip, topValuesDip[5])
 
    pairs <- outer(sample_IDs, sample_IDs, paste)
    plotData <- data.frame(values=gsm[row(gsm)>col(gsm)], pairs=paste0("  ",pairs[row(pairs)>col(pairs)]))
    minDip <- min(plotData$values)
    maxDip <- max(plotData$values)#ifelse(max(plotData$values)>bonferroni_cutoff_dip,max(plotData$values),NA)
    dipPlot <- ggplot(plotData, aes(values)) + 
        geom_histogram(color="blue",binwidth=.01,fill=I("blue")) + 
        ggtitle(paste0(subpop,collapse="_"))  + xlab("Similarity score") + 
        scale_x_continuous(expand=c(.2,0))+#limits=c(minDip, maxDip+1)) +
        theme_bw() +
        theme(plot.title = element_text(size=40), axis.title.x = element_text(size = 10), axis.title.y = element_blank(), axis.text.y=element_blank()) + 
        
        geom_vline(xintercept = bonferroni_cutoff_dip, color="red", linetype="dotted") + 
        geom_vline(data=subset(plotData, (values == maxDip & values>bonferroni_cutoff_dip)),aes(xintercept = values), color="blue", linetype="dotted") + 
        geom_vline(xintercept = median(gsm[row(gsm)!=col(gsm)]), color="black", linetype=1) + 
        
        geom_text(data=subset(plotData, values > label_cutoff), aes(values,label=pairs), y=0, angle = 80, hjust=0, size=3) +
        geom_text(data=subset(plotData, (values == maxDip & values>bonferroni_cutoff_dip)), x=maxDip, y=Inf, label=paste0("hat(phi)==", round(topValuesKinship[1],4),"  "),parse = TRUE, color="blue", angle = 0, size = 3, vjust = 1, hjust = 0) +
        
        annotate("text", x=bonferroni_cutoff_dip, y=Inf, label=paste0("alpha==",format(alphaCutoff/num_comparisons_dip, digits=1)),parse = TRUE, color="red", angle = 0, size = 3, vjust = 1, hjust = 1) +
        annotate("text", x=median(gsm[row(gsm)!=col(gsm)]), y=Inf, label=paste0("m=",round(median(gsm[row(gsm)!=col(gsm)]),3)," "), color="black", angle = 0, size = 3, vjust = 1, hjust = 1) 
    
    
    
    pdf(paste0("./plots/s_distributions/",outputDir,"/",paste0(subpop,collapse="_"), plotname, ".pdf"), width=4, height=4)
    print(dipPlot)    
    dev.off()
    
}
