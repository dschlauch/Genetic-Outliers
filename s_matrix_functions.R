calculateSMatrix <- function(subpop="CEU", filename="./data/combinedFiltered1000.gz", numberOfLines=40695, minVariants=5, alpha=.01){
    print(subpop)
    filename <- "./data/combinedFiltered1000.gz"
    
    con <- file(filename, "rt")
    system.time(genotypes <- apply(do.call(cbind, strsplit(readLines(con, numberOfLines)," ")), 1,as.numeric)[,hap.pop%in%subpop])
    close(con)
    
    hap.sampleIDs.subset <- hap.sampleIDs[hap.pop%in%subpop]
    sampleIDs.subset <- sampleIDs[pop%in%subpop]
    numSamples <- ncol(genotypes)
    sumVariants <- rowSums(genotypes)
    # reverse so that MAF<.5
    genotypes[sumVariants>(numSamples/2),] <- 1-genotypes[sumVariants>(numSamples/2),]
    sumVariants <- rowSums(genotypes)
    # remove < n variants
    genotypes <- genotypes[sumVariants>minVariants,]
    print("Number of used variants")
    print(nrow(genotypes))
    numFilteredVariants <- nrow(genotypes)
    sumFilteredVariants <- rowSums(genotypes)
    varcovMat <- cov(genotypes[,c(T,F)] + genotypes[,c(F,T)])
    
    totalPossiblePairs <- choose(numSamples,2)
    totalPairs <- choose(sumFilteredVariants,2)
    weights <- totalPossiblePairs/totalPairs
    p <- 1/weights
    var_s <- sum((1-p)/p)/(numFilteredVariants^2)
    print("variance of s")
    print(var_s)
    
    # Calculate expected values conditional on kinship
    meanOfWeights <- mean(weights)
    kinships <- seq(0,.25,.001)
    kinshipExpectation <- 1+kinships*(meanOfWeights-1)
        
    ggplot(data.frame(kinships, kinshipExpectation), aes(x=kinships, y=kinshipExpectation))+ geom_point()
    #     num_comparisons <- numSamples*(numSamples-1)/2
    #     bonferroni_cutoff <- qnorm((1-alpha/(2*num_comparisons)), sd=sqrt(var_s)) + 1
    
    s.i.j.numerator <- t(genotypes*weights)%*%genotypes
    s.i.j.denominator <- numFilteredVariants
    s.i.j <- s.i.j.numerator/s.i.j.denominator
    
    print(mean(s.i.j[row(s.i.j)!=col(s.i.j)]))
    print(median(s.i.j[row(s.i.j)!=col(s.i.j)]))
    
    estimatedKinship <- (s.i.j-1)/(meanOfWeights-1)
    plotFromGSM(subpop, s.i.j, var_s, estimatedKinship, hap.sampleIDs.subset, "haploid")
    
    # Collapse to diploid
    s.i.j.dip <- (s.i.j[c(T,F),c(T,F)] + s.i.j[c(F,T),c(T,F)] +s.i.j[c(T,F),c(F,T)] + s.i.j[c(F,T),c(F,T)])/4
    # very lazy variance estimate...
    varS <- var_s/4
    
    estimatedKinshipDip <- (s.i.j.dip-1)/(meanOfWeights-1)
    plotFromGSM(subpop,s.i.j.dip, varS, estimatedKinshipDip, sampleIDs.subset, "diploid")
    
    
    
    saveRDS(s.i.j, paste0("./plots/s_distributions/plotdata/",paste0(subpop,collapse="_"),"_sij_hap.rds"))
    saveRDS(s.i.j.dip, paste0("./plots/s_distributions/plotdata/",paste0(subpop,collapse="_"),"_sij_dip.rds"))
    saveRDS(varcovMat, paste0("./plots/s_distributions/plotdata/",paste0(subpop,collapse="_"),"_varcovMat.rds"))
    list(s_i_j = s.i.j.dip, varcovMat=varcovMat)
}
plotFromGSM <- function(subpop, gsm, varS, kinship, sample_IDs, plotname="", alpha=.01){
    print(mean(gsm[row(gsm)!=col(gsm)]))
    print(median(gsm[row(gsm)!=col(gsm)]))
    num_comparisons_dip <- choose(ncol(gsm),2)
    
    bonferroni_cutoff_dip <- qnorm((1-alpha)^(1/num_comparisons_dip), sd=sqrt(varS)) + 1
    
    topValuesDip <- sort(gsm[row(gsm)>col(gsm)], decreasing=T)
    topValuesKinship <- sort(kinship[row(kinship)>col(kinship)], decreasing=T)
    
    topDipIndices <- sapply(1:10, function(x){
        which(gsm==topValuesDip[x], arr.ind=T)[1,]
    })
    
    # This needed to prevent errors when there are no significant results
    if (sum(topValuesDip[1:10]>bonferroni_cutoff_dip)>0){
        topDipIndices <- topDipIndices[,topValuesDip[1:10]>bonferroni_cutoff_dip]
        mappedTopDipHits <- matrix(sample_IDs[topDipIndices],nrow=2)
        topHitsNamesDip <- apply(mappedTopDipHits,2,paste0,collapse="_")
        topValuesDip <- topValuesDip[topValuesDip>bonferroni_cutoff_dip]
        topValuesDip <- head(topValuesDip,10)
    } else {
        topHitsNamesDip <- NA
        topValuesDip <- 1
    }
    pairs <- outer(sample_IDs, sample_IDs, paste)
    plotData <- data.frame(values=gsm[row(gsm)>col(gsm)], pairs=paste0("  ",pairs[row(pairs)>col(pairs)]))
    minDip <- min(plotData$values)
    maxDip <- max(plotData$values)#ifelse(max(plotData$values)>bonferroni_cutoff_dip,max(plotData$values),NA)
    dipPlot <- ggplot(plotData, aes(values)) + 
        geom_histogram(color="blue",binwidth=.01,fill=I("blue")) + 
        ggtitle(paste0(subpop,collapse="_"))  + xlab("Similarity score") + 
        scale_x_continuous(limits=c(0, maxDip+1)) +
        theme_bw() +
        theme(plot.title = element_text(size=40), axis.title.x = element_text(size = 10), axis.title.y = element_blank(), axis.text.y=element_blank()) + 
        
        geom_vline(xintercept = bonferroni_cutoff_dip, color="red", linetype="dotted") + 
        geom_vline(data=subset(plotData, (values == maxDip & values>bonferroni_cutoff_dip)),aes(xintercept = values), color="blue", linetype="dotted") + 
        geom_vline(xintercept = median(gsm[row(gsm)!=col(gsm)]), color="red", linetype=1) + 
        
        geom_text(data=subset(plotData, values > bonferroni_cutoff_dip), aes(values,label=pairs), y=0, angle = 80, hjust=0, size=3) +
        geom_text(data=subset(plotData, (values == maxDip & values>bonferroni_cutoff_dip)), x=maxDip, y=Inf, label=paste0("hat(phi)==", round(topValuesKinship[1],3),"  "),parse = TRUE, color="blue", angle = 0, size = 3, vjust = 1, hjust = 0) +
        
        annotate("text", x=bonferroni_cutoff_dip, y=Inf, label=paste0("alpha==",format(alpha/num_comparisons_dip, digits=1)),parse = TRUE, color="red", angle = 0, size = 3, vjust = 1, hjust = 1) +
        annotate("text", x=median(gsm[row(gsm)!=col(gsm)]), y=Inf, label=paste0("m=",round(median(gsm[row(gsm)!=col(gsm)]),3)," "), color="red", angle = 0, size = 3, vjust = 1, hjust = 0) 
    
    
    
    pdf(paste0("./plots/s_distributions/",paste0(subpop,collapse="_"), plotname, ".pdf"), width=4, height=4)
    print(dipPlot)
    dev.off()
}