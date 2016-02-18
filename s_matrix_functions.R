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
    
    num_comparisons <- numSamples*(numSamples-1)/2
    bonferroni_cutoff <- qnorm((1-alpha/(2*num_comparisons)), sd=sqrt(var_s)) + 1

    s.i.j.numerator <- t(genotypes*weights)%*%genotypes
    s.i.j.denominator <- numFilteredVariants
    s.i.j <- s.i.j.numerator/s.i.j.denominator
    
    print(mean(s.i.j[row(s.i.j)!=col(s.i.j)]))
    print(median(s.i.j[row(s.i.j)!=col(s.i.j)]))
    
    plotFromGSM(s.i.j, var_s, hap.sampleIDs.subset, "haploid")
    
    # Collapse to diploid
    s.i.j.dip <- (s.i.j[c(T,F),c(T,F)] + s.i.j[c(F,T),c(T,F)] +s.i.j[c(T,F),c(F,T)] + s.i.j[c(F,T),c(F,T)])/4
    # very lazy variance estimate...
    varS <- var_s/4
 
    plotFromGSM(s.i.j.dip, varS, sampleIDs.subset, "diploid")
    
    
    
    saveRDS(s.i.j, paste0("./plots/s_distributions/plotdata/",paste0(subpop,collapse="_"),"_sij_hap.rds"))
    saveRDS(s.i.j.dip, paste0("./plots/s_distributions/plotdata/",paste0(subpop,collapse="_"),"_sij_dip.rds"))
    list(s_i_j = s.i.j.dip, varcovMat=varcovMat)
}
plotFromGSM <- function(gsm, varS, sample_IDs, plotname=""){
    print(mean(gsm[row(gsm)!=col(gsm)]))
    print(median(gsm[row(gsm)!=col(gsm)]))
    num_comparisons_dip <- choose(ncol(gsm),2)
    
    bonferroni_cutoff_dip <- qnorm((1-alpha)^(1/num_comparisons_dip), sd=sqrt(varS)) + 1
    
    topValuesDip <- sort(gsm[row(gsm)>col(gsm)], decreasing=T)
    
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
    
    plotData <- data.frame(values=gsm[row(gsm)>col(gsm)])
    minDip <- min(plotData$values)
    maxDip <- max(plotData$values)
    dipPlot <- ggplot(plotData, aes(values)) + 
        geom_histogram(color="red",binwidth=.01,fill=I("blue")) + 
        ggtitle(paste0(subpop,collapse="_"))  + 
        theme(plot.title = element_text(size=40), axis.title.x = element_text(size = 10)) + 
        xlab("s") + geom_vline(xintercept = bonferroni_cutoff_dip, color="red", linetype="dotted") +
        annotate("text", x=bonferroni_cutoff_dip -.015, y=200, label=paste0("Multiple testing cutoff, p=",format(1/num_comparisons_dip, digits=1)), color="red", angle = 90, size = 2, hjust = 0) +
        annotate("text", x=topValuesDip, y=10, label=topHitsNamesDip,angle = 80, hjust=0, size=2)
    
    pdf(paste0("./plots/s_distributions/",paste0(subpop,collapse="_"), plotname, ".pdf"), width=4, height=4)
    print(dipPlot)
    dev.off()
}