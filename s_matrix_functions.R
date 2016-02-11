calculateSMatrix <- function(subpop="CEU", filename="./data/combinedFiltered1000.gz", numberOfLines=10695, minVariants=5, alpha=.01){
    print(subpop)
    filename <- "./data/combinedFiltered1000.gz"
    con <- file(filename, "rt")
    #     system.time(genotypes <- apply(do.call(cbind, strsplit(readLines(con, numberOfLines)," ")), 1,as.numeric)[,hap.pop%in%subpop])
    system.time(genotypes <- apply(do.call(cbind, strsplit(readLines(con, numberOfLines)," ")), 1,as.numeric)[,hap.pop%in%subpop])
    hap.sampleIDs.subset <- hap.sampleIDs[hap.pop%in%subpop]
    sampleIDs.subset <- sampleIDs[pop%in%subpop]
    close(con)
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
    
    topValuesHap <- sort(s.i.j[row(s.i.j)>col(s.i.j)], decreasing=T)
   
    print(t(sapply(1:10, function(x){
        which(s.i.j==topValuesHap[x], arr.ind=T)[1,]
    })))
    topHapIndices <- sapply(1:10, function(x){
        which(s.i.j==topValuesHap[x], arr.ind=T)[1,]
    })
    
    # This needed to prevent errors when there are no significant results
    if (sum(topValuesHap[1:10]>bonferroni_cutoff)>0){
        topHapIndices <- topHapIndices[,topValuesHap[1:10]>bonferroni_cutoff]
        mappedTopHapHits <- matrix(hap.sampleIDs.subset[topHapIndices],nrow=2)
        topHitsNamesHap <- apply(mappedTopHapHits,2,paste0,collapse="_")
        topValuesHap <- topValuesHap[topValuesHap>bonferroni_cutoff]
        topValuesHap <- head(topValuesHap,10)
    } else {
        topHitsNamesHap <- NA
        topValuesHap <- 0
    }
    
    # Collapse to diploid
    s.i.j.dip <- (s.i.j[c(T,F),c(T,F)] + s.i.j[c(F,T),c(T,F)] +s.i.j[c(T,F),c(F,T)] + s.i.j[c(F,T),c(F,T)])/4
    print(mean(s.i.j.dip[row(s.i.j.dip)!=col(s.i.j.dip)]))
    print(median(s.i.j.dip[row(s.i.j.dip)!=col(s.i.j.dip)]))
    num_comparisons_dip <- choose(ncol(s.i.j.dip),2)
    
    # very lazy variance estimate...
    var_s_dip <- var_s/4
    bonferroni_cutoff_dip <- qnorm((1-alpha)^(1/num_comparisons_dip), sd=sqrt(var_s_dip)) + 1
    
    topValuesDip <- sort(s.i.j.dip[row(s.i.j.dip)>col(s.i.j.dip)], decreasing=T)
    
    topDipIndices <- sapply(1:10, function(x){
        which(s.i.j.dip==topValuesDip[x], arr.ind=T)[1,]
    })
    
    # This needed to prevent errors when there are no significant results
    if (sum(topValuesDip[1:10]>bonferroni_cutoff_dip)>0){
        topDipIndices <- topDipIndices[,topValuesDip[1:10]>bonferroni_cutoff_dip]
        mappedTopDipHits <- matrix(sampleIDs.subset[topDipIndices],nrow=2)
        topHitsNamesDip <- apply(mappedTopDipHits,2,paste0,collapse="_")
        topValuesDip <- topValuesDip[topValuesDip>bonferroni_cutoff_dip]
        topValuesDip <- head(topValuesDip,10)        
    } else {
        topHitsNamesDip <- NA
        topValuesDip <- 0
    }
    plotData <- data.frame(hap=(s.i.j[row(s.i.j)!=col(s.i.j)]))
    hapPlot <- ggplot(plotData, aes(hap, col="blue")) + geom_histogram(color="red",binwidth=.02,fill=I("blue")) + 
        ggtitle(paste0(subpop,collapse="_")) + 
        theme(plot.title = element_text(size=80), axis.title.x = element_text(size = 10))+ 
        xlab("s") + geom_vline(xintercept = bonferroni_cutoff, color="red", linetype="dotted") +
        annotate("text", x=bonferroni_cutoff -.06, y=400, label=paste0("Multiple testing cutoff, p=",format(1/num_comparisons, digits=1)), color="red",angle = 90, size = 10, hjust = 0) +
        annotate("text", x=topValuesHap, y=10, label=topHitsNamesHap,angle = 80, hjust=0, size = 10)
    
    
    tiff(paste0("./plots/s_distributions/",paste0(subpop,collapse="_"),"haploid.tiff"), width=960, height=960)
    print(hapPlot)
    dev.off()
    plotData <- data.frame(dip=s.i.j.dip[row(s.i.j.dip)!=col(s.i.j.dip)])
    dipPlot <- ggplot(plotData, aes(dip, col="blue")) + geom_histogram(color="red",binwidth=.02,fill=I("blue")) + 
        ggtitle(paste0(subpop,collapse="_"))  + 
        theme(plot.title = element_text(size=80), axis.title.x = element_text(size = 10)) + 
        xlab("s") + geom_vline(xintercept = bonferroni_cutoff_dip, color="red", linetype="dotted") +
        annotate("text", x=bonferroni_cutoff_dip -.02, y=200, label=paste0("Multiple testing cutoff, p=",format(1/num_comparisons_dip, digits=1)), color="red", angle = 90, size = 10, hjust = 0) +
        annotate("text", x=topValuesDip, y=10, label=topHitsNamesDip,angle = 80, hjust=0, size = 10)
    
    tiff(paste0("./plots/s_distributions/",paste0(subpop,collapse="_"),"diploid.tiff"), width=960, height=960)
    print(dipPlot)
    dev.off()
    
    saveRDS(s.i.j.dip, paste0("./plots/s_distributions/plotdata/",paste0(subpop,collapse="_"),"_sij.rds",collapse="_"))
    list(s_i_j = s.i.j.dip, varcovMat=varcovMat)
}