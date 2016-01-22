calculateSMatrix <- function(subpop="CEU", numberOfLines=10695, minVariants=5, alpha=.05){
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
    numFilteredVariants <- nrow(genotypes)
    sumFilteredVariants <- rowSums(genotypes)
    
    totalPossiblePairs <- choose(numSamples,2)
    totalPairs <- choose(sumFilteredVariants,2)
    weights <- totalPossiblePairs/totalPairs
    p <- 1/weights
    var_s <- sum((1-p)/p)/(numFilteredVariants^2)
    print("variance of s")
    print(var_s)
    
    alpha <- .05
    num_comparisons <- numSamples*(numSamples-1)/2
    bonferroni_cutoff <- qnorm((1-alpha)^(1/num_comparisons), sd=sqrt(var_s)) + 1

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
    topHapIndices <- topHapIndices[,topValuesHap[1:10]>bonferroni_cutoff]
    mappedTopHapHits <- matrix(hap.sampleIDs.subset[topHapIndices],nrow=2)
    topValuesHap <- topValuesHap[topValuesHap>bonferroni_cutoff]
    topValuesHap <- head(topValuesHap,10)
    
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
    topDipIndices <- topDipIndices[,topValuesDip[1:10]>bonferroni_cutoff_dip]
    mappedTopDipHits <- matrix(sampleIDs.subset[topDipIndices],nrow=2)
    topValuesDip <- topValuesDip[topValuesDip>bonferroni_cutoff_dip]
    topValuesDip <- head(topValuesDip,10)
    plotData <- data.frame(hap=(s.i.j[row(s.i.j)!=col(s.i.j)]))
    hapPlot <- ggplot(plotData, aes(hap)) + geom_histogram(color="blue",binwidth=.02) + 
        ggtitle(paste0("Distribution of haploid s, population: ",paste0(subpop,collapse="_"))) + xlab("s") +
        xlab("s") + geom_vline(xintercept = bonferroni_cutoff, color="red") +
        annotate("text", x=bonferroni_cutoff -.06, y=400, label=paste0("Multiple testing cutoff, p=",format(1/num_comparisons, digits=1)), color="red",angle = 90) +
        annotate("text", x=topValuesHap, y=10, label=apply(mappedTopHapHits,2,paste0,collapse="_"),angle = 80, hjust=0)
    
    
    tiff(paste0("./plots/s_distributions/",paste0(subpop,collapse="_"),"haploid.tiff"))
    print(hapPlot)
    dev.off()
    plotData <- data.frame(dip=s.i.j.dip[row(s.i.j.dip)!=col(s.i.j.dip)])
    dipPlot <- ggplot(plotData, aes(dip)) + geom_histogram(color="blue",binwidth=.02) + 
        ggtitle(paste0("Distribution of diploid s, population: ",paste0(subpop,collapse="_"))) + 
        xlab("s") + geom_vline(xintercept = bonferroni_cutoff_dip, color="red") +
        annotate("text", x=bonferroni_cutoff_dip -.02, y=200, label=paste0("Multiple testing cutoff, p=",format(1/num_comparisons_dip, digits=1)), color="red",angle = 90)+
        annotate("text", x=topValuesDip, y=10, label=apply(mappedTopDipHits,2,paste0,collapse="_"),angle = 45, hjust=0)
    
    tiff(paste0("./plots/s_distributions/",paste0(subpop,collapse="_"),"diploid.tiff"))
    print(dipPlot)
    dev.off()
    
    saveRDS(s.i.j, paste0("./plots/s_distributions/plotdata/",paste0(subpop,collapse="_"),"_sij.rds",collapse="_"))
    
}