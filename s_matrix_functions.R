calculateSMatrix <- function(subpop="CEU", filename="./data/combinedFiltered1000.gz", numberOfLines=40695, minVariants=5, alpha=.01){
    require(readr)
    print(subpop)
    
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
    
    # Fully simulated ---------------------------------------------------------
    numSamples <- 100
    numVariants<- 100000
    genotypes <- matrix(rbinom(numSamples*2*numVariants,1, .1), ncol=numSamples*2)
    numSamples <- ncol(genotypes)
    sumVariants <- rowSums(genotypes)
    genotypes <- genotypes[sumVariants>5,]
    subpop <- "Simulated"
    sum(rowSums(genotypes)<2)    
#     ################### Toggle this.  Adds a related individual for testing ################
    coefR <- .1
    probs <- rep((1-2*coefR)/ncol(genotypes),ncol(genotypes))
    probs[1:2] <- probs[1:2]+coefR
    indices <- replicate(nrow(genotypes), {sample(ncol(genotypes),2 , prob=probs, replace=F)})
#     indices1 <- sample(ncol(genotypes), nrow(genotypes), prob=probs, replace=T)
#     indices2 <- sample(ncol(genotypes), nrow(genotypes), prob=probs, replace=T)
    relatedHap1 <- mapply(function(x,y){genotypes[x,y]},1:nrow(genotypes),indices[1,])
    relatedHap2 <- mapply(function(x,y){genotypes[x,y]},1:nrow(genotypes),indices[2,])
    genotypes <- cbind(genotypes, relatedHap1, relatedHap2)
#     ################################################



    print("Number of used variants")
    print(nrow(genotypes))
    numFilteredVariants <- nrow(genotypes)
    sumFilteredVariants <- rowSums(genotypes)
    varcovMat <- cov(genotypes[,c(T,F)] + genotypes[,c(F,T)])
    
    totalPossiblePairs <- choose(numSamples,2)
    totalPairs <- choose(sumFilteredVariants,2)
    weights <- totalPossiblePairs/totalPairs
    p <- 1/weights
    var_s_hap <- sum((1-p)/p)/(numFilteredVariants^2)
    print("variance of s (haploid)")
    print(var_s_hap)
    
    # Calculate expected values conditional on kinship
    pkweightsMean <- mean(((sumFilteredVariants-2)/numSamples)*weights)
    kinships <- seq(0,.25,.001)
    kinshipExpectation <- 1+kinships*(pkweightsMean-1)
        
    s_matrix_numerator <- t(genotypes*weights)%*%genotypes
    s_matrix_denominator <- numFilteredVariants
    s_matrix_hap <- s_matrix_numerator/s_matrix_denominator
    
    print(mean(s_matrix_hap[row(s_matrix_hap)!=col(s_matrix_hap)]))
    print(median(s_matrix_hap[row(s_matrix_hap)!=col(s_matrix_hap)]))
    
    estimatedKinship <- (s_matrix_hap-1)/(pkweightsMean-1)
    
    # Collapse to diploid
    s_matrix_dip <- (s_matrix_hap[c(T,F),c(T,F)] + s_matrix_hap[c(F,T),c(T,F)] +s_matrix_hap[c(T,F),c(F,T)] + s_matrix_hap[c(F,T),c(F,T)])/4
    # very lazy variance estimate...
    var_s_dip <- var_s_hap/4
    

    list(s_matrix_dip=s_matrix_dip, s_matrix_hap=s_matrix_hap, pkweightsMean=pkweightsMean, var_s_dip=var_s_dip, var_s_hap=var_s_hap, varcovMat=varcovMat)
    

# For sim -----------------------------------------------------------------

var_s <- var_s_dip
gsm <- s_matrix_dip
sample_IDs <- paste0("Sample",1:ncol(gsm))
sample_IDs[length(sample_IDs)] <- "Related"

var_s <- var_s_hap
gsm <- s_matrix_hap
sample_IDs <- paste0("Sample",1:ncol(gsm))
sample_IDs[c(length(sample_IDs),length(sample_IDs)-1)] <- "Related"

alpha=.01

}
plotFromGSM <- function(subpop, gsm, var_s, pkweightsMean, sample_IDs, plotname="", outputDir=".",alpha=.01){
    print(mean(gsm[row(gsm)!=col(gsm)]))
    print(median(gsm[row(gsm)!=col(gsm)]))
    num_comparisons_dip <- choose(ncol(gsm),2)
    
    bonferroni_cutoff_dip <- qnorm((1-alpha)^(1/num_comparisons_dip), sd=sqrt(var_s)) + 1
    
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
        
        annotate("text", x=bonferroni_cutoff_dip, y=Inf, label=paste0("alpha==",format(alpha/num_comparisons_dip, digits=1)),parse = TRUE, color="red", angle = 0, size = 3, vjust = 1, hjust = 1) +
        annotate("text", x=median(gsm[row(gsm)!=col(gsm)]), y=Inf, label=paste0("m=",round(median(gsm[row(gsm)!=col(gsm)]),3)," "), color="black", angle = 0, size = 3, vjust = 1, hjust = 1) 
    
    
    
    pdf(paste0("./plots/s_distributions/",outputDir,"/",paste0(subpop,collapse="_"), plotname, ".pdf"), width=4, height=4)
    print(dipPlot)
    dev.off()
}
