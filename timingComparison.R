library(data.table)

##  A version of stego with all the pre-processing, data 
##  validation, post-analysis statistical tests steps removed. 
##  Additionally added the eigen 
##  step for comparison with princomp().
stego <- function(genotypesSubpop){
    numSamples <- ncol(genotypesSubpop)
    genotypesSubpop <- as.matrix(genotypesSubpop)
    
    numFilteredVariants <- nrow(genotypesSubpop)
    sumFilteredVariants <- rowSums(genotypesSubpop)

    totalPossiblePairs <- choose(numSamples,2)
    totalPairs <- choose(sumFilteredVariants,2)
    weights <- totalPossiblePairs/totalPairs
    p <- 1/weights
    
    var_s_hap <- sum(weights-1)/(numFilteredVariants^2)
    pkweightsMean <- mean(((sumFilteredVariants-2)/numSamples)*weights)
    
    s_matrix_numerator <- t(genotypesSubpop*weights)%*%genotypesSubpop
    s_matrix_denominator <- numFilteredVariants
    s_matrix_hap <- s_matrix_numerator/s_matrix_denominator
    eigen(s_matrix_hap)
}

numSamples <- seq(100,3600,100)
replicateIndex <- 32
numVariants <- 100000
timed_methods <- c(princomp, prcomp, stego)
method_names <- c("PCA","prcomp","STEGO")

sapply(numSamples, function(x){
    genotypes <- matrix(rbinom(x*numVariants,1,.2), ncol=x)
    sapply(seq_along(method_names), function(meth){
        timing <- c(system.time(timed_methods[meth][[1]](genotypes)), nsamp=x, replicateIndex=replicateIndex, method=method_names[meth])
        saveRDS(timing, paste0("./timing/",method_names[meth],"_",x,"_",replicateIndex, "_.rds"))
        cat(paste0(x, " "))
    })
})


# # # 1 time use
# apply(times,1,function(x){
#     saveRDS(x[2:9], paste0("./timing/",x[1]))
# })


library(ggplot2)
library(data.table)
# Plot all the times
times <- as.data.table(t(sapply(list.files("./timing/"), function(x) readRDS(paste0("./timing/",x)))),keep.rownames=TRUE)
times$elapsed <- as.numeric(times$elapsed)
times$nsamp <- as.numeric(times$nsamp)
times$replicateIndex <- as.numeric(times$replicateIndex)

times <- times[nsamp!=3700]
# times <- times[replicateIndex>=20]
times <- times[method!="STEGO_FULL"]

pdf("timing_comparison.pdf")
ggplot(times) + geom_point(aes(x=nsamp/2, y=elapsed, col=method), alpha=.5) + theme_bw() + xlab("Number of Samples") +ylab("Timing (seconds)")
dev.off()

# collapse the time replicates
average_times <- times[, mean(elapsed),by=.(method,nsamp)]
names(average_times)[3] <- "elapsed"

pdf("average_timing_comparison.pdf")
print(ggplot(average_times) + geom_point(aes(x=nsamp/2, y=elapsed, col=method)) + 
    theme_bw() + xlab("Number of Samples") +ylab("Average Timing (seconds)") +
    ggtitle("Computation time for 100k variants vs X samples in R") + 
    theme(plot.title = element_text(hjust = 0.5))+ 
    scale_colour_discrete(name="Method",
                          breaks=c("cor", "PCA", "prcomp","STEGO"),
                          labels=c("cor()", "princomp()", "prcomp()","stego()")))
dev.off()
