library(data.table)

speedTest <- function(genotypesSubpop){
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
}

numSamples <- seq(500,3000,100)
replicateIndex <- 2
numVariants <- 100000
timed_methods <- c(speedTest, princomp, cor)
method_names <- c("stego", "princomp")

sapply(numSamples, function(x){
    genotypes <- matrix(rbinom(x*numVariants,1,.2), ncol=x)
    sapply(seq_along(method_names), function(meth){
        timing <- c(system.time(timed_methods[meth][[1]](genotypes)), nsamp=x, replicateIndex=replicateIndex, method=method_names[meth])
        saveRDS(timing, paste0("./timing/",method_names[meth],"_",x,"_",replicateIndex, "_.rds"))
        cat(paste0(x, " "))
    })
})


# # 1 time use
# apply(times[replicateIndex==1],1,function(x){
#     saveRDS(x[c(3,4,5,6,7,2,8,9)], paste0("./timing/",x[1]))
#     }
# )


times <- as.data.table(t(sapply(list.files("./timing/"), function(x) readRDS(paste0("./timing/",x)))),keep.rownames=TRUE)
times$elapsed <- as.numeric(times$elapsed)
times$nsamp <- as.numeric(times$nsamp)

# collapse the times
times[,lapply(.SD, mean),by=c(method, nsamp)]

library(ggplot2)

pdf("timing_comparison.pdf")
ggplot(times[replicateIndex==1]) + geom_point(aes(x=nsamp/2, y=elapsed, col=method)) + theme_bw() + xlab("Number of Samples") +ylab("Timing (seconds)")
dev.off()
