library(ggplot2)
files <- list.files(path='~/1000GP/output_associations', pattern="*.csv")
allSNPs  <- do.call(rbind, lapply(file.path('~/1000GP/output_associations',files), read.csv, row.names=1))

allSNPs <- allSNPs[order(allSNPs$phenotype, allSNPs$MAF, allSNPs$Observed),]

# allSNPs$overallExp<-NA
# 
# sapply(levels(as.factor(allSNPs$phenotype)), function(x){
#     sapply(levels(as.factor(allSNPs$MAF)), function(y), x{
#     }     
# }))
fit <- princomp(mydata, cor=TRUE)

png(paste("~/1000GP/plots/chr", chr, "_differentialInflation.png",sep=""), width=800)
ggplot(allSNPs, aes(Expected, Observed, color=MAF))+ geom_abline(intercept = 0) + ggtitle("Phenotype Risk") + geom_point() + facet_wrap(~ phenotype)
dev.off()