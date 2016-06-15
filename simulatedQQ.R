
simResults <- readRDS(paste0("./plots/s_distributions/",outputDir,"/plotdata/Simulated_data.rds"))
simResults <- readRDS(paste0("./plots/s_distributions/",outputDir,"/plotdata/Simulated_data_cok0625.rds"))
simResults <- readRDS(paste0("./plots/s_distributions/",outputDir,"/plotdata/Simulated_0data.rds"))


simPopResults <- getPopResults(simResults)


hist(as.numeric(simPopResults$structureKSTest))
# hist(as.numeric(simPopResults$crypticPValue))

zScores <- getScoresFromSimResults(simResults, diploid=T, method="z-score")
pValues <- getScoresFromSimResults(simResults, diploid=T, method="p-value")

rowMediansZScores <- apply(zScores,1,median)
rowMediansPValues <- apply(pValues,1,median)


simQQPlot <- qplot(qnorm(1:nrow(zScores)/nrow(zScores)),rowMediansZScores) + geom_abline(intercept=0,slope=1) +
    ggtitle("QQ plot for simulated homogeneous population, 200 haplotypes, 19900 pairs") + xlab("Expected -log(p)") + ylab("Observed -log(p)") #+
#     annotate("text", x = -log(1/nrow(pValues)), y = -log(rowMedians[1]), hjust=1.1, label = "Related Pair, Phi=.0625")

pdf("./plots/simulatedQQ.pdf")
simQQPlot
dev.off()

pdf("./plots/simulatedNegLogQQ.pdf")
qplot(-log(1:nrow(pValues)/nrow(pValues)),-log(rowMediansPValues)) + geom_abline(intercept=0,slope=1) +
    ggtitle("QQ plot for simulated homogeneous population, 200 haplotypes, 19900 pairs") + xlab("Expected -log(p)") + ylab("Observed -log(p)") #+
#     annotate("text", x = -log(1/nrow(pValues)), y = -log(rowMedians[1]), hjust=1.1, label = "Related Pair, Phi=.0625")
dev.off()