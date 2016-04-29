
# Generate Plots ----------------------------------------------------------

results <- readRDS(paste0("./plots/s_distributions/",outputDir,"/plotdata/all_data.rds"))
sapply(names(results), function(pop_i){
    s_vector <- results[[pop_i]]$s_matrix_dip[row(results[[pop_i]]$s_matrix_dip)>col(results[[pop_i]]$s_matrix_dip)] 
    plotFromGSM(subpop=pop_i, gsm=results[[pop_i]]$s_matrix_dip, var_s=var(s_vector), pkweightsMean=results[[pop_i]]$pkweightsMean, "diploid", outputDir=outputDir)
})

# simulated
simResults <- readRDS(paste0("./plots/s_distributions/",outputDir,"/plotdata/Simulated_data.rds"))
simResults <- readRDS(paste0("./plots/s_distributions/",outputDir,"/plotdata/Simulated_data_cok0625.rds"))
# sapply(names(simResults), function(sim_i){
#     s_vector <- simResults[[sim_i]]$s_matrix_dip[row(simResults[[sim_i]]$s_matrix_dip)>col(simResults[[sim_i]]$s_matrix_dip)] 
#     c(sd(s_vector),sqrt(simResults[[sim_i]]$var_s_dip))
#     #     plotFromGSM(subpop=sim_i, gsm=simResults[[sim_i]]$s_matrix_dip, var_s=var(s_vector), pkweightsMean=simResults[[sim_i]]$pkweightsMean, "diploid", outputDir=outputDir)
})
# Calculate structure significance ----------------------------------------

getPopResults <- function(results){
    as.data.table(t(sapply(names(results), function(pop_i){
        s_vector <- sort(results[[pop_i]]$s_matrix_hap[row(results[[pop_i]]$s_matrix_hap)>col(results[[pop_i]]$s_matrix_hap)], decreasing=T)
        topKinship <- (s_vector[1]-1)/(results[[pop_i]]$pkweightsMean-1)
        btest <- binom.test(sum(s_vector>mean(s_vector)), length(s_vector), alternative="less")
        structureKSTest <- ks.test((s_vector-1)/sqrt(results[[pop_i]]$var_s_hap), "pnorm", alternative = c("less"))$p.value
        crypticSig <- ifelse((s_vector[1]-1)/sqrt(results[[pop_i]]$var_s_hap) > qnorm(1-.005/length(s_vector)), "YES+",
                                       ifelse((s_vector[1]-1)/sqrt(results[[pop_i]]$var_s_hap) > qnorm(1-.025/length(s_vector)),"YES","NO"))
        structureSig <- ifelse(structureKSTest<.01, "YES+",ifelse(structureKSTest<.05,"YES","NO"))
        c(structurePValue=btest$p.value, var_s=results[[pop_i]]$var_s_hap, sampleVariance=var(s_vector),
          structureKSTest=structureKSTest, closestRelatives=topKinship, crypticSig=crypticSig, structureSig=structureSig)
    })), keep.rownames=T)
}

popResults <- getPopResults(results)
simPopResults <- getPopResults(simResults)

hist(as.numeric(simPopResults$structureKSTest))
# hist(as.numeric(simPopResults$crypticPValue))
pValues <- sapply(simResults, function(res){
#     s_vector <- sort(res$s_matrix_hap[row(res$s_matrix_hap)>col(res$s_matrix_hap)], decreasing=T)
    s_vector <- res$s_matrix_hap[row(res$s_matrix_hap)>col(res$s_matrix_hap)]
    crypticPValue <- 1-pnorm((s_vector-1)/sqrt(res$var_s_hap))
    crypticPValue
})

# unrelated
pValues <- apply(pValues, 2, sort)

# related
related <-pValues[199,]
pValues <- rbind(related, apply(pValues[-199,], 2, sort))

rowMedians <- apply(pValues,1,mean)

pdf("./plots/simulatedQQ.pdf")
qplot(-log((1:nrow(pValues))/nrow(pValues)),-log(rowMedians)) + geom_abline(intercept=0,slope=1) +
    ggtitle("QQ plot for simulated homogeneous population, 200 haplotypes, 19900 pairs") + xlab("Expected -log(p)") + ylab("Observed -log(p)") +
    annotate("text", x = -log(1/nrow(pValues)), y = -log(rowMedians[1]), hjust=1.1, label = "Related Pair, Phi=.0625")
dev.off()

popResults <- popResults[order(rn)[rank(popGroup$pop)]]
popResults$super <- popGroup$group

write.table(popResults[,c("rn","super", "structureSig", "crypticSig"), with=F], 
            paste0("./plots/s_distributions/",outputDir,"/popTable.txt"),quote =F, row.names=F, sep="\t",
            col.names=c("Population","Super Population","Structure","Cryptic Relatedness"))

# tidy up and plot structure results --------------------------------------

names(popResults)[1] <- "pop"
popResults$group <- popGroup[popResults$pop,"group"] 
popResults <- popResults[order(group,pop)]
popResults$pop <- factor(popResults$pop, levels=popResults$pop)
maxYvalue <- 200
ggPVals <- ggplot(popResults, aes(y=-log(structureKSTest), x=pop)) +
#     geom_hline(yintercept=-log(.05)) + #geom_text(aes(13,-log(.05),label = "alpha = .05", vjust = -1),parse = T) + 
    geom_hline(yintercept=0) +
    geom_hline(yintercept=-log(.01), linetype='dashed', color="red") + #geom_text(aes(20,-log(.01),size=5,label = "alpha == .01", vjust = -1),parse = T) + 
    annotation_custom(
        grob = textGrob(label = expression(alpha==.01), hjust = 0, gp = gpar(col="red",fontsize=10, cex = 1.5)),
        ymin = -log(.01),      # Vertical position of the textGrob
        ymax = -log(.01),
        xmin = 27,         # Note: The grobs are positioned outside the plot area
        xmax = 27) +    
    ggtitle("p-value for structure in each population") + theme_bw() +
    ylab("-log(p-value)") + xlab("Population") + ggtitle("Structure Detected in 1000 Genomes Populations") +
    guides(color = guide_legend(title = "Super population")) + 
    ylim(0,maxYvalue) +
    geom_point(aes(color=group),size=3) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
gt <- ggplot_gtable(ggplot_build(ggPVals))
gt$layout$clip[gt$layout$name == "panel"] <- "off"

pdf(paste0("./plots/s_distributions/",outputDir,"/pValueForPop.pdf"), width=8, height=8)
grid.draw(gt)
dev.off()
# ggplot(subset(popResults,!pop%in%c("PEL","MXL","ASW","PUR")), aes(y=var_s, x=sampleVariance, label=pop))+ geom_point(col="red", size=5) + geom_text(aes(y=var_s, x=sampleVariance))+ geom_abline() +
#     xlim(0,.005)+ylim(0,.005)
cor(subset(popResults,!pop%in%c("PEL","MXL","ASW","PUR"))$sampleVariance,subset(popResults,!pop%in%c("PEL","MXL","ASW","PUR"))$var_s)
varianceRatio <- popResults$sampleVariance/popResults$var_s
names(varianceRatio)<-popResults$pop
cbind(sort(varianceRatio, decreasing = TRUE))
sum(popResults$structurePValue<.01)
popResults$structurePValue[popResults$structurePValue<.01]
popResults <- popResults[order(popResults$structurePValue),]
structuredPops <- popResults$pop[popResults$structurePValue<.01]


# generate summary table --------------------------------------------------

resTableAll <- data.table(do.call(rbind, lapply(unique(pop),function(pop_i){
    result <- readRDS(paste0("./plots/s_distributions/",outputDir,"/plotdata/",pop_i, "_data.rds"))
    
    result$s_matrix_dip[upper.tri(result$s_matrix_dip)] <- NA
    diag(result$s_matrix_dip) <- NA
    resTable <- melt(result$s_matrix_dip, na.rm=T)
    resTable$estimatedCoK <- (resTable$value-1)/(result$pkweightsMean-1)
    resTable$pop <- pop_i
    colnames(resTable) <- c("SampleID_1","SampleID_2","s","CoK","pop")
    resTable
})))

gazal_related_rev <- copy(gazal_related_orig)
setcolorder(gazal_related_rev, c(2,1,3,4,5))
names(gazal_related_rev)[1:2] <- names(gazal_related_rev)[2:1]
gazal_related <- rbind(gazal_related_orig,gazal_related_rev)
names(gazal_related)[1:2] <- names(resTableAll)[1:2]
resTableAll <- merge(resTableAll,gazal_related, by=names(gazal_related)[1:2],all.x=T)
resTableAll$group <- popGroup[resTableAll$pop,"group"]
sum(!is.na(resTableAll[,INFERED.RELATIONSHIP]))
resTableAll[SampleID_1=="HG00100"]
resTableAll[is.na(resTableAll[["INFERED.RELATIONSHIP"]]),"INFERED.RELATIONSHIP"]<- "Unrelated"
resTableAll$CombinedInfered <- "Unrelated"
resTableAll$CombinedInfered[resTableAll$INFERED.RELATIONSHIP=="CO"] <- "CO"
resTableAll$CombinedInfered[resTableAll$INFERED.RELATIONSHIP=="AV"|resTableAll$INFERED.RELATIONSHIP=="HS"] <- "AV or HS"
resTableAll$CombinedInfered[resTableAll$INFERED.RELATIONSHIP=="FS"|resTableAll$INFERED.RELATIONSHIP=="PO"] <- "FS or PO"
resTableAll[order(-CoK)]
# ggplot(resTableAll, aes(CombinedInfered, CoK)) +
#     geom_jitter(aes(alpha=(INFERED.RELATIONSHIP!="Unrelated")*.02)) + scale_x_discrete(limits=levels(resTableAll$INFERED.RELATIONSHIP)[c(6,2,1,4,3,5)])
# ggplot(resTableAll, aes(INFERED.RELATIONSHIP, CoK)) +
#     geom_boxplot() + scale_x_discrete(limits=levels(resTableAll$INFERED.RELATIONSHIP)[c(6,2,1,4,3,5)]) +xlab("Inferred Relationship (Gazal 2015)") + ylab("Estimated kinship")
ggCoKGazal <- ggplot(resTableAll, aes(CombinedInfered, CoK)) + theme_bw() +
    geom_jitter(aes(color=group), alpha=.5) + scale_x_discrete(limits=c("Unrelated","CO","AV or HS", "FS or PO")) + geom_violin(alpha=.4) +
    xlab("Inferred Relationship (Gazal 2015)") + ylab("Estimated kinship") + ggtitle("Estimated Kinship vs Inferred Relationship (Gazal 2015)") + 
    guides(color = guide_legend(title = "Super population")) +
    theme(plot.title = element_text(size=20))
tiff(paste0("./plots/s_distributions/",outputDir,"/EstimatedCoKvsGazal.tiff"), width=900, height=900)
print(ggCoKGazal)
dev.off()


# ggplot(subset(resTableAll, !pop%in%structuredPops), aes(CombinedInfered, CoK)) + theme_bw() +
#     geom_boxplot() + geom_jitter(aes(color=group), alpha=.5) + scale_x_discrete(limits=c("Unrelated","CO","AV or HS", "FS or PO")) +
#     xlab("Inferred Relationship (Gazal 2015)") + ylab("Estimated kinship") + ggtitle("Our estimated CoK vs Inferred Relationship (Gazal 2015)")
ggCoKGazalUnstructured <- ggplot(subset(resTableAll, !pop%in%structuredPops), aes(CombinedInfered, CoK)) + theme_bw() +
    geom_jitter(aes(color=group), alpha=.5) + scale_x_discrete(limits=c("Unrelated","CO","AV or HS", "FS or PO")) + geom_violin(alpha=.4) +
    xlab("Inferred Relationship (Gazal 2015)") + ylab("Estimated kinship") + ggtitle("Estimated Kinship vs Inferred Relationship (Gazal 2015)") +
    theme(plot.title = element_text(size=20)) +
    guides(color = guide_legend(title = "Super population"))
tiff(paste0("./plots/s_distributions/",outputDir,"/EstimatedCoKvsGazalUnstruct.tiff"), width=900, height=900)
print(ggCoKGazalUnstructured)
dev.off()
