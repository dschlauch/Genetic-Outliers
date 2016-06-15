source('~/1000GP/s_matrix_functions.R')
source('~/1000GP/read1000GPsupportFiles.R')

# outputDir <- 'filtered40_LD10'
# outputDir <- 'filtered40_TGP2261_LD10'

# Generate Plots ----------------------------------------------------------

results <- readRDS(paste0("./plots/s_distributions/",outputDir,"/plotdata/all_data.rds"))

# Calculate structure significance ----------------------------------------

popResults <- getPopResults(results, "var_s_dip", "s_matrix_dip")


popResults <- popResults[order(rn)[rank(popGroup$pop)]]
popResults$super <- popGroup$group

write.table(popResults[,c("rn","super", "structureSig", "crypticSig"), with=F], 
            paste0("./plots/s_distributions/",outputDir,"/popTable.txt"),quote =F, row.names=F, sep="\t",
            col.names=c("Population","Super Population","Structure","Cryptic Relatedness"))


# Get inbreeding values ---------------------------------------------------

# inbredStats <- lapply(results, function(x){
#     a <- x$s_matrix_hap
#     b <- diag(a[seq(1,nrow(a),2),seq(2,nrow(a),2)])
#     names(b) <- colnames(a)[c(T,F)]
#     b
# })
# 
# numInbred <- lapply(inbredStats, function(x){
#     sum(1-pnorm(scale(x))<.01)
# })
# sort(unlist(inbredStats))
# hist(unlist(b))

# tidy up and plot structure results --------------------------------------

names(popResults)[1] <- "pop"
popResults$group <- popGroup[popResults$pop,"group"] 
popResults <- popResults[order(group,pop)]
popResults$pop <- factor(popResults$pop, levels=popResults$pop)
maxYvalue <- 200

ggPVals <- ggplot(popResults, aes(y=-log(as.numeric(structureKSTest)), x=pop)) +
#     geom_hline(yintercept=-log(.05)) + #geom_text(aes(13,-log(.05),label = "alpha = .05", vjust = -1),parse = T) + 
    geom_hline(yintercept=0) +
    geom_hline(yintercept=-log(.01), linetype='dashed', color="red") + #geom_text(aes(20,-log(.01),size=5,label = "alpha == .01", vjust = -1),parse = T) + 
    annotation_custom(
        grob = textGrob(label = expression(alpha==.01), hjust = 0, gp = gpar(col="red",fontsize=10, cex = 1.5)),
        ymin = -log(.01),      # Vertical position of the textGrob
        ymax = -log(.01),
        xmin = 27,         # Note: The grobs are positioned outside the plot area
        xmax = 27) +    
    theme_bw() +
    ylab("-log(p-value)") + xlab("Population") + ggtitle("Population Substructure") +
    guides(color = guide_legend(title = "Super population")) + 
    ylim(0,maxYvalue) +
    geom_point(aes(color=group),size=3) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5, size=14),
          axis.title = element_text(size=30), plot.title=element_text(size=30))
gt <- ggplot_gtable(ggplot_build(ggPVals))
gt$layout$clip[gt$layout$name == "panel"] <- "off"

pdf(paste0("./plots/s_distributions/",outputDir,"/pValueForPop.pdf"), width=8, height=8)
grid.draw(gt)
dev.off()


# ggplot(subset(popResults,!pop%in%c("PEL","MXL","ASW","PUR")), aes(y=var_s, x=sampleVariance, label=pop))+ geom_point(col="red", size=5) + geom_text(aes(y=var_s, x=sampleVariance))+ geom_abline() +
#     xlim(0,.005)+ylim(0,.005)
cor(as.numeric(subset(popResults,!pop%in%c("PEL","MXL","ASW","PUR"))$sampleVariance),as.numeric(subset(popResults,!pop%in%c("PEL","MXL","ASW","PUR"))$var_s))
varianceRatio <- as.numeric(popResults$sampleVariance)/as.numeric(popResults$var_s)
names(varianceRatio)<-popResults$pop
cbind(sort(varianceRatio, decreasing = TRUE))
sum(popResults$structurePValue<.01)
popResults$structurePValue[as.numeric(popResults$structureKSTest)<.01]
popResults <- popResults[order(popResults$structurePValue),]
structuredPops <- popResults$pop[as.numeric(popResults$structureKSTest)<.01]


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
    theme(plot.title = element_text(size=20),
          axis.text = element_text(size=14),
          axis.title = element_text(size=20),
          legend.text = element_text(size=12))
tiff(paste0("./plots/s_distributions/",outputDir,"/EstimatedCoKvsGazal.tiff"), width=900, height=900)
print(ggCoKGazal)
dev.off()
pdf(paste0("./plots/s_distributions/",outputDir,"/EstimatedCoKvsGazal.pdf"), width=8, height=8)
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

