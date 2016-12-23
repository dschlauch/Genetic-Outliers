library(data.table)
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(cluster)
library(dendextend)

setwd("~/1000GP/")

sample <- read.table("./data/1000GP_Phase3.sample", sep=" ", header=T)

# Read in the U's and V's that were gathered as cumulative
cumulativeFiles <- list.files(pattern="output_20_")
V.list       <- lapply(file.path(cumulativeFiles, "all_V.csv"), read.csv, row.names=1)
U.list       <- lapply(file.path(cumulativeFiles, "all_U.csv"), read.csv, row.names=1)
V.list       <- lapply(V.list, as.matrix)
U.list       <- lapply(U.list, as.matrix)

alleleMax <- as.numeric(sub("output_20_","",cumulativeFiles))

names(V.list) <- alleleMax
names(U.list) <- alleleMax

V.list <- V.list[as.character(sort(alleleMax))]
U.list <- U.list[as.character(sort(alleleMax))]

U.interval.list <- lapply(2:length(U.list), function(index){
    U.list[[index]]-U.list[[index-1]]
})

V.interval.list <- lapply(2:length(V.list), function(index){
    V.list[[index]]-V.list[[index-1]]
})

names(U.interval.list) <- paste(names(V.list)[-length(V.list)],names(V.list)[-1],sep='-')
names(V.interval.list) <- paste(names(V.list)[-length(V.list)],names(V.list)[-1],sep='-')

U.interval.list[["20-50"]] <- U.list[["50"]]
V.interval.list[["20-50"]] <- V.list[["50"]]


# Read in the U's and V's that were gathered in an interval
interval.file.list <- list.files(pattern="output")[!(grepl("_20_",list.files(pattern="output"))|grepl("ass",list.files(pattern="output")))]

# Too much for my poor laptop to handle
interval.file.list <- interval.file.list[-c(2,10)]

V.list2       <- lapply(file.path(interval.file.list, "all_V.csv"), read.csv, row.names=1)
U.list2       <- lapply(file.path(interval.file.list, "all_U.csv"), read.csv, row.names=1)
V.list2       <- lapply(V.list2, as.matrix)
U.list2       <- lapply(U.list2, as.matrix)

allele.interval <- sub("_","-",sub("output_","",interval.file.list))

names(V.list2) <- allele.interval
names(U.list2) <- allele.interval

# Combine the two lists
U.interval.list <- c(U.interval.list, U.list2)
V.interval.list <- c(V.interval.list, V.list2)

print("Removing individual U and V matrices")
# Clean up
rm(V.list)
rm(U.list)
rm(V.list2)
rm(U.list2)
gc()
#reorder lists
listorder <- order(as.numeric(unlist(strsplit(names(U.interval.list),"-"))[c(T,F)]))

U.interval.list <- U.interval.list[listorder]
V.interval.list <- V.interval.list[listorder]

totalMA <- unlist(lapply(U.interval.list, function(x){sum(diag(x))}))

png("./plots/alleles_per_interval.png", width=800)
qplot(y=totalMA,x=factor(names(totalMA), levels=names(totalMA)), xlab="Allele Interval", ylab="Total alleles") + labs(title = "Total minor alleles per interval")
dev.off()

combinedJaccard <- function(jaccardMat){
   (jaccardMat[seq(1,nrow(jaccardMat),2),seq(1,nrow(jaccardMat),2)] + 
    jaccardMat[seq(2,nrow(jaccardMat),2),seq(1,nrow(jaccardMat),2)] +
    jaccardMat[seq(1,nrow(jaccardMat),2),seq(2,nrow(jaccardMat),2)] + 
    jaccardMat[seq(2,nrow(jaccardMat),2),seq(2,nrow(jaccardMat),2)])/4
}

U_tmp <- matrix(0, nrow(U.interval.list[[1]]), ncol(U.interval.list[[1]]))
V_tmp <- matrix(0, nrow(U.interval.list[[1]]), ncol(U.interval.list[[1]]))
jaccard.cumulative.list <- list()
jaccard.interval.list <- list()

# Must be run in for loop (Do not change to apply)
for (i in names(U.interval.list)[-1:-2]){
    print(i)
    U_tmp <- U_tmp + U.interval.list[[i]]
    V_tmp <- V_tmp + V.interval.list[[i]]
    jaccard.cumulative.list[[i]] <- combinedJaccard(V_tmp/U_tmp) 
    jaccard.interval.list[[i]] <- combinedJaccard(V.interval.list[[i]]/U.interval.list[[i]])
}

# Cleanup
rm(U.interval.list)
rm(V.interval.list)
gc()

getJaccardRatio <- function(jaccard, pop.1.name, pop.2.name, sample){
    pop.1.Samples <- sample$POP==pop.1.name
    pop.2.Samples <- sample$POP==pop.2.name
    
    within.1 <- jaccard[pop.1.Samples, pop.1.Samples]
    within.2 <- jaccard[pop.2.Samples, pop.2.Samples]
    between  <- jaccard[pop.2.Samples, pop.1.Samples]
    
    
    within <- (sum(within.1[row(within.1)!=col(within.1)])+sum(within.2[row(within.2)!=col(within.2)])) / 
        ((nrow(within.1)*(nrow(within.1)-1))+(nrow(within.2)*(nrow(within.2)-1)))
    within/mean(as.matrix(between))
}
populationPartition <- function(jaccard, popName1, popName2){
    subpop <-sample[,2]==popName1|sample[,2]==popName2
    heatmapData <- jaccard
    colnames(heatmapData) <- sample[,2]
    heatmapData <- heatmapData[subpop,subpop]
    diag(heatmapData) <- 0
    #my_palette <- colorRampPalette(c("green", "black", "red"))
    #heatmap.2(as.matrix(heatmapData), trace="none",col=bluered,breaks=c(0,0.0071,.008,.026),symbreaks=T)
    
    clusters <- pam(heatmapData,2)$clustering
    AintersectB1 <- sum(clusters[colnames(heatmapData)==popName1]==1)
    AunionB1 <- sum(clusters[colnames(heatmapData)==popName2]==1)+
        length(clusters[colnames(heatmapData)==popName1])
    AintersectB2 <- sum(clusters[colnames(heatmapData)==popName1]==2)
    AunionB2 <- sum(clusters[colnames(heatmapData)==popName2]==2)+
        length(clusters[colnames(heatmapData)==popName1])
    
    max(AintersectB1/AunionB1, AintersectB2/AunionB2)
}

populations <- unique(sample[,2:3])[order(unique(sample[,2])),]
ggData <- data.frame(ratio=numeric(0), clustering=numeric(0), pop=character(0), super1=character(0), super2=character(0), cumulative=character(0))
for(i in 1:(nrow(populations)-1)){
    for (j in (i+1):nrow(populations)){
        cat("!")
        ggData <- rbind(ggData, 
                        data.frame(ratio=unlist(lapply(jaccard.interval.list, getJaccardRatio, populations[i,1], populations[j,1], sample)), 
                                    clustering=unlist(lapply(jaccard.interval.list, populationPartition, populations[i,1], populations[j,1])),
                                    pop=paste(populations[i,1], populations[j,1], sep=" vs "),
                                    super1=populations[i,2],
                                    super2=populations[j,2],
                                    cumulative='interval')
        )
        ggData <- rbind(ggData, 
                        data.frame(ratio=unlist(lapply(jaccard.cumulative.list, getJaccardRatio, populations[i,1], populations[j,1], sample)), 
                                   clustering=unlist(lapply(jaccard.cumulative.list, populationPartition, populations[i,1], populations[j,1])),
                                   pop=paste(populations[i,1], populations[j,1], sep=" vs "),
                                   super1=populations[i,2],
                                   super2=populations[j,2],
                                   cumulative='cumulative')
        )
    }
}
ggData$interval <-factor(names(jaccard.interval.list), levels = names(jaccard.interval.list))
ggData$minAC <-unlist(strsplit(names(jaccard.interval.list),"-"))[seq(1,2*length(names(jaccard.interval.list)),2)]
ggData$maxAC <-unlist(strsplit(names(jaccard.interval.list),"-"))[seq(2,2*length(names(jaccard.interval.list)),2)]
interval_percent <- paste(as.numeric(ggData$minAC)/50,"%-", as.numeric(ggData$maxAC)/50,"%", sep="")
ggData$interval_percent <- factor(interval_percent, levels=unique(interval_percent))

plotPopComparison <- function(comparison, data){
    ggplot(data=subset(ggData,grepl(comparison, data$pop)&data$cumulative=="interval"), aes(x=interval_percent, y=ratio, group=1)) +
        geom_line(colour="red") + 
        geom_point(colour="red")  +
        ggtitle(paste("Jaccard Ratio vs. MAF interval for ", comparison,sep="")) + 
        xlab("MAF interval")+ 
        ylab("Jaccard ratio: (Mean Within)/(Mean Between)") +
        theme_bw() + guides(colour=guide_legend(title="") )
}

png("./plots/ITUvsSTU_JaccardRatio.png", width=800)
plotPopComparison("ITU vs STU", ggData)
dev.off()
png("./plots/IBSvsTSI_JaccardRatio.png", width=800)
plotPopComparison("IBS vs TSI", ggData)
dev.off()
png("./plots/CEUvsGBR_JaccardRatio.png", width=800)
plotPopComparison("CEU vs GBR", ggData)
dev.off()

#Color population vs YRI
png("./plots/YRIvsAll_pop_color_interval.png", width=800)
ggplot(data=subset(ggData,grepl("YRI",ggData$pop)&grepl("interval",ggData$cumulative)), aes(x=interval_percent, y=ratio, group=pop, colour=pop)) +
    geom_line() + geom_point() + ggtitle("Jaccard Ratio vs. MAF Interval for YRI vs X") + 
    xlab("MAF interval")+ 
    ylab("Jaccard ratio: (Mean Within)/(Mean Between)") +
    theme_bw() + guides(colour=guide_legend(title="Comparison") )
dev.off()

#Color continent vs YRI
png("./plots/YRIvsAll_continent_color_interval.png", width=800)
ggplot(data=subset(ggData,grepl("YRI",ggData$pop)&grepl("interval",ggData$cumulative)), aes(x=interval_percent, y=ratio, group=pop, colour=super1)) +
    geom_line() + geom_point() + ggtitle("Jaccard Ratio vs. MAF for YRI (Interval)") + 
    xlab("MAF interval")+ 
    ylab("Jaccard ratio: (Mean Within)/(Mean Between)")
dev.off()

#Color continent vs YRI (cumulative)
png("./plots/YRIvsAll_continent_color_cumulative.png", width=800)
ggplot(data=subset(ggData,grepl("YRI",ggData$pop)&grepl("cumulative",ggData$cumulative)), aes(x=interval_percent, y=ratio, group=pop, colour=super1)) +
    geom_line() + geom_point() + ggtitle("Jaccard Ratio vs. MAF for YRI (Cumulative)") + 
    xlab("MAF interval")+ 
    ylab("Jaccard ratio: (Mean Within)/(Mean Between)")
dev.off()

#Color continent vs YRI (interval)
png("./plots/YRIvsAll_continent_color_interval.png", width=1000)
ggplot(data=subset(ggData,grepl("YRI",ggData$pop)&grepl("interval",ggData$cumulative)), aes(x=interval_percent, y=ratio, group=pop, colour=super1)) +
    geom_line() + geom_point() + ggtitle("Jaccard Index vs. MAF Interval for YRI vs X") + 
    xlab("MAF interval")+ 
    ylab("Jaccard ratio: (Mean Within)/(Mean Between)")+
    theme_bw() + guides(colour=guide_legend(title="Continent of X") )
dev.off()
# k-medroid clustering

png("./plots/IBSvsTSI_clustering_efficiency.png", width=800)
ggplot(data=subset(ggData,grepl("IBS vs TSI",ggData$pop)), aes(x=interval_percent, y=clustering, group=cumulative, colour=cumulative)) +
    geom_line() + geom_point() + ggtitle("Clustering efficiency vs. MAF (IBS vs TSI)") + 
    xlab("MAF interval")+ 
    ylab("Jaccard index for clustering")
dev.off()

png("./plots/ITUvsSTU_clustering_efficiency.png", width=800)
ggplot(data=subset(ggData,grepl("ITU vs STU",ggData$pop)), aes(x=interval_percent, y=clustering, group=cumulative, colour=cumulative)) +
    geom_line() + geom_point() + ggtitle("Clustering efficiency vs. MAF (ITU vs STU)") + 
    xlab("MAF interval")+ 
    ylab("Jaccard index for clustering")
dev.off()

png("./plots/CEUvsGBR_clustering_efficiency.png", width=800)
ggplot(data=subset(ggData,grepl("CEU vs GBR",ggData$pop)), aes(x=interval_percent, y=clustering, group=cumulative, colour=cumulative)) +
    geom_line() + ggtitle("Clustering efficiency vs. MAF (CEU vs GBR)") + 
    xlab("MAF interval")+ 
    ylab("Jaccard index for clustering")
dev.off()

# subpop HCL
hclplot <- function(jaccardMat, pop1, pop2, name){
    rownames(jaccardMat) <- sample[,2]
    popsubset <- grepl(pop1, rownames(jaccardMat))|grepl(pop2, rownames(jaccardMat))
    jaccardMat <- jaccardMat[popsubset,popsubset]
    diag(jaccardMat) <- 0
    jaccardMat <- (max(jaccardMat)-jaccardMat)/max(jaccardMat)
    lowest.node <- max(jaccardMat)
    hclObj <- as.dendrogram(hclust(as.dist(jaccardMat),method="average"), hang=lowest.node*.1)
    labels_colors(hclObj) <- c('red','blue')[as.numeric(as.factor(rownames(jaccardMat)))[order.dendrogram(hclObj)]]
    plot(hclObj,ylim=c(1-lowest.node*1.25,1), main = paste0(pop1,', ',pop2,' Hierarchical Clustering\n', name), axes=FALSE)
}

createHCLDir <- function(jaccardList, pop1, pop2){
    dir <- paste('./plots/hcl_',pop1,'vs',pop2)
    dir.create(dir, showWarnings = FALSE, recursive = FALSE)
    invisible(lapply(names(jaccardList), function(name){
        png(paste(dir,'/', name,".png",sep=""), width=1600)
        hclplot(jaccard.interval.list[[name]], pop1, pop2, paste0("MAF: ",substring(name,2),"%"))
        dev.off()
    }))
    system(paste0("convert -delay 150 -loop 0 '", dir, "'/*.png '", dir, "'/animation.gif"))
}
names(jaccard.interval.list) <- paste0(LETTERS[16:1],gsub('%','',interval_percent[1:16]))
createHCLDir(jaccard.interval.list,'CEU','CHB')
createHCLDir(jaccard.interval.list,'CHB','JPT')
createHCLDir(jaccard.interval.list,'CHB','CHS')
createHCLDir(jaccard.interval.list,'IBS','TSI') 
createHCLDir(jaccard.interval.list,'ITU','STU')
createHCLDir(jaccard.interval.list,'LWK','YRI')
createHCLDir(jaccard.interval.list,'PJL','PUR')
createHCLDir(jaccard.interval.list,'GBR','FIN')
createHCLDir(jaccard.interval.list,'PJL','GBR')


