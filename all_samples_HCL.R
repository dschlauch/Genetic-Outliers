library(gplots)
colorCodes <- c(AFR="red", AMR="green", EUR="blue", SAS="yellow", EAS="black")
colorCodesPop <- c("dodgerblue2","#E31A1C", # red
                          "green4",
                          "#6A3D9A", # purple
                          "#FF7F00", # orange
                          "black","gold1",
                          "skyblue2","#FB9A99", # lt pink
                          "palegreen2",
                          "#CAB2D6", # lt purple
                          "#FDBF6F", # lt orange
                          "gray70", "khaki2",
                          "maroon","orchid1","deeppink1","blue1","steelblue4",
                          "darkturquoise","green1","yellow4","yellow3",
                          "darkorange4","brown","black")
names(colorCodesPop)<-unique(pop)
colorCodesPop <- colorCodesPop[order(names(colorCodesPop))]
labelCol <- function(x) {
  if (is.leaf(x)) {
    ## fetch label
    label <- attr(x, "label")
    #     code <- substr(label, 1, 1)
    ## use the following line to reset the label to one letter code
    # attr(x, "label") <- code
    attr(x, "nodePar") <- list(lab.col=colorCodes[label])
  }
  return(x)
}

# Load in data
# jaccardMatrix <- read.csv('~/1000GP/output_0_20/combined_jaccard.csv', row.names=1)
# results <- readRDS("~/1000GP/plots/s_distributions/filtered40_LD10_all/plotdata/Allsamples_data.rds")
results <- readRDS("~/1000GP/plots/s_distributions/filtered40_LD10_min100_all/plotdata/Allsamples_data.rds")
# results <- readRDS("~/1000GP/plots/s_distributions/plotdata/STU_ITU_data.rds")
# results <- readRDS("~/1000GP/plots/s_distributions/plotdata/CEU_YRI_data.rds")
# results <- readRDS("~/1000GP/plots/s_distributions/plotdata/IBS_TSI_data.rds")
jaccardMatrix <- results$s_matrix_dip
varcovMatrix <- results$varcovMat

# ### hierarchical clustering
# 
# subset <- sample(2504,100)
subset <- T
jaccardMatrix <- as.matrix(jaccardMatrix)[subset,subset]
varcovMatrix <- as.matrix(varcovMatrix)[subset,subset]
# diag(jaccardMatrix)<-0
whichPeople <- sampleIDs%in%rownames(jaccardMatrix) 
rownames(jaccardMatrix) <- group[whichPeople]
colnames(jaccardMatrix) <- pop[whichPeople]
jaccardMatrix <- jaccardMatrix[order(rownames(jaccardMatrix),colnames(jaccardMatrix)),order(rownames(jaccardMatrix),colnames(jaccardMatrix))]
rownames(varcovMatrix) <- group[whichPeople]
colnames(varcovMatrix) <- pop[whichPeople]
varcovMatrix <- varcovMatrix[order(rownames(varcovMatrix),colnames(varcovMatrix)),order(rownames(varcovMatrix),colnames(varcovMatrix))]

hc <- hclust(as.dist(max(jaccardMatrix)-jaccardMatrix),method="average")
d <- dendrapply(as.dendrogram(hc, hang=max(jaccardMatrix)*.1), labelCol)
plot(d, main="Hierarchical Clustering across superpopulations")
# 
# jaccardMatrix[jaccardMatrix>.005]<-.005
# dev.off()
# heatmap.2(jaccardMatrix, Rowv=d, Colv=d, dendrogram="none", trace="none",
#           labRow="",labCol="",key=FALSE, RowSideColors=colorCodes[rownames(jaccardMatrix)], ColSideColors=colorCodesPop[colnames(jaccardMatrix)],
# #           lwid = c(1,10),lhei = c(.01,5), margins = c(5,10))
# )
# legend(.1,.5, names(colorCodes), inset=c(-0.01,0),  lty=1, lwd=10, col=colorCodes, cex = 0.75)
# legend("topright", names(colorCodesPop[1:13]),  lty=1, lwd=10, col=colorCodesPop[1:13], cex = 0.5, horiz = T,  inset=c(0,.05))
# legend("topright", names(colorCodesPop[14:26]),  lty=1, lwd=10, col=colorCodesPop[14:26], cex = 0.5, horiz = T,  inset=c(0,.10))

plotHeatmap <-  function(x, subset=NA, titleText="GSM"){
    if(is.na(subset)){
        subset <- rep(T,nrow(x))
    }
    simMat <- as.matrix(x)[subset,subset]
    diag(simMat)<- NA
    
    
    simMat[] <- apply(scale(simMat),1,rank)
    simMat <- simMat + t(simMat)
    colorCodesPop <- colorCodesPop[unique(colnames(simMat))]
    breaks <- c(seq(quantile(simMat,.55, na.rm=T),quantile(simMat,.80, na.rm=T),length.out=100), # Blue breaks
                seq(quantile(simMat,.80, na.rm=T),quantile(simMat,.99, na.rm=T),length.out=100)) # Red breaks
    heatmap.2(simMat, Rowv=F, Colv=F, dendrogram="none", trace="none", main="", cex.main=12, col="bluered",
        labRow="",labCol="",key=FALSE, RowSideColors=colorCodes[rownames(simMat)], ColSideColors=colorCodesPop[colnames(simMat)],
        breaks=breaks)

    title(main=titleText,cex.main = 4, adj=.7)
#     legend(0.1,.5, names(colorCodes), inset=.1,  lty=1, lwd=30, col=colorCodes, cex = 2)
    text(.18,.63, "AFR - African", cex=1.5)
    text(.14,.47, "AMR - Ad Mixed American", cex=1.5)
    text(.17,.35, "EAS - East Asian", cex=1.5)
    text(.17,.20, "EUR - European", cex=1.5)
    text(.16,.05, "SAS - South Asian", cex=1.5)
    legend("bottomright", names(colorCodesPop[1:7]),  lty=1, lwd=30, col=colorCodesPop[1:7], cex = 1.5, horiz = F,  inset=c(.59,.77), y.intersp=1.1, title="AFR")
    legend("bottomright", names(colorCodesPop[8:11]),  lty=1, lwd=30, col=colorCodesPop[8:11], cex = 1.5, horiz = F,  inset=c(.45,.77), y.intersp=1.1, title="AMR")
    legend("bottomright", names(colorCodesPop[12:16]),  lty=1, lwd=30, col=colorCodesPop[12:16], cex = 1.5, horiz = F,  inset=c(.33,.77), y.intersp=1.1, title="EAS")
    legend("bottomright", names(colorCodesPop[17:21]),  lty=1, lwd=30, col=colorCodesPop[17:21], cex = 1.5, horiz = F,  inset=c(.18,.77), y.intersp=1.1, title="EUR")
    legend("bottomright", names(colorCodesPop[22:26]),  lty=1, lwd=30, col=colorCodesPop[22:26], cex = 1.5, horiz = F,  inset=c(.04,.77), y.intersp=1.1, title="SAS")
}

## Usage
## Commented so can source from crypticness.R
png(paste0("./plots/s_distributions/",outputDir,"/GSM.png"), width = 1400, height = 1200)
plotHeatmap(jaccardMatrix,title="Genetic Similarity Matrix")
dev.off()
png(paste0("./plots/s_distributions/",outputDir,"/varcovGSM.png"), width = 1400, height = 1200)
plotHeatmap(varcovMatrix,title="Varcov GSM")
dev.off()



# for removing possible related individuals
# related <- which(jaccardMatrix>.1,arr.ind=T)[,1]
