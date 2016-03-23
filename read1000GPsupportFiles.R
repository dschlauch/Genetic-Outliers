
#outputDir <- "filtered1000_fixed_phi"
#genotypeFile <- "./data/1000GP_Phase3_chr10.hap.gz"
genotypeFile <- "./data/combinedFiltered40.gz"
numberOfLines <- 5000
minVariants <- 10
numCores <- 4
args<-commandArgs(TRUE)
outputDir <- '.'
# gazalFilter <-"TGP2261"
gazalFilter <- "NA"
ldPrune <- 10
if(length(args)!=0){
    genotypeFile <- args[1]
    numberOfLines <- as.numeric(args[2])
    minVariants <- as.numeric(args[3])
    numCores <- as.numeric(args[4])
    outputDir <- args[5]
    gazalFilter <- args[6]
    ldPrune <- as.numeric(args[7])
    dir.create(paste0("~/1000GP/plots/s_distributions/",outputDir))
    dir.create(paste0("~/1000GP/plots/s_distributions/",outputDir,"/plotdata"))
}

# incorporating gazal filtering -------------------------------------------
if(gazalFilter=="NA"|is.na(gazalFilter)){
    qcFilter <- NULL
} else {
    gazal_filtered <- read.csv("./data/gazal_filtered.csv",stringsAsFactors=F)
    gazal_filtered <- gazal_filtered[order(gazal_filtered[,1]),]
    qcFilter <- gazal_filtered[,gazalFilter]=="YES"    
}

# Real data
sample <- read.table("~/1000GP/data/1000GP_Phase3.sample", sep=" ", header=T)
sampleIDs <- as.character(sample[,1])
pop <- as.character(sample[,2])
group <- as.character(sample[,3])
sex <- as.character(sample[,4])
hap.pop <- rep(pop,each=2)
hap.sampleIDs <- rep(as.character(sample[,1]),each=2)

popGroup <- melt(table(pop,group))
popGroup <- popGroup[popGroup$value>0,]
rownames(popGroup) <- popGroup$pop