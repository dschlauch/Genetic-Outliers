
# incorporating gazal filtering -------------------------------------------
if(!exists("gazalFilter")){
    gazalFilter <-NA
}
if(gazalFilter=="NA"|is.na(gazalFilter)){
    qcFilter <- NULL
} else {
    gazal_filtered <- read.csv("./data/gazal_filtered.csv",stringsAsFactors=F)
    gazal_filtered <- gazal_filtered[order(gazal_filtered[,1]),]
    qcFilter <- gazal_filtered[,gazalFilter]=="YES"    
}

gazal_related_orig <- data.table(read.csv("./data/gazal_related.csv"))
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