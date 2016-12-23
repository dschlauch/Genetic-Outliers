library(stego)

gazalFilter <-"TGP2261"
source('~/1000GP/read1000GPsupportFiles.R')

# genotypes <- readRDS("~/1000GP/data/combinedFiltered1000.rds")
# genotypes <- readRDS("~/1000GP/data/combinedFiltered100.rds")
genotypes <- readRDS("/n/regal/quackenbush_lab/dschlauch/genotypes19.rds")
# genotypes <- readRDS("~/1000GP/data/combinedFiltered200.rds")
genotypes <- readRDS("~/1000GP/data/genotypes19_all_itustu.rds")

# results <- run_stego(genotypes[,rep(qcFilter,each=2),with=F], phased=T, groups="all.together",sampleNames=sampleIDs[qcFilter], labels=pop[qcFilter], super=group[qcFilter], minVariants=5, blocksize=1, saveResult=NA, verbose=T)

directory <- "./allVariants_CHR14/"
dir.create(directory)
pairs <- do.call(cbind, lapply(unique(group), function(x){
    combn(unique(pop[group==x]),2)
}))

apply(pairs,2,function(x){
    popFilter <-  pop %in% x
    sampleFilter <- popFilter & qcFilter
    
    results <- run_stego(genotypes[,rep(sampleFilter,each=2),with=F], phased=T, groups="all.together",sampleNames=sampleIDs[sampleFilter], labels=pop[sampleFilter], super=group[sampleFilter], minVariants=2, blocksize=1, saveDir=NA, verbose=T, simFun=cor)
    saveFile <- paste(x, collapse="-")
    saveRDS(results, file=paste0(directory,saveFile,".rds"))
    
})

# results <- run_stego(genotypes, phased=T, groups="all.together", minVariants=2, blocksize=1, saveDir=NA, verbose=T, simFun=cov)
# saveRDS(results, file=paste0("all_ITU_STU.rds"))
# results <- run_stego(genotypes[,rep(qcFilter,each=2),with=F], phased=T, groups="each.separately",sampleNames=sampleIDs[qcFilter], labels=pop[qcFilter], super=group[qcFilter], minVariants=5, blocksize=1, saveResult=NA, verbose=T)

# pairwiseResults <- run_stego(genotypes[,rep(qcFilter,each=2),with=F], phased=T, groups="pairwise.within.superpop", sampleNames=sampleIDs[qcFilter], labels=pop[qcFilter], super=group[qcFilter], minVariants=2, blocksize=1, saveDir=NA, verbose=T, simFun=cov, cores=60)

# pairwiseResults <- readRDS(file="~/1000GP/data/pairwiseResults40.rds")
# pairwiseResults <- readRDS(file="~/1000GP/data/pairwiseResults200.rds")
# pairwiseResults <- readRDS(file="~/1000GP/data/pairwiseResults1000.rds")
# pairwiseResults <- readRDS(file="~/1000GP/data/pairwiseResults146.rds")
# saveRDS(pairwiseResults, file="~/1000GP/data/pairwiseResults19.rds")

library(stego)

gazalFilter <-"TGP2261"
source('~/1000GP/read1000GPsupportFiles.R')

# genotypes <- readRDS("~/1000GP/data/combinedFiltered1000.rds")
# genotypes <- readRDS("~/1000GP/data/combinedFiltered100.rds")
genotypes <- readRDS("/n/regal/quackenbush_lab/dschlauch/genotypes19_1.rds")
# genotypes <- readRDS("~/1000GP/data/combinedFiltered200.rds")

# results <- run_stego(genotypes[,rep(qcFilter,each=2),with=F], phased=T, groups="all.together",sampleNames=sampleIDs[qcFilter], labels=pop[qcFilter], super=group[qcFilter], minVariants=5, blocksize=1, saveResult=NA, verbose=T)

pairs <- do.call(cbind, lapply(unique(group), function(x){
    combn(unique(pop[group==x]),2)
}))
apply(pairs[,57:45],2,function(x){
    popFilter <-  pop %in% x
    sampleFilter <- popFilter & qcFilter
    
    results <- run_stego(genotypes[,rep(sampleFilter,each=2),with=F], phased=T, groups="all.together",sampleNames=sampleIDs[sampleFilter], labels=pop[sampleFilter], super=group[sampleFilter], minVariants=2, blocksize=1, saveDir=NA, verbose=T, simFun=cov)
    saveFile <- paste(x, collapse="-")
    saveRDS(results, file=paste0("./block19_1/",saveFile,".rds"))
    
})
