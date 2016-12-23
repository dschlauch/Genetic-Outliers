library(data.table)

gazalFilter <-"TGP2261"
source('~/1000GP/read1000GPsupportFiles.R')

# system.time(genotypes <- fread("zcat ~/1000GP/data/combinedFiltered19_1.gz", sep=" ", nrows=-1, header=F))
# saveRDS(genotypes[,pop%in%c("STU","ITU"), with=F], "/n/regal/quackenbush_lab/dschlauch/genotypes19_1_itustu.rds")
# saveRDS(genotypes, "/n/regal/quackenbush_lab/dschlauch/genotypes19_1.rds")
# rm(genotypes)
# gc()
# 
system.time(genotypes <- fread("zcat ~/1000GP/data/combinedFiltered19_2.gz", sep=" ", nrows=-1, header=F))
saveRDS(genotypes[,pop%in%c("STU","ITU"), with=F], "/n/regal/quackenbush_lab/dschlauch/genotypes19_2_itustu.rds")
saveRDS(genotypes, "/n/regal/quackenbush_lab/dschlauch/genotypes19_2.rds")
rm(genotypes)
gc()

# system.time(genotypes <- fread("zcat ~/1000GP/data/combinedFiltered19_3.gz", sep=" ", nrows=-1, header=F))
# saveRDS(genotypes[,pop%in%c("STU","ITU"), with=F], "/n/regal/quackenbush_lab/dschlauch/genotypes19_3_itustu.rds")
# saveRDS(genotypes, "/n/regal/quackenbush_lab/dschlauch/genotypes19_3.rds")
# rm(genotypes)
# gc()

# system.time(genotypes <- fread("zcat ~/1000GP/data/combinedFiltered19_4.gz", sep=" ", nrows=-1, header=F))
# saveRDS(genotypes[,pop%in%c("STU","ITU"), with=F], "/n/regal/quackenbush_lab/dschlauch/genotypes19_4_itustu.rds")
# saveRDS(genotypes, "/n/regal/quackenbush_lab/dschlauch/genotypes19_4.rds")
# rm(genotypes)
# gc()

# library(stego)
# 
# gazalFilter <-"TGP2261"
# source('~/1000GP/read1000GPsupportFiles.R')
# 
# 
# pairs <- do.call(cbind, lapply(unique(group), function(x){
#     combn(unique(pop[group==x]),2)
# }))
# 
# 
# genotypes <- readRDS("/n/regal/quackenbush_lab/dschlauch/genotypes19_3.rds")
# apply(pairs,2,function(x){
#     popFilter <-  pop %in% x
#     sampleFilter <- popFilter & qcFilter
#     
#     results <- run_stego(genotypes[,rep(sampleFilter,each=2),with=F], phased=T, groups="all.together",sampleNames=sampleIDs[sampleFilter], labels=pop[sampleFilter], super=group[sampleFilter], minVariants=2, blocksize=1, saveDir=NA, verbose=T, simFun=cov)
#     saveFile <- paste(x, collapse="-")
#     saveRDS(results, file=paste0("~/1000GP/block19_3/",saveFile,".rds"))
#     
# })
# 
# genotypes <- readRDS("/n/regal/quackenbush_lab/dschlauch/genotypes19_4.rds")
# apply(pairs,2,function(x){
#     popFilter <-  pop %in% x
#     sampleFilter <- popFilter & qcFilter
#     
#     results <- run_stego(genotypes[,rep(sampleFilter,each=2),with=F], phased=T, groups="all.together",sampleNames=sampleIDs[sampleFilter], labels=pop[sampleFilter], super=group[sampleFilter], minVariants=2, blocksize=1, saveDir=NA, verbose=T, simFun=cov)
#     saveFile <- paste(x, collapse="-")
#     saveRDS(results, file=paste0("~/1000GP/block19_4/",saveFile,".rds"))
#     
# })
# 
