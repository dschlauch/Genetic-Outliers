library(MASS)
library(stego)
library(data.table)

n_variants <- 20000
mean_af <- .05
range_af <- .003
n_group1 <- 500
n_group2 <- 500

set.seed(1)

sampleNames <- c(paste0("Pop1_",1:n_group1),paste0("Pop2_",1:n_group2))

af_ancestral <- rexp(n_variants,1/mean_af)
af_1 <- af_ancestral + runif(n_variants,-range_af,range_af)
af_2 <- af_ancestral + runif(n_variants,-range_af,range_af)

af_1[af_1<0]<-0
af_2[af_2<0]<-0

genotypes_group1 <- matrix(rbinom(n = n_group1*n_variants*2, size = 1, prob = af_1),ncol = n_group1*2)
genotypes_group2 <- matrix(rbinom(n = n_group2*n_variants*2, size = 1, prob = af_2),ncol = n_group2*2)
generateggplotDT <- function(stego.correction, varcov.correction, sampleNames, components=c(1,2,3)){
    components.df <- data.frame(samples=rep(sampleNames,2))
    components.df$Vector1 <- c(stego.correction[,components[1]],varcov.correction[,components[1]])
    components.df$Vector2 <- c(stego.correction[,components[2]],varcov.correction[,components[2]])
    components.df$Vector3 <- c(stego.correction[,components[3]],varcov.correction[,components[3]])
    components.df$pop <- c(rep("Group 1",n_group1),rep("Group 2",n_group1))
    components.df$method <- c(rep("STEGO",length(sampleNames)),rep("Variance-covariance",length(sampleNames)))
    data.table(components.df)
}

unphasedData <- cbind(genotypes_group1,genotypes_group2)[,c(T,F)]+cbind(genotypes_group1,genotypes_group2)[,c(F,T)]
stegoRes <-run_stego(unphasedData, phased = F)

# # Eigenstrat scaling method
scaledData <- t(apply(unphasedData,1,function(x){
    p <- (1+sum(x))/(2*(n_group1+n_group2)+2)
    (x-2*p)/sqrt(p*(1-p))
    }))

scaledData <- scaledData[!is.na(rowSums(scaledData)),]

varcov <- cov(scaledData)
# varcov <- cov(unphasedData)

stego.correction <- eigen(apply(stegoRes$s_matrix_dip,2,scale),symmetric=T)$vectors
varcov.correction <- eigen(apply(varcov,2,scale),symmetric=T)$vectors

components.dt <- generateggplotDT(stego.correction, varcov.correction, sampleNames, c(1,2,3))
pcpplots <- ggplot(components.dt, aes(x = Vector1, y = Vector2)) + 
    ggtitle("Principal Component Plots") + ylab("Component 2")+ xlab("Component 1") + theme_bw() + facet_grid(.~method, scales="free") +
    geom_point(aes(color=pop), size=2, alpha=.7) + 
    guides(colour=guide_legend(title="Populations:")) + 
    theme(plot.title = element_text(hjust = 0.5),legend.position="top")+ 
    scale_colour_brewer(palette="Set1") +
    coord_fixed(ratio = 10)
print(pcpplots)
pdf("./simulated_PCP_plots.pdf")
print(pcpplots)
dev.off()

withinVsBetween2 <-  function(components.dt){
    pair <- unique(components.dt$pop)
    totalmseSTEGO <- mean(components.dt[method=="STEGO"]$Vector1^2+components.dt[method=="STEGO"]$Vector2^2)
    Vector1DifStego <- components.dt[method=="STEGO",Vector1-mean(Vector1),by=pop]$V1
    Vector2DifStego <- components.dt[method=="STEGO",Vector2-mean(Vector2),by=pop]$V1
    stegoMSE <- mean(Vector1DifStego^2+Vector2DifStego^2)
    
    totalmseVARCOV <- mean(components.dt[method=="Variance-covariance"]$Vector1^2+components.dt[method=="Variance-covariance"]$Vector2^2)
    Vector1Difvarcov <- components.dt[method=="Variance-covariance",Vector1-mean(Vector1),by=pop]$V1
    Vector2Difvarcov <- components.dt[method=="Variance-covariance",Vector2-mean(Vector2),by=pop]$V1
    varcovMSE <- mean(Vector1Difvarcov^2+Vector2Difvarcov^2)
    c("stego"=stegoMSE/totalmseSTEGO, "Variance-covariance"=varcovMSE/totalmseVARCOV)  
}
withinVsBetween2(components.dt)

