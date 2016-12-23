library(ggplot2)
n <- 200
mac <- 1:n

w_k <- choose(n,2)/choose(mac,2)
w_k[1] <-0
df <- data.frame(w_k=w_k, mac=mac)


setEPS()
postscript("supplemental_figure_wk.eps")

ggplot(df[2:20,]) +geom_line(aes(x=mac,y=w_k)) +theme_bw() + 
    labs(y=expression(w[k]),x="Minor allele count") + ggtitle(expression(w[k]*" vs Minor allele count for 200 alleles")) +
    geom_point(data=df[1:20,], aes(x=mac,y=w_k))
dev.off()
