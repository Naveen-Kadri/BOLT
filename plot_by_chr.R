infile <- snakemake@input [["infile"]]
plot_file=snakemake@output[["plot_file"]]
maf_thresh <- as.numeric(snakemake@params[["maf_thresh"]])
#phenotype=snakemake@wildcards[["pheno"]]
#breed=snakemake@wildcards[["breed"]]
#sex=snakemake@wildcards[["sex"]]


#infile <- '/cluster/work/pausch/naveen/BOLT/BOLT/GRR/fv/female/forplot.txt'
#maf_thresh <- 0.005



cnames  <- c('chr', 'pos','freq', 'p')
res  <- data.frame (matrix(   scan (infile),ncol=4,byrow=T ))
colnames (res) <- cnames
res$maf <- pmin (res$freq , 1-res$freq)
res <- res [res$maf > maf_thresh, ]


##keep only below 1e-3
res <- res [res$p < 1e-3 ,]

##tiff (plot_file, height=12,width=20,units="in",res=300)
pdf (plot_file, height=12,width=20)
#par (mfrow = c (6,5), mar=c(2,2,1,1),oma=c(6,6,1,1))
chrs <- sort (unique (res$chr))
for (chr in 1:length(chrs)) {
    mychr <- res [res$chr == chr ,]
    mychr <- mychr [order (res$pos) ,]
    plot (mychr$pos/1e6, -log10(mychr$p), pch=20, xlab='',ylab='',bty='l',main=paste0('BTA',chr),las=1)
    cat (chr, "\n")
}



mtext (side=1, 'Position(Mb)', outer=TRUE, line=1,cex=1.5)
mtext (side=2, expression(-log[10]~p), outer=TRUE, line=1,cex=1.5)
dev.off ()
