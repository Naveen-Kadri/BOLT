infile <- snakemake@input [["infile"]]
plot_file=snakemake@output[["plot_file"]]
maf_thresh <- as.numeric(snakemake@params[["maf_thresh"]])
phenotype=snakemake@wildcards[["pheno"]]
breed=snakemake@wildcards[["breed"]]
sex=snakemake@wildcards[["sex"]]

mytitle=paste0(phenotype, "_", breed, "_", sex)
cat ("mytitle is ", mytitle, "\n")
cnames  <- c('chr', 'pos','freq', 'p')
res  <- data.frame (matrix(   scan (infile),ncol=4,byrow=T ))
colnames (res) <- cnames
res$maf <- pmin (res$freq , 1-res$freq)
res <- res [res$maf > maf_thresh, ]


tiff (plot_file, height=12,width=24,units="in",res=300)
chrcol <- 1
poscol <- 2
pcol <- 4
colors <- c("darkgray", "black")
fontsize=1.5  #font size for the ylab
mymar <- c (4, 6, 2, 2)
par (mar=mymar, cex.lab=1.5,cex.axis=1.5)
thresh <- -log10(0.05/1e6)

#alternate the color for consecutive chromosomes
res$mycol <- NA
chrs <- unique (res$chr) ; chrs
evens <- chrs [1:length(chrs)%%2==0]
odds <- chrs [1:length(chrs)%%2==1]
res [which (res[,chrcol] %in% evens), 'mycol' ]  <-  colors [1]
res [ which (res[,chrcol] %in% odds), 'mycol' ]  <-  colors [2]

#get chromosome sizes
sizes <- aggregate (res$pos, by=list(res$chr),max) [,2]

#continuous positions
toadd <- c (0, cumsum (sizes))
toadd <- toadd [-c(length(toadd))]
res$pos <- res[, poscol] + toadd [res [, chrcol]]

#pvalue range
cat ("myrange :", min (res[,pcol]), max(res[,pcol]), "\n" )
maxi <- max (-log10(res[, pcol]))
cat ('the most significant pvalue = ', maxi,'\n')
     

##quick plot
##res <- res [res$p < 1e-3 ,];nrow (res)
plot (res$pos, -log10(res[, pcol]), col=res$mycol, xaxt="n", ylab=expression (-log [10](italic("P"))), ylim =c(0, maxi*1.1), cex.main=1 , cex.lab=fontsize, cex.axis=1, xlab="", pch=20, cex=1.5,yaxt="n")
axis (side=2, las=2)

#chromosome labels
ats <- aggregate (res$pos, by=list(res$chr),mean) [,2]
text (ats, rep(maxi*1.1, length(ats)), labels=chrs, cex=1)
abline (v= sizes + toadd,col='gray', lty=2 )
abline (h =thresh,  col='darkseagreen', lty=2 )
dev.off ()
