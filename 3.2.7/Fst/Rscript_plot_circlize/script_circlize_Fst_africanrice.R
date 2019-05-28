# plus d'infos: https://cran.r-project.org/web/packages/circlize/vignettes/circlize.pdf
#install.packages("circlize")

# files
setwd("/media/thibaultleroy/Seagate Expansion Drive/AfricanRice/")
SNP_FstValues=read.table("out.weir.fst_SNPperSNP.withoutNA.chronly.withheader",header=TRUE)


### Colors
col_uniq =rep(c("darkgrey", "black"), 6)
col =rep(c("orangered", "red"), 6)
col2 = rep(c("steelblue3","steelblue4"),6)

###

library("circlize")
library(plotrix)

png(filename="test_circlize_png_Fst_africanrice_070519_allcolorscale.png",width=1000,height=1000,units ="px")

circos.clear()
par(mar=c(1,1,1,1), lwd=0.1, cex =0.7)
circos.par("track.height"=0.2,start.degree=90,gap.degree=4)
circos.initialize(factors = SNP_FstValues$CHROM, x = SNP_FstValues$POS)

### first track : Fst SNP/SNP 
circos.trackPlotRegion(factors=SNP_FstValues$CHROM, y=SNP_FstValues$WEIR_AND_COCKERHAM_FST, bg.border = NA, track.height=0.22, panel.fun=function(POS,WEIR_AND_COCKERHAM_FST){
  #circos.axis(h=1.5, labels.facing="inside", col="navyblue",labels.col="navyblue",lwd=2,labels.cex=0.9,major.at=c(100000,10000000,20000000,30000000,40000000,50000000),labels=c(0,"10Mb","20Mb","30Mb","40Mb","50Mb"))
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.index = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), max(ylim)+0.7, sector.index, facing = "inside", cex = 1.4,font=2)
  circos.yaxis("left", lwd=1,labels.cex=0.85,at=c(-0.2,0,0.2,0.4,0.6,0.8,1),labels=c(-0.2,0,0.2,0.4,0.6,0.8,1))
})

circos.trackPoints(SNP_FstValues$CHROM, SNP_FstValues$POS, SNP_FstValues$WEIR_AND_COCKERHAM_FST, pch = 20,col=color.scale(SNP_FstValues$WEIR_AND_COCKERHAM_FST,c(1,0.8,0.9),c(1,1,0),0),cex=0.1)



###  secondtrack: Fst values slid windows 10 kb (non-overlapping)

Wind_FstValues=read.table("out.windowed.weir.fst_slidwin10kb.chronly.withheader",header=TRUE)

circos.trackPlotRegion(factors=Wind_FstValues$CHROM, y=Wind_FstValues$WEIGHTED_FST, bg.border = NA, track.height=0.22, panel.fun=function(BINSTART,WEIGHTED_FST){
  circos.axis(h=-0.7, labels.facing="inside", col="navyblue",labels.col="navyblue",lwd=2,labels.cex=0.9,major.at=c(100000,10000000,20000000,30000000,40000000,50000000),labels=c(0,"10Mb","20Mb","30Mb","40Mb","50Mb"))
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.index = get.cell.meta.data("sector.index")
  #circos.text(mean(xlim), max(ylim), sector.index, facing = "inside", cex = 1)
  circos.yaxis("left", lwd=1,labels.cex=0.85,at=c(-0.2,0,0.2,0.4,0.6,0.8,1),labels=c(-0.2,0,0.2,0.4,0.6,0.8,1))
})
circos.trackLines(Wind_FstValues$CHROM, Wind_FstValues$BIN_START, Wind_FstValues$WEIGHTED_FST, type="h",baseline=0,lwd=ifelse(Wind_FstValues$WEIGHTED_FST>=0.8,0.8,ifelse(Wind_FstValues$WEIGHTED_FST>=0.5,0.5,0.12)),col=ifelse(Wind_FstValues$WEIGHTED_FST>=0.8,"red",ifelse(Wind_FstValues$WEIGHTED_FST>=0.5,"orange","black")))
#circos.trackLines(Wind_FstValues$CHROM, Wind_FstValues$BIN_START, Wind_FstValues$WEIGHTED_FST, type="h",baseline=0,lwd=ifelse(Wind_FstValues$WEIGHTED_FST>=0.8,0.8,ifelse(Wind_FstValues$WEIGHTED_FST>=0.5,0.5,0.12)),col=color.scale(SNP_FstValues$WEIR_AND_COCKERHAM_FST,c(1,0.8,0.9),c(1,1,0),0))


## add boxplot
#par(fig = c(0.05,0.25, 0.5, 0.5), mar = c(4,6,1,1), new = T)  
par(new=TRUE)
par(oma=c(1,1,0,1))
par(mfcol=c(3,3), mfg=c(2,2))
par(mar=c(1,4,0,1))

# Filled Density Plot
d <- density(SNP_FstValues$WEIR_AND_COCKERHAM_FST,bw = 0.015)
#plot(d, main="Distribution of FST values")
polygon(d, col="yellow", border="forestgreen") 
bidon=d[which(d$x < -0.0283749)]
polygon(bidon, col="blue", border="forestgreen") 

bidon=SNP_FstValues[which(SNP_FstValues$WEIR_AND_COCKERHAM_FST < -0.0283749),]
e <- density(SNP_FstValues[which(SNP_FstValues$WEIR_AND_COCKERHAM_FST < -0.0283749),]$WEIR_AND_COCKERHAM_FST,bw = 0.015)
#plot(density(SNP_FstValues$WEIR_AND_COCKERHAM_FST,bw = 0.01),main="")

f <- density(SNP_FstValues[which(SNP_FstValues$WEIR_AND_COCKERHAM_FST > 0.4058460),]$WEIR_AND_COCKERHAM_FST,bw = 0.015)
                           #plot(density(SNP_FstValues$WEIR_AND_COCKERHAM_FST,bw = 0.01),main="")

quantile(SNP_FstValues$WEIR_AND_COCKERHAM_FST, c(.05, 0.95)) 

dev.off()

### Distribution Fst (bw=0.015)
library(ggplot2)
dens <- density(SNP_FstValues$WEIR_AND_COCKERHAM_FST,bw=0.015)
df <- data.frame(x=dens$x, y=dens$y)
probs <- c(0.05, 0.90, 0.95, 0.99,0.999)
quantiles <- quantile(SNP_FstValues$WEIR_AND_COCKERHAM_FST, prob=probs)
df$quant <- factor(findInterval(df$x,quantiles))
ggplot(df, aes(x,y)) + geom_line() + geom_ribbon(aes(ymin=0, ymax=y, fill=quant)) + scale_x_continuous(breaks=quantiles) + scale_fill_brewer(guide="none",palette = "Spectral",direction=-1)+
  xlab(expression("F"[ST]))+ylab("density")+ theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))