# TL - 230519
library(ggplot2)

setwd("/media/thibaultleroy/Seagate Expansion Drive/projet_BNP-BHW/BayPass_100kSNPs")

### Generate a manhattan plot of XtX values
setwd("/media/thibaultleroy/Seagate Expansion Drive/projet_BNP-BHW/BayPass_100kSNPs")
realdata2=read.table("BayPass_merged_results_121017_Toc_Pluvio.txt.pseudoK.withheader",h=T)
head(realdata2)

### Generate a BF vs. XtX plot
# Mean annual Temperature (highlight SNP with BF values > 15 or 20 & XtX values > quantiles based on neutral simulations, see also ./3.3.6/Rscript_XtX/script_baypass_XtX_plots.R)
tiff("BFToc_XtX_250519.tiff", width = 8, height = 8, units = 'in', res = 600)
ggplot(realdata2, aes(x=xtx, y=BFMeanToc)) + geom_point(size=ifelse(realdata2$xtx >27.89115 & realdata2$BFMeanToc > 20,1,ifelse(realdata2$xtx>20.22704 & realdata2$BFMeanToc > 15,0.5,0.1)), colour=ifelse(realdata2$xtx >27.89115 & realdata2$BFMeanToc > 20,"red",ifelse(realdata2$xtx>20.22704 & realdata2$BFMeanToc > 15,"orange","grey")))+
xlab("XtX values")+ylab("Bayes Factors (covariable: Mean annual temperature)")+
geom_hline(yintercept = 15,linetype="dotted",color="black")+
geom_hline(yintercept = 20,linetype="dotted",color="black")+
geom_vline(xintercept = 27.89115, linetype="dotted",color="black")+
geom_vline(xintercept = 20.22704, linetype="dotted",color="black")+
theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ theme(legend.position = "none")+
theme(legend.title=element_text(size=12, face="bold"), legend.title.align=0.5, legend.text=element_text(size=10, face="italic"), axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
      axis.text.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.5,face="plain"),
      axis.text.y = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.5,face="plain"),  
      axis.title.x = element_text(colour="black",size=18,angle=0,hjust=.5,vjust=.2,face="italic"),
      axis.title.y = element_text(colour="black",size=18,angle=90,hjust=.5,vjust=.5,face="italic"))
dev.off()

# Precipitation sums
tiff("BFPrec_XtX_250519.tiff", width = 8, height = 8, units = 'in', res = 600)
ggplot(realdata2, aes(x=xtx, y=BFMeanPluvio)) + geom_point(size=ifelse(realdata2$xtx >27.89115 & realdata2$BFMeanPluvio > 20,1,ifelse(realdata2$xtx>20.22704 & realdata2$BFMeanPluvio > 15,0.5,0.1)), colour=ifelse(realdata2$xtx >27.89115 & realdata2$BFMeanPluvio > 20,"darkblue",ifelse(realdata2$xtx>20.22704 & realdata2$BFMeanPluvio > 15,"cornflowerblue","grey")))+
  xlab("XtX values")+ylab("Bayes Factors (covariable: Annual precipitation sums)")+
  geom_hline(yintercept = 15,linetype="dotted",color="black")+
  geom_hline(yintercept = 20,linetype="dotted",color="black")+
  geom_vline(xintercept = 27.89115, linetype="dotted",color="black")+
  geom_vline(xintercept = 20.22704, linetype="dotted",color="black")+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ theme(legend.position = "none")+
  theme(legend.title=element_text(size=12, face="bold"), legend.title.align=0.5, legend.text=element_text(size=10, face="italic"), axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=18,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=18,angle=90,hjust=.5,vjust=.5,face="italic"))
dev.off()

### Manhattan plots Bayes Factors
tiff("BF_Manhattanplots_chr1-6_250519.tiff", width = 7, height = 10, units = 'in', res = 600)
#png("BF_Manhattanplots_chr1-6_250519.png", width = 10, height = 6, units = 'in', res = 600)
par(mfrow = c(3, 2))
for (i in 1:6){
  # for each chromosome, keep only SNPs with XtX >= 15
  datacurrentchr=subset(realdata2,realdata2$chr==paste("Chr",i,sep="")&realdata2$xtx>=15&(realdata2$BFMeanToc>15|realdata2$BFMeanPluvio>15))
  # then plot these SNPs
  plot(datacurrentchr$BFMeanPluvio~datacurrentchr$pos.1,main=paste("Chromosome ",i),ylim=c(15,max(c(max(datacurrentchr$BFMeanToc),max(datacurrentchr$BFMeanPluvio)))),cex.main=1.5,cex.axis=1.2,cex.lab=1.4,pch=20,xlab=paste ("position - chr", i), ylab="BF",cex=1, col=ifelse(datacurrentchr$xtx >27.89115 & datacurrentchr$BFMeanPluvio > 20,"darkblue",ifelse(datacurrentchr$xtx>20.22704 & datacurrentchr$BFMeanPluvio > 15,"cornflowerblue","grey")))
  points(datacurrentchr$BFMeanToc~datacurrentchr$pos.1, pch=20, col=ifelse(datacurrentchr$xtx >27.89115 & datacurrentchr$BFMeanToc > 20,"red",ifelse(datacurrentchr$xtx>20.22704 & datacurrentchr$BFMeanToc > 15,"orange","grey")))
}
dev.off()

tiff("BF_Manhattanplots_chr7-12_250519.tiff", width = 7, height = 10, units = 'in', res = 600)
#png("BF_Manhattanplots_chr7-12_250519.png", width = 10, height = 6, units = 'in', res = 600)
par(mfrow = c(3, 2))
for (i in 7:12){
  # for each chromosome, keep only SNPs with XtX >= 15
  datacurrentchr=subset(realdata2,realdata2$chr==paste("Chr",i,sep="")&realdata2$xtx>=15&(realdata2$BFMeanToc>15|realdata2$BFMeanPluvio>15))
  # then plot these SNPs
  plot(datacurrentchr$BFMeanPluvio~datacurrentchr$pos.1,main=paste("Chromosome ",i),ylim=c(15,max(c(max(datacurrentchr$BFMeanToc),max(datacurrentchr$BFMeanPluvio)))),cex.main=1.5,cex.axis=1.2,cex.lab=1.4,pch=20,xlab=paste ("position - chr", i), ylab="BF",cex=1, col=ifelse(datacurrentchr$xtx >27.89115 & datacurrentchr$BFMeanPluvio > 20,"darkblue",ifelse(datacurrentchr$xtx>20.22704 & datacurrentchr$BFMeanPluvio > 15,"cornflowerblue","grey")))
  points(datacurrentchr$BFMeanToc~datacurrentchr$pos.1, pch=20, col=ifelse(datacurrentchr$xtx >27.89115 & datacurrentchr$BFMeanToc > 20,"red",ifelse(datacurrentchr$xtx>20.22704 & datacurrentchr$BFMeanToc > 15,"orange","grey")))
}
dev.off()
