library("poolfstat") # to install: install.packages("poolfstat")
library(reshape2)
library(ggplot2)

setwd("/home/thibaultleroy/Collab/Sessile_Reanalyse")
# read the sync file (Pileup > Sync: mpileup2sync.pl script from Popoolation2)
# subsample the sync file
pooldata=popsync2pooldata(sync.file="/media/thibaultleroy/Seagate Expansion Drive/projet_BNP-BHW/PileupSync/BNP-BHW_18pops_v2.3.pileup.sync.1",poolsizes=c(50,50,50,50,50,44,50,50,50,50,40,40,40,40,40,40,20,36),poolnames=c("9","97","124","204","217","218","219","233","253","256","L1","O1","L8","O8","O12","L12","O16","L16"),min.rc=10,min.cov.per.pool = 50, max.cov.per.pool = 250, min.maf=0.02,noindel=TRUE)

# Compute a pairwise Fst matrix
Fstmatrix=computePairwiseFSTmatrix(pooldata, method = "Anova",min.cov.per.pool = 50, max.cov.per.pool = 250, min.maf = 0.02,output.snp.values = FALSE)$PairwiseFSTmatrix

# reshape the matrix to be used by ggplot
melted_Fstmatrix<- melt(Fstmatrix)
head(melted_Fstmatrix)

tiff("Fst_pairwisematrix.tiff", width = 10, height = 6, units = 'in', res = 600)
#png("Fst_pairwisematrix.png", width = 10, height = 6, units = 'in', res = 600)

# generate plot
ggplot(data = melted_Fstmatrix, aes(x=Var1, y=Var2, fill=value))+ geom_tile()+
scale_fill_gradient2("Pairwise\nFST values", low = "navyblue", mid = "white",high = "red",  midpoint = 0)+
  xlab("")+ylab("")+
  scale_x_discrete(labels=c("Pool1" = "9", "Pool2" = "97", "Pool3" = "124","Pool4" = "204", "Pool5" = "217", "Pool6" = "218","Pool7" = "219", "Pool8" = "233", "Pool9" = "253","Pool10" = "256", "Pool11" = "L1", "Pool12" = "O1","Pool13" = "L8", "Pool14" = "O8", "Pool15" = "O12","Pool16" = "L12", "Pool17" = "O16", "Pool18" = "L16"))+
  scale_y_discrete(labels=c("Pool1" = "9", "Pool2" = "97", "Pool3" = "124","Pool4" = "204", "Pool5" = "217", "Pool6" = "218","Pool7" = "219", "Pool8" = "233", "Pool9" = "253","Pool10" = "256", "Pool11" = "L1", "Pool12" = "O1","Pool13" = "L8", "Pool14" = "O8", "Pool15" = "O12","Pool16" = "L12", "Pool17" = "O16", "Pool18" = "L16"))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ theme(legend.key = element_blank())+
  theme(legend.title=element_text(size=12, face="bold"), legend.title.align=0.5, legend.text=element_text(size=10, face="italic"), axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=12,angle=90,hjust=.5,vjust=.5,face="italic"))

 dev.off()

# Among-population FST (Global FST values computed among all populations for each SNP)
SNPvalues=computeFST(pooldata, method="Anova")$snp.FST
SNPvalues=as.data.frame(SNPvalues)
plot(density(SNPvalues))
ggplot(data=SNPvalues,aes(x=SNPvalues)) + geom_density() + xlim(0,1) + xlab("Among-population Fst values")+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ theme(legend.key = element_blank())+
  theme(legend.title=element_text(size=12, face="bold"), legend.title.align=0.5, legend.text=element_text(size=10, face="italic"), axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=12,angle=90,hjust=.5,vjust=.5,face="italic"))

# pairwise Fst for all locus
Fstlocus=computePairwiseFSTmatrix(pooldata, method = "Anova",min.cov.per.pool = 50, max.cov.per.pool = 250, min.maf = 0.02,output.snp.values = TRUE)$PairwiseSnpFST
# reshape the matrix to be used by ggplot
melted_Fstlocus <- melt(Fstlocus)
head(melted_Fstlocus)
ggplot(data=melted_Fstlocus,aes(x=value,color=Var2)) + geom_density(size=0.1) +  xlim(-0.1,1) + xlab("Pairwise Fst values")+
  scale_color_manual(values = c(rep("grey",153)))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ theme(legend.position = "none")+
  theme(legend.title=element_text(size=12, face="bold"), legend.title.align=0.5, legend.text=element_text(size=10, face="italic"), axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=12,angle=90,hjust=.5,vjust=.5,face="italic"))


### 
tiff("Fst_distribution_pairwise_among-pop.tiff", width = 10, height = 6, units = 'in', res = 600)
#png("Fst_distribution_pairwise_among-pop.png", width = 10, height = 6, units = 'in', res = 600)

ggplot(data=melted_Fstlocus,aes(x=value,color=Var2)) + geom_density(size=0.1) + geom_density(data=SNPvalues,aes(x=SNPvalues),col="black") + xlim(-0.1,1) + xlab("Pairwise (grey) and among-population (black) Fst values")+
  scale_color_manual(values = c(rep("grey",153)))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ theme(legend.position = "none")+
  theme(legend.title=element_text(size=12, face="bold"), legend.title.align=0.5, legend.text=element_text(size=10, face="italic"), axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=12,angle=90,hjust=.5,vjust=.5,face="italic"))

dev.off()