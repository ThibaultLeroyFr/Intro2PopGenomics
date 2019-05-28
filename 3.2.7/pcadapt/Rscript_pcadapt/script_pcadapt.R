library("pcadapt")
setwd("/media/thibaultleroy/Seagate Expansion Drive/AfricanRice/")

# use Plink to convert vcf to bed
# e.g. ./plink --allow-extra-chr --make-bed --noweb --out /media/thibaultleroy/Seagate\ Expansion\ Drive/AfricanRice/merged_joint_bwa_mem_mdup_raw.filtered.vcf.PASS.withheader.bed --vcf /media/thibaultleroy/Seagate\ Expansion\ Drive/AfricanRice/merged_joint_bwa_mem_mdup_raw.filtered.vcf.PASS.withheader.vcf

filename <- read.pcadapt("merged_joint_bwa_mem_mdup_raw.filtered.vcf.PASS.withheader.bed.bed",type="bed")

# define populations
poplist.int <- c(rep(1, 18), rep(2, 25), rep(1, 5)) # 1 = O.barthii individuals, 2 = O. glaberrima individuals
poplist.names <- c(rep("O.barthii", 18),rep("O.glaberrima", 25),rep("O.barthii", 5))


### number of components
x <- pcadapt(input = filename, K = 20) 
# screeplot
plot(x, option = "screeplot")
# here, most of the variation is explained by the first four components

### PCA
# show axis 1 vs 2
plot(x, option = "scores", pop = poplist.names)
# show axis 1 vs 3 (feel free to change i and j)
plot(x, option = "scores", i = 1, j = 3, pop = poplist.names)

### Genome scan considering K=2
x <- pcadapt(filename, K = 2)
# manhattan plot
plot(x , option = "manhattan")
plot(x, option = "qqplot")
hist(x$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
plot(x, option = "stat.distribution")

### Based on the scree plot (first four components), we can perform a genome scan considering K=4
x <- pcadapt(filename, K = 4)
par(mfrow = c(2, 2))
for (i in 1:4){
  plot(x$loadings[, i], pch = 19, cex = .3, ylab = paste0("Loadings PC", i))
}
# Observed: some clear clusters of SNPs that are associated => Effect of linkage disequilibrium => require SNP pruning

### SNP prunning (here 2k SNP)
res <- pcadapt(filename, K = 20, LD.clumping = list(size = 2000, thr = 0.1))
plot(res, option = "screeplot")
### check if the plot of loading looks better
for (i in 1:4){
  plot(x$loadings[, i], pch = 19, cex = .3, ylab = paste0("Loadings PC", i))
}
# Observed: now it works better (neighboring SNPs cluster less)

# analysis for K=3 after LD thinning
x <- pcadapt(filename, K = 3, LD.clumping = list(size = 2000, thr = 0.1))
plot(x, option = "qqplot")
par(mfrow = c(1, 1))
# generate a manhattan plot, highlighting SNPs with pvalue < 0.0001 (red) or 0.01 (orange)
plot(x , option = "manhattan",col=ifelse(x$pvalues<0.00001,"red",ifelse(x$pvalues<0.01,"orange","darkgrey")))
write.table(x$pvalues,file="pvalues_SNP.txt",na="NA",quote=FALSE,sep="\t")


###### 3 main plots (PCA + screeplot)
# PCA axis 1 vs 2
png("PCA_pcadapt_axis1-2.png", width = 5, height = 4, units = 'in', res = 600)
#tiff("PCA_pcadapt_axis1-2", width = 5, height = 4, units = 'in', res = 600)
plot(res, option = "scores", pop = poplist.names)
dev.off()

# PCA axis 1 vs 3 (feel free to change i and j)
png("PCA_pcadapt_axis1-3.png", width = 5, height = 4, units = 'in', res = 600)
#tiff("PCA_pcadapt_axis1-3", width = 5, height = 4, units = 'in', res = 600)
plot(res, option = "scores", i = 1, j = 3, pop = poplist.names)
dev.off()

png("screeplot_pcadapt.png", width = 5, height = 4, units = 'in', res = 600)
#tiff("screeplot_pcadapt", width = 8, height = 4, units = 'in', res = 600)
plot(res, option = "screeplot")
dev.off()


### Manhattan plots

# merge the file containing pvalues (pcadapt) and genome position of each SNP using the initial vcf (the two files need to have the same number of lines)
#awk '{print $1"  "$2}' merged_joint_bwa_mem_mdup_raw.filtered.vcf.PASS > chrinfo_SNP
# paste chrinfo_SNP pvalues_SNP.txt > pvalues_withchrinfo.txt

myfile=read.table("pvalues_withchrinfo.txt",h=F,sep="\t")
#
#
#

#tiff("PCAdapt_scan_chr1-6.tiff", width = 6, height = 8, units = 'in', res = 600)
png("PCAdapt_scan_chr1-6.png", width = 6, height = 8, units = 'in', res = 600)
par(mfrow = c(3, 2))
for (i in 1:6){
  # substract data for each chromosome
  datacurrentchr=subset(myfile,myfile$V1==i)
  plot(-log10(datacurrentchr$V4)~datacurrentchr$V2,main=paste("Chromosome ",i),pch=20,cex=0.1,xlab=paste ("position - chr", i), ylab="-log10(p-values)",col=ifelse(datacurrentchr$V4<0.00001,"red",ifelse(datacurrentchr$V4<0.01,"orange","darkgrey")))
}
dev.off()

#tiff("PCAdapt_scan_chr7-12.tiff", width = 6, height = 8, units = 'in', res = 600)
png("PCAdapt_scan_chr7-12.png", width = 6, height = 8, units = 'in', res = 600)
par(mfrow = c(3, 2))
for (i in 7:12){
  datacurrentchr=subset(myfile,myfile$V1==i)
  plot(-log10(datacurrentchr$V4)~datacurrentchr$V2,main=paste("Chromosome ",i),pch=20,cex=0.1,xlab=paste ("position - chr", i), ylab="-log10(p-values)",col=ifelse(datacurrentchr$V4<0.00001,"red",ifelse(datacurrentchr$V4<0.01,"orange","darkgrey")))
}
dev.off()
