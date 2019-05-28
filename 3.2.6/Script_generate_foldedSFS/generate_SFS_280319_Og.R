# TL - 270319
library(ggplot2)
setwd("/home/tleroy/work2/Collab/AfricanRice/")
SFS=read.table("merged_joint_bwa_mem_mdup_raw.filtered.vcf.PASS.withheader.geno.Og.maf",h=F)

pdf(file="SFS_AfricanRice_Oglaberrima_individuals.pdf",height=8,width=6)

SFSpoly=subset(SFS,V8>0)
ggplot(SFSpoly, aes(x=V8)) + geom_histogram(color="black", fill="darkgreen", binwidth = 1/50) + xlab("Folded SFS (MAF - O. glaberrima - 25 individuals)")

ggplot(SFS, aes(x=V7)) + geom_histogram(color="black", fill="white")
ggplot(SFS, aes(x=V8)) + geom_histogram(color="black", fill="white")

dev.off()
