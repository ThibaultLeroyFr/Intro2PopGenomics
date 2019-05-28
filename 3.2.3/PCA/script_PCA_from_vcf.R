# TL - 140219

setwd("~/Collab/AfricanRice/")

# DANS CONSOLE
#wget --no-check-certificate https://github.com/zhengxwen/gdsfmt/tarball/master -O gdsfmt_latest.tar.gz
#wget --no-check-certificate https://github.com/zhengxwen/SNPRelate/tarball/master -O SNPRelate_latest.tar.gz
#R CMD INSTALL gdsfmt_latest.tar.gz
#R CMD INSTALL SNPRelate_latest.tar.gz
library(gdsfmt)
library(SNPRelate)



########## SET 1 #####""
# vcf
#vcf.fn <- "/home/thibaultleroy/Capture/Zobo/merged_joint_bwa_mem_mdup_raw.filtered.vcf.2000SNPs.withheader"
vcf.fn <- "merged_joint_bwa_mem_mdup_raw.filtered.vcf.PASS.withoutNA.3pc.withheader"


# VCF => GDS
snpgdsVCF2GDS(vcf.fn, "merged_joint_bwa_mem_mdup_raw.filtered.vcf.PASS.withoutNA.3pc.withheader.gds", method="biallelic.only")

snpgdsSummary("merged_joint_bwa_mem_mdup_raw.filtered.vcf.PASS.withoutNA.3pc.withheader.gds")

# Open the GDS file
genofile <- snpgdsOpen("merged_joint_bwa_mem_mdup_raw.filtered.vcf.PASS.withoutNA.3pc.withheader.gds")

#pop_code <- scan("bidon", what=character())
#pop_code <- read.gdsn(index.gdsn(genofile, path="test.gds"))

#table(pop_code)

#snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2) #0.2
#snpset.id <- unlist(snpset)

#pca <- snpgdsPCA(genofile, snp.id=genofile$snp.id, num.thread=2)
pca <- snpgdsPCA(genofile,autosome.only=FALSE)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))

tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  EV3 = pca$eigenvect[,3],    # the second eigenvector
                  EV4 = pca$eigenvect[,4],    # the second eigenvector
                  stringsAsFactors = FALSE)

print(tab)

pdf(paste("merged_joint_bwa_mem_mdup_raw.filtered.vcf.PASS.withoutNA.3pc.withheader.","PCA.pdf",sep=""),12,10)
par(mfrow=c(2,2))


plot(tab$EV1,tab$EV2, cex=1.4,xlab=paste("set1 PC1 - ",round(pc.percent[1],3),"%",sep=""), ylab=paste("set1 PC2 - ",round(pc.percent[2],3),"%",sep=""),pch=20,col=ifelse(substr(tab$sample.id,1,2)=="Ob","red",ifelse(substr(tab$sample.id,1,2)=="Og","navyblue","forestgreen")))
#text(tab$EV1, tab$EV2, labels = tab$sample.id, pos = 4,cex=0.8)

plot(tab$EV1,tab$EV3, cex=1.4,xlab=paste("set1 PC1 - ",round(pc.percent[1],3),"%",sep=""), ylab=paste("set1 PC3 - ",round(pc.percent[3],3),"%",sep=""),pch=20,col=ifelse(substr(tab$sample.id,1,2)=="Ob","red",ifelse(substr(tab$sample.id,1,2)=="Og","navyblue","forestgreen")))
#text(tab$EV1, tab$EV3, labels = tab$sample.id, pos = 4,cex=0.8)

plot(tab$EV1,tab$EV4, cex=1.4,xlab=paste("set1 PC1 - ",round(pc.percent[1],4),"%",sep=""), ylab=paste("set1 PC4 - ",round(pc.percent[4],3),"%",sep=""),pch=20,col=ifelse(substr(tab$sample.id,1,2)=="Ob","red",ifelse(substr(tab$sample.id,1,2)=="Og","navyblue","forestgreen")))
#text(tab$EV1, tab$EV4, labels = tab$sample.id, pos = 4,cex=0.8)

plot(tab$EV2,tab$EV3, cex=1.4,xlab=paste("set1 PC2 - ",round(pc.percent[2],4),"%",sep=""), ylab=paste("set1 PC3 - ",round(pc.percent[3],3),"%",sep=""),pch=20,col=ifelse(substr(tab$sample.id,1,2)=="Ob","red",ifelse(substr(tab$sample.id,1,2)=="Og","navyblue","forestgreen")))
#text(tab$EV2, tab$EV3, labels = tab$sample.id, pos = 4,cex=0.8)


## with annot
plot(tab$EV1,tab$EV2, xlab=paste("set1 PC1 - ",round(pc.percent[1],3),"%",sep=""), ylab=paste("set1 PC2 - ",round(pc.percent[2],3),"%",sep=""),pch=20)
text(tab$EV1, tab$EV2, labels = tab$sample.id, pos = 4,cex=0.8)

plot(tab$EV1,tab$EV3, xlab=paste("set1 PC1 - ",round(pc.percent[1],3),"%",sep=""), ylab=paste("set1 PC3 - ",round(pc.percent[3],3),"%",sep=""),pch=20)
text(tab$EV1, tab$EV3, labels = tab$sample.id, pos = 4,cex=0.8)

plot(tab$EV1,tab$EV4, xlab=paste("set1 PC1 - ",round(pc.percent[1],4),"%",sep=""), ylab=paste("set1 PC4 - ",round(pc.percent[4],3),"%",sep=""),pch=20)
text(tab$EV1, tab$EV4, labels = tab$sample.id, pos = 4,cex=0.8)

plot(tab$EV2,tab$EV3, xlab=paste("set1 PC2 - ",round(pc.percent[2],4),"%",sep=""), ylab=paste("set1 PC3 - ",round(pc.percent[3],3),"%",sep=""),pch=20)
text(tab$EV2, tab$EV3, labels = tab$sample.id, pos = 4,cex=0.8)

#dev.off()
# tab without Zobo22, 14, 7 et 8

#dev.off()
write.table(tab,file="merged_joint_bwa_mem_mdup_raw.filtered.vcf.PASS.withoutNA.3pc.withheader.PCAcoord",sep="\t",quote=FALSE)
dev.off()
