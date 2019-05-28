library(ggplot2)
library(cowplot)
library("gridExtra")
setwd("/home/thibault/Collab/Blé/dataset_140617/Structure_2k_10k")

#png("barplots_K2",1200,900)
par(mfrow = c(3, 1),     # pour 2 lignes, c(3,1) pour 3 lignes
    oma = c(2, 2, 2, 2), #, # two rows of text at the outer left and bottom margin #  c(bottom, left, top, right) 
    mar = c(2, 5, 1, 2)+0.1) # A numerical vector of the form c(bottom, left, top, right) which gives the number of lines of margin to be 

admixK2=read.table("merged_joint_bwa_mem_mdup_raw.filtered.vcf.PASS.withoutNA.3pc.withheader.K.2.meanQ.withheader.sed", header=FALSE)
admixK3=read.table("merged_joint_bwa_mem_mdup_raw.filtered.vcf.PASS.withoutNA.3pc.withheader.K.3.meanQ.withheader.sed", header=FALSE)
admixK4=read.table("merged_joint_bwa_mem_mdup_raw.filtered.vcf.PASS.withoutNA.3pc.withheader.K.4.meanQ.withheader.sed", header=FALSE)
admixK5=read.table("merged_joint_bwa_mem_mdup_raw.filtered.vcf.PASS.withoutNA.3pc.withheader.K.5.meanQ.withheader.sed", header=FALSE)

a1 <- ggplot(admixK2, aes(x=V3,y=V7, fill=V6)) +  geom_bar(width=1, stat="identity") + 
  scale_fill_manual(values = c("#ffe458","#cb745b"))+ 
  xlab("")+ ylab("Ancestry [K=2]")+theme(legend.position="none")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 0), 
        axis.text.x = element_text(colour="black",size=0,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))+
#  annotate("text", x = 540, y = 0.5, label="K=2",size=8, fontface =2)+
  # O. barthii
  geom_segment(x = 23.5, y = -0.1, xend = 23.5, yend = 1, colour = "black")+
  geom_segment(x = 20.5, y = 0, xend = 20.5, yend = 1, colour = "black",lty=2)+
  # O. glaberrima
  geom_segment(x = 28.5, y = 0, xend = 28.5, yend = 1, colour = "black",lty=2)+
  geom_segment(x = 38.5, y = 0, xend = 38.5, yend = 1, colour = "black",lty=2)+
  # annotations
  annotate("text", x = 10, y = -0.05, label="O.barthii",size=6, fontface =3,colour = " black")+
  annotate("text", x = 37, y = -0.05, label="O.glaberrima",size=6, fontface =3,colour = " black")+
  annotate("text", x = 10, y = 0.05, label="Mali",size=4, fontface =1,colour = "white")+
  annotate("text", x = 22, y = 0.05, label="Nigeria",size=4, fontface =1,colour = "white")+
  annotate("text", x = 26, y = 0.05, label="Ghana",size=4, fontface =1,colour = "white")+
  annotate("text", x = 34, y = 0.05, label="Mali",size=4, fontface =1,colour = "white")+
  annotate("text", x = 44, y = 0.05, label="Nigeria",size=4, fontface =1,colour = "white")

a2 <- ggplot(admixK3, aes(x=V3,y=V7, fill=V6)) +  geom_bar(width=1, stat="identity") + 
  scale_fill_manual(values = c("darkgreen","#cb745b","#ffe458"))+
  xlab("")+ ylab("Ancestry [K=3]")+theme(legend.position="none")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 0), 
        axis.text.x = element_text(colour="black",size=0,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))+
  #  annotate("text", x = 540, y = 0.5, label="K=2",size=8, fontface =2)+
  # O. barthii
  geom_segment(x = 23.5, y = -0.1, xend = 23.5, yend = 1, colour = "black")+
  geom_segment(x = 20.5, y = 0, xend = 20.5, yend = 1, colour = "black",lty=2)+
  # O. glaberrima
  geom_segment(x = 28.5, y = 0, xend = 28.5, yend = 1, colour = "black",lty=2)+
  geom_segment(x = 38.5, y = 0, xend = 38.5, yend = 1, colour = "black",lty=2)+
  # annotations
  annotate("text", x = 10, y = -0.05, label="O.barthii",size=6, fontface =3,colour = " black")+
  annotate("text", x = 37, y = -0.05, label="O.glaberrima",size=6, fontface =3,colour = " black")+
  annotate("text", x = 10, y = 0.05, label="Mali",size=4, fontface =1,colour = "white")+
  annotate("text", x = 22, y = 0.05, label="Nigeria",size=4, fontface =1,colour = "white")+
  annotate("text", x = 26, y = 0.05, label="Ghana",size=4, fontface =1,colour = "white")+
  annotate("text", x = 34, y = 0.05, label="Mali",size=4, fontface =1,colour = "white")+
  annotate("text", x = 44, y = 0.05, label="Nigeria",size=4, fontface =1,colour = "white")
  
  
  

a3 <- ggplot(admixK4, aes(x=V3,y=V7, fill=V6)) +  geom_bar(width=1, stat="identity") + 
  scale_fill_manual(values = c("#ffe458","cornflowerblue","#cb745b","darkgreen"))+
  xlab("")+ ylab("Ancestry [K=4]")+theme(legend.position="none")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 0), 
        axis.text.x = element_text(colour="black",size=0,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))+
  #  annotate("text", x = 540, y = 0.5, label="K=2",size=8, fontface =2)+
  # O. barthii
  geom_segment(x = 23.5, y = -0.1, xend = 23.5, yend = 1, colour = "black")+
  geom_segment(x = 20.5, y = 0, xend = 20.5, yend = 1, colour = "black",lty=2)+
  # O. glaberrima
  geom_segment(x = 28.5, y = 0, xend = 28.5, yend = 1, colour = "black",lty=2)+
  geom_segment(x = 38.5, y = 0, xend = 38.5, yend = 1, colour = "black",lty=2)+
  # annotations
  annotate("text", x = 10, y = -0.05, label="O.barthii",size=6, fontface =3,colour = " black")+
  annotate("text", x = 37, y = -0.05, label="O.glaberrima",size=6, fontface =3,colour = " black")+
  annotate("text", x = 10, y = 0.05, label="Mali",size=4, fontface =1,colour = "white")+
  annotate("text", x = 22, y = 0.05, label="Nigeria",size=4, fontface =1,colour = "white")+
  annotate("text", x = 26, y = 0.05, label="Ghana",size=4, fontface =1,colour = "white")+
  annotate("text", x = 34, y = 0.05, label="Mali",size=4, fontface =1,colour = "white")+
  annotate("text", x = 44, y = 0.05, label="Nigeria",size=4, fontface =1,colour = "white")
  
  
a4 <- ggplot(admixK5, aes(x=V3,y=V7, fill=V6)) +  geom_bar(width=1, stat="identity") + 
  scale_fill_manual(values = c("deeppink4","cornflowerblue","#ffe458","#cb745b","darkgreen"))+
  xlab("")+ ylab("Ancestry [K=5]")+theme(legend.position="none")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 0), 
        axis.text.x = element_text(colour="black",size=0,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))+
  #  annotate("text", x = 540, y = 0.5, label="K=2",size=8, fontface =2)+
  # O. barthii
  geom_segment(x = 23.5, y = -0.1, xend = 23.5, yend = 1, colour = "black")+
  geom_segment(x = 20.5, y = 0, xend = 20.5, yend = 1, colour = "black",lty=2)+
  # O. glaberrima
  geom_segment(x = 28.5, y = 0, xend = 28.5, yend = 1, colour = "black",lty=2)+
  geom_segment(x = 38.5, y = 0, xend = 38.5, yend = 1, colour = "black",lty=2)+
  # annotations
  annotate("text", x = 10, y = -0.05, label="O.barthii",size=6, fontface =3,colour = " black")+
  annotate("text", x = 37, y = -0.05, label="O.glaberrima",size=6, fontface =3,colour = " black")+
  annotate("text", x = 10, y = 0.05, label="Mali",size=4, fontface =1,colour = "white")+
  annotate("text", x = 22, y = 0.05, label="Nigeria",size=4, fontface =1,colour = "white")+
  annotate("text", x = 26, y = 0.05, label="Ghana",size=4, fontface =1,colour = "white")+
  annotate("text", x = 34, y = 0.05, label="Mali",size=4, fontface =1,colour = "white")+
  annotate("text", x = 44, y = 0.05, label="Nigeria",size=4, fontface =1,colour = "white")
  
grid.arrange(a1,a2, a3,a4,ncol=1)

grid.arrange(a1,a2,ncol=1)