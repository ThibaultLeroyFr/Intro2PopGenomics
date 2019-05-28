# TL - 110519
library(ggplot2)
setwd("/home/thibaultleroy/Collab/Sessile_Reanalyse")
ratios=read.table("Poolseq_vs_individual_100xpoolvs20xNbindividuals.sed",sep="\t",h=T)

#tiff("Poolseq_vs_individual_sequencing_110519.tiff", width = 9, height = 6, units = 'in', res = 600)
png("Poolseq_vs_individual_sequencing_110519.png", width = 9, height = 6, units = 'in', res = 600)

ggplot(ratios, aes(x=Number_of_individuals,y=Ratio_No_ExpError,color=Biases)) + geom_line(size=1.5) +  ylim(0,2) + xlim(0,100)+
  xlab("Number of samples (individually sequenced at 20X)") +  ylab("SD pool-seq / SD individally-sequenced")+
  scale_color_manual(values=c("#E69F00", "#56B4E9")) +
  geom_hline(yintercept = 1,colour="black",linetype="solid") +
  geom_segment(aes(x= 26.3,xend=26.3, y=0,yend=1),colour="#E69F00",linetype="dashed") +
  geom_segment(aes(x= 23.3,xend=23.3, y=0,yend=1),colour="#56B4E9",linetype="dashed") +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  annotate("text", x = 70, y = 0.95, label = "Below this line: the pool-seq strategy performs better")+
  annotate("text", x = 70, y = 1.05, label = "Above this line: the individual strategy performs better")+
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))

dev.off()