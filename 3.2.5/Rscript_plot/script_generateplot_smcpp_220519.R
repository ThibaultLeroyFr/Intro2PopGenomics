# TL - 220519
library("ggplot2")
library("scales") # install.packages("scales")
setwd("/home/thibaultleroy/Collab/AfricanRice/smcpp_quentin")
popsize=read.table("plot_generation.csv",sep=",",h=T)
head(popsize)

### generate plot
tiff("smcpp_barthii_glaberrima.tiff", width = 10, height = 6, units = 'in', res = 600)
#png("smcpp_barthii_glaberrima.png", width = 10, height = 6, units = 'in', res = 600)
#jpeg("smcpp_barthii_glaberrima.jpeg", width = 10, height = 6, units = 'in', res = 600)

ggplot(data=popsize,aes(x=x,y=y,color=label)) + geom_line(size=1.4)+ xlab("Time (in years before present)") + ylab("Ne")+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
  theme_bw()+ annotation_logticks() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ 
  theme(legend.title=element_text(size=0, face="bold"), legend.title.align=0.5, legend.text=element_text(size=10, face="italic"), axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=16,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=16,angle=90,hjust=.5,vjust=.5,face="italic"))

dev.off()