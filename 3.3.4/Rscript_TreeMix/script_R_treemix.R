# TL - 140519
library(ape) # install.packages("ape") if you need to install it 

source('/home/thibaultleroy/Software/treemix-1.13/src/plotting_funcs.R')
setwd("/media/thibaultleroy/Seagate Expansion Drive/projet_BNP-BHW/")

### First tree (assuming no migration tree)
# variance explained

##### Explained variance over bootstraps
varianceexp = get_f("BNP-BHW_18pops_v2.3_freqSNP_mincount10_filterd_biallelic_sites_counts.outfiletreemixK1000m.0")
print(paste("The inferred tree explained", round(varianceexp,3)*100, "% of the total variance based on the variance-covariance matrix"))
# 0.8980568 in my case

###### 1/ draw a tree (without migration noded)
## To generate classic Treemix visualization, do:
plot_tree("BNP-BHW_18pops_v2.3_freqSNP_mincount10_filterd_biallelic_sites_counts.outfiletreemixK1000m.0") #filename before the ".treeout.gz"
## Show the residuals (remainding variance that is not explained by the tree)
#write a file with the pop names
pop_order=c("9","97","124","204","217","218","219","233","253","256","L1","O1","L8","O8","L12","O12","L16","O16")
write(pop_order,file="pop_order.txt")
# show the residuals for this tree, with pop in the desired order
plot_resid("BNP-BHW_18pops_v2.3_freqSNP_mincount10_filterd_biallelic_sites_counts.outfiletreemixK1000m.0","pop_order.txt")

### To redo the first treemix figure in the chapter: 
# first copy the tree in a text file
# bash CMD: zmore BNP-BHW_18pops_v2.3_freqSNP_mincount10_filterd_biallelic_sites_counts.outfiletreemixK1000m.0.treeout.gz > tree_M0.txt
# read the file
MyTree <- read.tree("tree_M0.txt")
MyTree<-ladderize(MyTree) 
# draw an unrooted tree + the residuals in a file
#tiff("~/Treemix_withoutmigrationnodes.tiff", width = 10, height = 6, units = 'in', res = 600)
png("~/Treemix_withoutmigrationnodes.png", width = 10, height = 6, units = 'in', res = 600)
par(mfrow=c(1,2))
plot(MyTree, edge.color="darkgrey",cex=0.8,edge.width=3,type="unrooted")
plot_resid("BNP-BHW_18pops_v2.3_freqSNP_mincount10_filterd_biallelic_sites_counts.outfiletreemixK1000m.0","pop_order.txt")
dev.off()


###### 2/ variance explained for simulations from m=0 to m=10 (migration edges)
varianceexp=NULL
for (i in 0:10){
  out <- paste("m",i, sep="")
  get_f(paste('BNP-BHW_18pops_v2.3_freqSNP_mincount10_filterd_biallelic_sites_counts.outfiletreemixK1000m.',i,sep="")) -> out
  varianceexp = rbind(varianceexp, as.numeric(cbind(i,out)))
}
print(varianceexp)
varianceexp=as.data.frame(varianceexp)

tiff("~/Treemix_varianceexplained.tiff", width = 10, height = 6, units = 'in', res = 600)
#png("~/Treemix_varianceexplained.png", width = 10, height = 6, units = 'in', res = 600)


# plot
ggplot(varianceexp, aes(x=V1, y=V2)) + geom_point(pch=20,col="navyblue",size=2) + geom_line(linetype="dashed") + xlab("migration nodes")+ ylab("% explained variance")+
  scale_x_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9,10))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ theme(legend.key = element_blank())+
  theme(legend.title=element_text(size=16, face="bold"), legend.title.align=0.5, legend.text=element_text(size=14, face="italic"), axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), axis.text.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=16,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=16,angle=90,hjust=.5,vjust=.5,face="italic"))

dev.off()

## tree m=1
par(mfrow=c(1,2))
plot_tree("BNP-BHW_18pops_v2.3_freqSNP_mincount10_filterd_biallelic_sites_counts.outfiletreemixK1000m.1") #filename before the ".treeout.gz"
plot_resid("BNP-BHW_18pops_v2.3_freqSNP_mincount10_filterd_biallelic_sites_counts.outfiletreemixK1000m.1","pop_order.txt")

## tree m=2
par(mfrow=c(1,2))
plot_tree("BNP-BHW_18pops_v2.3_freqSNP_mincount10_filterd_biallelic_sites_counts.outfiletreemixK1000m.2") #filename before the ".treeout.gz"
plot_resid("BNP-BHW_18pops_v2.3_freqSNP_mincount10_filterd_biallelic_sites_counts.outfiletreemixK1000m.2","pop_order.txt")
# and so on

### 

### To redo the second treemix figure in the chapter: 
# first copy the tree in a text file
# bash CMD: zmore BNP-BHW_18pops_v2.3_freqSNP_mincount10_filterd_biallelic_sites_counts.outfiletreemixK1000m.1.treeout.gz > tree_M1.txt
# read the file
MyTree <- read.tree("tree_M1.txt")
MyTree<-ladderize(MyTree) 
# draw an unrooted tree + the residuals in a file
#tiff("~/Treemix_withoutmigrationnodes.tiff", width = 10, height = 6, units = 'in', res = 600)
#png("~/Treemix_withoutmigrationnodes.png", width = 10, height = 6, units = 'in', res = 600)
par(mfrow=c(1,2))
plot(MyTree, edge.color="darkgrey",cex=0.8,edge.width=3,type="unrooted")
add.scale.bar(x=0,y=0,cex = 0.75, font = 2, lcol="darkgrey",col = "darkgrey",lwd=3)
text(0.001,0.0003,"Treemix drift parameter",cex=0.75,font=2,col="darkgrey")
#nodelabels(cex=0.65,bg="yellow")
#node O8-O12 = 30 ; node 233-ç7: 35

# add a scale for the migration weight
col2func <- colorRampPalette(c("yellow","red"))
col2function1000 <-col2func(1000)
for (i in 1:1000){
  dec=i*0.0000015
  rect(0, 0.004+dec, 0.0004, 0.00401+dec, lwd=0,col=col2func(1000)[i])
}
text(0.0002,0.006,"Migration\nweight",font=2, cex=0.75)
text(0.00058,0.00545,"0.5",cex=0.6)
text(0.00054,0.0041,"0",cex=0.6)

# add the corresponding migration node
MyMigNodes<-read.table("tree_M1.txt",skip=1,sep=(" ")) # skip=1 to only keep information related to migration nodes
migcolor=round(1000*2*(as.numeric(MyMigNodes[1,1])),0) # si 0 à 0.5 alors faire *2
# here the arrows was manually set because the tree is unrooted, but it is possible to use the matrix generated by "dist.nodes(MyTree)" to generate all distances from a root, which can be very useful to draw more easily the migration nodes (arrows)
#0.393964 NA NA NA (97:0.00170865,233:0.00202388):0.00030483 (O12:0.00267841,O8:0.00347784):0.000529587

arrows(x0=0.0041,y0=0.00458,x1=0.00292,y1=0.0027,col=col2function1000[migcolor], length=0.1) # xo = V5 # x1 (pointe fleche) = V6
plot_resid("BNP-BHW_18pops_v2.3_freqSNP_mincount10_filterd_biallelic_sites_counts.outfiletreemixK1000m.1","pop_order.txt")
#dev.off()