# plus d'infos: https://cran.r-project.org/web/packages/circlize/vignettes/circlize.pdf
#install.packages("circlize")

# files
setwd("/home/thibaultleroy/Collab/AfricanRice/")
FileWithpiValues=read.table("Genomic_pi_ROD_010419.withoutNA.txt",header=TRUE)


### Colors
col_uniq =rep(c("darkgrey", "black"), 6)
col =rep(c("orangered", "red"), 6)
col2 = rep(c("steelblue3","steelblue4"),6)

###

library("circlize")
circos.clear()
par(mar=c(1,1,1,1), lwd=0.1, cex =0.7)
circos.par("track.height"=0.15,start.degree=90,gap.degree=4)
circos.initialize(factors = FileWithpiValues$chr, x = FileWithpiValues$medianpos)

### first track : GC 
circos.trackPlotRegion(factors=FileWithpiValues$chr, y=FileWithpiValues$GC_Oryza_barthii, bg.border = NA, track.height=0.13, panel.fun=function(medianpos,GC_Oryza_barthii){
  circos.axis(h=0.52, labels.facing="inside", col="navyblue",labels.col="navyblue",lwd=2,labels.cex=0.9,major.at=c(100000,10000000,20000000,30000000,40000000,50000000),labels=c(0,"10Mb","20Mb","30Mb","40Mb","50Mb"))
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.index = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), max(ylim)+0.15, sector.index, facing = "inside", cex = 1.4,font=2)
  circos.yaxis("left", lwd=1,labels.cex=0.8)
})

circos.trackLines(FileWithpiValues$chr, FileWithpiValues$medianpos, FileWithpiValues$GC_Oryza_barthii, col=col_uniq, border = "black",lwd=0.9)

###  secondtrack: pi values
circos.trackPlotRegion(factors=FileWithpiValues$chr, y=FileWithpiValues$pi_site_Oryza_barthii, bg.border = NA, track.height=0.13, panel.fun=function(medianpos,pi_site_Oryza_barthii){
  #circos.axis(h=-0.005, labels.facing="inside", lwd=2,labels.cex=0.6)
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.index = get.cell.meta.data("sector.index")
  #circos.text(mean(xlim), max(ylim), sector.index, facing = "inside", cex = 1)
  circos.yaxis("left", lwd=1,labels.cex=0.8,at=c(0.0001,0.005,0.01,0.015),labels=c(0,0.005,0.01,0.015))
})
circos.trackLines(FileWithpiValues$chr, FileWithpiValues$medianpos, FileWithpiValues$pi_site_Oryza_barthii, col=col, border = "black",lwd=0.9)
circos.trackLines(FileWithpiValues$chr, FileWithpiValues$medianpos, FileWithpiValues$pi_site_Oryza_glaberrima, col=col2, border = "black",lwd=1.1)



### third track: ROD = Reduction of diversity [= 1 - (pi Oryza barthii / pî Oryza glaberrima)]
circos.trackPlotRegion(factors=FileWithpiValues$chr, y=FileWithpiValues$RODmorethanMinus1, bg.border = NA, track.height=0.15, panel.fun=function(medianpos,RODmorethanMinus1){
  #circos.axis(h=-1.2, labels.facing="inside", lwd=2,labels.cex=0.6)
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.index = get.cell.meta.data("sector.index")
  #circos.text(mean(xlim), max(ylim)+0.4, sector.index, facing = "inside", cex = 1.2,font=2)
  circos.yaxis("left", lwd=1,labels.cex=0.8,at=c(-0.9,0,0.9),labels=c(-0.9,0,0.9))
})
circos.trackLines(FileWithpiValues$chr, FileWithpiValues$medianpos, FileWithpiValues$RODmorethanMinus1, type = "h",baseline=0,col=ifelse(FileWithpiValues$RODmorethanMinus1>0.8,"red",ifelse(FileWithpiValues$RODmorethanMinus1< 0,'green3','darkgoldenrod2')),lwd=ifelse(FileWithpiValues$RODmorethanMinus1>0.8,0.6,0.35))



### forth track: Tajima's D Oryza glaberrima
# median(FileWithpiValues$D_Pop_Oryza_glaberrima) = 0.170747
# variation +2 / -2 around median value => 2.170747 / -1.829253

circos.trackPlotRegion(factors=FileWithpiValues$chr, y=FileWithpiValues$D_Pop_Oryza_glaberrima, bg.border = NA, track.height=0.15, panel.fun=function(medianpos,D_Pop_Oryza_glaberrima){
  #circos.axis(h=-4, labels.facing="inside", lwd=2,labels.cex=0.8)
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.index = get.cell.meta.data("sector.index")
  #circos.text(mean(xlim), max(ylim)+0.4, sector.index, facing = "inside", cex = 1.2,font=2)
  circos.yaxis("left", lwd=1,labels.cex=0.8) #,at=c(-0.9,0,0.9),labels=c(-0.9,0,0.9))
})
circos.trackLines(FileWithpiValues$chr, FileWithpiValues$medianpos, FileWithpiValues$D_Pop_Oryza_glaberrima, type = "h",baseline=0.170747, col=ifelse(FileWithpiValues$D_Pop_Oryza_glaberrima < -1.829253,"red",ifelse(FileWithpiValues$D_Pop_Oryza_glaberrima> 2.170747,"purple","grey")),lwd=ifelse(FileWithpiValues$D_Pop_Oryza_glaberrima < -1.829253,0.6,0.35))

