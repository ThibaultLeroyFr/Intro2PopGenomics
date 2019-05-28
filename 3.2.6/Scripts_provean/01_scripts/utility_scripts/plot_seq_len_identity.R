#!/usr/bin/Rscript

data <- read.table("04_blasts/data_to_plot")

pdf(file="distribution_plot.pdf", 20,15)
plot(data,pch=16, col=rgb(0,0,1, 0.2), xlab="sequence identity", ylab="seq length")
dev.off()
