#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)

pdf(args[2], width = 11.69, height = 8.27)

data = read.table(args[1])

plot(data,
     type = "l",
     lwd  = 3,
     col  = rgb(255,0,0,150,maxColorValue=255),
     main = "Insert size distribution",    
     xlab = "Insert size",
     ylab = "Fraction of properly-paired pairs",
     xlim = c(0, quantile(rep(data$V1,data$V2*1000),.98)))

dev.off()
