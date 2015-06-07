#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)

png(args[2])

data = read.table(args[1])

ylimit = quantile(data$V5, 0.98)
if ( length(args) == 3 )
    ylimit = as.numeric(args[3])


x = data$V4
y = data$V5
plot(x, y, 
     main = "GC content vs Coverage",    
     xlab = "GC content",
     ylab = "Coverage",
     col=rgb(0,100,0,50,maxColorValue=255), 
     pch=16,
     ylim = c(0, ylimit))

dev.off()
