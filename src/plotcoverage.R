#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)

pdf(args[2], width = 11.69, height = 8.27)

data = read.table(args[1])
datax = data[-1,]
xlimit = quantile(rep(datax$V1, datax$V2), 0.98)
bases = sum(as.numeric(data$V2))
data$V2 = data$V2 / bases
data = data[-1,]

if ( length(args) == 3 )
    xlimit = as.numeric(args[3])

plot(data,
     type = "l",
     lwd  = 3,
     col  = rgb(255,0,0,150,maxColorValue=255),
     main = "Depth of coverage",    
     xlab = "Coverage",
     ylab = "Fraction of genome covered",
     xlim = c(0, xlimit))

dev.off()
