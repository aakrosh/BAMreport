#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)

pdf(args[2], width = 11.69, height = 8.27)

par(mfrow = c(1,2), oma = c(0, 0, 3, 0))

# create the boxplot for the quality values
data = read.table(args[1])

# for the first read
plotdata = data[which(data$V1 == "1"),]
plotdata = plotdata[,3:7]

conf = matrix(0, nrow = 2, ncol = nrow(plotdata))
conf[1,] = 0
conf[2,] = 60

z = list(stats = t(as.matrix(plotdata)), 
         n     = rep(100, nrow(plotdata)), 
         conf  = conf, 
         out   = vector(), 
         group = vector(), 
         names = seq(1, nrow(plotdata)))

bxp(z, 
    outline = F, 
    xlab = "Position on the read",
    ylab = "Quality value",
    axes = F,
    boxfill = rgb(255,0,0,150,maxColorValue=255), 
    whisklty = 3,
    ylim = c(0, 60),
    main = "Read 1"
    )

axis(1, at = seq(1, nrow(plotdata), 5), labels = seq(1, nrow(plotdata), 5))
axis(2, at = seq(0, 60, 10), labels = seq(0, 60, 10))
#axis(4, at = seq(0, 60, 10), labels = seq(0, 60, 10))

# for the second read
plotdata = data[which(data$V1 == "2"),]
plotdata = plotdata[,3:7]

conf = matrix(0, nrow = 2, ncol = nrow(plotdata))
conf[1,] = 0
conf[2,] = 60

z = list(stats = t(as.matrix(plotdata)), 
         n     = rep(100, nrow(plotdata)), 
         conf  = conf, 
         out   = vector(), 
         group = vector(), 
         names = seq(1, nrow(plotdata)))

bxp(z, 
    outline = F, 
    xlab = "Position on the read",
    ylab = "Quality value",
    axes = F,
    boxfill = rgb(255,0,0,150,maxColorValue=255), 
    whisklty = 3,
    ylim = c(0, 60),
    main = "Read 2"
    )

axis(1, at = seq(1, nrow(plotdata), 5), labels = seq(1, nrow(plotdata), 5))
axis(2, at = seq(0, 60, 10), labels = seq(0, 60, 10))
#axis(4, at = seq(0, 60, 10), labels = seq(0, 60, 10))


mtext("Quality value variation", outer = TRUE, cex = 1.5)

dev.off()
