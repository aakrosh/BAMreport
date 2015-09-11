#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)

pdf(args[2], width = 11.69, height = 8.27)

par(mfrow = c(1,1), oma = c(0, 0, 3, 0))

data = read.table(args[1])

# for the first read
plotdata = data[which(data$V1 == "1"),]
mp = barplot(plotdata$V3,
        xlab = "5' Clip position",
        ylab = "Fraction of read1s",
        axes = F,
        main = "Read 1")

index = seq(0, max(plotdata$V2), 10)
index[1] = 1

axis(1, at = mp[index], labels = index)
axis(2)

mtext("Clipping position distribution", outer = TRUE, cex = 1.5)

dev.off()
