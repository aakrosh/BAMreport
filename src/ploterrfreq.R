#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)

pdf(args[2], width = 11.69, height = 8.27)

par(mfrow = c(1,2), oma = c(0, 0, 3, 0))

data = read.table(args[1])

# for the first read
plotdata = data[which(data$V1 == "1"),]
mp = barplot(plotdata$V3,
        xlab = "Position on the read",
        ylab = "Fraction of read1s with errors",
        axes = F,
        main = "Read 1")

index = seq(0, max(plotdata$V2), 10)
index[1] = 1

axis(1, at = mp[index], labels = index)
axis(2)

# for the second read
plotdata = data[which(data$V1 == "2"),]
mp = barplot(plotdata$V3,
        xlab = "Position on the read",
        ylab = "Fraction of read2s with errors",
        axes = F,
        main = "Read 2")

index = seq(0, max(plotdata$V2), 10)
index[1] = 1

axis(1, at = mp[index], labels = index)
axis(2)

mtext("Error rate with read position", outer = TRUE, cex = 1.5)

dev.off()
