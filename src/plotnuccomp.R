#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)

pdf(args[2], width = 11.69, height = 8.27)

par(mfcol = c(5,2), oma = c(0, 0, 3, 0), mar=c(4,4,1,0.5))

# create the boxplot for the quality values
data = read.table(args[1])

colors = c("red", "blue", "yellow", "green", "purple")

# for the first read
plotdata = data[which(data$V1 == "1" & data$V8 != 0),]
rA = plotdata$V3/plotdata$V8
rC = plotdata$V4/plotdata$V8
rG = plotdata$V5/plotdata$V8
rT = plotdata$V6/plotdata$V8
rN = plotdata$V7/plotdata$V8
yl = max(c(rA,rC,rG,rT,rN))

barplot(rA, ylim = c(0,yl),main = "A",ylab = "",xlab = "",col = colors[1]) 
barplot(rC, ylim = c(0,yl),main = "C",ylab = "",xlab = "",col = colors[2]) 
barplot(rG, ylim = c(0,yl),main = "G",ylab = "",xlab = "",col = colors[3]) 
barplot(rT, ylim = c(0,yl),main = "T",ylab = "",xlab = "",col = colors[4]) 
barplot(rN, ylim = c(0,yl),main = "N",ylab = "",xlab = "",col = colors[5]) 

# for the second read
plotdata = data[which(data$V1 == "2" & data$V8 != 0),]
rA = plotdata$V3/plotdata$V8
rC = plotdata$V4/plotdata$V8
rG = plotdata$V5/plotdata$V8
rT = plotdata$V6/plotdata$V8
rN = plotdata$V7/plotdata$V8
yl = max(c(rA,rC,rG,rT,rN))

plotdata = data[which(data$V1 == "1"),]
barplot(rA, ylim = c(0,yl),main = "A",ylab = "",xlab = "",col = colors[1]) 
barplot(rC, ylim = c(0,yl),main = "C",ylab = "",xlab = "",col = colors[2]) 
barplot(rG, ylim = c(0,yl),main = "G",ylab = "",xlab = "",col = colors[3]) 
barplot(rT, ylim = c(0,yl),main = "T",ylab = "",xlab = "",col = colors[4]) 
barplot(rN, ylim = c(0,yl),main = "N",ylab = "",xlab = "",col = colors[5]) 

mtext("Nucleotide composition variation (Read1, Read2)", outer=TRUE, cex = 1.5)

dev.off()
