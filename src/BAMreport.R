#!/usr/bin/Rscript

library(MASS)

## Default when nothing is passed
args = commandArgs(trailingOnly = TRUE)
if(length(args) == 0){
    args <- c("--help")    
}

## Help section
if("--help" %in% args) {
  cat("
      Prepare a report for the input BAM file.
 
      Arguments:
      --stats=someValue   - char, name of file with the stats
      --output=someValue  - char, name of the output PDF file
      --rlens=someValue   - char, name of the file with the read length distr.
      --rnucs=someValue   - char, name of the file with nucleotide distr.
      --rqual=someValue   - char, name of the file with quality distr.
      --insrt=someValue   - char, name of the file with insert-length distr.
      --rcovs=someValue   - char, name of the file with coverage distr.
      --gccov=someValue   - char, name of the file with the GC coverage info.
      --fclip=someValue   - char, name of the file with 5' clipping info.
      --mm=someValue      - char, name of the file with mismatch info.
      --indel=someValue   - char, name of the file with indel info.
      --help              - print this text
 
      Example:
      ./BAMreport.R --stats=stats.txt --output=report.pdf \n\n")
 
  q(save="no")
}

## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1
 
# set defaults
if(is.null(argsL$stats)) {
    argsL$stats = NULL
}
if(is.null(argsL$output)) {
    argsL$output = "report.pdf"
}
if(is.null(argsL$rlens)) {
    argsL$rlens = NULL
}
if(is.null(argsL$rnucs)) {
    argsL$rnucss = NULL
}
if(is.null(argsL$rqual)) {
    argsL$rqual = NULL
}
if(is.null(argsL$insrt)) {
    argsL$insrt = NULL
}
if(is.null(argsL$rcovs)) {
    argsL$rcovs = NULL
}
if(is.null(argsL$gccov)) {
    argsL$gccov = NULL
}
if(is.null(argsL$fclip)) {
    argsL$fclip = NULL
}
if(is.null(argsL$mm)) {
    argsL$mm = NULL
}
if(is.null(argsL$indel)) {
    argsL$indel = NULL
}

# the output file
pdf(argsL$output, width = 14, height = 10)
 
# create the stats page
# ---------------------
if (!is.null(argsL$stats)) {
    stats = argsL$stats
    fileConn = file(stats)
    data = readLines(fileConn)
    close(fileConn)

    plot.new()
    indx = 0
    for (l in data) {
        mtext(l, side = 3, line = indx, adj = 0)
        indx = indx - 1
    }
}

# create the read length distribution
# -----------------------------------
if (!is.null(argsL$rlens)) {
    data = read.table(argsL$rlens)
    tmp  = data[which(data$V1 == "2"),]

    par(cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5)
    if (max(tmp$V3) == 0) {
        par(mfrow = c(1,1), oma = c(0, 0, 3, 0))
    } else {         
        par(mfrow = c(1,2), oma = c(0, 0, 3, 0))
    }

    # for the first read
    plotdata = data[which(data$V1 == "1"),]
    total = sum(as.numeric(plotdata$V3))
    mp = barplot(plotdata$V3/total,
            xlab = "Length of read",
            ylab = "Fraction of read1's",
            axes = F,
            main = "Read 1")
    
    index = seq(0, max(plotdata$V2), 10)
    index[1] = 1
    
    axis(1, at = mp[index], labels = index)
    axis(2)
    
    if (max(tmp$V3) > 0) {
        # for the second read
        plotdata = data[which(data$V1 == "2"),]
        total = sum(as.numeric(plotdata$V3))
        mp = barplot(plotdata$V3/total,
                xlab = "Length of read",
                ylab = "Fraction of read2's",
                axes = F,
                main = "Read 2")
        
        index = seq(0, max(plotdata$V2), 10)
        index[1] = 1
        
        axis(1, at = mp[index], labels = index)
        axis(2)
    }
    mtext("Read length distribution", outer = TRUE, cex = 1.5)
}

# create the nucleotide distribution
# -----------------------------------
if (!is.null(argsL$rnucs)) {
    colors = c("red", "blue", "yellow", "green", "purple")
    data = read.table(argsL$rnucs)
    tmp  = data[which(data$V1 == "2"),]
    
    par(cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5)
    if (max(tmp$V8) == 0) {
        par(mfcol = c(5,1), oma = c(0, 0, 3, 0))
    } else {
        par(mfcol = c(5,2), oma = c(0, 0, 3, 0))
    }
    
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
    
    if (max(tmp$V8) != 0) {
        # for the second read
        plotdata = data[which(data$V1 == "2" & data$V8 != 0),]
        rA = plotdata$V3/plotdata$V8
        rC = plotdata$V4/plotdata$V8
        rG = plotdata$V5/plotdata$V8
        rT = plotdata$V6/plotdata$V8
        rN = plotdata$V7/plotdata$V8
        yl = max(c(rA,rC,rG,rT,rN))
        
        barplot(rA,ylim = c(0,yl),main= "A",ylab = "",xlab = "",col = colors[1]) 
        barplot(rC,ylim = c(0,yl),main= "C",ylab = "",xlab = "",col = colors[2]) 
        barplot(rG,ylim = c(0,yl),main= "G",ylab = "",xlab = "",col = colors[3]) 
        barplot(rT,ylim = c(0,yl),main= "T",ylab = "",xlab = "",col = colors[4]) 
        barplot(rN,ylim = c(0,yl),main= "N",ylab = "",xlab = "",col = colors[5]) 
    }

    mtext("Nucleotide composition variation (Read1, Read2)", outer=TRUE, cex = 1.5)
}

# create the quality distribution
# -------------------------------
if (!is.null(argsL$rqual)) {
    data = read.table(argsL$rqual)
    tmp  = data[which(data$V1 == "2"),]
    par(cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5)
    if (max(tmp$V7) == 0) {
        par(mfrow = c(1,1), oma = c(0, 0, 3, 0))
    } else {
        par(mfrow = c(1,2), oma = c(0, 0, 3, 0))
    }

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
    
    if (max(tmp$V7) != 0) {
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
        
        axis(1,at=seq(1, nrow(plotdata), 5), labels = seq(1, nrow(plotdata), 5))
        axis(2, at = seq(0, 60, 10), labels = seq(0, 60, 10))
    }

    mtext("Quality value variation (Read1,Read2)", outer = TRUE, cex = 1.5)
}

# plot the insert length distribution
# -----------------------------------
if (!is.null(argsL$insrt)) {
    par(cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5)
    par(mfrow = c(1,1), oma = c(0, 0, 3, 0))
    colors = c(rgb(255,0,0,150,maxColorValue=255),
               rgb(0,0,255,100,maxColorValue=255))
    data = read.table(argsL$insrt)
    total = sum(as.numeric(data$V2))
    data$V2 = data$V2 / total
    xmax = quantile(rep(data$V1,data$V2*1000),.99)
    
    # find the best bit normal distribution
    smpl = rep(data$V1,data$V2*1000)
    fit = fitdistr(smpl, "normal")
    x = seq(0,xmax,length=10*xmax)*fit$estimate[2]
    hx = dnorm(x, fit$estimate[1], fit$estimate[2])

    ymax = max(max(hx),max(data$V2))    

    plot(data,
         type = "l",
         lwd  = 3,
         col  = colors[1],
         main = "Insert length distribution",    
         xlab = "Insert length",
         ylab = "Fraction of properly-paired pairs",
         xlim = c(0, xmax),
         ylim = c(0, ymax))
    lines(x, hx, col = colors[2], lwd = 3)
    
    legend("topright",
           legend = c("observed",paste("best normal fit (",round(fit$estimate[1],2),",",round(fit$estimate[2],2), ")")),
           fill = colors,
           cex = 1.5)
}

# plot the coverage distribution
# ------------------------------
if (!is.null(argsL$rcovs)) {
    par(cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5)
    par(mfrow = c(1,1), oma = c(0, 0, 3, 0))
    colors = c(rgb(255,0,0,150,maxColorValue=255),
               rgb(0,0,255,100,maxColorValue=255))

    data = read.table(argsL$rcovs)

    datax = data[-1,]
    xlimit = quantile(rep(datax$V1, datax$V2), 0.98)

    bases = sum(as.numeric(data$V2))
    data$V2 = data$V2 / bases
    data = data[-1,]
    repfrac = sum(as.numeric(data$V2))
 
    mfac = 100 / data$V2[1]
    smpl = rep(data$V1,data$V2 * mfac)
    fit = fitdistr(smpl, "normal")
    x = seq(0,xlimit,length=100*xlimit)*fit$estimate[2]
    hx = dnorm(x, fit$estimate[1], fit$estimate[2]) * repfrac
    ymax = max(max(hx),max(data$V2))    

    plot(data,
     type = "l",
     lwd  = 3,
     col  = colors[1],
     main = "Depth of coverage",    
     xlab = "Coverage",
     ylab = "Fraction of genome covered",
     xlim = c(0, xlimit),
     ylim = c(0, ymax))

    lines(x, hx, col = colors[2], lwd = 3)
    legend("topright",
       legend = c("observed",paste("best normal fit (",round(fit$estimate[1],2),",",round(fit$estimate[2],2), ")")),
       fill = colors, 
       cex = 1.5)
}

# plot the normalized coverage distribution (more informative with low cov)
# -----------------------------------------
if (!is.null(argsL$gccov)) {
    data = read.table(argsL$gccov)
    coverage = data$V5
    avgcoverage = mean(coverage)
    relativecov = round(coverage / avgcoverage, 3)
    
    h = hist(relativecov, xlim = c(0,2), breaks = 1000)
    
    rel  = subset(relativecov,relativecov <= 2)
    xfit = seq(0,2,length=1000)
    yfit = dnorm(xfit, mean = mean(rel), sd = sd(rel)) 
    yfit = yfit * max(h$counts) / max(yfit)
    lines(xfit, yfit, col="blue", lwd=2)
}

# plot the GC coverage distribution
# -------------------------------
if (!is.null(argsL$gccov)) {
    data = read.table(argsL$gccov)
    ylimit = quantile(subset(data$V5, data$V5 != 0.00), 0.98)
    par(cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5)
    
    x = data$V4
    y = data$V5
    plot(x, y, 
         main = "GC content vs Coverage from aligned reads",    
         xlab = "GC content",
         ylab = "Coverage",
         col=rgb(0,100,0,50,maxColorValue=255), 
         pch=16,
         ylim = c(0, ylimit))
}

# plot the 5' clipping positions
# ---------------------------------
if (!is.null(argsL$fclip)) {
    data = read.table(argsL$fclip)
    tmp  = data[which(data$V1 == "2"),]
    par(cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5)
    if (max(tmp$V4) == 0) {
        par(mfrow = c(1,1), oma = c(0, 0, 3, 0))
    } else {
        par(mfrow = c(1,2), oma = c(0, 0, 3, 0))
    }    

    # for the first read
    plotdata = data[which(data$V1 == "1"),]
    mp = barplot(plotdata$V3/plotdata$V4,
            xlab = "5' Clip position",
            ylab = "Fraction of read1's",
            axes = F,
            main = "Read 1")
    
    index = seq(0, max(plotdata$V2), 10)
    index[1] = 1
    
    axis(1, at = mp[index], labels = index)
    axis(2)
    
    if (max(tmp$V4) != 0) {
        # for the second read
        plotdata = data[which(data$V1 == "2"),]
        mp = barplot(plotdata$V3/plotdata$V4,
                xlab = "5' Clip position",
                ylab = "Fraction of read2's",
                axes = F,
                main = "Read 2")
        
        index = seq(0, max(plotdata$V2), 10)
        index[1] = 1
        
        axis(1, at = mp[index], labels = index)
        axis(2)
    }

    mtext("5' Clipping locations on aligned reads", outer = TRUE, cex = 1.5)
}

# plot the position of mismatches
# -------------------------------
if (!is.null(argsL$mm)) {
    data = read.table(argsL$mm)
    tmp  = data[which(data$V1 == "2"),]
    par(cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5)
    if (max(tmp$V4) == 0) {
        par(mfrow = c(1,1), oma = c(0, 0, 3, 0))
    } else {
        par(mfrow = c(1,2), oma = c(0, 0, 3, 0))
    }

       
    # for the first read
    plotdata = data[which(data$V1 == "1"),]
    mp = barplot(plotdata$V3 / plotdata$V4,
            xlab = "Position on the read",
            ylab = "Fraction of read1s with mismatches",
            axes = F,
            main = "Read 1")
    
    index = seq(0, max(plotdata$V2), 10)
    index[1] = 1
    
    axis(1, at = mp[index], labels = index)
    axis(2)
    
    if (max(tmp$V4) != 0) {
        # for the second read
        plotdata = data[which(data$V1 == "2"),]
        mp = barplot(plotdata$V3 / plotdata$V4,
                xlab = "Position on the read",
                ylab = "Fraction of read2s with mismatches",
                axes = F,
                main = "Read 2")
        
        index = seq(0, max(plotdata$V2), 10)
        index[1] = 1
        
        axis(1, at = mp[index], labels = index)
        axis(2)
    }

    mtext("Mismatch positions (vs reference) on aligned reads", outer = TRUE, cex = 1.5)
}

# plot the position of indels
# ---------------------------------
if (!is.null(argsL$indel)) {
    data = read.table(argsL$indel)
    tmp  = data[which(data$V1 == "2"),]
    par(cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5)
    if (max(tmp$V4) == 0) {
        par(mfrow = c(1,1), oma = c(0, 0, 3, 0))
    } else {
        par(mfrow = c(1,2), oma = c(0, 0, 3, 0))
    }

   
    # for the first read
    plotdata = data[which(data$V1 == "1"),]
    mp = barplot(plotdata$V3 / plotdata$V4,
            xlab = "Position on the read",
            ylab = "Fraction of read1s with indels",
            axes = F,
            main = "Read 1")
    
    index = seq(0, max(plotdata$V2), 10)
    index[1] = 1
    
    axis(1, at = mp[index], labels = index)
    axis(2)
    
    if (max(tmp$V4) != 0) {
        # for the second read
        plotdata = data[which(data$V1 == "2"),]
        mp = barplot(plotdata$V3 / plotdata$V4,
                xlab = "Position on the read",
                ylab = "Fraction of read2s with indels",
                axes = F,
                main = "Read 2")
        
        index = seq(0, max(plotdata$V2), 10)
        index[1] = 1
        
        axis(1, at = mp[index], labels = index)
        axis(2)
    }    
    mtext("Indel positions (vs reference) on the read", outer = TRUE, cex = 1.5)
}

dev.off()
