#!/usr/bin/Rscript

library(MASS)
library(ggplot2)
library(RColorBrewer)

## Default when nothing is passed
args = commandArgs(trailingOnly = TRUE)
if(length(args) == 0){
    args <- c("--help")    
}

## Help section
if("--help" %in% args) {
  cat("
      Prepare a comparison report for the different prep methods.
 
      Arguments:
      --stats=someValue   - char, name of file with the stats
      --windw=someValue   - char, comma-separated files with the window stats
      --insrt=someValue   - char, comma-separated files with insert length freq
      --output=someValue  - char, name of the output PDF file
      --help              - print this text
 
      Example:
      ./BAMreport.R --stats=stats.txt --output=report.pdf \n\n")
 
  q(save="no")
}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot
# objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
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
if(is.null(argsL$windw)) {
    argsL$windw = NULL
}
if(is.null(argsL$insrt)) {
    argsL$insrt = NULL
}

# the output file
pdf(argsL$output, width = 14, height = 10)
 
# create the stats page
# ---------------------
if (!is.null(argsL$stats)) {
    stats = argsL$stats

    data = read.table(stats, header = TRUE)
    data$aligned = data$aligned / data$generated
    data$mapped = data$mapped / data$generated
    data$dups = data$dups / data$generated

    p1 = ggplot(data, aes(x = prep, y = aligned)) + geom_boxplot()
    p1 = p1 + labs(title = "Fraction aligned", y = "fraction of aligned reads")  
    p1 = p1 + ylim(0,1)

    p2 = ggplot(data, aes(x = prep, y = mapped)) + geom_boxplot()
    p2 = p2 + labs(title = "Fraction mapped", y = "fraction of mapped reads")  
    p2 = p2 + ylim(0,1)

    p3 = ggplot(data, aes(x = prep, y = dups)) + geom_boxplot()
    p3 = p3 + labs(title = "Fraction duplicates", y = "fraction of duplicates")  
    p3 = p3 + ylim(0,1)

    multiplot(p1,p2,p3, cols = 3)
}

# create the window report
if (!is.null(argsL$windw)) {
    files = strsplit(argsL$windw, ",")

    # print the coverage distribution
    ndata = read.table(files[[1]][1], header = TRUE)
    names = c(strsplit(names(ndata)[1], "_")[[1]][1])
    coverage = ndata$cov
    avgcoverage = mean(coverage)
    relativecov = round(coverage / avgcoverage, 3)
    data = as.data.frame(relativecov)
    times = c(nrow(data))

    for (f in files[[1]][2:length(files[[1]])]) {
        ndata = read.table(f, header = TRUE)
        names = cbind(names,strsplit(names(ndata)[1], "_")[[1]][1])
        coverage = ndata$cov
        avgcoverage = mean(coverage)
        relativecov = round(coverage / avgcoverage, 3)
        tmp = as.data.frame(relativecov)
        times = c(times, nrow(tmp))
        data = rbind(data, tmp)
    }
    data$grp = rep(factor(names), times = times)
    p = ggplot(data = data, aes(x = relativecov, fill = grp)) 
    p = p + geom_density(alpha = 0.2) 
    p = p + xlim(0,2) 
    p = p + labs(x = "Window coverage / Median window coverage")
    p = p + labs(y = "Fraction of windows")
    p = p + labs(title = "Relative depth of coverage distribution")
    print(p)

    # print the boxplot for the number of start points
    ndata = read.table(files[[1]][1], header = TRUE)
    names = c(strsplit(names(ndata)[1], "_")[[1]][1])
    numstarts = ndata$numstarts
    numstarts = numstarts / median(numstarts)
    data = diff(numstarts)
    data = as.data.frame(data)
    names(data) = c("numstarts")
    times = c(nrow(data))

    for (f in files[[1]][2:length(files[[1]])]) {
        ndata = read.table(f, header = TRUE)
        names = cbind(names,strsplit(names(ndata)[1], "_")[[1]][1])
        numstarts = ndata$numstarts
        numstarts = numstarts / median(numstarts)
        tmp = diff(numstarts)
        tmp = as.data.frame(tmp)
        times = c(times, nrow(tmp))
        names(tmp) = c("numstarts")
        data = rbind(data, tmp)
    }
    data$grp = rep(factor(names), times = times)
    p = ggplot(data = data, aes(x = grp, y = numstarts, fill = grp))
    p = p + geom_boxplot(outlier.shape = NA)
    p = p + scale_y_continuous(limits = quantile(data$numstarts, c(0.1, 0.9)))
    p = p + labs(x = "Run")
    p = p + labs(y = "Normalized difference between neighboring windows")
    p = p + labs(title = "Variation in difference between neighboring windows")
    print(p)


    # print the boxplot for the number of fragments
}

# create the insert length report
if (!is.null(argsL$insrt)) {
    files = strsplit(argsL$insrt, ",")
    
    data = read.table(files[[1]][1], header = TRUE)
    times = c(nrow(data))
    names = c(names(data)[2])
    colnames(data) = c("V1","V2")
    data$V2 = data$V2 / sum(data$V2)

    for (f in files[[1]][2:length(files[[1]])]) {
        ndata = read.table(f, header = TRUE)
        times = cbind(times, nrow(ndata))
        names = cbind(names, names(ndata)[2])
        colnames(ndata) = c("V1","V2")
        ndata$V2 = ndata$V2 / sum(ndata$V2)

        data = rbind(data, ndata)
    }
    
    data$grp = rep(factor(names), times = times)
    p = ggplot(data = data, aes(x = V1, y = V2, colour = grp)) + geom_line()
    p = p + labs(x = "Insert distance")
    p = p + labs(y = "Fraction of pairs")
    p = p + labs(title = "Insert length distribution")
    print(p)
}

dev.off()
