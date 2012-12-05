#!/usr/bin/env Rscript
usage.string <- 
"%prog [options] genomeCoverageBedoutput [genomeCoverageBedoutput2,...]
Prints plots from genome/contig coverage file computed with genomeCoverageBed.

Results in 3 output files:
genomeCoverageBed.contig.png
genomeCoverageBed.base.png
genomeCoverageBed.length.png
"

# Parse arguments
parse_arguments <- function() {
  suppressMessages({
    library(optparse)
  })
  option_list <- list()
  option_parser <- OptionParser(option_list=option_list,
                     usage=usage.string)
  opt <- parse_args(option_parser, positional_arguments=TRUE)
  return(opt)
}
opt <- parse_arguments()
args <- opt$args
if (length(args) == 0) {
    stop(paste(sep="\n","Zero arguments",usage.string))
}

for (arg in args) {
    # Read BEDTools GenomeCoverageBed histogram
    cov <- read.table(arg)

    # Determine per contig mean coverage
    meancovcon <- aggregate(cov$V2[cov$V1 != "genome"] 
                            * cov$V3[cov$V1 != "genome"] 
                            / cov$V4[cov$V1 != "genome"],
                            list(cov$V1[cov$V1 != "genome"]),
                            sum)
    png(filename=paste(sep="",arg,".contig.png"),width=1000,height=1000)
    hist(meancovcon$x,
         breaks=1000,
         xlab="mean coverage",
         ylab="# of contigs",
         main=paste("Per contig mean coverage of ",
                    arg))
    suppressMessages(dev.off())

    # Determine per contig 1K mean coverage
    meancovcon <- aggregate(cov$V2[cov$V1 != "genome" & cov$V4 > 1000] 
                            * cov$V3[cov$V1 != "genome" & cov$V4 > 1000] 
                            / cov$V4[cov$V1 != "genome" & cov$V4 > 1000],
                            list(cov$V1[cov$V1 != "genome" & cov$V4 > 1000]),
                            sum)
    png(filename=paste(sep="",arg,".contig.1K.png"),width=1000,height=1000)
    hist(meancovcon$x,
         breaks=1000,
         xlab="mean coverage",
         ylab="# of contigs",
         main=paste("Per contig mean coverage of ",
                    arg))
    suppressMessages(dev.off())

    # Determine per base coverage
    png(filename=paste(sep="",arg,".base.png"),width=1000,height=1000)
    plot(cov$V2[cov$V1 == "genome"],
         cov$V3[cov$V1 == "genome"],
         xlab="coverage",
         ylab="Number of bases",
         main=paste("Per base coverage of ",
                    arg))
    # Vertical line at the median, a legend, and the text value
    median.cov <- median(rep(cov$V2,cov$V3))
    abline(v=median.cov, col = "red")
    legend("topright", inset=.05, c("Median coverage"), fill="red")
    text(x=median.cov, y=0,
         labels=paste(as.character(median.cov), "x"), 
         col="red")
    suppressMessages(dev.off())

    # Determine contig size histogram
    png(filename=paste(sep="",arg,".length.png"),width=1000,height=1000)
    contig.lengths <- aggregate(cov$V4[cov$V1 != "genome"], 
                                list(cov$V1[cov$V1 != "genome"]), 
                                min)
    hist(contig.lengths$x,
         xlab="Contig length",
         ylab="# of contigs",
         main=paste("Contig length distribution ",
                    arg))
    sorted.contig.lengths <- sort(contig.lengths$x)
    n50 <- floor(length(sorted.contig.lengths) / 2)
    l50 <- sorted.contig.lengths[[n50]]
    points(x=l50,y=n50,col="red")
    text(x=l50,
         y=n50, 
         labels=paste("(L50 ",as.character(l50),
                      "; N50 ", as.character(n50),
                      ")",sep=""),
         col="red",
         pos=4)
    legend("topright",c(paste("Total contigs",
                              as.character(length(unique(cov$V1)) - 1)),
                        paste("Total bases",
                              as.character(cov$V4[cov$V1 ==
                              "genome"][[1]]))))
    suppressMessages(dev.off())
}
