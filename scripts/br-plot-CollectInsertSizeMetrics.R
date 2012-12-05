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
  opt <- parse_args(option_parser,positional_arguments=TRUE)
  return(opt)
}
opt <- parse_arguments()
args <- opt$args
if (length(args) == 0) {
    stop(paste(sep="\n","Zero arguments",usage.string))
}

for (arg in args) {
    # Read BEDTools CollectInsertSizeMetrics histogram
    insert.sizes <- read.table(arg,skip=10,header=TRUE)

    png(filename=paste(sep="",arg,".hist.png"),width=1000,height=1000)
    plot(insert.sizes$insert_size,insert.sizes$All_Reads.fr_count,
         xlab="insert size",ylab="number of pairs",
         main=paste("Insert size distribution of",arg))
    # Vertical line at the median, a legend, and the text value
    median.insert.size = median(rep(insert.sizes$insert_size,
                                    insert.sizes$All_Reads.fr_count))
    abline(v=median.insert.size, col = "red")
    legend("topright", inset=.05, c("Median insert size"), fill="red")
    text(x=median.insert.size, y=0,
         labels=as.character(median.insert.size), 
         col="red")
    suppressMessages(dev.off())
}
