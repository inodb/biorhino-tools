#!/usr/bin/env Rscript
usage.string <- 
"%prog [options] table [table2,...]
Plots R summary of given tables
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
    #stop(paste(sep="\n","Zero arguments",usage.string))
    args <- c("stdin")
}

for (arg in args) {
    t <- read.table(arg)
    print(summary(t))
}
