#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

# test arguments
if (length(args) < 3) {
  stop("Three files required: 1. tsv file with barcodes; 2. tsv file with sequences; 3. file for output", call. = FALSE)
} else if (length(args) == 3) {
  args[4] <- '1'
}

library(readr)
library(dplyr)

barcodes <- args[1]
sequences <- args[2]
outfile <- args[3]

min.reads <- as.numeric(args[4])

joined <- inner_join(read_tsv(barcodes), read_tsv(sequences), by = 'read') %>%
  group_by(across(-read)) %>%
  tally(name = 'assembly.count') %>%
  filter(assembly.count >= min.reads)

write_tsv(joined, outfile)
