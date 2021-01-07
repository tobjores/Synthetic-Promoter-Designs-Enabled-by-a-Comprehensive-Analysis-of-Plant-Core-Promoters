#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

# test arguments
if (length(args) < 2) {
  stop("Two files required: 1. tsv file with subassembly; 2. fasta file with WT sequences", call. = FALSE)
}

library(readr)
library(dplyr)
library(stringr)

# set input file
infile <- args[1]

# set output file
if (length(args) == 2) {
  args[3] <- paste0(str_extract(infile, '^.*(?=\\.)'), '_variant.tsv')
}

outfile <- args[3]

# read wild type sequences
WT.seqs <- read_lines(args[2])

WT <- tibble(
  gene = str_replace(WT.seqs[seq(1, length(WT.seqs), 2)], fixed('>'), ''),
  sequence = WT.seqs[seq(2, length(WT.seqs), 2)]
)

rm(WT.seqs)

# calculate stop coordinate
get_stop <- function(cigar, start) {
  stop = start + sum(as.numeric(unlist(str_split(str_replace_all(cigar, '[0-9]+I', ''), '[=XD]'))), na.rm = TRUE) - 1
}

# extract mutations from cigar string
get_mut <- function(gene, cigar, start, seq) {

  if (grepl('^[0-9]+=$', cigar)) {

    mutations <- 'WT'

  } else {

    bases <- unlist(str_split(cigar, '[=XID]'))
    type <- unlist(str_split(cigar, '[0-9]+'))

    mutations <- tibble(
        bases = as.numeric(bases[bases != '']),
        type = type[type != '']
      ) %>%
      mutate(
        pos.wt = cumsum(if_else(type == 'I', 0, bases)) + start - 1,
        pos.mut = cumsum(if_else(type == 'D', 0, bases)),
      ) %>%
      filter(type != '=') %>%
      mutate(
        mutation = case_when(
          type == 'I' ~ paste(pos.wt, '_', pos.wt + 1, 'ins', str_sub(seq, start = pos.mut - bases + 1, end = pos.mut), sep = ''),
          type == 'D' & bases == 1 ~ paste(pos.wt - bases + 1, 'del', sep = ''),
          type == 'D' & bases > 1 ~ paste(pos.wt - bases + 1, '_', pos.wt, 'del', sep = ''),
          type == 'X' & bases == 1 ~ paste(pos.wt, str_sub(WT[WT$gene == gene, 'sequence'], start = pos.wt, end = pos.wt), '>', str_sub(seq, start = pos.mut, end = pos.mut), sep = ''),
          type == 'X' & bases > 1 ~ paste(pos.wt - bases + 1, '_', pos.wt, 'delins', str_sub(seq, start = pos.mut - bases + 1, end = pos.mut), sep = '')
        )
      ) %>%
      pull(mutation)

    mutations <- paste(mutations, collapse = ';')

  }

  return(mutations)

}

# annotate subassembly with variant info
subassembly <- read_tsv(infile) %>%
  rowwise() %>%
  mutate(
    stop = get_stop(cigar, start),
    variant = get_mut(gene, cigar, start, sequence)
  ) %>%
  select(-cigar, -sequence) %>%
  select(barcode, gene, start, stop, variant, everything())

write_tsv(subassembly, outfile)
