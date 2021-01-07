#!/usr/bin/env Rscript

library(readr)
library(dplyr)

for (species in c('Arabidopsis', 'Maize', 'Sorghum')) {
  noUTR <- read_lines(paste0(species, '_noUTR.txt'))

  promoters <- read_tsv(paste0(species, '_all_promoters.tsv')) %>%
    mutate(
      UTR = (! gene %in% noUTR)
    ) %>%
    group_by(sequence, type, strand, mutations) %>%
    summarise(
      gene = paste0(gene, collapse = ';'),
      UTR = all(UTR)
    ) %>%
    ungroup() %>%
    group_by(sequence) %>%
    summarise(
      across(-UTR, ~if_else(n_distinct(.x) == 2, paste0(.x, collapse = '/'), first(.x))),
      UTR = all(UTR)
    ) %>%
    ungroup() %>%
    select(gene, type, sequence, strand, UTR, mutations) %>%
    arrange(gene)
  
  write_tsv(promoters, paste0(species, '_all_promoters_unique.tsv'), na = '')

  fasta <- promoters %>%
    select(gene, sequence) %>%
    bind_rows(c(gene = '35Spr', sequence = 'GCAAGACCCTTCCTCTATATAAGGAAGTTCATTTCATTTGGAGAGGACACG')) %>%
    mutate(
      gene = paste0('>', gene)
    ) %>%
    select(gene, sequence)

  write_delim(fasta, paste0(species, '_all_promoters_unique.fa'), delim = '\n', col_names = FALSE)

  assign(species, promoters)
}

all.promoters <- bind_rows(Arabidopsis, Maize, Sorghum) %>%
  select(gene, sequence) %>%
  bind_rows(c(gene = '35Spr', sequence = 'GCAAGACCCTTCCTCTATATAAGGAAGTTCATTTCATTTGGAGAGGACACG')) %>%
  mutate(
    gene = paste0('>', gene)
  ) %>%
  select(gene, sequence)

write_delim(all.promoters, 'all_promoters_unique.fa', delim = '\n', col_names = FALSE)
