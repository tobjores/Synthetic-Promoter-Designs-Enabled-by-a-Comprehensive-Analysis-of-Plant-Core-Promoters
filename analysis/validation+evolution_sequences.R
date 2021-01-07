library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(universalmotif)


### load motifs
CPEs <- read_meme('../data/misc/CPEs.meme')
TFs <- read_meme('../data/misc/TF-clusters.meme')


### load main data
load('RData/main_data.Rdata')


### load promoter sequences
all.sequences <- bind_rows(
  'At' = read_tsv('../data/promoter_annotation/Arabidopsis_all_promoters_unique.tsv'),
  'Zm' = read_tsv('../data/promoter_annotation/Maize_all_promoters_unique.tsv'),
  'Sb' = read_tsv('../data/promoter_annotation/Sorghum_all_promoters_unique.tsv'),
  .id = 'sp'
)

all.seqs <- Biostrings::DNAStringSet(deframe(select(all.sequences, gene, sequence)))


### select sequences for promoters detected in both expression systems in the dark with enhancer
detected.promoters <- promoter.strength %>%
  filter(! light & enhancer) %>%
  select(sys, gene) %>%
  group_by(gene) %>%
  filter(any(sys == 'leaf') & any(sys == 'proto')) %>%
  ungroup() %>%
  distinct(gene) %>%
  pull()

detected.seqs <- all.seqs[detected.promoters]


### load export functions
source('export_functions.R')


### pick sequences for validation (first validation set)

## TATA
TATA.scores <- as_tibble(scan_sequences(filter_motifs(CPEs, altname = 'TATA'), detected.seqs, threshold = 1, nthreads = 0)) %>%
  mutate(
    rel.score = (score - min.score) / (max.score - min.score)
  ) %>%
  rename('gene' = sequence) %>%
  inner_join(select(all.sequences, -strand), by = 'gene')

TATA.high.score <- TATA.scores %>%
  filter(rel.score >= 0.75) %>%
  group_by(gene) %>%
  filter(n() == 1) %>%
  ungroup() %>%
  filter(rel.score >= 0.9 & between(start, 125, 140)) %>%
  filter(type == 'protein_coding' & is.na(mutations) & ! grepl('[/;]', gene)) %>%
  filter(substr(match, 5, 5) == 'T' & substr(match, 7, 7) == 'T') %>%
  group_by(sp) %>%
  slice_sample(n = 10) %>%
  ungroup()

# mutations: start + 4 T -> G; start + 6 T -> G
TATA.high.sequences <- TATA.high.score %>%
  mutate(
    motif = 'TATA',
    WT = sequence,
    mutA = `str_sub<-`(sequence, start + 4, start + 4, value = 'G'),
    mutB = `str_sub<-`(sequence, start + 6, start + 6, value = 'G'),
    mutAB = `str_sub<-`(mutA, start + 6, start + 6, value = 'G')
  )

# export motif logos for WT and mutated sequences
TATA.variants <- TATA.high.sequences %>%
  mutate(
    across(c(WT, mutA, mutB, mutAB), str_sub, start = start + 4, end = start + 10)
  ) %>%
  select(WT, mutA, mutB, mutAB)

TATA.motifs <- lapply(TATA.variants, create_motif)

lapply(seq_along(TATA.motifs), function(x) PWM.to.LaTeX(TATA.motifs[[x]], paste0('../figures/rawData/motif_mutTATA_', names(TATA.motifs)[[x]], '.tsv')))


TATA.high.sequences <- TATA.high.sequences %>%
  select(sp, gene, motif, WT, starts_with('mut'), -mutations) %>%
  pivot_longer(
    c(WT, starts_with('mut')),
    names_to = 'variant',
    values_to = 'sequence'
  )

TATA.high.to.order <- TATA.high.score %>%
  mutate(
    motif = 'TATA',
    sequence = `str_sub<-`(sequence, start + 4, start + 4, value = 'K'),
    sequence = `str_sub<-`(sequence, start + 6, start + 6, value = 'K')
  ) %>%
  select(sp, motif, gene, sequence)


TATA.low.score <- TATA.scores %>%
  filter(rel.score >= 0.7) %>%
  group_by(gene) %>%
  filter(n() == 1) %>%
  ungroup() %>%
  filter(rel.score <= 0.75 & between(start, 125, 140)) %>%
  filter(type == 'protein_coding' & is.na(mutations) & ! grepl('[/;]', gene)) %>%
  group_by(sp) %>%
  slice_sample(n = 10) %>%
  ungroup()

# mutations: start +4 NNNNNNN -> TATAAAT (+TATA) or TAGAAAT (+mutTATA)
TATA.low.sequences <- TATA.low.score %>%
  mutate(
    motif = 'noTATA',
    WT = sequence,
    `+TATA` = `str_sub<-`(sequence, start + 4, start + 10, value = 'TATAAAT'),
    `+mutTATA` = `str_sub<-`(sequence, start + 4, start + 10, value = 'TAGAAAT')
  ) 

# export motif logos for WT and mutated sequences
TATA.low.variants <- TATA.low.sequences %>%
  mutate(
    across(c(WT, `+TATA`, `+mutTATA`), str_sub, start = start + 4, end = start + 10)
  ) %>%
  select(WT, `+TATA`, `+mutTATA`)

TATA.low.motifs <- lapply(TATA.low.variants, create_motif)

lapply(seq_along(TATA.low.motifs), function(x) PWM.to.LaTeX(TATA.low.motifs[[x]], paste0('../figures/rawData/motif_insTATA_', names(TATA.low.motifs)[[x]], '.tsv')))


TATA.low.sequences <- TATA.low.sequences %>%
  select(sp, motif, gene, WT, ends_with('TATA')) %>%
  pivot_longer(
    c(WT, ends_with('TATA')),
    names_to = 'variant',
    values_to = 'sequence'
  )

TATA.low.to.order <- TATA.low.score %>%
  mutate(
    motif = '+TATA',
    sequence = `str_sub<-`(sequence, start + 4, start + 10, value = 'TAKAAAT')
  ) %>%
  bind_rows(mutate(TATA.low.score, motif = 'noTATA')) %>%
  select(sp, motif, gene, sequence)


## TCP
TCP.scores <- as_tibble(scan_sequences(filter_motifs(TFs, altname = 'TCP'), detected.seqs, threshold = 1, RC = TRUE, nthreads = 0)) %>%
  mutate(
    rel.score = (score - min.score) / (max.score - min.score)
  ) %>%
  rename('gene' = sequence) %>%
  inner_join(select(all.sequences, -strand), by = 'gene')

TCP.high.score <- TCP.scores %>%
  filter(rel.score >= 0.75) %>%
  group_by(gene, motif) %>%
  filter(n() == 1) %>%
  ungroup() %>%
  filter(rel.score >= 0.9) %>%
  filter(type == 'protein_coding' & is.na(mutations) & ! grepl('[/;]', gene)) %>%
  mutate(
    base1 = case_when(
      strand == '+' & motif == 'TF-cluster_15' ~ start + 2,
      strand == '+' & motif == 'TF-cluster_22' ~ start + 1,
      strand == '-' ~ stop + 2
    )
  ) %>%
  filter(substr(sequence, base1, base1 + 1) == 'GG') %>%
  group_by(sp, motif) %>%
  slice_sample(n = 5) %>%
  ungroup()

# mutations: TF-cluster 15, forward: start + 2 GG -> TT; cluster 15, reverse: stop + 2 GG -> TT
#            TF-cluster 22, forward: start + 1 GG -> TT; cluster 22, reverse: stop + 2 GG -> TT
TCP.high.sequences <- TCP.high.score %>%
  mutate(
    motif = paste0('TCP(', if_else(motif == 'TF-cluster_15', '15', '22'), ')'),
    WT = sequence,
    mutA = if_else(strand == '+', `str_sub<-`(sequence, base1, base1, value = 'T'), NA_character_),
    mutB = if_else(strand == '+', `str_sub<-`(sequence, base1 + 1, base1 + 1, value = 'T'), NA_character_),
    mutAB = if_else(strand == '+', `str_sub<-`(mutA, base1 + 1, base1 + 1, value = 'T'), NA_character_),
    mutC = if_else(strand == '-', `str_sub<-`(sequence, base1, base1, value = 'T'), NA_character_),
    mutD = if_else(strand == '-', `str_sub<-`(sequence, base1 + 1, base1 + 1, value = 'T'), NA_character_),
    mutCD = if_else(strand == '-', `str_sub<-`(mutC, base1 + 1, base1 + 1, value = 'T'), NA_character_)
  ) %>%
  select(sp, motif, gene, WT, starts_with('mut'), -mutations) %>%
  pivot_longer(
    c(WT, starts_with('mut')),
    names_to = 'variant',
    values_to = 'sequence'
  ) %>%
  filter(! is.na(sequence))

TCP.high.to.order <- TCP.high.score %>%
  mutate(
    motif = 'TCP',
    sequence = `str_sub<-`(sequence, base1, base1 + 1, value = 'KK')
  ) %>%
  select(sp, motif, gene, sequence)


## HSF
HSF.scores <- as_tibble(scan_sequences(filter_motifs(TFs, altname = 'HSF/S1Fa-like'), detected.seqs, threshold = 1, RC = TRUE, nthreads = 0)) %>%
  mutate(
    rel.score = (score - min.score) / (max.score - min.score)
  ) %>%
  rename('gene' = sequence) %>%
  inner_join(select(all.sequences, -strand), by = 'gene')

HSF.high.score <- HSF.scores %>%
  filter(rel.score >= 0.76) %>%
  group_by(gene) %>%
  filter(n() == 1) %>%
  ungroup() %>%
  filter(rel.score >= 0.9) %>%
  filter(type == 'protein_coding' & is.na(mutations) & ! grepl('[/;]', gene)) %>%
  mutate(
    base1 = if_else(strand == '+', start + 6, stop + 1),
    base2 = if_else(strand == '+', start + 10, stop + 5),
  ) %>%
  filter(substr(sequence, base1, base1) == 'T' & substr(sequence, base2, base2) == 'G') %>%
  group_by(sp) %>%
  slice_sample(n = 10) %>%
  ungroup()

# mutations: forward: start + 6 T -> G; start + 10 G -> T; reverse: stop + 1 T -> G; stop + 5 G -> T
HSF.high.sequences <- HSF.high.score %>%
  mutate(
    motif = 'HSF',
    WT = sequence,
    mutA = if_else(strand == '+', `str_sub<-`(sequence, base1, base1, value = 'G'), NA_character_),
    mutB = if_else(strand == '+', `str_sub<-`(sequence, base2, base2, value = 'T'), NA_character_),
    mutAB = if_else(strand == '+', `str_sub<-`(mutA, base2, base2, value = 'T'), NA_character_),
    mutC = if_else(strand == '-', `str_sub<-`(sequence, base1, base1, value = 'G'), NA_character_),
    mutD = if_else(strand == '-', `str_sub<-`(sequence, base2, base2, value = 'T'), NA_character_),
    mutCD = if_else(strand == '-', `str_sub<-`(mutC, base2, base2, value = 'T'), NA_character_),
  ) %>%
  select(sp, motif, gene, WT, starts_with('mut'), -mutations) %>%
  pivot_longer(
    c(WT, starts_with('mut')),
    names_to = 'variant',
    values_to = 'sequence'
  ) %>%
  filter(! is.na(sequence))

HSF.high.to.order <- HSF.high.score %>%
  mutate(
    motif = 'HSF',
    sequence = `str_sub<-`(sequence, base1, base1, value = 'K'),
    sequence = `str_sub<-`(sequence, base2, base2, value = 'K')
  ) %>%
  select(sp, motif, gene, sequence)


## WRKY
WRKY.scores <- as_tibble(scan_sequences(filter_motifs(TFs, altname = 'WRKY/C3H'), detected.seqs, threshold = 1, RC = TRUE, nthreads = 0)) %>%
  mutate(
    rel.score = (score - min.score) / (max.score - min.score)
  ) %>%
  rename('gene' = sequence) %>%
  inner_join(select(all.sequences, -strand), by = 'gene')

WRKY.high.score <- WRKY.scores %>%
  filter(rel.score >= 0.75) %>%
  group_by(gene) %>%
  filter(n() == 1) %>%
  ungroup() %>%
  filter(rel.score >= 0.9) %>%
  filter(type == 'protein_coding' & is.na(mutations) & ! grepl('[/;]', gene)) %>%
  mutate(
    base1 = if_else(strand == '+', start + 2, stop + 1),
    base2 = if_else(strand == '+', start + 1, stop + 3)
  ) %>%
  filter(substr(sequence, base1, base1) == 'T' & substr(sequence, base2, base2) == 'G') %>%
  group_by(sp) %>%
  slice_sample(n = 10) %>%
  ungroup()

# mutations: forward: start + 1 GT -> TG; reverse: stop + 1 T -> G; stop + 3 G -> T
WRKY.high.sequences <- WRKY.high.score %>%
  mutate(
    motif = 'WRKY',
    WT = sequence,
    mutA = if_else(strand == '+', `str_sub<-`(sequence, base1, base1, value = 'G'), NA_character_),
    mutB = if_else(strand == '+', `str_sub<-`(sequence, base2, base2, value = 'T'), NA_character_),
    mutAB = if_else(strand == '+', `str_sub<-`(mutA, base2, base2, value = 'T'), NA_character_),
    mutC = if_else(strand == '-', `str_sub<-`(sequence, base1, base1, value = 'G'), NA_character_),
    mutD = if_else(strand == '-', `str_sub<-`(sequence, base2, base2, value = 'T'), NA_character_),
    mutCD = if_else(strand == '-', `str_sub<-`(mutC, base2, base2, value = 'T'), NA_character_)
  ) %>%
  select(sp, motif, gene, WT, starts_with('mut'), -mutations) %>%
  pivot_longer(
    c(WT, starts_with('mut')),
    names_to = 'variant',
    values_to = 'sequence'
  ) %>%
  filter(! is.na(sequence))

WRKY.high.to.order <- WRKY.high.score %>%
  mutate(
    motif = 'WRKY',
    sequence = `str_sub<-`(sequence, base1, base1, value = 'K'),
    sequence = `str_sub<-`(sequence, base2, base2, value = 'K')
  ) %>%
  select(sp, motif, gene, sequence)


## PIF
PIF.scores <- as_tibble(scan_sequences(filter_motifs(TFs, name = 'TF-cluster_6'), detected.seqs, threshold = 1, RC = FALSE, nthreads = 0)) %>%
  mutate(
    rel.score = (score - min.score) / (max.score - min.score)
  ) %>%
  rename('gene' = sequence) %>%
  inner_join(select(all.sequences, -strand), by = 'gene')

PIF.high.score <- PIF.scores %>%
  filter(rel.score >= 0.7) %>%
  group_by(gene) %>%
  filter(n() == 1) %>%
  ungroup() %>%
  filter(rel.score >= 0.9) %>%
  filter(type == 'protein_coding' & is.na(mutations) & ! grepl('[/;]', gene)) %>%
  mutate(
    base1 = start + 3,
    base2 = start + 5
  ) %>%
  filter(substr(sequence, base1, base1) == 'G' & substr(sequence, base2, base2) == 'G') %>%
  group_by(sp) %>%
  slice_sample(n = 10) %>%
  ungroup()

# mutations: forward: start + 3 G -> T; start + 5 G -> T
PIF.high.sequences <- PIF.high.score %>%
  mutate(
    motif = 'PIF',
    WT = sequence,
    mutA = if_else(strand == '+', `str_sub<-`(sequence, base1, base1, value = 'T'), NA_character_),
    mutB = if_else(strand == '+', `str_sub<-`(sequence, base2, base2, value = 'T'), NA_character_),
    mutAB = if_else(strand == '+', `str_sub<-`(mutA, base2, base2, value = 'T'), NA_character_)
  ) %>%
  select(sp, motif, gene, WT, starts_with('mut'), -mutations) %>%
  pivot_longer(
    c(WT, starts_with('mut')),
    names_to = 'variant',
    values_to = 'sequence'
  )

PIF.high.to.order <- PIF.high.score %>%
  mutate(
    motif = 'PIF',
    sequence = `str_sub<-`(sequence, base1, base1, value = 'K'),
    sequence = `str_sub<-`(sequence, base2, base2, value = 'K')
  ) %>%
  select(sp, motif, gene, sequence)


### combine sequences for first validation set
validation.seqs <- bind_rows(TATA.high.sequences, TATA.low.sequences, TCP.high.sequences, HSF.high.sequences, WRKY.high.sequences, PIF.high.sequences) %>%
  mutate(
    name = paste(gene, motif, variant, sep = '.')
  ) %>%
  select('gene' = name, sp, sequence)

validation.seqs.to.order <- bind_rows(TATA.high.to.order, TATA.low.to.order, TCP.high.to.order, HSF.high.to.order, WRKY.high.to.order, PIF.high.to.order) %>%
  mutate(
    sequence = paste0('agggaagactcgctg', sequence, 'cccgtcgtcttcagg')
  ) %>%
  pull(sequence)

write_tsv(validation.seqs, 'validation_sequences/validation_seqs.tsv')

write_lines(validation.seqs.to.order, 'validation_sequences/validation_sequences.txt')

validation.refseq <- validation.seqs %>%
  mutate(
    line = paste0('>', gene, '\n', sequence)
  ) %>%
  pull(line)

write_lines(validation.refseq, 'validation_sequences/validation_refseq.fa')


### pick sequences for validation (second validation set)

## BREu/d
TATA.scores <- as_tibble(scan_sequences(filter_motifs(CPEs, altname = 'TATA'), detected.seqs, threshold = 1, nthreads = 0)) %>%
  mutate(
    rel.score = (score - min.score) / (max.score - min.score)
  ) %>%
  rename('gene' = sequence) %>%
  filter(rel.score >= 0.85) %>%
  group_by(gene) %>%
  filter(n() == 1) %>%
  ungroup() %>%
  filter(between(0.5 * (start + stop), 107, 143)) %>%
  inner_join(
    select(all.sequences, -strand),
    by = 'gene'
  ) %>%
  mutate(
    upstream = str_sub(paste0('GCTG', sequence), start - 4 + 4, start + 3 + 4),
    downstream = str_sub(paste0(sequence, case_when(sp == 'At' ~ 'TCTC', sp == 'Zm' ~ 'CCCG', sp == 'Sb' ~ 'TACC')), stop - 3, stop + 3)
  )

up.seqs <- Biostrings::DNAStringSet(deframe(select(TATA.scores, gene, upstream)))

BREu.scores <- as_tibble(scan_sequences(filter_motifs(CPEs, altname = 'BREu'), up.seqs, threshold = 1, nthreads = 0)) %>%
  mutate(
    BREu.score = (score - min.score) / (max.score - min.score)
  ) %>%
  rename('gene' = sequence) %>%
  inner_join(
    select(all.sequences, gene, sp),
    by = 'gene'
  )

down.seqs <- Biostrings::DNAStringSet(deframe(select(TATA.scores, gene, downstream)))

BREd.scores <- as_tibble(scan_sequences(filter_motifs(CPEs, altname = 'BREd'), down.seqs, threshold = 1, nthreads = 0)) %>%
  mutate(
    BREd.score = (score - min.score) / (max.score - min.score)
  ) %>%
  rename('gene' = sequence) %>%
  inner_join(
    select(all.sequences, gene, sp),
    by = 'gene'
  )

# mutate BREu
BREu.consensus <- 'AGCGCGCC'
BREu.defective <- 'AGTGCAAC'

remove.BREu <- function(seq, pos) {
  str_sub(seq, pos + 2, pos + 6) <- paste0('T', str_sub(seq, pos + 3, pos + 4), 'AA')
  return(seq)
}

BREu.high.genes <- BREu.scores %>%
  filter(BREu.score >= 0.85) %>%
  group_by(sp) %>%
  slice_sample(n = 10) %>%
  pull(gene)

if (length(BREu.high.genes) < 30) {
  additional.genes <- BREu.scores %>%
    filter(BREu.score >= 0.85 & ! gene %in% BREu.high.genes) %>%
    slice_sample(n = 30 - length(BREu.high.genes)) %>%
    pull(gene)
  
  BREu.high.genes <- c(BREu.high.genes, additional.genes)
}

BREu.high.seqs <- tibble(
  gene = rep(BREu.high.genes, each = 2),
  variant = rep(c('WT', 'mut'), times = length(BREu.high.genes))
) %>%
  inner_join(
    select(all.sequences, gene, sequence, sp),
    by = 'gene'
  ) %>%
  inner_join(
    select(TATA.scores, gene, start),
    by = 'gene'
  ) %>%
  mutate(
    start = start - 4,
    sequence = if_else(variant == 'mut', remove.BREu(sequence, start), sequence),
    motif = 'BREu'
  )

# export motifs
BREu.high.variants <- BREu.high.seqs %>%
  mutate(
    sequence = str_sub(sequence, start, start + 22)
  ) %>%
  select(gene, variant, sequence) %>%
  pivot_wider(
    names_from = variant,
    values_from = sequence
  ) %>%
  select(-gene)

BREu.high.motifs <- lapply(BREu.high.variants, create_motif)

lapply(seq_along(BREu.high.motifs), function(x) PWM.to.LaTeX(BREu.high.motifs[[x]], paste0('../figures/rawData/motif_mutBREu_', names(BREu.high.motifs)[[x]], '.tsv')))


BREu.high.seqs <- BREu.high.seqs %>%
  select(sp, gene, motif, variant, sequence)

# mutate BREu
BREd.consensus <- 'GTTTGTT'
BREd.defective <- 'GATAGAT'

remove.BREd <- function(seq, pos) {
  str_sub(seq, pos + 1, pos + 5) <- paste0('A', str_sub(seq, pos + 2, pos + 2), 'A', str_sub(seq, pos + 4, pos + 4), 'A')
  return(seq)
}

BREd.high.genes <- BREd.scores %>%
  filter(BREd.score >= 0.85) %>%
  group_by(sp) %>%
  slice_sample(n = 10) %>%
  pull(gene)

if (length(BREd.high.genes) < 30) {
  additional.genes <- BREd.scores %>%
    filter(BREd.score >= 0.85 & ! gene %in% BREd.high.genes) %>%
    slice_sample(n = 30 - length(BREd.high.genes)) %>%
    pull(gene)
  
  BREd.high.genes <- c(BREd.high.genes, additional.genes)
}

BREd.high.seqs <- tibble(
  gene = rep(BREd.high.genes, each = 2),
  variant = rep(c('WT', 'mut'), times = length(BREd.high.genes))
) %>%
  inner_join(
    select(all.sequences, gene, sequence, sp),
    by = 'gene'
  ) %>%
  inner_join(
    select(TATA.scores, gene, stop),
    by = 'gene'
  ) %>%
  mutate(
    start = stop - 3,
    sequence = if_else(variant == 'mut', remove.BREd(sequence, start), sequence),
    motif = 'BREd'
  ) 

# export motifs
BREd.high.variants <- BREd.high.seqs %>%
  mutate(
    sequence = str_sub(sequence, start - 16, start + 6)
  ) %>%
  select(gene, variant, sequence) %>%
  pivot_wider(
    names_from = variant,
    values_from = sequence
  ) %>%
  select(-gene)

BREd.high.motifs <- lapply(BREd.high.variants, create_motif)

lapply(seq_along(BREd.high.motifs), function(x) PWM.to.LaTeX(BREd.high.motifs[[x]], paste0('../figures/rawData/motif_mutBREd_', names(BREd.high.motifs)[[x]], '.tsv')))


BREd.high.seqs <- BREd.high.seqs %>%
  select(sp, gene, motif, variant, sequence)


# insert BREu/d
add.BREu <- function(seq, start) {
  str_sub(seq, start, start + 7) <- BREu.consensus
  return(seq)
}

add.BREd <- function(seq, start) {
  str_sub(seq, start, start + 6) <- BREd.consensus
  return(seq)
}

no.BREu <- BREu.scores %>%
  filter(BREu.score < 0.5) %>%
  pull(gene)

no.BREd <- BREd.scores %>%
  filter(BREd.score < 0.5) %>%
  pull(gene)

no.BRE.seqs <- all.sequences %>%
  select(gene, sp, sequence) %>%
  filter(gene %in% no.BREu & gene %in% no.BREd) %>%
  group_by(sp) %>%
  slice_sample(n = 10) %>%
  ungroup() %>%
  inner_join(
    select(TATA.scores, gene, start, stop),
    by = 'gene'
  ) %>%
  mutate(
    BREu.start = start - 4,
    BREd.start = stop - 3,
    variant = list(c('WT', 'BREu', 'BREd'))
  ) %>%
  unnest_longer(variant) %>%
  mutate(
    sequence = case_when(
      variant == 'BREu' ~ add.BREu(sequence, BREu.start),
      variant == 'BREd' ~ add.BREd(sequence, BREd.start),
      variant == 'WT' ~ sequence
    ),
    motif = 'noBRE'
  )

# export motifs
no.BRE.variants <- no.BRE.seqs %>%
  mutate(
    sequence = str_sub(sequence, start - 4, stop + 3),
    variant = if_else(variant == 'WT', variant, paste0('+', variant))
  ) %>%
  select(gene, variant, sequence) %>%
  pivot_wider(
    names_from = variant,
    values_from = sequence
  ) %>%
  select(-gene)

no.BRE.motifs <- lapply(no.BRE.variants, create_motif)

lapply(seq_along(no.BRE.motifs), function(x) PWM.to.LaTeX(no.BRE.motifs[[x]], paste0('../figures/rawData/motif_insBRE_', names(no.BRE.motifs)[[x]], '.tsv')))


no.BRE.seqs <- no.BRE.seqs %>%
  select(sp, gene, motif, variant, sequence)

all.BRE.seqs <- bind_rows(BREu.high.seqs, BREd.high.seqs, no.BRE.seqs)


## pick strong/intermediate/weak promoters for evolution
promoter.strength.mod <- promoter.strength %>%
  filter(is.na(mutations)) %>%
  filter(enhancer & ! light) %>%
  filter(n.experiments == 2) %>%
  select(sys, sp, gene, enrichment) %>%
  group_by(sys) %>%
  mutate(
    strength = case_when(
      between(enrichment, quantile(enrichment, 0.8), quantile(enrichment, 1)) ~ 'strong',
      between(enrichment, quantile(enrichment, 0.4), quantile(enrichment, 0.6)) ~ 'intermediate',
      between(enrichment, quantile(enrichment, 0), quantile(enrichment, 0.2)) ~ 'weak',
      TRUE ~ 'in between'
    )
  ) %>%
  filter(strength != 'in between') %>%
  select(-enrichment) %>%
  pivot_wider(
    names_from = sys,
    values_from = strength
  ) %>%
  mutate(
    level = if_else(leaf == proto, leaf, paste(leaf, proto, sep = '/'))
  ) %>%
  filter(
    ! is.na(level)
  )

strong.genes <- promoter.strength.mod %>%
  filter(level == 'strong') %>%
  group_by(sp) %>%
  slice_sample(n = 10) %>%
  pull(gene)

intermediate.genes <- promoter.strength.mod %>%
  filter(level == 'intermediate') %>%
  group_by(sp) %>%
  slice_sample(n = 10) %>%
  pull(gene)

weak.genes <- promoter.strength.mod %>%
  filter(level == 'weak') %>%
  group_by(sp) %>%
  slice_sample(n = 10) %>%
  pull(gene)

weak.strong.genes <- promoter.strength.mod %>%
  filter(level == 'weak/strong') %>%
  group_by(sp) %>%
  slice_sample(n = 10) %>%
  pull(gene)

if (length(weak.strong.genes) < 30) {
  additional.genes <- promoter.strength.mod %>%
    filter(level == 'weak/strong' & ! gene %in% weak.strong.genes) %>%
    slice_sample(n = 30 - length(weak.strong.genes)) %>%
    pull(gene)
  
  weak.strong.genes <- c(weak.strong.genes, additional.genes)
}

strong.weak.genes <- promoter.strength.mod %>%
  filter(level == 'strong/weak') %>%
  group_by(sp) %>%
  slice_sample(n = 10) %>%
  pull(gene)

if (length(strong.weak.genes) < 30) {
  additional.genes <- promoter.strength.mod %>%
    filter(level == 'strong/weak' & ! gene %in% strong.weak.genes) %>%
    slice_sample(n = 30 - length(strong.weak.genes)) %>%
    pull(gene)
  
  strong.weak.genes <- c(strong.weak.genes, additional.genes)
}

native.promoters <- c(strong.genes, intermediate.genes, weak.genes, strong.weak.genes, weak.strong.genes)

native.promoter.seqs <- all.sequences %>%
  filter(gene %in% native.promoters) %>%
  left_join(select(promoter.strength.mod, gene, level), by = 'gene') %>%
  mutate(
    motif = paste0('native(', level, ')'),
    variant = 'WT'
  ) %>%
  select(sp, gene, motif, variant, sequence)


## generate synthetic promoters
freqs.At <- all.sequences %>%
  filter(sp == 'At') %>%
  mutate(
    A = 1 - (nchar(gsub('A', '', sequence)) / nchar(sequence)),
    C = 1 - (nchar(gsub('C', '', sequence)) / nchar(sequence)),
    G = 1 - (nchar(gsub('G', '', sequence)) / nchar(sequence)),
    T = 1 - (nchar(gsub('T', '', sequence)) / nchar(sequence))
  ) %>%
  summarise(
    across(c(A, C, G, T), mean)
  ) %>%
  pivot_longer(
    everything()
  ) %>%
  deframe()

freqs.Zm <- all.sequences %>%
  filter(sp == 'Zm') %>%
  mutate(
    A = 1 - (nchar(gsub('A', '', sequence)) / nchar(sequence)),
    C = 1 - (nchar(gsub('C', '', sequence)) / nchar(sequence)),
    G = 1 - (nchar(gsub('G', '', sequence)) / nchar(sequence)),
    T = 1 - (nchar(gsub('T', '', sequence)) / nchar(sequence))
  ) %>%
  summarise(
    across(c(A, C, G, T), mean)
  ) %>%
  pivot_longer(
    everything()
  ) %>%
  deframe()

backgrounds.At <- create_sequences(seqnum = 100000, seqlen = 170, freqs = freqs.At)
backgrounds.Zm <- create_sequences(seqnum = 100000, seqlen = 170, freqs = freqs.Zm)

# filter out sequences with a BsaI or BbsI site
RE.sites <- c('(G|^)GTCT(C|$)', '(G|^)AGA(CC|C$|$)', '(G|^)AAGA(C|$)', '(G|^)TCTT(C|$)')

backgrounds.At <- backgrounds.At[! grepl(paste0('(', paste0(RE.sites, collapse = ')|('), ')'), backgrounds.At)]
backgrounds.Zm <- backgrounds.Zm[! grepl(paste0('(', paste0(RE.sites, collapse = ')|('), ')'), backgrounds.Zm)]

# filter out sequences with high scores for TATA, Inr, Ypatch, or some TFs (NAC, TCP, and HSF)
TATA <- filter_motifs(CPEs, altname = 'TATA')[[1]]
Inr <- filter_motifs(CPEs, altname = 'Inr')[[1]]
Ypatch <- filter_motifs(CPEs, altname = 'Ypatch')[[1]]

thresh <- 0.75

for (species in c('At', 'Zm')) {
  sequences <- get(paste0('backgrounds.', species))
  
  for (motif in c('TATA', 'Inr', 'Ypatch')) {
    if (motif == 'Inr') {
      test.seqs <- Biostrings::DNAStringSet(sequences, start = 151, end = 170)
    } else {
      test.seqs <- sequences
    }
    
    no.motif <- as_tibble(scan_sequences(get(motif), test.seqs, threshold = 1, nthreads = 0)) %>%
      group_by('gene' = sequence) %>%
      filter(score == max(score)) %>%
      ungroup() %>%
      mutate(
        score = (score - min.score) / (max.score - min.score)
      ) %>%
      filter(score < thresh) %>%
      pull(sequence) %>%
      as.numeric()
    
    sequences <- sequences[no.motif]
  }
  
  for (TF in c('TF-cluster_1', 'TF-cluster_15', 'TF-cluster_16', 'TF-cluster_22')) {
    test.seqs <- sequences
    
    no.motif <- as_tibble(scan_sequences(filter_motifs(TFs, name = TF), test.seqs, threshold = 1, nthreads = 0)) %>%
      group_by('gene' = sequence) %>%
      filter(score == max(score)) %>%
      ungroup() %>%
      mutate(
        score = (score - min.score) / (max.score - min.score)
      ) %>%
      filter(score < thresh) %>%
      pull(sequence) %>%
      as.numeric()
    
    sequences <- sequences[no.motif]
  }
  
  assign(paste0('backgrounds.filt.', species), sequences)
}

backgrounds.At <- sample(as.character(backgrounds.filt.At), 10)
backgrounds.Zm <- sample(as.character(backgrounds.filt.Zm), 10)

# functions to add core promoter elements
add.initiator <- function(seq) {
  str_sub(seq, 160, 169) <- gsub('G', 'C', str_sub(seq, 160, 169))
  str_sub(seq, 160, 169) <- gsub('A', 'T', str_sub(seq, 160, 169))
  str_sub(seq, 164, 166) <- 'TCA'
  return(seq)
}

add.TATA <- function(seq) {
  str_sub(seq, 133, 140) <- 'TATAAATA'
  return(seq)
}

add.Ypatch <- function(seq) {
  str_sub(seq, 147, 154) <- gsub('A|G', 'C', str_sub(seq, 147, 154))
  return(seq)
}

# create synthetic promoters
synthetic.promoters <- tibble(
  'At' = backgrounds.At,
  'Zm' = backgrounds.Zm,
  'bkg.id' = seq_along(backgrounds.At)
) %>%
  pivot_longer(
    -bkg.id,
    names_to = 'bkg.sp',
    values_to = 'WT'
  ) %>%
  mutate(
    TATA = add.TATA(WT),
    Inr = add.initiator(WT),
    Ypatch = add.Ypatch(WT),
    TATA.Inr = add.initiator(TATA),
    TATA.Ypatch = add.Ypatch(TATA),
    Inr.Ypatch = add.Ypatch(Inr),
    TATA.Inr.Ypatch = add.Ypatch(TATA.Inr)
  ) %>% pivot_longer(
    -starts_with('bkg.'),
    names_to = 'variant',
    values_to = 'sequence'
  ) %>%
  mutate(
    gene = paste0('syn', str_sub(bkg.sp, 1, 1), bkg.id),
    variant = gsub('.', '+', variant, fixed = TRUE),
    variant = if_else(variant == 'WT', variant, paste0('+', variant)),
    motif = 'synthetic'
  ) %>%
  select('sp' = bkg.sp, gene, motif, variant, sequence)

# export motifs
synthetic.variants <- synthetic.promoters %>%
  filter(variant %in% c('+Inr', '+Ypatch', '+TATA')) %>%
  mutate(
    variant = sub('+', '', variant, fixed = TRUE),
    sequence = case_when(
      variant == 'Inr' ~ str_sub(sequence, 160, 169),
      variant == 'Ypatch' ~ str_sub(sequence, 147, 154),
      variant == 'TATA' ~ str_sub(sequence, 133, 140)
    )
  ) %>%
  select(gene, variant, sequence) %>%
  pivot_wider(
    names_from = variant,
    values_from = sequence
  ) %>%
  select(-gene)

synthetic.motifs <- lapply(synthetic.variants, create_motif)

lapply(seq_along(synthetic.motifs), function(x) PWM.to.LaTeX(synthetic.motifs[[x]], paste0('../figures/rawData/motif_synPRO_', names(synthetic.motifs)[[x]], '.tsv')))


## export promoters for evolution
promoters.evolution <- synthetic.promoters %>%
  mutate(
    name = paste0(gene, '(', variant, ')')
  ) %>%
  select(name, sequence) %>%
  bind_rows(select(native.promoter.seqs, 'name' = gene, sequence))

write_tsv(promoters.evolution, 'validation_sequences/promoters_for_evolution.tsv')


### synthetic promoters with added TF binding sites
# TF binding sites will be inserted into synthetic promoters containing a TATA-box

## functions to add TF binding sites
# NAC (TF-cluster 1)
add.NAC <- function(seq, pos) {
  str_sub(seq, pos - 7, pos - 1) <- 'TTACGTG'
  str_sub(seq, pos + 4, pos + 8) <- 'ACAAG'
  return(seq)
}
# TCP 1 (TF-cluster 15)
add.TCP1 <- function(seq, pos) {
  str_sub(seq, pos - 4, pos + 5) <- 'TGGGGCCCAC'
  return(seq)
}
# HSF (TF-cluster 16)
add.HSF <- function(seq, pos) {
  str_sub(seq, pos - 6, pos + 6) <- 'GAAGCTTCTAGAA'
  return(seq)
}
# TCP 2 (TF-cluster 22)
add.TCP2 <- function(seq, pos) {
  str_sub(seq, pos - 3, pos + 4) <- 'GGGACCAC'
  return(seq)
}
# wrapper function
add.TF <- function(seq, TF, pos) {
  if (TF != 'none') {
    add.fun <- get(paste0('add.', TF))
    seq <- add.fun(seq, pos)
  }
  return(seq)
}


## scan TF positons
TF.pos.scan <- synthetic.promoters %>%
  filter(variant == '+TATA') %>%
  mutate(
    motif = 'TFscan',
    pos = list(c(10, 30, 40, 50, 60, 70, 80, 90, 100, 105, 110, 115, 120, 150, 155, 160))
  ) %>%
  unnest_longer(
    pos
  ) %>%
  mutate(
    NAC = add.NAC(sequence, pos),
    TCP1 = add.TCP1(sequence, pos),
    TCP2 = add.TCP2(sequence, pos),
    HSF = add.HSF(sequence, pos)
  ) %>%
  select(-sequence) %>%
  pivot_longer(
    c(NAC, TCP1, TCP2, HSF),
    names_to = 'TF',
    values_to = 'sequence'
  ) %>%
  mutate(
    variant = paste(TF, pos, sep = '_')
  ) %>%
  select(sp, gene, motif, variant, sequence)


## TF combinations
test.TFs <- c('none', 'NAC', 'TCP1', 'TCP2', 'HSF')

TF.combos <- expand_grid(TF1 = test.TFs, TF2 = test.TFs, TF3 = test.TFs) %>%
  filter(! (TF1 == 'none' & TF2 == 'none' & TF3 == 'none')) %>%
  mutate(
    combo = paste(TF1, TF2, TF3, sep = '-')
  ) %>%
  pull (combo)

TF.combo.seqs <- synthetic.promoters %>%
  filter(variant == '+TATA') %>%
  mutate(
    motif = 'TFcombo',
    variant = list(TF.combos)
  ) %>%
  unnest_longer(
    variant
  ) %>%
  rowwise() %>%
  mutate(
    sequence = add.TF(sequence, unlist(str_split(variant, '-'))[c(TRUE, FALSE, FALSE)], 35),
    sequence = add.TF(sequence, unlist(str_split(variant, '-'))[c(FALSE, TRUE, FALSE)], 65),
    sequence = add.TF(sequence, unlist(str_split(variant, '-'))[c(FALSE, FALSE, TRUE)], 95)
  ) %>%
  ungroup()


## fix restriction enzyme sites in TF sequences
RE.sites <- paste0('(', paste0(c('(G|^)GTCT(C|$)', '(G|^)AGA(CC|C$|$)', '(G|^)AAGA(C|$)', '(G|^)TCTT(C|$)'), collapse = ')|('), ')')

# manual inspection revealed that for the original sequences, only the BbsI site GAAGAC was introduced and in all cases it can be resolved by mutating a G in the background to C
all.TF.seqs <- bind_rows(TF.pos.scan, TF.combo.seqs) %>%
  mutate(
    sequence = if_else(grepl(RE.sites, sequence), gsub('GAAGAC', 'GAACAC', sequence), sequence)
  )


### combine sequences for first validation set
validation.seqs.two <- bind_rows(all.BRE.seqs, native.promoter.seqs, synthetic.promoters, all.TF.seqs)%>%
  mutate(
    name = paste(gene, motif, variant, sep = '.')
  )

write_tsv(validation.seqs.two, 'validation_sequences/validation_seqs_two.tsv')


### combine with evolved promoters and save sequences to order
# the in silico evolution has to be performed before this code chunk
# this can be performed independent of the above code (after loading the required libraries) to prevent a change in the selected sequences
validation.seqs.two <- read_tsv('validation_sequences/validation_seqs_two.tsv')

native.sp <- validation.seqs.two %>%
  filter(grepl('native', motif, fixed = TRUE)) %>%
  select(gene, sp) %>%
  deframe()

validation.seqs.two.to.order <- read_tsv('../CNN/evolution_data.tsv') %>%
  filter(opt_for != 'start') %>%
  group_by(origin, opt_for) %>%
  filter(round %in% c(3, max(round))) %>%
  ungroup() %>%
  bind_rows(validation.seqs.two) %>%
  distinct(sequence) %>%
  mutate(
    sequence = paste0('agggaagactcgctg', sequence, 'cccgtcgtcttcagg')
  ) %>%
  pull()

write_lines(validation.seqs.two.to.order, 'validation_sequences/validation_sequences_two.txt')

validation.seqs.evo <- read_tsv('../CNN/evolution_data.tsv') %>%
  filter(opt_for != 'start') %>%
  group_by(origin, opt_for) %>%
  filter(round %in% c(3, max(round))) %>%
  ungroup() %>%
  rename('gene' = origin) %>%
  group_by(round, gene, sequence, pred_leaf, pred_proto) %>%
  summarise(
    opt_for = if_else(n() == 1, first(opt_for), paste0(opt_for, collapse = '/'))
  ) %>%
  ungroup() %>%
  mutate(
    variant = paste(opt_for, round, sep = '-'),
    motif = 'evolution',
    name = paste(gene, motif, variant, sep = '.')
  ) %>%
  bind_rows(validation.seqs.two) %>%
  mutate(
    sp = case_when(
      ! grepl('evolution', motif, fixed = TRUE) ~ sp,
      grepl('synA', gene, fixed = TRUE) ~ 'At',
      grepl('synZ', gene, fixed = TRUE) ~ 'Zm',
      TRUE ~ native.sp[gene]
    )
  ) %>%
  select('gene' = name, sp, sequence)

write_tsv(validation.seqs.evo, 'validation_sequences/validation_seqs_evo.tsv')

validation.refseq.two <- validation.seqs.evo %>%
  mutate(
    line = paste0('>', gene, '\n', sequence)
  ) %>%
  pull(line)

write_lines(validation.refseq.two, 'validation_sequences/validation_refseq_two.fa')

write_lines(c(validation.refseq, validation.refseq.two), 'validation_sequences/validation_refseq_all.fa')
