library(dplyr)
library(tidyr)
library(readr)
library(tibble)

# load files (the 'eglab_*' files were obtained from Aimer Gutierrez Diaz from the Grotewold lab; Mejía-Guerra et al., 2015, Plant Cell, doi: 10.1105/tpc.15.00630)
shoot <- read_tsv('eglab_dtss_agpv4_shoot.pc.gff', comment = '#', col_names = c('chromosome', 'source', 'type', 'start', 'stop', 'score', 'strand', 'phase', 'annotation')) %>%
  mutate(
    gene = substr(annotation, 19, 32)
  )

root <- read_tsv('eglab_dtss_agpv4_root.pc.gff', comment = '#', col_names = c('chromosome', 'source', 'type', 'start', 'stop', 'score', 'strand', 'phase', 'annotation')) %>%
  mutate(
    gene = substr(annotation, 19, 32)
  )

reference <- read_tsv('Maize_protein_coding_genes.gff3', comment = '#', col_names = c('chromosome', 'source', 'type', 'start', 'stop', 'score', 'strand', 'phase', 'annotation')) %>%
  mutate(
    gene = substr(annotation, 9, 22)
  )

# merge shoot and root dTSSs
correctedTSS <- union(shoot, root)

# get reference TSSs
referenceTSS <- reference %>%
  mutate(
    TSS = if_else(strand == '+', start, stop)
  ) %>%
  select(gene, TSS) %>%
  deframe()

# extract dTSSs linked to 2+ genes and link them to the closest reference TSS
multipleTSS <- correctedTSS %>%
  filter(grepl(',', annotation, fixed = TRUE)) %>%
  rename('gene1' = gene) %>%
  mutate(
    gene2 = substr(annotation, 34, 47),
    gene3 = substr(annotation, 49, 62),
    gene4 = substr(annotation, 64, 77)
  ) %>%
  pivot_longer(
    cols = starts_with('gene'),
    values_to = 'gene',
    names_to = 'gene_id',
    names_prefix = 'gene'
  ) %>%
  filter(grepl('Zm00001d', gene, fixed = TRUE)) %>%
  mutate(
    refTSS = referenceTSS[gene]
  ) %>%
  group_by(across(c(-gene, -gene_id, -refTSS))) %>%
  summarise(
    gene = first(gene[abs(stop - refTSS) == min(abs(stop - refTSS))])
  ) %>%
  ungroup() %>%
  select(annotation, gene) %>%
  deframe()

correctedTSS <- correctedTSS %>%
  mutate(
    gene = if_else(annotation %in% names(multipleTSS), multipleTSS[annotation], gene)
  ) %>%
  select(gene, stop) %>%
  deframe()
  
# export names of genes with a dTSS in the Grotewold data
write_lines(names(correctedTSS), 'Maize_corrected_TSS.txt')

# correct TSS annotation in B73 reference annotation
TSS <- reference %>%
  mutate(
    start = if_else(strand == '+' & gene %in% names(correctedTSS), correctedTSS[gene], start),
    stop = if_else(strand == '-' & gene %in% names(correctedTSS), correctedTSS[gene], stop)
  ) %>% 
  select(-gene)

# save corrected annotation to file
write_tsv(TSS, 'Maize_protein_coding_genes.gff3', col_names = FALSE)
