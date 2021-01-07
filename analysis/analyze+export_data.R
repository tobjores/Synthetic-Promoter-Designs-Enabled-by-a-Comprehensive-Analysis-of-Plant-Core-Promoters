library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(stringr)
library(gprofiler2)
library(universalmotif)
library(openxlsx)

### create folder for data export
if (! file.exists('../figures/rawData')){
  dir.create('../figures/rawData')
}

### set species and other names
all.species <- tibble(
  common = c('Maize', 'Sorghum', 'Arabidopsis'),
  short = c('Zm', 'Sb', 'At'),
  latin = c('zmays', 'sbicolor', 'athaliana')
)

short.common <- deframe(select(all.species, short, common))
short.latin <- deframe(select(all.species, short, latin))

all.systems <- tibble(
  short = c('leaf', 'proto'),
  long = c('tobacco leaves', 'maize protoplasts')
)

system.name <- deframe(all.systems)
enhancer.name <- c(`FALSE` = 'noENH', `TRUE` = 'withENH')
light.name <- c(`FALSE` = 'dark', `TRUE` = 'light')

sys.sp <- expand_grid(species = all.species$short, system = all.systems$short)

# threshold for motif hits
thresh <- 0.85


### load main data and group promoters by GC content
load('RData/main_data.Rdata')

promoter.strength <- promoter.strength %>%
  group_by(sp) %>%
  mutate(
    GC.bins = as.character(cut_number(GC, 5, dig.lab = 2)),
    GC.low = as.numeric(unlist(strsplit(as.character(GC.bins), '\\(|\\[|,|\\)|\\]'))[seq(2, 3 * length(GC.bins), 3)]),
    GC.high = as.numeric(unlist(strsplit(as.character(GC.bins), '\\(|\\[|,|\\)|\\]'))[seq(3, 3 * length(GC.bins), 3)]),
    GC.label = paste(GC.low, GC.high, sep = '\nto\n')
  ) %>%
  ungroup()

enhancer.effect <- enhancer.effect %>%
  group_by(sp) %>%
  mutate(
    GC.bins = as.character(cut_number(GC, 5, dig.lab = 2)),
    GC.low = as.numeric(unlist(strsplit(as.character(GC.bins), '\\(|\\[|,|\\)|\\]'))[seq(2, 3 * length(GC.bins), 3)]),
    GC.high = as.numeric(unlist(strsplit(as.character(GC.bins), '\\(|\\[|,|\\)|\\]'))[seq(3, 3 * length(GC.bins), 3)]),
    GC.label = paste(GC.low, GC.high, sep = '\nto\n')
  ) %>%
  ungroup()

light.dependency <- light.dependency %>%
  group_by(sp) %>%
  mutate(
    GC.bins = as.character(cut_number(GC, 5, dig.lab = 2)),
    GC.low = as.numeric(unlist(strsplit(as.character(GC.bins), '\\(|\\[|,|\\)|\\]'))[seq(2, 3 * length(GC.bins), 3)]),
    GC.high = as.numeric(unlist(strsplit(as.character(GC.bins), '\\(|\\[|,|\\)|\\]'))[seq(3, 3 * length(GC.bins), 3)]),
    GC.label = paste(GC.low, GC.high, sep = '\nto\n')
  ) %>%
  ungroup()

load('RData/main_data_validation.Rdata')


### load promoter sequences
all.sequences <- bind_rows(
  'At' = read_tsv('../data/promoter_annotation/Arabidopsis_all_promoters_unique.tsv'),
  'Zm' = read_tsv('../data/promoter_annotation/Maize_all_promoters_unique.tsv'),
  'Sb' = read_tsv('../data/promoter_annotation/Sorghum_all_promoters_unique.tsv'),
  .id = 'sp'
) %>%
  mutate(
    GC = nchar(gsub('A|T', '', sequence)) / nchar(sequence)
  )

all.seqs <- Biostrings::DNAStringSet(deframe(select(all.sequences, gene, sequence)))


### load CPE and TF motifs and scores
CPEs <- read_meme('../data/misc/CPEs.meme')
TFs <- read_meme('../data/misc/TF-clusters.meme')

motif.scores <- read_tsv('motifs/motif-scores.tsv.gz')


### load export functions
source('export_functions.R')


### data export

## internal controls (Figure 1b)
load('RData/controls.Rdata')

controls.filtered <- controls %>%
  filter(grepl('withPRO', variant, fixed = TRUE)) %>%
  filter(enhancer & ! light) %>%
  mutate(
    control = factor(
      variant,
      levels = c('control-withPRO-noENH', 'control-withPRO-withENH'),
      labels = c('noENH', 'withENH'),
      ordered = TRUE
    )
  ) %>%
  group_by(across(c(-count.out, -count.inp, -rep, -enrichment))) %>%
  summarise(
    enrichment = mean(enrichment)
  ) %>%
  ungroup()

for (system in all.systems$short) {
  df <- controls.filtered %>%
    filter(sys == system) %>%
    mutate(
      sample.name = factor(
        paste(sp, control, sep = '_'),
        levels = c(paste(all.species$short, 'noENH', sep = '_'), paste(all.species$short, 'withENH', sep = '_')),
        ordered = TRUE
      )
    ) %>%
    select(enrichment, sample.name) %>%
    group_by(sample.name) %>%
    mutate(
     id = seq_len(n())
    ) %>%
    pivot_wider(
      names_from = sample.name,
      values_from = enrichment
    ) %>%
    select(-id)

  write_tsv(df, paste0('../figures/rawData/enrichment_controls_', system, '.tsv'), na = 'NaN')
}


## replicate correlation (Figure 1c,d and Supp Figure 1)
load('RData/ag_experiments.Rdata')

ag.experiments.mod <- ag.experiments %>%
  select(sys, sp, enhancer, light, gene, rep, enrichment) %>%
  mutate(
    enhancer = enhancer.name[as.character(enhancer)],
    light = light.name[as.character(light)]
  )

for (i in 1:nrow(sys.sp)) {
  system <- sys.sp[i, 'system', drop = TRUE]
  species <- sys.sp[i, 'species', drop = TRUE]

  for (condition in light.name) {

    if (condition == 'light' & system == 'proto') {
      next
    }

    df <- ag.experiments.mod %>%
      filter(sys == system & sp == species & light == condition) %>%
      mutate(
        sample.name = enhancer
      ) %>%
      select(gene, sample.name, rep, enrichment) %>%
      pivot_wider(
        names_from = rep,
        names_prefix = 'rep',
        values_from = enrichment
      ) %>%
      slice_sample(prop = 1)

    correlation <- df %>%
      group_by(sample.name) %>%
      summarise(
        spearman = sprintf('$\\rho=%.2f$', cor(rep1, rep2, method = 'spearman', use = 'complete.obs')),
        rsquare = sprintf('$R^2=%.2f$', cor(rep1, rep2, method = 'pearson', use = 'complete.obs')^2)
      )

    write_tsv(select(df, -gene), paste0('../figures/rawData/enrichment_correlation_', species, '_', if_else(system == 'leaf', condition, system), '.tsv'), na = 'NaN')
    write_tsv(correlation, paste0('../figures/rawData/enrichment_correlation_', species, '_', if_else(system == 'leaf', condition, system), '_stats.tsv'), na = 'NaN')
  }
}


## comparison promoter strength in tobacco leaves vs maize protoplasts (Figure 1e and Supp Figure 1)
promoter.strength.filtered <- promoter.strength %>%
  filter(! light) %>%
  select(sys, sp, enhancer, gene, enrichment) %>%
  mutate(
    enhancer = enhancer.name[as.character(enhancer)]
  )

for (species in all.species$short) {
  df <- promoter.strength.filtered %>%
    filter(sp == species) %>%
    mutate(
      sample.name = enhancer
    ) %>%
    select(gene, sample.name, sys, enrichment) %>%
    pivot_wider(
      names_from = sys,
      values_from = enrichment
    ) %>%
    slice_sample(prop = 1)

  correlation <- df %>%
    group_by(sample.name) %>%
    summarise(
      spearman = sprintf('$\\rho=%.2f$', cor(leaf, proto, method = 'spearman', use = 'complete.obs')),
      rsquare = sprintf('$R^2=%.2f$', cor(leaf, proto, method = 'pearson', use = 'complete.obs')^2)
    ) %>%
    ungroup()

  write_tsv(select(df, -gene), paste0('../figures/rawData/enrichment_leaf-vs-proto_', species, '.tsv'), na = 'NaN')
  write_tsv(correlation, paste0('../figures/rawData/enrichment_leaf-vs-proto_', species, '_stats.tsv'), na = 'NaN')
}


## promoter strength without and with enhancer in dark (Figure 1f,g)
promoter.strength.filtered <- promoter.strength %>%
  filter(! light) %>%
  select(sys, sp, enhancer, gene, enrichment) %>%
  mutate(
    enhancer = enhancer.name[as.character(enhancer)]
  )

for (i in 1:nrow(sys.sp)) {
  system <- sys.sp[i, 'system', drop = TRUE]
  species <- sys.sp[i, 'species', drop = TRUE]

  df <- promoter.strength.filtered %>%
    filter(sys == system & sp == species) %>%
    mutate(
      sample.name = factor(
        enhancer,
        levels = c('noENH', 'withENH'),
        ordered = TRUE
      ),
      LaTeX.label = if_else(enhancer == 'withENH', '$+$', '$-$')
    )

  LaTeX.violinplot(
    data = df,
    samples_from = sample.name,
    values_from = enrichment,
    file = paste('../figures/rawData/enrichment_enh', species, system, sep = '_'),
    LaTeX.label,
    outliers = FALSE,
    p_values = FALSE,
    half = FALSE,
    data_range = range(reps.mean.filt %>% filter(! light) %>% pull(enrichment))
  )
}


## GO-term enrichments (Figure 1h and 5c)
gmt.IDs <- read_tsv('GO-terms/gprofiler_gmt_IDs.tsv')
n.genes <- 1000

# GO-terms promoter strength without enhancer in dark (Figure 1h)
for (i in 1:nrow(sys.sp)) {
  system <- sys.sp[i, 'system', drop = TRUE]
  species <- sys.sp[i, 'species', drop = TRUE]
  
  gmt.ID <- gmt.IDs %>%
    filter(sp == short.latin[species]) %>%
    pull(ID)
  
  promoter.strength.filtered <- promoter.strength %>%
    filter(sys == system & sp == species & ! enhancer & ! light)
  
  all.genes <- promoter.strength.filtered %>%
    distinct(gene) %>%
    pull()
  all.genes <- unlist(strsplit(all.genes, ';|/'))
  
  high.strength <- promoter.strength.filtered %>%
    arrange(desc(enrichment)) %>%
    distinct(gene) %>%
    pull()
  high.strength <- unlist(strsplit(high.strength, ';|/'))[1:n.genes]
  
  high.gost <- gost(high.strength, organism = gmt.ID, custom_bg = all.genes, significant = FALSE)
  
  pvalues <- tibble(
      term_id = high.gost$result$term_id,
      term_name = high.gost$result$term_name,
      p_value = high.gost$result$p_value,
      significant = high.gost$result$significant
    ) %>%
    mutate(
      sp = species,
      sys = system
    )
  
  assign(paste('pvalues', species, system, sep = '.'), pvalues)
}

pvalues.LaTeX <- bind_rows(mget(paste('pvalues', rep(all.species$short, each = 2), rep(all.systems$short, times = 3), sep = '.'))) %>%
  filter(term_id %in% c('GO:0005576', 'GO:0000786', 'GO:0006950', 'GO:0016491', 'GO:0051082')) %>%
  mutate(
    p_value = -log10(p_value),
    p = if_else(significant, p_value, NA_real_),
    ns = if_else(significant, NA_real_, p_value),
    term_id = factor(term_id, levels = c('GO:0005576', 'GO:0000786', 'GO:0006950', 'GO:0016491', 'GO:0051082'), ordered = TRUE)
  ) %>%
  select(-p_value, -significant) %>%
  pivot_wider(
    names_from = sp,
    values_from = c(p, ns)
  ) %>%
  arrange(
    term_id
  )

write_tsv(pvalues.LaTeX %>% filter(sys == 'leaf') %>% select(-sys), '../figures/rawData/GO_enrichment_leaf.tsv', na = 'NaN')
write_tsv(pvalues.LaTeX %>% filter(sys == 'proto') %>% select(-sys), '../figures/rawData/GO_enrichment_proto.tsv', na = 'NaN')

# GO-terms light-dependency without enhancer (Figure 5c)
for (species in all.species$short) {
  gmt.ID <- gmt.IDs %>%
    filter(sp == short.latin[species]) %>%
    pull(ID)
  
  light.dependency.filtered <- light.dependency %>%
    filter(sp == species & ! enhancer)
  
  all.genes <- light.dependency.filtered %>%
    distinct(gene) %>%
    pull()
  all.genes <- unlist(strsplit(all.genes, ';|/'))
  
  high.light <- light.dependency.filtered %>%
    arrange(desc(light.dependency)) %>%
    distinct(gene) %>%
    pull()
  high.light <- unlist(strsplit(high.light, ';|/'))[1:n.genes]
  
  high.gost <- gost(high.light, organism = gmt.ID, custom_bg = all.genes, significant = FALSE)
  
  pvalues <- tibble(
      term_id = high.gost$result$term_id,
      term_name = high.gost$result$term_name,
      p_value = high.gost$result$p_value,
      significant = high.gost$result$significant
    ) %>%
    mutate(
      sp = species
    )
  
  assign(paste('pvalues', species, sep = '.'), pvalues)
}

pvalues.LaTeX <- bind_rows(mget(paste('pvalues', all.species$short, sep = '.'))) %>%
  filter(term_id %in% c('GO:0009579', 'GO:0009536', 'GO:0015979')) %>%
  mutate(
    p_value = -log10(p_value),
    p = if_else(significant, p_value, NA_real_),
    ns = if_else(significant, NA_real_, p_value),
    term_id = factor(term_id, levels = c('GO:0009579', 'GO:0009536', 'GO:0015979'), ordered = TRUE)
  ) %>%
  select(-p_value, -significant) %>%
  pivot_wider(
    names_from = sp,
    values_from = c(p, ns)
  ) %>%
  arrange(
    term_id
  )

write_tsv(pvalues.LaTeX, '../figures/rawData/GO_light-dep.tsv', na = 'NaN')


## promoter strength by type without enhancer in dark (Figure 1i,j)
promoter.strength.filtered <- promoter.strength %>%
  filter(! light & ! enhancer)

for (i in 1:nrow(sys.sp)) {
  system <- sys.sp[i, 'system', drop = TRUE]
  species <- sys.sp[i, 'species', drop = TRUE]

  df <- promoter.strength.filtered %>%
    filter(sys == system & sp == species) %>%
    mutate(
      sample.name = factor(
        case_when(
          type == 'miRNA' ~ 'miRNA',
          type != 'miRNA' & UTR ~ 'withUTR',
          TRUE ~ 'noUTR'
        ),
        levels = c('miRNA', 'withUTR', 'noUTR'),
        ordered = TRUE
      ),
      LaTeX.label = case_when(
        sample.name == 'miRNA' ~ 'miRNA',
        sample.name == 'withUTR' ~ 'protein coding\\\\with 5\\textprime\\ UTR',
        TRUE ~ 'protein coding\\\\without 5\\textprime\\ UTR'
      )
    )

  LaTeX.violinplot(
    data = df,
    samples_from = sample.name,
    values_from = enrichment,
    file = paste('../figures/rawData/enrichment_type', species, system, sep = '_'),
    LaTeX.label,
    outliers = FALSE,
    p_values = TRUE,
    data_range = range(promoter.strength.filtered %>% filter(sys == system) %>% pull(enrichment))
  )
}


## density GC content (Figure 2a)
GC.density <- all.sequences %>%
  group_by(sp) %>%
  summarise(
    x = list(density(GC, from = min(GC), to = max(GC))$x),
    y = list(density(GC, from = min(GC), to = max(GC))$y)
  ) %>%
  pivot_wider(
    names_from = sp,
    values_from = c(x, y),
    names_sep = '.'
  ) %>%
  unnest(everything())


GC.average <- all.sequences %>%
  group_by(sp) %>%
  summarise(
    promoters = mean(GC)
  ) %>%
  left_join(
    tibble(
      sp = c('At', 'Zm', 'Sb'),
      genome = c(0.3606, 0.4686, 0.4384),
      y = list(c(0, 10))
    ),
    by = 'sp'
  ) %>%
  pivot_wider(
    names_from = sp,
    values_from = c(promoters, genome),
    names_sep = '.'
  ) %>%
  unnest(y)

write_tsv(GC.density, '../figures/rawData/GC_density.tsv', na = 'NaN')
write_tsv(GC.average, '../figures/rawData/GC_average.tsv', na = 'NaN')


## correlation between GC content and promoter strength by position (Figure 2c)
window.size <- 10
n.windows <- nchar(all.sequences$sequence[1]) - window.size + 1


all.sequences.mod <- all.sequences %>%
  rowwise() %>%
  mutate(
    GC = list(
      nchar(gsub('A|T', '', str_sub(sequence, seq(from = 1, length.out = n.windows), seq(from = window.size, length.out = n.windows)))) / window.size
    ),
    pos = list(seq(from = 0.5 * (window.size + 1), length.out = n.windows))
  )

correlation <- promoter.strength %>%
  filter(! enhancer & ! light) %>%
  select(-starts_with('GC')) %>%
  inner_join(select(all.sequences.mod, gene, GC, pos), by = 'gene') %>%
  unnest(c(GC, pos)) %>%
  group_by(sp, sys, pos) %>%
  summarise(
    pearson = cor(enrichment, GC)
  ) %>%
  arrange(pos) %>%
  pivot_wider(
    values_from = pearson,
    names_from = c(sp, sys)
  )

write_tsv(correlation, '../figures/rawData/GC_by_pos.tsv', na = 'NaN')


## promoter strength by GC content in dark without enhancer (Figure 2b,d)
promoter.strength.filtered <- promoter.strength %>%
  filter(! light & ! enhancer)

for (i in 1:nrow(sys.sp)) {
  system <- sys.sp[i, 'system', drop = TRUE]
  species <- sys.sp[i, 'species', drop = TRUE]

  df <- promoter.strength.filtered %>%
    filter(sys == system & sp == species) %>%
    mutate(
      sample.name = factor(
        GC.low,
        levels = sort(unique(GC.low)),
        labels = paste0('bin', seq_along(unique(GC.low))),
        ordered = TRUE
      ),
      LaTeX.label = sprintf('%.2f\\GCto %.2f', GC.low, GC.high)
    )


  # export data
  LaTeX.violinplot(
    data = df,
    samples_from = sample.name,
    values_from = enrichment,
    file = paste('../figures/rawData/enrichment_GC', species, system, sep = '_'),
    LaTeX.label,
    outliers = FALSE,
    p_values = TRUE,
    data_range = range(promoter.strength.filtered %>% filter(sys == system) %>% pull(enrichment))
  )
}


## TATA-box promoters (Figure 3a-c, 4c,d, and 5d and Supp Figure 3)

# histogram of TATA-box positions (Figure 3a)
TATA <- filter_motifs(CPEs, altname = 'TATA')[[1]]

TATA.scores <- as_tibble(scan_sequences(TATA, all.seqs, threshold = 1, nthreads = 0))

TATA.pos <- TATA.scores %>%
  mutate(
    score = (score - min.score) / (max.score - min.score)
  ) %>%
  filter(score >= thresh) %>%
  group_by(sequence) %>%
  filter(score == max(score)) %>%
  mutate(
    pos = 0.5 * (start + stop),
    weight = 1 / n()
  ) %>%
  ungroup()

TATA.pos.mod <- TATA.pos %>%
  group_by('gene' = sequence) %>%
  summarise(
    TATA.pos = if_else(any(between(pos, 107, 143)), 'prefPos', 'nonPref')
  )

data <- TATA.pos %>%
  inner_join(select(all.sequences, 'sequence' = gene, sp)) %>%
  group_by(sp, pos) %>%
  summarise(
    percentage = sum(weight) / sum(all.sequences$sp == unique(sp)) * 100
  ) %>%
  ungroup() %>%
  pivot_wider(
    names_from = sp,
    values_from = percentage
  ) %>%
  arrange(pos)

write_tsv(data, '../figures/rawData/TATA_position.tsv', na = 'NaN')


# promoter strength by TATA position in dark without enhancer (Figure 3b,c)
promoter.strength.filtered <- left_join(
  promoter.strength,
    TATA.pos.mod,
    by = 'gene'
  ) %>%
  filter(! enhancer & ! light)

for (i in 1:nrow(sys.sp)) {
  system <- sys.sp[i, 'system', drop = TRUE]
  species <- sys.sp[i, 'species', drop = TRUE]

  df <- promoter.strength.filtered %>%
    filter(sys == system, sp == species) %>%
    mutate(
      sample.name = factor(
        if_else(is.na(TATA.pos), 'noTATA', TATA.pos),
        levels = c('noTATA', 'nonPref', 'prefPos'),
        ordered = TRUE
      ),
      LaTeX.label = if_else(sample.name == 'noTATA', '$-$', '$+$')
    )

  LaTeX.violinplot(
    data = df,
    samples_from = sample.name,
    values_from = enrichment,
    file = paste('../figures/rawData/enrichment_TATA-pos', species, system, sep = '_'),
    LaTeX.label,
    outliers = FALSE,
    p_values = TRUE,
    data_range = range(promoter.strength.filtered %>% filter(sys == system) %>% pull(enrichment))
  )
}


# promoter strength by TATA and GC in dark without enhancer (Supp Figure 3)
promoter.strength.filtered <- left_join(
    promoter.strength,
    TATA.pos.mod,
    by = 'gene'
  ) %>%
  filter(! enhancer & ! light)

for (i in 1:nrow(sys.sp)) {
  system <- sys.sp[i, 'system', drop = TRUE]
  species <- sys.sp[i, 'species', drop = TRUE]

  df <- promoter.strength.filtered %>%
    filter(sys == system & sp == species) %>%
    mutate(
      sample = paste0(GC.low, if_else(! is.na(TATA.pos), '+', '-')),
      sample.name = factor(
        sample,
        levels = sort(unique(sample)),
        labels = paste0('bin', rep(seq_along(unique(GC.low)), each = 2), rep(c('-', '+'), times = length(unique(GC.low)))),
        ordered = TRUE
      ),
      LaTeX.label = sprintf('%.2f\\GCto %.2f', GC.low, GC.high)
    )

  LaTeX.violinplot(
    data = df,
    samples_from = sample.name,
    values_from = enrichment,
    file = paste('../figures/rawData/enrichment_TATA+GC', species, system, sep = '_'),
    LaTeX.label,
    outliers = FALSE,
    p_values = TRUE,
    half = TRUE,
    data_range = range(promoter.strength.filtered %>% filter(sys == system) %>% pull(enrichment))
  )
}


# enhancer responsiveness  by TATA position in dark (Figure 4c,d)
enhancer.effect.filtered <- left_join(
    enhancer.effect,
    TATA.pos.mod,
    by = 'gene'
  ) %>%
  filter(! light)

for (i in 1:nrow(sys.sp)) {
  system <- sys.sp[i, 'system', drop = TRUE]
  species <- sys.sp[i, 'species', drop = TRUE]

  df <- enhancer.effect.filtered %>%
    filter(sys == system, sp == species) %>%
    mutate(
      sample.name = factor(
        if_else(is.na(TATA.pos), 'noTATA', TATA.pos),
        levels = c('noTATA', 'nonPref', 'prefPos'),
        ordered = TRUE
      ),
      LaTeX.label = if_else(sample.name == 'noTATA', '$-$', '$+$')
    )

  LaTeX.violinplot(
    data = df,
    samples_from = sample.name,
    values_from = enhancer.effect,
    file = paste('../figures/rawData/enh-effect_TATA-pos', species, system, sep = '_'),
    LaTeX.label,
    outliers = FALSE,
    p_values = TRUE,
    data_range = range(enhancer.effect.filtered %>% filter(sys == system) %>% pull(enhancer.effect))
  )
}


# light-dependency by TATA and GC without enhancer (Figure 5d)
light.dependency.filtered <- left_join(
    light.dependency,
    TATA.pos.mod,
    by = 'gene'
  ) %>%
  filter(! enhancer)

for (species in all.species$short) {
  df <- light.dependency.filtered %>%
    filter(sp == species) %>%
    mutate(
      sample = paste0(GC.low, if_else(! is.na(TATA.pos), '+', '-')),
      sample.name = factor(
        sample,
        levels = sort(unique(sample)),
        labels = paste0('bin', rep(seq_along(unique(GC.low)), each = 2), rep(c('-', '+'), times = length(unique(GC.low)))),
        ordered = TRUE
      ),
      LaTeX.label = sprintf('%.2f\\GCto %.2f', GC.low, GC.high)
    )

  LaTeX.violinplot(
    data = df,
    samples_from = sample.name,
    values_from = light.dependency,
    file = paste('../figures/rawData/light-dep_TATA+GC', species, sep = '_'),
    LaTeX.label,
    outliers = FALSE,
    p_values = TRUE,
    half = TRUE,
    data_range = range(light.dependency.filtered$light.dependency)
  )
}


## promoter strength by core promoter elements in dark without enhancer (Supp Figure 4a-d, 6b,c, and 7)
promoter.strength.filtered <- promoter.strength %>%
  filter(! enhancer & ! light) %>%
  inner_join(motif.scores)

for (CPE in c('Inr', 'TCT', 'Ypatch')) {
  for (i in 1:nrow(sys.sp)) {
    system <- sys.sp[i, 'system', drop = TRUE]
    species <- sys.sp[i, 'species', drop = TRUE]
    
    df <- promoter.strength.filtered %>%
      filter(sys == system & sp == species) %>%
      mutate(
        sample = paste0(GC.low, if_else(!! sym(CPE) >= thresh, '+', '-')),
        sample.name = factor(
          sample,
          levels = sort(unique(sample)),
          labels = paste0('bin', rep(seq_along(unique(GC.low)), each = 2), rep(c('-', '+'), times = length(unique(GC.low)))),
          ordered = TRUE
        ),
        LaTeX.label = sprintf('%.2f\\GCto %.2f', GC.low, GC.high)
      )
    
    LaTeX.violinplot(
      data = df,
      samples_from = sample.name,
      values_from = enrichment,
      file = paste('../figures/rawData/enrichment', paste0(CPE, '+GC'), species, system, sep = '_'),
      LaTeX.label,
      outliers = FALSE,
      p_values = TRUE,
      half = TRUE,
      data_range = range(promoter.strength.filtered %>% filter(sys == system) %>% pull(enrichment))
    )
  }
}

for (CPE in c('BREu', 'BREd')) {
  for (system in all.systems$short) {
    
    df <- promoter.strength.filtered %>%
      filter(sys == system) %>%
      filter(TATA >= 0.7) %>%
      mutate(
        sample = paste0(sp, if_else(!! sym(CPE) >= thresh, '+', '-')),
        sample.name = factor(
          sample,
          levels = paste0(rep(c('At', 'Zm', 'Sb'), each = 2), c('-', '+')),
          ordered = TRUE
        ),
        LaTeX.label = short.common[as.character(sp)]
      )
    
    LaTeX.violinplot(
      data = df,
      samples_from = sample.name,
      values_from = enrichment,
      file = paste('../figures/rawData/enrichment', CPE, system, sep = '_'),
      LaTeX.label,
      outliers = FALSE,
      p_values = TRUE,
      half = TRUE,
      data_range = range(promoter.strength.filtered$enrichment)
    )
  }
}


## histogram of Y patch positions (Supp Figure 6a)
Ypatch <- filter_motifs(CPEs, altname = 'Ypatch')[[1]]

Ypatch.scores <- as_tibble(scan_sequences(Ypatch, all.seqs, RC = TRUE, threshold = 1, nthreads = 0))

Ypatch.pos <- Ypatch.scores %>%
  mutate(
    score = (score - min.score) / (max.score - min.score)
  ) %>%
  filter(score >= thresh) %>%
  group_by(sequence) %>%
  filter(score == max(score)) %>%
  mutate(
    pos = 0.5 * (start + stop),
    weight = 1 / n()
  ) %>%
  ungroup()

data <- Ypatch.pos %>%
  inner_join(select(all.sequences, 'sequence' = gene, sp)) %>%
  group_by(sp, pos) %>%
  summarise(
    percentage = sum(weight) / dim(all.sequences)[1] * 100
  ) %>%
  ungroup() %>%
  pivot_wider(
    names_from = sp,
    values_from = percentage
  ) %>%
  arrange(pos)

write_tsv(data, '../figures/rawData/Ypatch_position.tsv', na = 'NaN')


## promoter strength by TF binding site in dark without enhancer (Supp Figure 8)

# TCP
TF.scores <- motif.scores %>%
  select(gene, 'score1' = TF_15, 'score2' = TF_22)

promoter.strength.filtered <- promoter.strength %>%
  mutate(
    binding.site = gene %in% pull(filter(TF.scores, score1 >= thresh | score2 >= thresh), gene)
  ) %>%
  filter(! enhancer & ! light)

for (i in 1:nrow(sys.sp)) {
  system <- sys.sp[i, 'system', drop = TRUE]
  species <- sys.sp[i, 'species', drop = TRUE]

  df <- promoter.strength.filtered %>%
    filter(sys == system & sp == species) %>%
    mutate(
      sample = paste0(GC.low, if_else(binding.site, '+', '-')),
      sample.name = factor(
        sample,
        levels = sort(unique(sample)),
        labels = paste0('bin', rep(seq_along(unique(GC.low)), each = 2), rep(c('-', '+'), times = length(unique(GC.low)))),
        ordered = TRUE
      ),
      LaTeX.label = sprintf('%.2f\\GCto %.2f', GC.low, GC.high)
    )

  LaTeX.violinplot(
    data = df,
    samples_from = sample.name,
    values_from = enrichment,
    file = paste('../figures/rawData/enrichment_TCP', species, system, sep = '_'),
    LaTeX.label,
    outliers = FALSE,
    p_values = TRUE,
    half = TRUE,
    data_range = range(promoter.strength.filtered %>% filter(sys == system) %>% pull(enrichment))
  )
}


# HSF
TF.scores <- motif.scores %>%
  select(gene, 'score' = TF_16)

promoter.strength.filtered <- promoter.strength %>%
  mutate(
    binding.site = gene %in% pull(filter(TF.scores, score >= thresh), gene)
  ) %>%
  filter(! enhancer & ! light)

for (i in 1:nrow(sys.sp)) {
  system <- sys.sp[i, 'system', drop = TRUE]
  species <- sys.sp[i, 'species', drop = TRUE]

  df <- promoter.strength.filtered %>%
    filter(sys == system & sp == species) %>%
    mutate(
      sample = paste0(GC.low, if_else(binding.site, '+', '-')),
      sample.name = factor(
        sample,
        levels = sort(unique(sample)),
        labels = paste0('bin', rep(seq_along(unique(GC.low)), each = 2), rep(c('-', '+'), times = length(unique(GC.low)))),
        ordered = TRUE
      ),
      LaTeX.label = sprintf('%.2f\\GCto %.2f', GC.low, GC.high)
    )

  LaTeX.violinplot(
    data = df,
    samples_from = sample.name,
    values_from = enrichment,
    file = paste('../figures/rawData/enrichment_HSF', species, system, sep = '_'),
    LaTeX.label,
    outliers = FALSE,
    p_values = TRUE,
    half = TRUE,
    data_range = range(promoter.strength.filtered %>% filter(sys == system) %>% pull(enrichment))
  )
}


## TF position histograms and effects in dark without enhancer (Supp Figure 9)

# select sequences with a single strong TATA-box
TATA.genes <- motif.scores %>%
  filter(TATA >= thresh) %>%
  pull(gene)

TATA.seqs <- all.seqs[TATA.genes]

TATA.pos <- as_tibble(scan_sequences(filter_motifs(CPEs, name = 'TATA'), TATA.seqs, threshold = 1, nthreads = 0)) %>%
  group_by(sequence) %>%
  filter(score == max(score)) %>%
  filter(n() == 1) %>%
  ungroup() %>%
  select(sequence, 'TATA.start' = start, 'TATA.stop' = stop)

# TCP
TCP.pos <- as_tibble(scan_sequences(filter_motifs(TFs, altname = 'TCP'), all.seqs, RC = TRUE, threshold = 1, nthreads = 0)) %>%
  mutate(
    score = (score - min.score) / (max.score - min.score)
  ) %>%
  filter(score >= thresh) %>%
  group_by(sequence) %>%
  filter(score == max(score)) %>%
  mutate(
    pos = 0.5 * (start + stop)
  ) %>%
  ungroup()

# histogram
data <- TCP.pos %>%
  inner_join(select(all.sequences, 'sequence' = gene, sp)) %>%
  group_by(sp, pos) %>%
  summarise(
    percentage = n()
  ) %>%
  ungroup() %>%
  pivot_wider(
    names_from = sp,
    values_from = percentage
  ) %>%
  arrange(pos)

write_tsv(data, '../figures/rawData/TCP_position.tsv', na = 'NaN')

# promoter strength (rel. to TATA position)
TCP.pos.TATA <- TCP.pos %>%
  inner_join(TATA.pos, by = 'sequence') %>%
  mutate(
    rel.pos = case_when(
      stop < TATA.start ~ 'upstream',
      start > TATA.stop ~ 'downstream'
    )
  ) %>%
  filter(! is.na(rel.pos))

promoter.strength.filtered <- promoter.strength %>%
  filter(! enhancer & ! light & sys == 'leaf') %>%
  inner_join(select(TCP.pos.TATA, 'gene' = sequence, rel.pos), by = 'gene')

for (species in all.species$short) {
  df <- promoter.strength.filtered %>%
    filter(sp == species) %>%
    mutate(
      sample.name = ordered(
        rel.pos,
        levels = c('upstream', 'downstream')
      ),
      LaTeX.label = rel.pos
    )
  
  LaTeX.violinplot(
    data = df,
    samples_from = sample.name,
    values_from = enrichment,
    file = paste('../figures/rawData/enrichment_TCPpos', species, 'leaf', sep = '_'),
    LaTeX.label,
    outliers = FALSE,
    p_values = TRUE,
    half = FALSE,
    data_range = range(promoter.strength.filtered$enrichment)
  )
}


# HSF
HSF.pos <- as_tibble(scan_sequences(filter_motifs(TFs, name = 'TF-cluster_16'), all.seqs, RC = TRUE, threshold = 1, nthreads = 0)) %>%
  mutate(
    score = (score - min.score) / (max.score - min.score)
  ) %>%
  filter(score >= thresh) %>%
  group_by(sequence) %>%
  filter(score == max(score)) %>%
  mutate(
    pos = 0.5 * (start + stop)
  ) %>%
  ungroup()

# histogram
data <- HSF.pos %>%
  inner_join(select(all.sequences, 'sequence' = gene, sp)) %>%
  group_by(sp, pos) %>%
  summarise(
    percentage = n()
  ) %>%
  ungroup() %>%
  pivot_wider(
    names_from = sp,
    values_from = percentage
  ) %>%
  arrange(pos)

write_tsv(data, '../figures/rawData/HSF_position.tsv', na = 'NaN')

# promoter strength (rel. to TATA position)
HSF.pos.TATA <- HSF.pos %>%
  inner_join(TATA.pos, by = 'sequence') %>%
  mutate(
    rel.pos = case_when(
      stop < TATA.start ~ 'upstream',
      start > TATA.stop ~ 'downstream'
    )
  ) %>%
  filter(! is.na(rel.pos))

promoter.strength.filtered <- promoter.strength %>%
  filter(! enhancer & ! light & sys == 'proto') %>%
  inner_join(select(HSF.pos.TATA, 'gene' = sequence, rel.pos), by = 'gene')

for (species in all.species$short) {
  df <- promoter.strength.filtered %>%
    filter(sp == species) %>%
    mutate(
      sample.name = ordered(
        rel.pos,
        levels = c('upstream', 'downstream')
      ),
      LaTeX.label = rel.pos
    )
  
  LaTeX.violinplot(
    data = df,
    samples_from = sample.name,
    values_from = enrichment,
    file = paste('../figures/rawData/enrichment_HSFpos', species, 'proto', sep = '_'),
    LaTeX.label,
    outliers = FALSE,
    p_values = TRUE,
    half = FALSE,
    data_range = range(promoter.strength.filtered$enrichment)
  )
}


## enhancer responsiveness by promoter type in dark (Supp Figure 10)
enhancer.effect.filtered <- enhancer.effect %>%
  filter(! light)

for (i in 1:nrow(sys.sp)) {
  system <- sys.sp[i, 'system', drop = TRUE]
  species <- sys.sp[i, 'species', drop = TRUE]
  
  df <- enhancer.effect.filtered %>%
    filter(sys == system & sp == species) %>%
    mutate(
      sample.name = factor(type, levels = c('miRNA', 'protein_coding'), ordered = TRUE),
      LaTeX.label = if_else(sample.name == 'miRNA', 'miRNA', 'protein\\\\coding')
    )
  
  LaTeX.violinplot(
    data = df,
    samples_from = sample.name,
    values_from = enhancer.effect,
    file = paste('../figures/rawData/enh-effect_type', species, system, sep = '_'),
    LaTeX.label,
    outliers = FALSE,
    p_values = TRUE,
    data_range = range(enhancer.effect.filtered %>% filter(sys == system) %>% pull(enhancer.effect))
  )
}


## enhancer responsiveness by GC content in dark (Figure 4e,f)
enhancer.effect.filtered <- enhancer.effect %>%
  filter(! light)

for (i in 1:nrow(sys.sp)) {
  system <- sys.sp[i, 'system', drop = TRUE]
  species <- sys.sp[i, 'species', drop = TRUE]
  
  df <- enhancer.effect.filtered %>%
    filter(sys == system & sp == species) %>%
    mutate(
      sample.name = factor(
        GC.low,
        levels = sort(unique(GC.low)),
        labels = paste0('bin', seq_along(unique(GC.low))),
        ordered = TRUE
      ),
      LaTeX.label = sprintf('%.2f\\GCto %.2f', GC.low, GC.high)
    )

  LaTeX.violinplot(
    data = df,
    samples_from = sample.name,
    values_from = enhancer.effect,
    file = paste('../figures/rawData/enh-effect_GC', species, system, sep = '_'),
    LaTeX.label,
    outliers = FALSE,
    p_values = TRUE,
    data_range = range(enhancer.effect.filtered %>% filter(sys == system) %>% pull(enhancer.effect))
  )
}


## enhancer responsiveness by tissue specificity in dark (Figure 4a,b)

# load expression data and calculate specificity
for (species in all.species$short) {
  if (! file.exists(paste0('expression_data/', short.common[species], '_FPKM.tsv'))) {
    if (species == 'At') {
      experiment <- 'E-MTAB-7978'
    } else if (species == 'Zm') {
      experiment <- 'E-GEOD-50191'
    } else if (species == 'Sb') {
      experiment <- 'E-MTAB-5956'
    } else {
      stop(paste('Unknown species:', species))
    }
    download.file(
      paste0('https://www.ebi.ac.uk/gxa/experiments-content/', experiment, '/resources/ExperimentDownloadSupplier.RnaSeqBaseline/fpkms.tsv'),
      paste0('expression_data/', short.common[species], '_FPKM.tsv')
    )
  }
  
  tissue.specificity <- read_tsv(paste0('expression_data/', short.common[species], '_FPKM.tsv'), comment = '#') %>%
    mutate_if(is.numeric, replace_na, 0.1) %>%
    filter_if(is.numeric, any_vars(. >= 1)) %>%
    mutate_if(is.numeric, log2) %>%
    select(-`Gene Name`) %>%
    rename('gene' = `Gene ID`) %>%
    pivot_longer(
      -gene,
      names_to = 'sample',
      values_to = 'expression'
    ) %>%
    group_by(gene) %>%
    summarise(
      Tau = sum(1 - (expression / max(expression))) / (n() - 1)
    ) %>%
    filter(! is.na(Tau)) %>%
    mutate(
      specificity = cut_number(Tau, n = 3, ordered = TRUE, labels = c('low', 'medium', 'high'))
    )
  
  assign(paste0('tissue.specificity.', species), tissue.specificity)
}


tissue.specificity <- bind_rows(
  tissue.specificity.At,
  tissue.specificity.Zm,
  tissue.specificity.Sb,
)

# merge tissue specificity with enhancer responsiveness and export data
enhancer.effect.filtered <- enhancer.effect %>%
  filter(! light) %>%
  inner_join(tissue.specificity, by = 'gene')


for (i in 1:nrow(sys.sp)) {
  system <- sys.sp[i, 'system', drop = TRUE]
  species <- sys.sp[i, 'species', drop = TRUE]
  
  df <- enhancer.effect.filtered %>%
    filter(sys == system & sp == species) %>%
    mutate(
      sample.name = specificity,
      LaTeX.label = sample.name
    )
  
  LaTeX.violinplot(
    data = df,
    samples_from = sample.name,
    values_from = enhancer.effect,
    file = paste('../figures/rawData/enh-effect_specificity', species, system, sep = '_'),
    LaTeX.label,
    outliers = FALSE,
    p_values = TRUE,
    data_range = range(enhancer.effect.filtered %>% filter(sys == system) %>% pull(enhancer.effect))
  )
}


## enhancer responsiveness by TF binding site in dark (Supp Figure 11)
enhancer.effect.filtered <- enhancer.effect %>%
  filter(! light)

# TCPs
TF.scores <- motif.scores %>%
  select(gene, 'score1' = TF_15, 'score2' = TF_22)

enhancer.effect.mod <- enhancer.effect.filtered %>%
  mutate(
    binding.site = gene %in% pull(filter(TF.scores, score1 >= thresh | score2 >= thresh), gene)
  ) %>%
  filter(! light)

for (i in 1:nrow(sys.sp)) {
  system <- sys.sp[i, 'system', drop = TRUE]
  species <- sys.sp[i, 'species', drop = TRUE]
  
  df <- enhancer.effect.mod %>%
    filter(sys == system & sp == species) %>%
    mutate(
      sample = paste0(GC.low, if_else(binding.site, '+', '-')),
      sample.name = factor(
        sample,
        levels = sort(unique(sample)),
        labels = paste0('bin', rep(seq_along(unique(GC.low)), each = 2), rep(c('-', '+'), times = length(unique(GC.low)))),
        ordered = TRUE
      ),
      LaTeX.label = sprintf('%.2f\\GCto %.2f', GC.low, GC.high)
    )
  
  LaTeX.violinplot(
    data = df,
    samples_from = sample.name,
    values_from = enhancer.effect,
    file = paste('../figures/rawData/enh-effect_TCP', species, system, sep = '_'),
    LaTeX.label,
    outliers = FALSE,
    p_values = TRUE,
    half = TRUE,
    data_range = range(enhancer.effect.mod %>% filter(sys == system) %>% pull(enhancer.effect))
  )
}

# B3
TF.scores <- motif.scores %>%
  select(gene, 'score' = TF_35)

enhancer.effect.mod <- enhancer.effect.filtered %>%
  mutate(
    binding.site = gene %in% pull(filter(TF.scores, score >= thresh), gene)
  ) %>%
  filter(! light)

for (i in 1:nrow(sys.sp)) {
  system <- sys.sp[i, 'system', drop = TRUE]
  species <- sys.sp[i, 'species', drop = TRUE]
  
  df <- enhancer.effect.mod %>%
    filter(sys == system & sp == species) %>%
    mutate(
      sample = paste0(GC.low, if_else(binding.site, '+', '-')),
      sample.name = factor(
        sample,
        levels = sort(unique(sample)),
        labels = paste0('bin', rep(seq_along(unique(GC.low)), each = 2), rep(c('-', '+'), times = length(unique(GC.low)))),
        ordered = TRUE
      ),
      LaTeX.label = sprintf('%.2f\\GCto %.2f', GC.low, GC.high)
    )
  
  LaTeX.violinplot(
    data = df,
    samples_from = sample.name,
    values_from = enhancer.effect,
    file = paste('../figures/rawData/enh-effect_B3', species, system, sep = '_'),
    LaTeX.label,
    outliers = FALSE,
    p_values = TRUE,
    half = TRUE,
    data_range = range(enhancer.effect.mod %>% filter(sys == system) %>% pull(enhancer.effect))
  )
}

# WRKY
TF.scores <- motif.scores %>%
  select(gene, 'score' = TF_2)

enhancer.effect.mod <- enhancer.effect.filtered %>%
  mutate(
    binding.site = gene %in% pull(filter(TF.scores, score >= thresh), gene)
  ) %>%
  filter(! light)

for (i in 1:nrow(sys.sp)) {
  system <- sys.sp[i, 'system', drop = TRUE]
  species <- sys.sp[i, 'species', drop = TRUE]
  
  df <- enhancer.effect.mod %>%
    filter(sys == system & sp == species) %>%
    mutate(
      sample = paste0(GC.low, if_else(binding.site, '+', '-')),
      sample.name = factor(
        sample,
        levels = sort(unique(sample)),
        labels = paste0('bin', rep(seq_along(unique(GC.low)), each = 2), rep(c('-', '+'), times = length(unique(GC.low)))),
        ordered = TRUE
      ),
      LaTeX.label = sprintf('%.2f\\GCto %.2f', GC.low, GC.high)
    )
  
  LaTeX.violinplot(
    data = df,
    samples_from = sample.name,
    values_from = enhancer.effect,
    file = paste('../figures/rawData/enh-effect_WRKY', species, system, sep = '_'),
    LaTeX.label,
    outliers = FALSE,
    p_values = TRUE,
    half = TRUE,
    data_range = range(enhancer.effect.mod %>% filter(sys == system) %>% pull(enhancer.effect))
  )
}


## light-dependency without and with enhancer (Figure 5b)
light.dependency.mod <- light.dependency %>%
  mutate(
    enhancer = enhancer.name[as.character(enhancer)]
  )

for (species in all.species$short) {
  df <- light.dependency.mod %>%
    filter(sp == species) %>%
    mutate(
      sample.name = factor(
        enhancer,
        levels = c('noENH', 'withENH'),
        ordered = TRUE
      ),
      LaTeX.label = if_else(enhancer == 'withENH', '$+$', '$-$')
    )
  
  LaTeX.violinplot(
    data = df,
    samples_from = sample.name,
    values_from = light.dependency,
    file = paste('../figures/rawData/light-dep_enh', species, sep = '_'),
    LaTeX.label,
    outliers = FALSE,
    p_values = FALSE,
    half = FALSE
  )
}


## light-dependency by TF binding site without enhancer (Figure 5e,f)

# WRKY
TF.scores <- motif.scores %>%
  select(gene, 'score' = TF_2)

light.dependency.filtered <- light.dependency %>%
  mutate(
    binding.site = gene %in% pull(filter(TF.scores, score >= thresh), gene)
  ) %>%
  filter(! enhancer)

for (species in all.species$short) {
  df <- light.dependency.filtered %>%
    filter(sp == species) %>%
    mutate(
      sample = paste0(GC.low, if_else(binding.site, '+', '-')),
      sample.name = factor(
        sample,
        levels = sort(unique(sample)),
        labels = paste0('bin', rep(seq_along(unique(GC.low)), each = 2), rep(c('-', '+'), times = length(unique(GC.low)))),
        ordered = TRUE
      ),
      LaTeX.label = sprintf('%.2f\\GCto %.2f', GC.low, GC.high)
    )
  
  LaTeX.violinplot(
    data = df,
    samples_from = sample.name,
    values_from = light.dependency,
    file = paste('../figures/rawData/light-dep_WRKY', species, sep = '_'),
    LaTeX.label,
    outliers = FALSE,
    p_values = TRUE,
    half = TRUE,
    data_range = range(light.dependency.filtered$light.dependency)
  )
}

# TCPs
TF.scores <- motif.scores %>%
  select(gene, 'score1' = TF_15, 'score2' = TF_22)

light.dependency.filtered <- light.dependency %>%
  mutate(
    binding.site = gene %in% pull(filter(TF.scores, score1 >= thresh | score2 >= thresh), gene)
  ) %>%
  filter(! enhancer)

for (species in all.species$short) {
  df <- light.dependency.filtered %>%
    filter(sp == species) %>%
    mutate(
      sample = paste0(GC.low, if_else(binding.site, '+', '-')),
      sample.name = factor(
        sample,
        levels = sort(unique(sample)),
        labels = paste0('bin', rep(seq_along(unique(GC.low)), each = 2), rep(c('-', '+'), times = length(unique(GC.low)))),
        ordered = TRUE
      ),
      LaTeX.label = sprintf('%.2f\\GCto %.2f', GC.low, GC.high)
    )
  
  LaTeX.violinplot(
    data = df,
    samples_from = sample.name,
    values_from = light.dependency,
    file = paste('../figures/rawData/light-dep_TCP', species, sep = '_'),
    LaTeX.label,
    outliers = FALSE,
    p_values = TRUE,
    half = TRUE,
    data_range = range(light.dependency.filtered$light.dependency)
  )
}


## correlation with validation libraries (Supp Figure 2)
promoters <- promoter.strength %>%
  select(sys, enhancer, light, gene, enrichment)

validation <- promoter.strength.val %>%
  filter(variant == 'WT' & motif != 'synthetic') %>%
  select(set, sys, enhancer, light, gene, enrichment)

both <- inner_join(
    promoters,
    validation,
    by = c('sys', 'enhancer', 'light', 'gene'),
    suffix = c('.pro', '.val')
  ) %>%
  mutate(
    enhancer = enhancer.name[as.character(enhancer)],
    light = light.name[as.character(light)]
  )

helper <- tibble(
  set = rep(c('PROval', 'PROevo'), each = 3),
  system = rep(c('leaf', 'proto', 'leaf'), times = 2),
  condition = rep(c('dark', 'dark', 'light'), times = 2)
)

for (i in 1:nrow(helper)) {
  val.set <- helper[i, 'set', drop = TRUE]
  system <- helper[i, 'system', drop = TRUE]
  condition <- helper[i, 'condition', drop = TRUE]
  
  df <- both %>%
    filter(set == val.set & sys == system & light == condition) %>%
    mutate(
      sample.name = enhancer
    ) %>%
    select(sample.name, enrichment.pro, enrichment.val) %>%
    slice_sample(prop = 1)
  
  correlation <- df %>%
    group_by(sample.name) %>%
    summarise(
      spearman = sprintf('$\\rho=%.2f$', cor(enrichment.pro, enrichment.val, method = 'spearman', use = 'complete.obs')),
      rsquare = sprintf('$R^2=%.2f$', cor(enrichment.pro, enrichment.val, method = 'pearson', use = 'complete.obs')^2)
    )
  
  write_tsv(df, paste0('../figures/rawData/enrichment_correlation_', val.set, '_', if_else(system == 'leaf', condition, system), '.tsv'), na = 'NaN')
  write_tsv(correlation, paste0('../figures/rawData/enrichment_correlation_', val.set, '_', if_else(system == 'leaf', condition, system), '_stats.tsv'), na = 'NaN')
}


## validation of TATA mutants without enhancer in dark (Figure 3e,g)
TATA.val <- promoter.strength.val %>%
  filter(set == 'PROval' & grepl('TATA', motif, fixed = TRUE)) %>%
  filter(! enhancer & ! light)

# set WT to 0
TATA.norm <- TATA.val %>%
  group_by(sys, enhancer, light, gene, motif) %>%
  mutate(
    enrichment = enrichment - median(enrichment[variant == 'WT'])
  ) %>%
  ungroup() %>%
  filter(variant != 'WT' & ! is.na(enrichment))

# mutations in TATA-box
for (system in all.systems$short) {
  df <- TATA.norm %>%
    filter(sys == system) %>%
    mutate(
      sample.name = ordered(variant, levels = c('mutA', 'mutB', 'mutAB'))
    ) %>%
    filter(! is.na(sample.name)) %>%
    group_by(sample.name) %>%
    mutate(
      p.value = wilcox.test(enrichment)$p.value
    ) %>%
    ungroup()
  
  LaTeX.boxplot(
    data = df,
    samples_from = sample.name,
    values_from = enrichment,
    file = paste('../figures/rawData/enrichment_mutateTATA', system, sep = '_'),
    p.value,
    outliers = TRUE,
    p_values = FALSE
  )
}

# insert TATA-box
for (system in all.systems$short) {
  df <- TATA.norm %>%
    filter(sys == system) %>%
    mutate(
      sample.name = ordered(variant, levels = c('+TATA', '+mutTATA'))
    ) %>%
    filter(! is.na(sample.name)) %>%
    group_by(sample.name) %>%
    mutate(
      p.value = wilcox.test(enrichment)$p.value
    ) %>%
    ungroup()
  
  LaTeX.boxplot(
    data = df,
    samples_from = sample.name,
    values_from = enrichment,
    file = paste('../figures/rawData/enrichment_insertTATA', system, sep = '_'),
    p.value,
    outliers = TRUE,
    p_values = FALSE
  )
}


## validation of BRE mutants without enhancer in dark (Supp Figure 4h)
BRE.val <- promoter.strength.val %>%
  filter(grepl('BRE', motif, fixed = TRUE)) %>%
  filter(! enhancer & ! light)

# set WT to 0
BRE.norm <- BRE.val %>%
  group_by(sys, enhancer, light, gene, motif) %>%
  mutate(
    enrichment = enrichment - median(enrichment[variant == 'WT'])
  ) %>%
  ungroup() %>%
  filter(variant != 'WT' & ! is.na(enrichment))

for (system in all.systems$short) {
  df <- BRE.norm %>%
    filter(sys == system) %>%
    mutate(
      sample.name = ordered(
        paste(motif, variant, sep = '.'),
        levels = c('BREu.mut', 'BREd.mut', 'noBRE.+BREu', 'noBRE.+BREd')
      )
    ) %>%
    filter(! is.na(sample.name)) %>%
    group_by(sample.name) %>%
    mutate(
      p.value = wilcox.test(enrichment)$p.value
    ) %>%
    ungroup() %>%
    mutate(
      LaTeX.label = sub('noBRE', 'no BRE', sample.name, fixed = TRUE),
      LaTeX.label = sub('+BRE', '+ BRE', LaTeX.label, fixed = TRUE),
      LaTeX.label = str_replace_all(LaTeX.label, 'BRE([ud])', 'BRE\\\\\\textsuperscript{\\1}'),
      LaTeX.label = sub('.', '\\\\[-.25\\baselineskip](', LaTeX.label, fixed = TRUE),
      LaTeX.label = paste0(LaTeX.label, ')')
    )
  
  LaTeX.boxplot(
    data = df,
    samples_from = sample.name,
    values_from = enrichment,
    file = paste('../figures/rawData/enrichment_mutateBRE', system, sep = '_'),
    p.value,
    LaTeX.label,
    outliers = TRUE,
    p_values = FALSE
  )
}


## synthetic promoters without enhancer in dark (Figure 6b,c)
synthetic <- promoter.strength.val %>%
  filter(motif == 'synthetic') %>%
  filter(! enhancer & ! light)

for (i in 1:nrow(sys.sp)) {
  system <- sys.sp[i, 'system', drop = TRUE]
  species <- sys.sp[i, 'species', drop = TRUE]
  
  if (species == 'Sb') {
    next
  }
  
  df <- synthetic %>%
    filter(sys == system & sp == species) %>%
    mutate(
      sample.name = ordered(
        variant,
        levels = c('WT', '+Inr', '+Ypatch', '+TATA', '+Inr+Ypatch', '+TATA+Inr', '+TATA+Ypatch', '+TATA+Inr+Ypatch')
      ),
      TATA = if_else(grepl('TATA', variant, fixed = TRUE), '$+$', '$-$'),
      Inr = if_else(grepl('Inr', variant, fixed = TRUE), '$+$', '$-$'),
      Ypatch = if_else(grepl('Ypatch', variant, fixed = TRUE), '$+$', '$-$')
    ) %>%
    filter(! is.na(sample.name))
  
  LaTeX.boxplot(
    data = df,
    samples_from = sample.name,
    values_from = enrichment,
    file = paste('../figures/rawData/enrichment_syntheticPRO', species, system, sep = '_'),
    TATA,
    Inr,
    Ypatch,
    outliers = TRUE,
    p_values = FALSE
  )
}


## TF combinations in synthetic promoters without enhancer in dark (Figure 6e,f and Supp Figure 13)
TFcombo <- promoter.strength.val %>%
  filter(motif == 'TFcombo' | (grepl('syn', gene, fixed = TRUE) & variant == '+TATA')) %>%
  filter(! enhancer & ! light)

# set WT (synthetic promoter +TATA) to 0
TFcombo.norm <- TFcombo %>%
  group_by(sys, enhancer, light, gene) %>%
  mutate(
    enrichment = enrichment - median(enrichment[variant == '+TATA'])
  ) %>%
  ungroup() %>%
  filter(! is.na(enrichment) & variant != '+TATA')

# count number of TFs
TFcombo.norm <- TFcombo.norm %>%
  mutate(
    TFs = 3 - str_count(variant, 'none')
  )

# combine TCPs (individual plots are very similar to each other)
TFcombo.norm <- TFcombo.norm %>%
  mutate(
    variant = gsub('TCP[12]', 'TCP', variant)
  )

# promoter strength by number of TFs
for (system in all.systems$short) {
  df <- TFcombo.norm %>%
    filter(sys == system) %>%
    mutate(
      sample.name = ordered(TFs)
    ) %>%
    filter(! is.na(sample.name)) %>%
    group_by(sample.name) %>%
    mutate(
      p.value = wilcox.test(enrichment)$p.value
    ) %>%
    ungroup()
  
  LaTeX.boxplot(
    data = df,
    samples_from = sample.name,
    values_from = enrichment,
    file = paste('../figures/rawData/enrichment_TFcombo_number', system, sep = '_'),
    p.value,
    outliers = TRUE,
    p_values = FALSE
  )
}

# single TFs
singleTF <- TFcombo.norm %>%
  filter(TFs == 1) %>%
  mutate(
    variant = gsub('-?none-?', '', variant)
  )

for (system in all.systems$short) {
  df <- singleTF %>%
    filter(sys == system) %>%
    mutate(
      sample.name = ordered(
        variant,
        levels = c('TCP', 'NAC', 'HSF')
      )
    ) %>%
    filter(! is.na(sample.name)) %>%
    group_by(sample.name) %>%
    mutate(
      p.value = wilcox.test(enrichment)$p.value
    ) %>%
    ungroup()

  LaTeX.boxplot(
    data = df,
    samples_from = sample.name,
    values_from = enrichment,
    file = paste('../figures/rawData/enrichment_TFcombo_single', system, sep = '_'),
    p.value,
    outliers = TRUE,
    p_values = FALSE
  )
}

# double TF combinations
doubleTF <- TFcombo.norm %>%
  filter(TFs == 2) %>%
  separate(
    variant,
    into = c('TF1', 'TF2', 'TF3'),
    sep = '-'
  ) %>%
  mutate(
    TF1 = if_else(TF1 == 'none', TF3, TF1),
    TF2 = if_else(TF2 == 'none', TF3, TF2)
  ) %>%
  rowwise() %>%
  mutate(
    combo = paste0(sort(c(TF1, TF2)), collapse = '+'),
    TCP = str_count(combo, 'TCP'),
    NAC = str_count(combo, 'NAC'),
    HSF = str_count(combo, 'HSF')
  ) %>%
  ungroup()

for (system in all.systems$short) {
  df <- doubleTF %>%
    filter(sys == system) %>%
    mutate(
      sample.name = ordered(
        combo,
        levels = sort(unique(combo))
      )
    ) %>%
    group_by(sample.name) %>%
    mutate(
      p.value = wilcox.test(enrichment)$p.value
    ) %>%
    ungroup()
  
  LaTeX.boxplot(
    data = df,
    samples_from = sample.name,
    values_from = enrichment,
    file = paste('../figures/rawData/enrichment_TFcombo_double', system, sep = '_'),
    TCP,
    NAC,
    HSF,
    p.value,
    outliers = TRUE,
    p_values = FALSE
  )
}

# triple TF combinations
tripleTF <- TFcombo.norm %>%
  filter(TFs == 3) %>%
  separate(
    variant,
    into = c('TF1', 'TF2', 'TF3'),
    sep = '-'
  ) %>%
  rowwise() %>%
  mutate(
    combo = paste0(sort(c(TF1, TF2, TF3)), collapse = '+'),
    TCP = str_count(combo, 'TCP'),
    NAC = str_count(combo, 'NAC'),
    HSF = str_count(combo, 'HSF')
  ) %>%
  ungroup()

for (system in all.systems$short) {
  df <- tripleTF %>%
    filter(sys == system) %>%
    mutate(
      sample.name = ordered(
        combo,
        levels = sort(unique(combo))
      )
    ) %>%
    group_by(sample.name) %>%
    mutate(
      p.value = wilcox.test(enrichment)$p.value
    ) %>%
    ungroup()
  
  LaTeX.boxplot(
    data = df,
    samples_from = sample.name,
    values_from = enrichment,
    file = paste('../figures/rawData/enrichment_TFcombo_triple', system, sep = '_'),
    TCP,
    NAC,
    HSF,
    p.value,
    outliers = TRUE,
    p_values = FALSE
  )
}

# export PWMs
PWM.to.LaTeX(filter_motifs(TFs, name = 'TF-cluster_1')[[1]], '../figures/rawData/motif_NAC.tsv')
PWM.to.LaTeX(filter_motifs(TFs, name = 'TF-cluster_15')[[1]], '../figures/rawData/motif_TCP15.tsv')
PWM.to.LaTeX(filter_motifs(TFs, name = 'TF-cluster_22')[[1]], '../figures/rawData/motif_TCP22.tsv')
PWM.to.LaTeX(filter_motifs(TFs, name = 'TF-cluster_16')[[1]], '../figures/rawData/motif_HSF.tsv')


## promoter evolution in dark (Figure 7c-f)

# with enhancer
evolution.val <- promoter.strength.val %>%
  filter(enhancer & ! light) %>%
  filter(motif == 'evolution' | motif == 'synthetic' | grepl('native', motif, fixed = 'TRUE')) %>%
  mutate(
    gene = if_else(motif == 'synthetic', paste0(gene, '(', variant, ')'), gene),
    variant = if_else(motif == 'evolution', variant, 'leaf/proto/both-0')
  ) %>%
  separate(
    variant,
    into = c('opt_for', 'round'),
    sep = '-'
  ) %>%
  mutate(
    round = as.numeric(round)
  ) %>%
  separate_rows(
    opt_for,
    sep = '/'
  )

sys.opt_for <- expand_grid(optimization = c('leaf', 'proto', 'both'), system = all.systems$short)

for (i in 1:nrow(sys.opt_for)) {
  system <- sys.opt_for[i, 'system', drop = TRUE]
  optimization <- sys.opt_for[i, 'optimization', drop = TRUE]
  
  df <- evolution.val %>%
    filter(sys == system & opt_for == optimization) %>%
    mutate(
      sample.name = ordered(round)
    )
  
  LaTeX.violinplot(
    data = df,
    samples_from = sample.name,
    values_from = enrichment,
    file = paste('../figures/rawData/enrichment_evolution', optimization, system, sep = '_'),
    outliers = FALSE,
    p_values = TRUE,
    data_range = range(evolution.val %>% filter(sys == system) %>% pull(enrichment))
  )
}

# without enhancer
evolution.val <- promoter.strength.val %>%
  filter(! enhancer & ! light) %>%
  filter(motif == 'evolution' | motif == 'synthetic' | grepl('native', motif, fixed = 'TRUE')) %>%
  mutate(
    gene = if_else(motif == 'synthetic', paste0(gene, '(', variant, ')'), gene),
    variant = if_else(motif == 'evolution', variant, 'leaf/proto/both-0')
  ) %>%
  separate(
    variant,
    into = c('opt_for', 'round'),
    sep = '-'
  ) %>%
  mutate(
    round = as.numeric(round)
  ) %>%
  separate_rows(
    opt_for,
    sep = '/'
  )

for (i in 1:nrow(sys.opt_for)) {
  system <- sys.opt_for[i, 'system', drop = TRUE]
  optimization <- sys.opt_for[i, 'optimization', drop = TRUE]
  
  df <- evolution.val %>%
    filter(sys == system & opt_for == optimization) %>%
    mutate(
      sample.name = ordered(round)
    )
  
  LaTeX.violinplot(
    data = df,
    samples_from = sample.name,
    values_from = enrichment,
    file = paste('../figures/rawData/enrichment_evolution_noEnh', optimization, system, sep = '_'),
    outliers = FALSE,
    p_values = TRUE,
    data_range = range(evolution.val %>% filter(sys == system) %>% pull(enrichment))
  )
}


## TF position-dependency without enhancer in dark (Figure 6g,h)
TFscan <- promoter.strength.val %>%
  filter(motif == 'TFscan' | (grepl('syn', gene, fixed = TRUE) & variant == '+TATA') | str_count(variant, 'none') == 2) %>%
  filter(! enhancer & ! light)

# set WT (synthetic promoter +TATA) to 0
TFscan.norm <- TFscan %>%
  group_by(sys, enhancer, light, gene) %>%
  mutate(
    enrichment = enrichment - median(enrichment[variant == '+TATA'])
  ) %>%
  ungroup() %>%
  filter(! is.na(enrichment) & variant != '+TATA')

# get TF and position
TFscan.norm <- TFscan.norm %>%
  mutate(
    variant = gsub('-none-none', '_35', variant, fixed = TRUE),
    variant = case_when(
      grepl('none-.*-none', variant) ~ paste(gsub('-?none-?', '', variant), 65, sep = '_'),
      grepl('none-none-', variant, fixed = TRUE) ~ paste(gsub('-?none-?', '', variant), 95, sep = '_'),
      TRUE ~ variant
    )
  ) %>%
  separate(
    variant,
    into = c('motif', 'pos'),
    sep = '_'
  ) %>%
  mutate(
    pos = as.numeric(pos)
  )

# combine TCPs (individual plots are very similar to each other)
TFscan.norm <- TFscan.norm %>%
  mutate(
    motif = if_else(grepl('TCP', motif, fixed = TRUE), 'TCP', motif)
  )

# export data: TCP for leaf system, HSF for protoplasts
for (system in all.systems$short) {
  df <- TFscan.norm %>%
    filter(sys == system & motif == if_else(system == 'leaf', 'TCP', 'HSF')) %>%
    mutate(
      sample.name = ordered(pos)
    ) %>%
    filter(! is.na(sample.name)) %>%
    group_by(sample.name) %>%
    mutate(
      p.value = wilcox.test(enrichment)$p.value
    ) %>%
    ungroup()
  
  LaTeX.boxplot(
    data = df,
    samples_from = sample.name,
    values_from = enrichment,
    file = paste('../figures/rawData/enrichment_TFscan', if_else(system == 'leaf', 'TCP', 'HSF'), system, sep = '_'),
    p.value,
    outliers = TRUE,
    p_values = FALSE
  )
}


### Supplementary tables

## promoter annotation (Supplementary Table 1)

# load subassemblies
subassemblies <- bind_rows(
    'At' = read_tsv('subassembly/subassembly_At.tsv'),
    'Zm' = read_tsv('subassembly/subassembly_Zm.tsv'),
    'Sb' = read_tsv('subassembly/subassembly_Sb.tsv'),
    .id = 'sp'
  ) %>%
  filter(variant == 'WT' & FL) %>%
  group_by(sp, gene, type, GC, UTR, mutations) %>%
  summarise(
    barcodes = n()
  ) %>%
  ungroup() %>%
  mutate(
    name = gene
  ) %>%
  separate_rows(
    gene,
    sep = '[/;]'
  )

# load coordinates
bed.cols <- c('chr', 'start', 'end', 'gene', 'score', 'strand')

promoter.coords <- bind_rows(
    read_tsv('../data/promoter_annotation/Arabidopsis_miRNA_promoters.bed', col_names = bed.cols),
    read_tsv('../data/promoter_annotation/Arabidopsis_protein_coding_promoters.bed', col_names = bed.cols),
    read_tsv('../data/promoter_annotation/Maize_miRNA_promoters.bed', col_names = bed.cols),
    read_tsv('../data/promoter_annotation/Maize_protein_coding_promoters.bed', col_names = bed.cols),
    read_tsv('../data/promoter_annotation/Sorghum_miRNA_promoters.bed', col_names = bed.cols),
    read_tsv('../data/promoter_annotation/Sorghum_protein_coding_promoters.bed', col_names = bed.cols)
  ) %>%
  mutate(
    # correct start coordinate (start is 0-based and stop is 1-based in .bed format)
    start = start + 1
  ) %>%
  select(-score)

# combine subassembly with coordinates and sequence
subassemblies.mod <- subassemblies %>%
  left_join(promoter.coords, by = 'gene') %>%
  group_by(sp, type, strand, GC, UTR, mutations, barcodes, name) %>%
  summarise(
    across(c(chr, start, end), paste0, collapse = ';')
  ) %>%
  ungroup() %>%
  group_by(sp, type, GC, UTR, mutations, barcodes, 'gene' = name) %>%
  summarise(
    across(c(chr, start, end, strand),  paste0, collapse = '/'),
    UTR = all(UTR)
  ) %>%
  ungroup() %>%
  left_join(
    select(all.sequences, gene, sequence),
    by = 'gene'
  ) %>%
  select(sp, gene, barcodes, type, chr, start, end, strand, GC, UTR, mutations, sequence)

# create an excel file with the promoter annotation
wb <- createWorkbook()

modifyBaseFont(wb, fontSize = 10, fontName = 'Arial')

bold <- createStyle(textDecoration = 'bold')
wrap <- createStyle(wrapText = TRUE)
two.digit <- createStyle(numFmt = '0.00')
seq.font <- createStyle(fontName = 'Courier New')

for (species in all.species$common) {
  addWorksheet(wb, sheetName = species)
}

for (species in all.species$short) {
  writeData(
    wb,
    sheet = short.common[species],
    'Supplementary Table 1 | Promoters analyzed in this study.',
    startCol = 1,
    startRow = 1
  )
  
  mergeCells(wb, sheet = short.common[species], rows = 1, cols = 1:11)
  addStyle(wb, sheet = short.common[species], style = bold, rows = 1, cols = 1:11)
  
  writeData(
    wb,
    sheet = short.common[species],
    "Promoters for the genes listed in this table were present in our library and linked to the indicated number of barcodes. The type of the gene and whether it contains an annotated 5' UTR is indicated. The chromosomal coordinates, orientation, and GC content of the promoters is listed. To remove recognition sites for the restriction enzymes used in library creation, the indicated mutations were introduced and the resulting sequences were array-synthesized. The promoter sequences for some genes were identical and the corresponding annotations were merged (separated by ';' for genes in the same orientation, or '/' for genes in opposite orientation).",
    startCol = 1,
    startRow = 2
  )
  
  mergeCells(wb, sheet = short.common[species], rows = 2, cols = 1:11)
  setRowHeights(wb, sheet = short.common[species], rows = 2, heights = 65)
  addStyle(wb, sheet = short.common[species], style = wrap, rows = 2, cols = 1:11)

  df <- subassemblies.mod %>%
    filter(sp == species) %>%
    select(-sp)
  
  writeData(
    wb,
    sheet = short.common[species],
    df,
    startCol = 1,
    startRow = 4,
    headerStyle = bold
  )
  
  addStyle(wb, sheet = short.common[species], style = two.digit, cols = 8, rows = 5:(dim(df)[1] + 4))
  addStyle(wb, sheet = short.common[species], style = seq.font, cols = 11, rows = 5:(dim(df)[1] + 4))
  
  setColWidths(wb, sheet = short.common[species], cols = 1, widths = 24)
  setColWidths(wb, sheet = short.common[species], cols = c(3, 5, 6, 10, 11), widths = 12)

  freezePane(wb, sheet = short.common[species], firstActiveRow = 5, firstActiveCol = 2)
}

saveWorkbook(wb, '../data/supp_tables/Supplementary Table 1.xlsx', overwrite = TRUE)


## promoter strengths (Supplementary Table 2)
promoter.strength.mod <- promoter.strength %>%
  select(sys, sp, enhancer, light, gene, enrichment) %>%
  mutate(
    enhancer = enhancer.name[as.character(enhancer)],
    light = light.name[as.character(light)]
  ) %>%
  pivot_wider(
    names_from = c(sys, enhancer, light),
    values_from = enrichment
  )

# create an excel file with the promoter strength
wb <- createWorkbook()

modifyBaseFont(wb, fontSize = 10, fontName = 'Arial')

options('openxlsx.numFmt' = '0.00')

bold <- createStyle(textDecoration = 'bold')
wrap <- createStyle(wrapText = TRUE)
boldwrap <- createStyle(textDecoration = 'bold', wrapText = TRUE, halign = 'center')

for (species in all.species$common) {
  addWorksheet(wb, sheetName = species)
}

for (species in all.species$short) {
  writeData(
    wb,
    sheet = short.common[species],
    'Supplementary Table 2 | Strength of plant promoters.',
    startCol = 1,
    startRow = 1
  )
  
  mergeCells(wb, sheet = short.common[species], rows = 1, cols = 1:7)
  addStyle(wb, sheet = short.common[species], style = bold, rows = 1, cols = 1:7)
  
  writeData(
    wb,
    sheet = short.common[species],
    "Promoter strength (log2; normalized to the 35S minimal promoter) was determined by STARR-seq in the indicated condition and assay system. #N/A indicates missing data.",
    startCol = 1,
    startRow = 2
  )
  
  mergeCells(wb, sheet = short.common[species], rows = 2, cols = 1:7)
  setRowHeights(wb, sheet = short.common[species], rows = 2, heights = 26)
  addStyle(wb, sheet = short.common[species], style = wrap, rows = 2, cols = 1:7)
  
  writeData(
    wb,
    sheet = short.common[species],
    list(
      'promoter',
      'no enhancer, dark,\ntobacco leaves',
      'no enhancer, dark,\nmaize protoplasts',
      'with enhancer, dark,\ntobacco leaves',
      'with enhancer, dark,\nmaize protoplasts',
      'no enhancer, light,\ntobacco leaves',
      'with enhancer, light,\ntobacco leaves'
    ),
    startCol = 1,
    startRow = 4
  )
  
  addStyle(wb, sheet = short.common[species], style = boldwrap, rows = 4, cols = 1:7)
  setRowHeights(wb, sheet = short.common[species], rows = 4, heights = 26)  

  df <- promoter.strength.mod %>%
    filter(sp == species) %>%
    select(gene, leaf_noENH_dark, proto_noENH_dark, leaf_withENH_dark, proto_withENH_dark, leaf_noENH_light, leaf_withENH_light)
  
  writeData(
    wb,
    sheet = short.common[species],
    df,
    startCol = 1,
    startRow = 5,
    colNames = FALSE,
    keepNA = TRUE
  )
  
  setColWidths(wb, sheet = short.common[species], cols = 1, widths = 24)
  setColWidths(wb, sheet = short.common[species], cols = 2:7, widths = 18)
  
  freezePane(wb, sheet = short.common[species], firstActiveRow = 5, firstActiveCol = 2)
}

saveWorkbook(wb, '../data/supp_tables/Supplementary Table 2.xlsx', overwrite = TRUE)


## summary of CPE and TF effects on promoter strength (Supplementary Table 4)
difference <- promoter.strength %>%
  inner_join(motif.scores, by = 'gene') %>%
  group_by(sys, sp, enhancer, light, GC.bins) %>%
  summarise(
    across(all_of(colnames(motif.scores)[-1]), function(x) median(enrichment[x >= thresh]) - median(enrichment[x < thresh]))
  ) %>%
  ungroup() %>%
  pivot_longer(
    c(-sys, -sp, -enhancer, -light, -GC.bins),
    names_to = 'motif',
    values_to = 'difference'
  ) %>%
  group_by(sys, enhancer, light, motif) %>%
  summarise(
    mean = mean(difference, na.rm = TRUE),
    sd = sd(difference, na.rm = TRUE),
    signif = wilcox.test(difference)$p.value
  ) %>%
  ungroup() %>%
  mutate(
    enhancer = enhancer.name[as.character(enhancer)],
    light = light.name[as.character(light)]
  ) %>%
  pivot_wider(
    names_from = c(sys, enhancer, light),
    values_from = c(mean, sd, signif)
  )

difference.mod <- difference %>%
  mutate(
    motif = ordered(
      motif,
      levels = c('TATA', 'Ypatch', 'Inr', 'BREu', 'BREd', 'TCT', paste('TF', seq_along(TFs), sep = '_'))
    )
  ) %>%
  arrange(motif) %>%
  select(
    motif,
    mean_leaf_noENH_dark, sd_leaf_noENH_dark, signif_leaf_noENH_dark,
    mean_proto_noENH_dark, sd_proto_noENH_dark, signif_proto_noENH_dark,
    mean_leaf_withENH_dark, sd_leaf_withENH_dark, signif_leaf_withENH_dark,
    mean_proto_withENH_dark, sd_proto_withENH_dark, signif_proto_withENH_dark,
    mean_leaf_noENH_light, sd_leaf_noENH_light, signif_leaf_noENH_light,
    mean_leaf_withENH_light, sd_leaf_withENH_light, signif_leaf_withENH_light,
  )

# create an excel file with the motif summary
wb <- createWorkbook()

options('openxlsx.numFmt' = '0.00')
modifyBaseFont(wb, fontSize = 10, fontName = 'Arial')

sciFmt <- createStyle(numFmt = '0.0E+00')
bold <- createStyle(textDecoration = 'bold')
wrap <- createStyle(wrapText = TRUE)
boldwrap <- createStyle(textDecoration = 'bold', wrapText = TRUE, halign = 'center')
sigStyle <- createStyle(bgFill = '#92D050')

addWorksheet(wb, sheetName = 'effect summary')

writeData(
  wb,
  sheet = 1,
  'Supplementary Table 4 | Summary of the effects of core promoter elements and TF binding sites on promoter strength.',
  startCol = 1,
  startRow = 1
)

mergeCells(wb, sheet = 1, rows = 1, cols = 1:19)
addStyle(wb, sheet = 1, style = bold, rows = 1, cols = 1:19)

writeData(
  wb,
  sheet = 1,
  'For each species, the promoters were grouped by GC content to yield 5 groups of approximately similar size each. The difference in promoter strength between promotors with and without the indicated core promoter element or TF binding site was calculated for each library and assay system and the mean difference (mean) and its standard deviation (sd) across all 15 groups (3 species x 5 GC content groups) are listed here. The Wilcoxon rank-sum test was used to compare the differences to a null distribution and the corresponding p-value is indicated. Significant (p-value < 0.0005) effects are highlighted in green. TF_x, binding site for TFs from cluster x.',
  startCol = 1,
  startRow = 2
)

mergeCells(wb, sheet = 1, rows = 2, cols = 1:19)
setRowHeights(wb, sheet = 1, rows = 2, heights = 39)
addStyle(wb, sheet = 1, style = wrap, rows = 2, cols = 1:19)

writeData(
  wb,
  sheet = 1,
  as.list(
    c(
      'element/TF binding site',
      rep(
        c(
          'no enhancer, dark,\ntobacco leaves',
          'no enhancer, dark,\nmaize protoplasts',
          'with enhancer, dark,\ntobacco leaves',
          'with enhancer, dark,\nmaize protoplasts',
          'no enhancer, light,\ntobacco leaves',
          'with enhancer, light,\ntobacco leaves'
        ),
        each = 3
      )
    )
  ),
  startCol = 1,
  startRow = 4
)
mergeCells(wb, sheet = 1, rows = 4:5, cols = 1)
mergeCells(wb, sheet = 1, rows = 4, cols = 2:4)
mergeCells(wb, sheet = 1, rows = 4, cols = 5:7)
mergeCells(wb, sheet = 1, rows = 4, cols = 8:10)
mergeCells(wb, sheet = 1, rows = 4, cols = 11:13)
mergeCells(wb, sheet = 1, rows = 4, cols = 14:16)
mergeCells(wb, sheet = 1, rows = 4, cols = 17:19)

addStyle(wb, sheet = 1, style = boldwrap, rows = 4, cols = 1:19)
setRowHeights(wb, sheet = 1, rows = 4, heights = 26)
setColWidths(wb, sheet = 1, cols = 1, widths = 12)

writeData(
  wb,
  sheet = 1,
  as.list(rep(c('mean', 'sd', 'p-value'), times = 6)),
  startCol = 2,
  startRow = 5
)

addStyle(wb, sheet = 1, style = bold, rows = 5, cols = 1:19)

writeData(
  wb,
  sheet = 1,
  difference.mod,
  startCol = 1,
  startRow = 6,
  colNames = FALSE
)

for (col in seq(4, 19, 3)) {
  addStyle(wb, sheet = 1, style = sciFmt, rows = 6:83, cols = col)
  conditionalFormatting(wb, 1, cols = c(col - 2, col), rows = 6:83, rule = paste0('$', LETTERS[col],'6 < 5e-4'), style = sigStyle)
}

freezePane(wb, sheet = 1, firstActiveRow = 6, firstActiveCol = 2)

saveWorkbook(wb, '../data/supp_tables/Supplementary Table 4.xlsx', overwrite = TRUE)


## Evolved promoters (Supplemetary Table 5)
evolved.promoters <- promoter.strength.val %>%
  filter(motif == 'evolution' | motif == 'synthetic' | grepl('native', motif, fixed = 'TRUE')) %>%
  mutate(
    gene = if_else(motif == 'synthetic', paste0(gene, '(', variant, ')'), gene),
    variant = if_else(motif == 'evolution', variant, '-0')
  ) %>%
  separate(
    variant,
    into = c('opt_for', 'round'),
    sep = '-'
  ) %>%
  mutate(
    round = as.numeric(round)
  ) %>%
  select(sys, enhancer, light, gene, enrichment, opt_for, round) %>%
  mutate(
    enhancer = enhancer.name[as.character(enhancer)],
    light = light.name[as.character(light)]
  ) %>%
  pivot_wider(
    names_from = c(sys, enhancer, light),
    values_from = enrichment
  )

evolved.seqs <- read_tsv('../CNN/evolution_data.tsv') %>%
  select('gene' = origin, opt_for, round, sequence) %>%
  group_by(gene, round, sequence) %>%
  summarise(
    opt_for = paste0(opt_for, collapse = '/')
  ) %>%
  ungroup() %>%
  mutate(
    opt_for = if_else(opt_for == 'start', '', opt_for)
  )

evolved.promoters.mod <- evolved.promoters %>%
  left_join(evolved.seqs, by = c('gene', 'opt_for', 'round')) %>%
  mutate(
    opt_for = gsub('leaf', 'tobacco', opt_for, fixed = TRUE),
    opt_for = gsub('proto', 'maize', opt_for, fixed = TRUE)
  )

# 35S promoter strength
strength.35S <- promoter.strength.val %>%
  filter(set == 'PROevo' & ((enhancer & variant == 'withPRO-withENH') | (! enhancer & variant == 'withPRO-noENH'))) %>%
  mutate(
    enhancer = enhancer.name[as.character(enhancer)],
    light = light.name[as.character(light)]
  ) %>%
  mutate(
    experiment = paste(sys, enhancer, light, sep = '_')
  ) %>%
  select(experiment, enrichment) %>%
  deframe()

# create an excel file with strength of evolved promoters
wb <- createWorkbook()

modifyBaseFont(wb, fontSize = 10, fontName = 'Arial')

options('openxlsx.numFmt' = '0.00')

bold <- createStyle(textDecoration = 'bold')
wrap <- createStyle(wrapText = TRUE)
boldwrap <- createStyle(textDecoration = 'bold', wrapText = TRUE, halign = 'center')
no.decimal <- createStyle(numFmt = '0')
seq.font <- createStyle(fontName = 'Courier New')
strongStyle <- createStyle(bgFill = '#92D050')

addWorksheet(wb, sheetName = 'evolved promoters')

writeData(
  wb,
  sheet = 1,
  'Supplementary Table 5 | Strength of evolved promoters.',
  startCol = 1,
  startRow = 1
)

mergeCells(wb, sheet = 1, rows = 1, cols = 1:10)
addStyle(wb, sheet = 1, style = bold, rows = 1, cols = 1:10)

writeData(
  wb,
  sheet = 1,
  "Promoter strength (log2; normalized to the 35S minimal promoter) was determined by STARR-seq in the indicated condition and asssay system. In silico evolution was performed for the indicated number of rounds using a CNN model trained on data from the tobacco leaf system (tobacco), a CNN model trained on data from the maize protoplast system (maize), or both models together (both). The names for synthetic promoters start with synA for synthetic promoters with an Arabidopsis-like nucleotide frequency and with synZ for promoters with a maize-like base composition. Promoters stronger than the 35S minimal promoter are highlighted in green. #N/A indicates missing data.",
  startCol = 1,
  startRow = 2
)

mergeCells(wb, sheet = 1, rows = 2, cols = 1:10)
setRowHeights(wb, sheet = 1, rows = 2, heights = 52)
addStyle(wb, sheet = 1, style = wrap, rows = 2, cols = 1:10)

writeData(
  wb,
  sheet = 1,
  list(
    'promoter',
    'model',
    'round',
    'no enhancer, dark,\ntobacco leaves',
    'no enhancer, dark,\nmaize protoplasts',
    'with enhancer, dark,\ntobacco leaves',
    'with enhancer, dark,\nmaize protoplasts',
    'no enhancer, light,\ntobacco leaves',
    'with enhancer, light,\ntobacco leaves',
    'sequence'
  ),
  startCol = 1,
  startRow = 4
)

addStyle(wb, sheet = 1, style = boldwrap, rows = 4, cols = 1:10)
setRowHeights(wb, sheet = 1, rows = 4, heights = 26)  

df <- evolved.promoters.mod %>%
  select(gene, opt_for, round, leaf_noENH_dark, proto_noENH_dark, leaf_withENH_dark, proto_withENH_dark, leaf_noENH_light, leaf_withENH_light, sequence)

writeData(
  wb,
  sheet = 1,
  df,
  startCol = 1,
  startRow = 5,
  colNames = FALSE,
  keepNA = TRUE
)

addStyle(wb, sheet = 1, style = no.decimal, cols = 3, rows = 5:(dim(df)[1] + 4))
addStyle(wb, sheet = 1, style = seq.font, cols = 10, rows = 5:(dim(df)[1] + 4))

for (i in seq_len(6)) {
  experiment <- c('leaf_noENH_dark', 'proto_noENH_dark', 'leaf_withENH_dark', 'proto_withENH_dark', 'leaf_noENH_light', 'leaf_withENH_light')[i]
  conditionalFormatting(wb, 1, cols = 3 + i, rows = 5:(dim(df)[1] + 4), rule = paste0('>', strength.35S[experiment]), style = strongStyle)
}


setColWidths(wb, sheet = 1, cols = 1, widths = 24)
setColWidths(wb, sheet = 1, cols = 4:9, widths = 18)
setColWidths(wb, sheet = 1, cols = 10, widths = 12)

freezePane(wb, sheet = 1, firstActiveRow = 5, firstActiveCol = 2)

saveWorkbook(wb, '../data/supp_tables/Supplementary Table 5.xlsx', overwrite = TRUE)
