library(readr)
library(dplyr)
library(tidyr)

### load and preprocess experiment data (barcode count files)
# input:        character or connection; data file with input read counts
# dark:         character or connection; data file with output read counts from dark condition
# light:        optional; character or connection; data file with output read counts from light condition (only available for tobacco leaf system)
# subassembly:  R object; subassembly data frame
# rcc:          numeric vector with two or three elements; read count cutoffs to be applied to the input and output;
#               first element: cutoff for input; second element: cutoff for output 'dark'; third element: cutoff for output 'light' (uses second element if no third element given)
load.experiment.bc <- function(input, dark, light = NA, subassembly, rcc = c(1, 1, 1)) {
  
  # configure read count cutoff
  if (length(rcc) == 2) {
    rcc[3] <- rcc[2]
  } else if (length(rcc) != 3 | ! is.numeric(rcc)) {
    stop('rcc must contain two to three numeric values')
  }
  
  # load and merge files (only keep barcodes detected in both input and output)
  data.inp <- read_table(input, col_names = c('count', 'barcode'))
  
  data.both <- read_table(dark, col_names = c('count', 'barcode')) %>%
    inner_join(data.inp, by = 'barcode', suffix = c('.out', '.inp')) %>%
    mutate(
      light = FALSE
    )
  
  # add 'light' data if supplied
  if (! is.na(light)) {
    data.light <- read_table(light, col_names = c('count', 'barcode')) %>%
      inner_join(data.inp, by = 'barcode', suffix = c('.out', '.inp')) %>%
      mutate(
        light = TRUE
      )
    
    data.both <- bind_rows(data.both, data.light)
  }
  
  # merge with subassembly
  data.both <- inner_join(data.both, subassembly, by = 'barcode')
  
  # read count cutoff
  data.both <- data.both %>%
    filter(
      (count.inp >= rcc[1]) & (count.out >= if_else(light, rcc[3], rcc[2]))
    )
  
  # calculate enrichment
  data.both <- data.both %>%
    group_by(light) %>%
    mutate(
      enrichment = log2((count.out / sum(count.out)) / (count.inp / sum(count.inp)))
    ) %>%
    ungroup()
  
  # return result
  return(data.both)
}



### native promoter libraries

## load subassemblies
subassembly.At <- read_tsv('subassembly/subassembly_At.tsv')
subassembly.Zm <- read_tsv('subassembly/subassembly_Zm.tsv')
subassembly.Sb <- read_tsv('subassembly/subassembly_Sb.tsv')


## load barcode data

# create helper table with file names (requires systematic file names)
sys <- c('leaf', 'proto')
sp <- c('At', 'Zm', 'Sb')
enhancer <- c(TRUE, FALSE)
rep <- c(1, 2)

# load data for genomic promoter libraries
all.experiments <- expand_grid(sys, sp, enhancer, rep) %>%
  nest_by(across(everything())) %>%
  summarise(
    basename = paste0(
      '../data/barcode_counts/', sys, '/', sp, '_Rep', rep, '/barcodes_pPSup_', sp, 'PRO_', if_else(enhancer, '35SEnh', 'noEnh'), '_', as.roman(rep)
    ),
    input_name = if_else(sys == 'proto' & sp == 'At' & rep == 2, gsub('Rep2', 'Rep1', gsub('II', 'I', basename, fixed = TRUE), fixed = TRUE), basename),
    load.experiment.bc(
      input = paste0(input_name, '_inp.count.gz'),
      dark = paste0(basename, '_dark.count.gz'),
      light = if_else(sys == 'leaf', paste0(basename, '_light.count.gz'), NA_character_),
      subassembly = get(paste0('subassembly.', sp)),
      rcc = c(5, 5, 5)
    )
  ) %>%
  select(-basename, -input_name) %>%
  ungroup()


## filter for WT sequences only
all.experiments <- all.experiments %>%
  filter(grepl('control', variant, fixed = TRUE) | (variant == 'WT' & FL) & gene != '35Spr')


## normalize replicates to no enhancer control
all.experiments <- all.experiments %>%
  group_by(sys, sp, enhancer, rep, light) %>%
  mutate(
    enrichment = enrichment - median(enrichment[variant == 'control-withPRO-noENH'])
  ) %>%
  ungroup()


## aggregate by promoter
ag.experiments <- all.experiments %>%
  filter(! grepl('control', variant, fixed = TRUE)) %>%
  group_by(across(c(-count.inp, -count.out, -barcode, -enrichment))) %>%
  summarise(
    n.bc = n(),
    sd = sd(enrichment),
    enrichment = median(enrichment)
  ) %>%
  ungroup()


## calculate mean across replicates
promoter.strength <- ag.experiments %>%
  group_by(across(c(-n.bc, -sd, -enrichment, -rep))) %>%
  summarise(
    n.experiments = n(),
    min.bc = min(n.bc),
    enrichment = mean(enrichment)
  ) %>%
  ungroup()


## calculate light dependency
light.dependency <- promoter.strength %>%
  filter(! is.na(type)) %>%
  group_by(across(c(-light, -enrichment, -n.experiments, -min.bc))) %>%
  filter(n() == 2) %>%
  summarise(
    min.bc = min(min.bc),
    min.experiments = min(n.experiments),
    dark = enrichment[! light],
    light = enrichment[light],
    light.dependency = light - dark
  ) %>%
  ungroup()


## calculate 35S enhancer effect
enhancer.effect <- promoter.strength %>%
  filter(! is.na(type)) %>%
  group_by(across(c(-enhancer, -enrichment, -n.experiments, -min.bc))) %>%
  filter(n() == 2) %>%
  summarise(
    min.bc = min(min.bc),
    min.experiments = min(n.experiments),
    withENH = enrichment[enhancer],
    noENH = enrichment[! enhancer],
    enhancer.effect = withENH - noENH
  ) %>%
  ungroup()


## save R objects
if (! file.exists('RData')){
  dir.create('RData')
}

save(promoter.strength, light.dependency, enhancer.effect, file = 'RData/main_data.Rdata')
save(ag.experiments, file = 'RData/ag_experiments.Rdata')

controls <- all.experiments %>%
  filter(grepl('control', variant))

save(controls, file = 'RData/controls.Rdata')



### validation libraries

## load subassemblies
subassembly.PROval <- read_tsv('subassembly/subassembly_PROval.tsv')
subassembly.PROevo <- read_tsv('subassembly/subassembly_PROevo.tsv')


## load barcode data

# create helper table with file names (requires systematic file names)
sys <- c('leaf', 'proto')
set <- c('PROval', 'PROevo')
enhancer <- c(TRUE, FALSE)
rep <- c(1, 2)

# load data for genomic promoter libraries
all.experiments.val <- expand_grid(sys, set, enhancer, rep) %>%
  nest_by(across(everything())) %>%
  summarise(
    basename = paste0(
      '../data/barcode_counts/', sys, '/', set, '_Rep', rep, '/barcodes_pPSup_', set, '_', if_else(enhancer, '35SEnh', 'noEnh'), '_', as.roman(rep)
    ),
    input_name = if_else(sys == 'proto' & rep == 2, gsub('Rep2', 'Rep1', gsub('II', 'I', basename, fixed = TRUE), fixed = TRUE), basename),
    load.experiment.bc(
      input = paste0(input_name, '_inp.count.gz'),
      dark = paste0(basename, '_dark.count.gz'),
      light = if_else(sys == 'leaf', paste0(basename, '_light.count.gz'), NA_character_),
      subassembly = get(paste0('subassembly.', set)),
      rcc = c(5, 5, 5)
    )
  ) %>%
  select(-basename, -input_name) %>%
  ungroup()


## filter for WT sequences only
all.experiments.val <- all.experiments.val %>%
  filter(grepl('control', variant, fixed = TRUE) | (variant == 'WT' & FL)) %>%
  mutate(
    gene = if_else(grepl('control', variant, fixed = TRUE), paste(variant, 'control', sub('control-', '', variant, fixed = TRUE), sep = '.'), gene)
  ) %>%
  select(-variant, -start, -stop, -FL)


## split "gene" into original annotation
all.experiments.val <- all.experiments.val %>%
  separate(
    gene,
    into = c('gene', 'motif', 'variant'),
    sep = '\\.'
  )


## normalize replicates to no enhancer control
all.experiments.val <- all.experiments.val %>%
  group_by(sys, set, enhancer, rep, light) %>%
  mutate(
    enrichment = enrichment - median(enrichment[gene == 'control-withPRO-noENH'])
  ) %>%
  ungroup()


## aggregate by promoter
ag.experiments.val <- all.experiments.val %>%
  filter(! grepl('control', variant, fixed = TRUE)) %>%
  group_by(across(c(-count.inp, -count.out, -barcode, -enrichment))) %>%
  summarise(
    n.bc = n(),
    sd = sd(enrichment),
    enrichment = median(enrichment)
  ) %>%
  ungroup()


## calculate mean across replicates
promoter.strength.val <- ag.experiments.val %>%
  group_by(across(c(-n.bc, -sd, -enrichment, -rep))) %>%
  summarise(
    n.experiments = n(),
    min.bc = min(n.bc),
    enrichment = mean(enrichment)
  ) %>%
  ungroup()


## calculate light dependency
light.dependency.val <- promoter.strength.val %>%
  filter(motif != 'control') %>%
  group_by(across(c(-light, -enrichment, -n.experiments, -min.bc))) %>%
  filter(n() == 2) %>%
  summarise(
    min.bc = min(min.bc),
    min.experiments = min(n.experiments),
    dark = enrichment[! light],
    light = enrichment[light],
    light.dependency = light - dark
  ) %>%
  ungroup()


## calculate 35S enhancer effect
enhancer.effect.val <- promoter.strength.val %>%
  filter(motif != 'control') %>%
  group_by(across(c(-enhancer, -enrichment, -n.experiments, -min.bc))) %>%
  filter(n() == 2) %>%
  summarise(
    min.bc = min(min.bc),
    min.experiments = min(n.experiments),
    withENH = enrichment[enhancer],
    noENH = enrichment[! enhancer],
    enhancer.effect = withENH - noENH
  ) %>%
  ungroup()


## save R objects
if (! file.exists('RData')){
  dir.create('RData')
}

save(promoter.strength.val, light.dependency.val, enhancer.effect.val, file = 'RData/main_data_validation.Rdata')
save(ag.experiments.val, file = 'RData/ag_experiments_validation.Rdata')

controls.val <- all.experiments.val %>%
  filter(motif == 'control')

save(controls.val, file = 'RData/controls_validation.Rdata')
