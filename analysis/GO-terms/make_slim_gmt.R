library(readr)
library(dplyr)
library(tidyr)
library(gprofiler2)


### set species names
all.species <- tibble(
  common = c('Maize', 'Sorghum', 'Arabidopsis'),
  short = c('Zm', 'Sb', 'At'),
  latin = c('zmays', 'sbicolor', 'athaliana')
)


### get source files
# files to map GO terms to GOslim terms from the AgBase GOSlimViewer (https://agbase.arizona.edu/cgi-bin/tools/goslimviewer_select.pl)
if (! file.exists('GO-terms/term2slim_plant.tsv') | ! file.exists('GO-terms/term2slim_generic.tsv') | ! file.exists('GO-terms/godef.tsv')) {
  download.file('http://agbase.arizona.edu/tools/goslimviewer/goslimviewer_standalone.zip', 'GO-terms/goslimviewer_standalone.zip')
  unzip('GO-terms/goslimviewer_standalone.zip', c('term2slim_plant.tsv', 'term2slim_generic.tsv', 'godef.tsv'), exdir = 'GO-terms')
}

# gmt files to map GO terms to Arabidopsis, maize, and sorghum gene names from g:Profiler (https://biit.cs.ut.ee/gprofiler/gost)
for (species in all.species$latin) {
  if (! file.exists(paste0('GO-terms/gprofiler_full_', species, '.ENSG.gmt'))) {
    download.file(paste0('https://biit.cs.ut.ee/gprofiler/static/gprofiler_full_', species, '.ENSG.gmt'), paste0('GO-terms/gprofiler_full_', species, '.ENSG.gmt'))
  }
}


### load mapping files
go.def <- read_tsv('GO-terms/godef.tsv', col_types = 'c-c', col_names = c('GO', 'name'))
go2slimP <- read_tsv('GO-terms/term2slim_plant.tsv', col_names = c('GO', 'GOslim'))
go2slimG <- read_tsv('GO-terms/term2slim_generic.tsv', col_names = c('GO', 'GOslim'))
go2slim <- union(go2slimP, go2slimG) %>%
  bind_rows(
    # add a GOslim term for nucleosomes
    tibble( # nucleosome
      'GO' = c('GO:0000786', 'GO:0000787', 'GO:0000788', 'GO:0043505'),
      'GOslim' = 'GO:0000786'
    )
  )


### map GO terms to GO slim terms
for (species in all.species$latin) {
  gaf <- read_fwf(paste0('GO-terms/gprofiler_full_', species,'.ENSG.gmt'), fwf_positions(c(1, 12), c(11, NA), col_names = c('GO', 'gene'))) %>%
    mutate(
      gene = strsplit(gene, '\t')
    ) %>%
    unnest(gene) %>%
    group_by(GO) %>%
    filter(gene != first(gene)) %>%
    ungroup()
  
  gaf.slim <- gaf %>%
    inner_join(go2slim, by = 'GO')
  
  gmt.slim <- gaf.slim %>%
    select(- GO) %>%
    rename('GO' = GOslim) %>%
    group_by(GO) %>%
    summarise(
      gene = paste0(unique(unlist(gene)), collapse = '\t')
    ) %>%
    inner_join(go.def, by = 'GO') %>%
    select(GO, name, gene)
  
  write.table(gmt.slim, paste0('GO-terms/gprofiler_slim_', species,'.ENSG.gmt'), quote = FALSE, sep = '\t', col.names = FALSE, row.names = FALSE)
  
  # upload GMT files to gprofiler
  # gmt.ID <- upload_GMT_file(paste0('GO-terms/gprofiler_slim_', species,'.ENSG.gmt'))
  
  # assign(paste0('gmt.ID.', species), gmt.ID)
}


### save IDs for the custom slim GMT files at gprofiler
custom.gmt <- tibble(
    sp = all.species$latin
  ) %>%
  rowwise() %>%
  mutate(
    ID = get(paste0('gmt.ID.', sp))
  ) %>%
  ungroup()

write_tsv(custom.gmt, 'GO-terms/gprofiler_gmt_IDs.tsv')