library(readr)
library(dplyr)
library(tidyr)
library(tibble)


### load main data
load('RData/main_data.Rdata')


### load motif scores
motif.scores <- read_tsv('motifs/motif-scores.tsv.gz')


### load promoter sequences
all.sequences <- bind_rows(
  'At' = read_tsv('../data/promoter_annotation/Arabidopsis_all_promoters_unique.tsv'),
  'Zm' = read_tsv('../data/promoter_annotation/Maize_all_promoters_unique.tsv'),
  'Sb' = read_tsv('../data/promoter_annotation/Sorghum_all_promoters_unique.tsv'),
  .id = 'sp'
)


### split data into train and test set (libraries with enhancer in dark)
set.seed(1)

promoter.strength.mod <- promoter.strength %>%
  filter(enhancer & ! light) %>%
  inner_join(motif.scores, by = 'gene') %>%
  group_by(sys, sp, enhancer, light) %>%
  mutate(
    set = factor(sample(c('train', 'test'), n(), replace = TRUE, prob = c(.9, .1)))
  ) %>%
  ungroup()


# ### save data to train a CNN
# CNN.data <- promoter.strength.mod %>%
#   inner_join(select(all.sequences, gene, sequence), by = 'gene') %>%
#   select(sys, set, sp, gene, sequence, enrichment)
# 
# for (system in unique(CNN.data$sys)) {
#   write_tsv(CNN.data %>% filter(sys == system & set == 'train') %>% select(-sys), paste0('../CNN/CNN_train_', system, '.tsv'))
#   write_tsv(CNN.data %>% filter(sys == system & set == 'test') %>% select(-sys), paste0('../CNN/CNN_test_', system, '.tsv'))
# }


### define terms for linear model
CPEs <- c('TATA', 'Inr', 'TCT', 'Ypatch', 'BREu', 'BREd')
TFs <- motif.scores %>% select(starts_with('TF')) %>% names()

terms <- c('GC', CPEs, TFs)

model.formula <- as.formula(paste('enrichment', paste(terms, collapse = ' + '), sep = ' ~ '))


### train and evaluate linear model for each system
# build the model
models <- promoter.strength.mod %>%
  filter(set == 'train') %>%
  nest_by(sys) %>% 
  summarise(
    model = list(lm(
      model.formula,
      data = data
    )),
    .groups = 'drop'
  )

# predict values
predictions <- promoter.strength.mod %>%
  nest_by(sys) %>%
  inner_join(models, by = 'sys') %>%
  mutate(
    set = list(data$set),
    sp = list(data$sp),
    enrichment = list(data$enrichment),
    prediction = list(predict(model, data))
  ) %>%
  select(sys, set, sp, enrichment, prediction) %>%
  unnest(c(set, sp, enrichment, prediction)) %>%
  ungroup()

# calculate correlation
correlation <- predictions %>%
  group_by(sys, set) %>%
  summarise(
    sp = 'all',
    pearson = cor(enrichment, prediction, use = 'complete.obs'),
    spearman =  cor(enrichment, prediction, use = 'complete.obs', method = 'spearman'),
    rsquare = pearson^2
  ) %>%
  ungroup()

correlation.sp <- predictions %>%
  group_by(sys, set, sp) %>%
  summarise(
    pearson = cor(enrichment, prediction, use = 'complete.obs'),
    spearman =  cor(enrichment, prediction, use = 'complete.obs', method = 'spearman'),
    rsquare = pearson^2
  ) %>%
  ungroup()

correlation.all <- bind_rows(correlation, correlation.sp)


### export results for plotting in LaTeX
for (system in unique(predictions$sys)) {
  df <- predictions %>%
    filter(set == 'test' & sys == system) %>%
    select('sample.name' = sp, enrichment, prediction) %>%
    slice_sample(prop = 1)
  
  cor.df <- correlation.all %>%
    filter(set == 'test' & sys == system) %>%
    mutate(
      spearman = sprintf('$\\rho=%.2f$', spearman),
      rsquare = sprintf('$R^2=%.2f$', rsquare)
    ) %>%
    select('sample.name' = sp, spearman, rsquare)
  
  write_tsv(df, paste0('../figures/rawData/linear-model_', system, '_pred.tsv'), na = 'NaN')
  write_tsv(cor.df, paste0('../figures/rawData/linear-model_', system, '_stats.tsv'), na = 'NaN') 
}
