### define functions to export data for plotting in LaTeX with pgfplots

## function to export data for boxplots
# data:             R object; data frame to export
# samples_from:     unquoted column specification; column containing the sample names by which to summarise
# values_from:      unquoted column specification; column containing the values to be summarised and exported
# file:             character; base name of the export file without extension (full name is: "<file>_boxplot.tsv")
# ...:              unquoted column specifications; additional columns to include in the export; should contain only a single value per sample!
# outliers:         boolean; whether to create a file for the outlier points (will be saved as: "<file>_boxplot_outliers.tsv")
# p_values:         boolean; whether to calculate all possible pairwise p values using a Mann-Whitney-Wilcox test (will be saved as: "<file>_pvalues.tsv")
# use_sample_name:  boolean; whether to use the sample name instead of a numeric index as column name for the outliers
# data_range:       numeric vector, length 2; min and max coordinates to be considered for pvalue coordinates (the range of <data> is used if NULL)

LaTeX.boxplot <- function(data, samples_from, values_from, file, ..., outliers = TRUE, p_values = TRUE, use_sample_name = FALSE, data_range = NULL) {
  sample_col <- enquo(samples_from)
  value_col <- enquo(values_from)
  additional_cols <- enquos(...)
  
  box <- data %>%
    select(sample = !! sample_col, value =  !! value_col, !!! additional_cols) %>%
    group_by(sample, !!! additional_cols) %>%
    summarise(
      med = median(value, na.rm = TRUE),
      lq = quantile(value, .25, na.rm = TRUE),
      uq = quantile(value, .75, na.rm = TRUE),
      iqr = 1.5 * IQR(value, na.rm = TRUE),
      lw = max(lq - iqr, min(value)),
      uw = min(uq + iqr, max(value)),
      n = n()
    ) %>%
    select(-iqr) %>%
    ungroup() %>%
    arrange(sample) %>%
    mutate(
      id = seq_len(n())
    )
  
  if (length(box$sample) != length(unique(box$sample))) {
    warning('The additional columns contain more than one value per sample. The results might be unexpected!')
  }
  
  write_tsv(box, paste0(file, '_boxplot.tsv'), na = 'NaN')
  
  if (outliers) {
    outlier.data <- data %>%
      select(value =  !! value_col, sample = !! sample_col) %>%
      group_by(sample) %>%
      summarise(
        outlier = list(value[outliers(value) != 0]),
        id = list(seq_along(unlist(outlier)))
      ) %>%
      unnest(c(outlier, id), keep_empty = TRUE) %>%
      ungroup() %>%
      arrange(sample) %>%
      pivot_wider(
        names_from = sample,
        values_from = outlier,
        id_cols = id
      ) %>%
      select(-id)
    
    if (use_sample_name) {
      colnames(outlier.data) <- paste0('outlier.', colnames(outlier.data))
    } else {
      colnames(outlier.data) <- paste0('outlier.', seq_len(dim(outlier.data)[2])) 
    }
    
    write_tsv(outlier.data, paste0(file, '_boxplot_outliers.tsv'), na = 'NaN')
  }
  
  if (p_values) {
    wilcox.all(
      data = data,
      samples_from = !! sample_col,
      values_from = !! value_col,
      file = file,
      data_range = data_range,
      use_sample_name = use_sample_name
    )
  }
}


## function to export data for violin plots
# data:             R object; data frame to export
# samples_from:     unquoted column specification; column containing the sample names by which to summarise
# values_from:      unquoted column specification; column containing the values to be summarised and exported
# file:             character; base name of the export file without extension (full name is: "<file>_violin.tsv")
# boxpllot:         boolean; whether to create a file for a boxplot (will be saved as: "<file>_boxplot.tsv")
# ...:              unquoted column specifications; additional columns to include in the boxplot export; should contain only a single value per sample!
# outliers:         boolean; whether to create a file for the outlier points (will be saved as: "<file>_boxplot_outliers.tsv")
# p_values:         boolean; whether to calculate all possible pairwise p values using a Mann-Whitney-Wilcox test (will be saved as: "<file>_pvalues.tsv")
# half:             boolean; whether to export only one half of the violinplot for split violins
# use_sample_name:  boolean; whether to use the sample name instead of a numeric index as column name
# data_range:       numeric vector, length 2; min and max coordinates to be considered for pvalue coordinates (the range of <data> is used if NULL)

LaTeX.violinplot <- function(
  data, samples_from, values_from, file, ..., boxplot = TRUE, outliers = FALSE, p_values = TRUE, half = FALSE, use_sample_name = FALSE, data_range = NULL
) {
  sample_col <- enquo(samples_from)
  value_col <- enquo(values_from)
  
  if (boxplot) {
    LaTeX.boxplot(
      data = data,
      samples_from = !! sample_col,
      values_from = !! value_col,
      file = file,
      ...,
      p_values = FALSE,
      outliers = outliers,
      use_sample_name = use_sample_name,
      data_range = data_range
    )
  }
  
  violin <- data %>%
    select(value =  !! value_col, sample = !! sample_col) %>%
    group_by(sample) %>%
    summarise(
      x = list(density(value, from = min(value), to = max(value), n = 256)$y),
      y = list(density(value, from = min(value), to = max(value), n = 256)$x)
    ) %>%
    arrange(sample) %>%
    ungroup()
  
  if (! use_sample_name) {
    violin <- violin %>%
      mutate(
        sample = ordered(seq_len(n()))
      )
  }
  
  violin <- violin %>%
    pivot_wider(names_from = sample, values_from = c(x, y), names_sep = '.') %>%
    unnest(everything()) %>%
    ungroup %>%
    mutate(
      id = row_number()
    )
  
  # normalize width (width of widest violin is set to 1)
  max.x <- max(select_at(violin, vars(starts_with('x.'))))
  
  violin <- violin %>%
    mutate_at(
      vars(starts_with('x.')),
      function(x) {x / max.x * .5}
    )
  
  if (half) {
    max.y <- sapply(colnames(violin), function(x) max(pull(violin, x)))
    max.y[grepl('x.', names(max.y), fixed = TRUE)] <- 0
    min.y <- sapply(colnames(violin), function(x) min(pull(violin, x)))
    min.y[grepl('x.', names(min.y), fixed = TRUE)] <- 0
    
    violin <- bind_rows(min.y, violin, max.y) %>%
      select(-id)
  } else {
    violin <- bind_rows(violin, arrange(mutate_at(violin, vars(starts_with('x.')), function(x) {-x}), desc(id))) %>%
      select(-id)
  }
  
  write_tsv(violin, paste0(file, '_violin.tsv'), na = 'NaN')
  
  if (p_values) {
    wilcox.all(
      data = data,
      samples_from = !! sample_col,
      values_from = !! value_col,
      file = file,
      data_range = data_range,
      use_sample_name = use_sample_name,
      half = half
    )
  }
}


## calculate and export all pairwise p values using a Mann-Whitney-Wilcox test
# data:             R object; data frame to export
# samples_from:     unquoted column specification; column containing the sample names by which to summarise
# values_from:      unquoted column specification; column containing the values to be summarised and exported
# file:             character; base name of the export file without extension (full name is: "<file>_pvalues.tsv")
# use_sample_name:  boolean; whether to use the sample name instead of a numeric index as column name
# data_range:       numeric vector, length 2; min and max coordinates to be considered for pvalue coordinates (the range of <data> is used if NULL)
# half:             boolean; whether to export a simplified version of the pvalues for half violin plots

wilcox.all <- function(data, samples_from, values_from, file, data_range = NULL, use_sample_name = FALSE, half = FALSE) {
  sample_col <- enquo(samples_from)
  value_col <- enquo(values_from)
  
  if (is.null(data_range)) {
    data_range <- range(pull(data, !! value_col))
  } else if (length(data_range) != 2 | ! is.numeric(data_range)) {
    stop('"data_range" must be a two component numerical vector')
  }
  
  samples <- data %>% distinct(!! sample_col) %>% arrange(!! sample_col) %>% pull(!! sample_col)
  
  if (! is.ordered(samples)) {
    samples <- ordered(samples) 
  }
  
  
  pvalues <- tibble(
    exp1 = combn(samples, 2)[1, ],
    exp2 = combn(samples, 2)[2, ]
  ) %>%
    rowwise() %>%
    mutate(
      p.value = wilcox.test(
        pull(filter(data, !! sample_col == exp1), !! value_col),
        pull(filter(data, !! sample_col == exp2), !! value_col)
      )$p.value,
      max1 = max(pull(filter(data, !! sample_col == exp1), !! value_col)),
      max2 = max(pull(filter(data, !! sample_col == exp2), !! value_col)),
      x1 = as.numeric(exp1),
      x2 = as.numeric(exp2)
    )
  
  if (half) {
    pvalues <- pvalues %>%
      filter(as.numeric(exp1) + 1 == as.numeric(exp2) & as.numeric(exp1) %% 2 == 1) %>%
      mutate(
        max1 = max(max1, max2) - 0.1 * diff(data_range),
        max2 = max1
      )
  }
  
  if (! use_sample_name) {
    pvalues <- pvalues %>%
      mutate(
        exp1 = x1,
        exp2 = x2
      )
  }
  
  pvalues <- pvalues %>%
    mutate(
      x = list(c(rep(x1, 2), (x1 + x2) / 2, rep(x2, 2))),
      y = list(c(max1, rep(data_range[2] + 0.1 * diff(data_range), 3), max2)),
      p.value = list(c(rep(NA_real_, 2), p.value, rep(NA_real_, 2))),
      comparison = paste(exp1, exp2, sep = '_')
    ) %>%
    select(comparison, x, y, p.value) %>%
    pivot_wider(
      names_from = comparison,
      values_from = c(x, y, p.value),
      names_sep = '.'
    ) %>%
    unnest(everything())
  
  write_tsv(pvalues, paste0(file, '_pvalues.tsv'), na = 'NaN')
}


## function to export a PWM logo
# PWM:             R object; the PWM as a universalmotif class object
# file:            character; name of the export file

PWM.to.LaTeX <- function(motif, file) {
  icm <- as_tibble(t(convert_type(motif, 'ICM')['motif'])) %>%
    mutate(
      pos = seq_len(n())
    ) %>%
    pivot_longer(
      -pos,
      names_to = 'base',
      values_to = 'IC'
    ) %>%
    group_by(pos) %>%
    arrange(IC, .by_group = TRUE) %>%
    mutate(
      plot = seq_len(n())
    ) %>%
    pivot_wider(
      names_from = plot,
      values_from = c(IC, base)
    )
  
  write_tsv(icm, file, na = 'NaN')
}


## annotate outliers (high outliers -> 1, low outliers -> -1, no outlier -> 0)
# data: numeric vector; data for which to annotate outliers
outliers <- function(data) {
  case_when(
    data > quantile(data, .75, na.rm = TRUE) + 1.5 * IQR(data, na.rm = TRUE) ~ 1,
    data < quantile(data, .25, na.rm = TRUE) - 1.5 * IQR(data, na.rm = TRUE) ~ -1,
    TRUE ~ 0
  )
}