library(dplyr)
library(magrittr)
source("constants.R")

# Set seed
seed <- 452016632
set.seed(seed, sample.kind = "Rounding")

# Read sample metadata
metadata_df <- read.csv(metadata_csv, stringsAsFactors = F)

# Restrict to Definite/No Evidence samples for generating cross-validation folds
metadata_df %>%
  dplyr::filter(LRTI_adjudication %in% c("Definite", "No Evidence")) ->
  cv_samples

# Generate folds for 5-fold cross-validation; ensure each fold has >=9 No Evidence samples
success <- FALSE
while (!success) {

  cv_samples %>%
    dplyr::select(sample_name, LRTI_adjudication) %>%
    dplyr::mutate(fold=sample(rep(1:5, length.out=nrow(.)))) ->
    cv_folds

  cv_folds %>%
    dplyr::select(LRTI_adjudication, fold) %>%
    table() ->
    cv_fold_table

  print("Generated CV folds with following counts:")
  print(cv_fold_table)

  cv_fold_table %>%
    .["No Evidence", ] %>%
    min() %>%
    {. >= 9} ->
    success

  if (!success) {
    print("At least one fold has too few No Evidence samples. Regenerating folds...")
  }
}

# Add folds to sample metadata and output the table
metadata_df %>% 
  dplyr::left_join(cv_folds %>% select(sample_name, fold), by = "sample_name") %>% 
  write.csv(classifier_meta_csv, row.names=F)