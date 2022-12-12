library(dplyr)
library(magrittr)
source("constants.R")

## Read all patient sample metadata
metadata_df <- read.csv(metadata_csv, stringsAsFactors = F) %>% 
  filter(!is_water)

## Viral scores

# Read viral hits
virus_results <- read.csv(viral_hits_csv, stringsAsFactors = F)

# Calculate viral score = sum(nt_rpm) of all pathogenic viruses present in a sample after background filtering
virus_results %>%
  dplyr::group_by(sample_name) %>%
  dplyr::summarize(viral_score = sum(nt_rpm)) %>%
  # add samples without virus
  dplyr::right_join(metadata_df %>% select(sample_name, LRTI_adjudication), by = "sample_name") %>%
  # add 0s for the scores of samples without virus
  tidyr::replace_na(list(viral_score=0)) %>%
  dplyr::select(sample_name, LRTI_adjudication, viral_score) ->
  viral_score

## Bacterial/fungal (bf) scores

# Read RBM results, limiting to respiratory pathogens and adding the total nonhost read counts from the sample metadata 
rbm_results <- read.csv(rbm_csv, stringsAsFactors = F) %>%
  filter(pathogen) %>% 
  dplyr::left_join(., y = metadata_df %>% select(sample_name, nonhost_reads), by = "sample_name")
  
# Calculate bf score = proportion of all potentially pathogenic RBM hits out of the total nonhost counts
rbm_results %>%
  dplyr::group_by(sample_name) %>%
  dplyr::summarize(bf_score=sum(nt_count / nonhost_reads)) %>%
  # add back the samples without bacterial/fungal hits and give them a score of 0
  dplyr::right_join(metadata_df %>% select(sample_name)) %>%
  tidyr::replace_na(list(bf_score=0)) ->
  bf_score

# Write out the combined microbial scores
viral_score %>% 
  left_join(., y = bf_score, by = "sample_name") %>% 
  write.csv(microbe_scores_csv, row.names = F)