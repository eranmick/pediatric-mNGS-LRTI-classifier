library(dplyr)
library(magrittr)
library(pROC)
source("constants.R")

# Read classifier metadata table, limiting to Definite/No Evidence samples used for CV
classifier_meta <- read.csv(classifier_meta_csv, stringsAsFactors = F) %>% 
  filter(LRTI_adjudication %in% c("Definite", "No Evidence"))

# Function to output the per-fold AUC table and ROC curve from a 
# classifier's out-of-fold LRTI probabilities. The filenames 
# will end with output_suffix to distinguish different classifier versions.
output_auc_roc <- function(preds_df, output_suffix){
  
  # Get dataframe of per-fold AUC values
  preds_df %>%
    dplyr::inner_join(classifier_meta %>% select(sample_name, fold), by = "sample_name") %>%
    dplyr::group_by(fold) %>%
    dplyr::summarize(auc=as.numeric(pROC::roc(LRTI_adjudication=='Definite', pred)$auc)) %>%
    dplyr::ungroup() ->
    auc_df
  
  # Output it
  write.csv(auc_df, paste0(cv_auc_prefix, output_suffix, ".csv"), row.names = F)

  # Generate ROC plot
  for (i in 1:max(auc_df$fold)) {
    preds_df %>%
      dplyr::inner_join(classifier_meta) %>%
      dplyr::filter(fold==i) %>%
      with(pROC::roc(LRTI_adjudication=='Definite', pred)) %>%
      plot(add=(i != 1), col=i)
  }
  
  # Output it
  dev.copy(pdf, paste0(cv_roc_prefix, output_suffix, ".pdf"))
  dev.off()
}

# Read out-of-fold probabilities of the lasso-RF host classifier
read.csv(paste0(lassoRF_results_prefix, "preds.csv"), stringsAsFactors = F) ->
  host_preds_df
# Output AUC and ROC
output_auc_roc(host_preds_df, "host_lassoRF")

# Read out-of-fold probabilities of the integrated host+microbe classifier
read.csv(host_microbe_cv_csv, stringsAsFactors = F) ->
  integ_preds_df
# Output AUC and ROC
output_auc_roc(integ_preds_df, "integrated")