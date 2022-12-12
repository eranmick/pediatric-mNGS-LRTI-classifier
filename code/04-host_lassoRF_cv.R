library(dplyr)
library(magrittr)
library(randomForest)
source("constants.R")

seed <- 241493047
set.seed(seed)

# Read classifier metadata table, limiting to Definite/No Evidence samples that are used for CV
classifier_meta <- read.csv(classifier_meta_csv, stringsAsFactors = F) %>% 
  filter(LRTI_adjudication %in% c("Definite", "No Evidence"))
# Read host genes VST values
vsd_mat <- read.csv(cv_vst_csv, row.names=1)
# Output variable to predict
y <- as.factor(classifier_meta$LRTI_adjudication)
# Table of features selected by lasso on each train-test split
lasso_features <- read.csv(paste0(lasso_results_prefix, "fold_coefs.csv"),
                           stringsAsFactors = F)

# Function to run random forest on lasso-selected genes for a single
# train-test split
lassoRF_cv_outer_fold <- function(test_fold) {
  # Extract the lasso features for this split
  lasso_features %>%
    dplyr::filter(fold==test_fold, gene !='(Intercept)') %>%
    .$gene ->
    keep

  # Subset to selected features
  X <- t(vsd_mat)[,keep]
  # Boolean for the training samples
  train <- classifier_meta$fold != test_fold
  # Fit the random forest
  rf <- randomForest(X[train,], y[train], ntree=n_rf_trees)
  # Return random forest and predictions on the test set
  list(test_fold=test_fold, mod=rf,
       pred=predict(rf, newdata=X[!train,], type='prob')[,"Definite"])
}

# Run on all the train-test splits
cv_list <- lapply(1:max(classifier_meta$fold),
                  function(i) lassoRF_cv_outer_fold(i))

# Save the out-of-fold predictions
cv_list %>%
  lapply(function(x) x$pred) %>%
  do.call(what=c) %>%
  {data.frame(pred=.)} %>%
  tibble::rownames_to_column("sample_name") %>%
  left_join(., y=classifier_meta %>% select(sample_name, LRTI_adjudication), by="sample_name") %>% 
  dplyr::relocate(LRTI_adjudication, .after=sample_name) %>% 
  write.csv(paste0(lassoRF_results_prefix, "preds.csv"), row.names=F)

# Save the training set's out-of-bag probabilities for each train-test
# split. We'll use these for training the integrated host+microbe classifier 
# in the context of cross-validation to avoid needing a third level of nesting
cv_list %>%
  lapply(function(x) x$mod$votes %>%
                       as.data.frame() %>%
                       tibble::rownames_to_column("sample_name") %>%
                       dplyr::mutate(test_fold=x$test_fold) %>%
                       dplyr::rename(pred=Definite) %>%
                       dplyr::select(test_fold, sample_name, pred)) %>%
  do.call(what=rbind) %>%
  left_join(., y=classifier_meta %>% select(sample_name, LRTI_adjudication), by="sample_name") %>% 
  dplyr::relocate(LRTI_adjudication, .after=sample_name) %>% 
  write.csv(paste0(lassoRF_results_prefix, "votes.csv"), row.names=F)