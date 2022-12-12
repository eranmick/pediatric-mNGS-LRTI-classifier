library(dplyr)
library(magrittr)
source("constants.R")

# Read classifier metadata table, limiting to Definite/No Evidence samples that are used for CV
classifier_meta <- read.csv(classifier_meta_csv, stringsAsFactors = F) %>% 
  filter(LRTI_adjudication %in% c("Definite", "No Evidence"))
# Read host genes VST values
vsd_mat <- read.csv(cv_vst_csv, row.names=1)

# Vector mapping ensembl gene IDs to gene symbols
read.csv(host_counts_csv, stringsAsFactors = F) %>%
  {`names<-`(.$gene_symbol, .$X)} ->
  ensg2gene

# Regression variables
X <- t(vsd_mat)
y <- classifier_meta$LRTI_adjudication == 'Definite'

# Function to extract coefficients from fitted lasso model
lasso_coef_df <- function(mod) {
  # Use the most regularized value of the tuning parameter that is
  # within 1 standard error of the optimum
  coef(mod, s='lambda.1se', gamma=c("gamma.1se"))[,1] %>%
    # Filter to nonzero coefficients
    .[. != 0] %>%
    {data.frame(coef=.)} %>%
    # Convert rowname to column, then grab the gene symbol for the
    # ensembl gene name
    tibble::rownames_to_column("gene") %>%
    dplyr::mutate(gene_symbol=ensg2gene[gene])
}

# Function to run a single train-test split. Fits a lasso model on the
# training set, outputs predictions on the test set and the selected
# coefficients for this split.
lasso_cv_outer_fold <- function(test_fold, ...) {
  # Boolean of the training samples
  train <- classifier_meta$fold != test_fold
  # For the nested folds within the training set, make them range from
  # 1 to 4 so glmnet is happy
  foldid <- classifier_meta$fold[train]
  foldid <- as.integer(as.factor(foldid))
  # Fit lasso logistic regression on the training set, with a nested
  # cross-validation within the folds of the training set
  mod <- glmnet::cv.glmnet(X[train,], y[train], family='binomial',
                           foldid=foldid, ...)

  # Generate predictions on the test set
  predict(mod, X[!train,], type='response',
          s='lambda.1se', gamma=c("gamma.1se"))[,1] %>%
    {data.frame(pred=.)} %>%
    tibble::rownames_to_column("sample_name") ->
    pred

  # Get the nonzero coefficients, also storing which train-test split
  # they were from
  lasso_coef_df(mod) %>%
    dplyr::mutate(fold=test_fold) ->
    nonzero_coefs

  list(pred=pred, mod=mod, coef=nonzero_coefs)
}

# Function to run the model on all train-test splits and output results to the out_prefix path. 
# The ellipsis... are extra arguments passed to glmnet
lasso_cv <- function(out_prefix, ...) {
  cv_list <- lapply(1:max(classifier_meta$fold),
                    function(x) lasso_cv_outer_fold(x, ...))

  lapply(cv_list, function(x) x$pred) %>%
    do.call(what=rbind) %>%
    write.csv(paste0(out_prefix, "preds.csv"), row.names=F)

  lapply(cv_list, function(x) x$coef) %>%
    do.call(what=rbind) %>%
    write.csv(paste0(out_prefix, "fold_coefs.csv"), row.names=F)

  glmnet::cv.glmnet(X, y, family='binomial',
                    foldid=classifier_meta$fold, ...) %>%
    lasso_coef_df() ->
    coefs

  coefs %>%
    write.csv(paste0(out_prefix, "coefs.csv"), row.names=F)
}

# Run the lasso
lasso_cv(lasso_results_prefix)