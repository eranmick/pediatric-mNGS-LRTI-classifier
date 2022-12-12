library(dplyr)
library(magrittr)
library(randomForest)
source("constants.R")

seed <- 326582
set.seed(seed)

# Read classifier metadata table, limiting to Definite/No Evidence samples that are used for CV
classifier_meta <- read.csv(classifier_meta_csv, stringsAsFactors = F) %>% 
  filter(LRTI_adjudication %in% c("Definite", "No Evidence"))

## Read CV host classifier results
# Out-of-fold RF probabilities for test set
lasso_rf_preds <- read.csv(paste0(lassoRF_results_prefix, "preds.csv"),
                           stringsAsFactors = F)
# Per-fold out-of-bag RF probabilities within training set
lasso_rf_votes <- read.csv(paste0(lassoRF_results_prefix, "votes.csv"),
                           stringsAsFactors = F)

## Read microbial scores, limit to Definite/No Evidence, and split into viral and bf
microbe_scores <- read.csv(microbe_scores_csv, stringsAsFactors = F) %>% 
  filter(sample_name %in% classifier_meta$sample_name)

microbe_scores %>%
  select(sample_name, viral_score) ->
  viral_score_df

microbe_scores %>%
  select(sample_name, bf_score) ->
  bf_score_df

## Function to fit the integrated host+microbe model for a single train-test split
integ_host_microbe_fold <- function(curr_test_fold) {
        # host scores for current split.
        rbind(
                # out-of-fold probabilities on the test set
                lasso_rf_preds %>%
                        dplyr::inner_join(classifier_meta) %>%
                        dplyr::filter(fold==curr_test_fold) %>%
                        dplyr::select(sample_name, pred),
                # out-of-bag probabilities on the train set
                lasso_rf_votes %>%
                        dplyr::filter(test_fold==curr_test_fold) %>%
                        dplyr::select(sample_name, pred)
        ) %>%
                dplyr::rename(host_prob=pred) %>%
                dplyr::mutate(
                        # regularize so we can pass thru logit function
                        host_prob=host_prob*(n_rf_trees)/(n_rf_trees+1) +
                                1 / (n_rf_trees+1),
                        # numerically stable way to compute 1 - regularized score
                        one_minus_prob=(1-host_prob)*(n_rf_trees)/(n_rf_trees+1) +
                                1 / (n_rf_trees+1)
                ) %>%
                # apply logit transformation
                dplyr::mutate(host_score = log(host_prob / (one_minus_prob))) %>%
                dplyr::select(sample_name, host_score) ->
                curr_host_scores
        
        # Compute bacterial scores for current split. add epsilon to scores
        # so we can log them. epsilon=smallest nonzero score in the train
        # set, divided by 10
        
        curr_bf_scores <- bf_score_df
        
        # compute epsilon=smallest nonzero score in train set, divided by 10
        curr_bf_scores %>%
                dplyr::inner_join(classifier_meta) %>%
                dplyr::filter(fold != curr_test_fold, bf_score > 0) %>%
                with(min(bf_score) / 10) ->
                bf_epsilon
        
        # log-transform the bacterial score
        curr_bf_scores$bf_score <- log10(curr_bf_scores$bf_score + bf_epsilon)
        
        
        # calculated epsilon for viral scores in a similar manner
        viral_score_df %>%
                dplyr::inner_join(classifier_meta) %>%
                dplyr::filter(fold != curr_test_fold, viral_score > 0) %>%
                with(min(viral_score) / 10) ->
                viral_epsilon
        
        # log-transform the viral score
        viral_score_df %>%
                dplyr::mutate(viral_score=log10(viral_score + viral_epsilon)) ->
                curr_viral_scores
        
        # combined dataframe with host, bacterial, viral scores
        classifier_meta %>%
                dplyr::inner_join(curr_host_scores) %>%
                dplyr::inner_join(curr_bf_scores) %>%
                dplyr::inner_join(curr_viral_scores) %>%
                dplyr::mutate(LRTI_adjudication=factor(
                        LRTI_adjudication, levels=c("No Evidence", "Definite"))) ->
                curr_df
        
        # fit logistic regression model
        logistic_mod <- glm(LRTI_adjudication ~ host_score + bf_score + viral_score,
                            family="binomial", data=dplyr::filter(
                                    curr_df, fold != curr_test_fold))
        
        curr_df %>%
                # filter to test set
                dplyr::filter(fold == curr_test_fold) %>%
                # add logistic prediction
                dplyr::mutate(logit_pred=predict(logistic_mod,
                                                 newdata=., type="response")) %>%
                # return the important variables
                dplyr::select(sample_name, logit_pred)
}

## Run on all 5 train-test splits
lapply(1:max(classifier_meta$fold), integ_host_microbe_fold) %>%
  do.call(what=rbind) %>%
  left_join(., y = classifier_meta %>% select(sample_name, LRTI_adjudication), by = "sample_name") %>% 
  dplyr::relocate(LRTI_adjudication, .after=sample_name) %>% 
  write.csv(host_microbe_cv_csv, row.names=F)