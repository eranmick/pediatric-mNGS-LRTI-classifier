library(dplyr)
library(magrittr)
library(hablar)
source("constants.R")

seed <- 241493047
set.seed(seed)

# Read microbial scores of all samples
microbe_scores <- read.csv(microbe_scores_csv, stringsAsFactors = F)

# Read host out-of-bag probabilities of training samples (Def/NoEvi)
train_rf_votes <- read.csv(host_final_def_noevi_rf_votes_csv, stringsAsFactors = F) %>% 
  dplyr::rename("rf_pred"=rf_votes_pred)

# Read host probabilities of test samples (Sus/Ind)
pred_rf_prob <- read.csv(host_sus_ind_rf_preds_csv, stringsAsFactors = F) %>% 
  select(-LRTI_adjudication)

# Bind all host probs
train_rf_votes %>% 
  bind_rows(pred_rf_prob) ->
  host_prob

# Combine all features
microbe_scores %>% 
  left_join(host_prob, by = "sample_name") ->
  all_features

## Apply feature transformations

# Calculate minimal non-zero viral score in training set
all_features %>% 
  filter(LRTI_adjudication %in% c("Definite", "No Evidence"),
         viral_score > 0) %>%
  pull(viral_score) %>% 
  min -> 
  min_train_viral_score

# Calculate minimal non-zero bf score in training set
all_features %>% 
  filter(LRTI_adjudication %in% c("Definite", "No Evidence"),
         bf_score > 0) %>%
  pull(bf_score) %>% 
  min -> 
  min_train_bf_score

# log10 transform viral and bf scores after adding the respective minimal non-zero score divided by 10;
# logit transform host probabilities after regularization
all_features %>% 
  mutate(trans_viral_score = log10(viral_score + min_train_viral_score/10),
         trans_bf_score = log10(bf_score + min_train_bf_score/10),
         reg_host_prob = rf_pred*(n_rf_trees)/(n_rf_trees+1) + 1/(n_rf_trees+1),
         one_minus_reg_prob=(1-reg_host_prob)*(n_rf_trees)/(n_rf_trees+1) + 1/(n_rf_trees+1),
         trans_host_score = log(reg_host_prob / one_minus_reg_prob)) ->
  all_features

# Fit integrated host+microbe logistic regression model on the training set
logistic_mod <- glm(LRTI_adjudication ~ trans_host_score + trans_bf_score + trans_viral_score,
                    family="binomial", 
                    data=all_features %>% filter(LRTI_adjudication %in% c("Definite", "No Evidence")) %>%
                      mutate("LRTI_adjudication" = factor(LRTI_adjudication, levels=c("No Evidence", "Definite"))))

# Predict using logistic regression model and output all features and predictions
all_features %>% 
  filter(LRTI_adjudication %in% c("Suspected", "Indeterminate")) %>% 
  mutate("integrated_prob" = predict(logistic_mod, newdata = ., type="response")) %>% 
  dplyr::rename("raw_host_prob" = rf_pred) %>% 
  write.csv(integ_sus_ind_preds_csv, row.names = F)