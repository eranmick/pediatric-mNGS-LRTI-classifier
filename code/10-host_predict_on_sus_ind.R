library(dplyr)
library(magrittr)
library(tibble)
library(DESeq2)
library(randomForest)
source("constants.R")

seed <- 241493047
set.seed(seed)

# Read classifier metadata table, limiting to patient samples
classifier_meta <- read.csv(classifier_meta_csv, stringsAsFactors = F) %>% 
  filter(!is_water)

# Split into training (Definite/No Evidence) and prediction (Suspected/Indeterminate) sets
classifier_meta %>% 
  filter(LRTI_adjudication %in% c("Definite", "No Evidence")) ->
  training_samples

classifier_meta %>% 
  filter(LRTI_adjudication %in% c("Suspected", "Indeterminate")) ->
  prediction_samples

# Read host counts (remove gene symbol column)
host_counts <- read.csv(host_counts_csv, row.names = 1)[,-1]

# Subset counts to the the training samples and to genes that have >10 counts in >20% of samples
# (this is the same as what was used in the host classifier CV)
train_counts <- host_counts[, training_samples$sample_name]
keep <- rowSums(train_counts > 10) > 0.2*ncol(train_counts)
train_counts <- train_counts[keep, ]

# Create DESeq2 object with the training samples only
dds_train <- DESeq2::DESeqDataSetFromMatrix(countData = train_counts,
                                            colData = training_samples,
                                            design = ~1)
dds_train <- estimateSizeFactors(dds_train)
dds_train <- estimateDispersions(dds_train)

# Apply VST to training sample counts only 
# (this is the same as what was used in the host classifier CV)
vsd_train <- varianceStabilizingTransformation(dds_train) %>% 
  assay %>% 
  round(., digits=2)

# These are the counts of all the samples (training+prediction), 
# limited to the genes present in the training data by the above criteria
predict_counts <- host_counts[keep, classifier_meta$sample_name]

# Create DESeq2 object
dds_predict <- DESeq2::DESeqDataSetFromMatrix(countData = predict_counts,
                                              colData = classifier_meta,
                                              design = ~1)
dds_predict <- estimateSizeFactors(dds_predict)
dispersionFunction(dds_predict) <- dispersionFunction(dds_train) # assign the dispersion function from the training data alone

# Apply VST to training+prediction sample counts using the training-only dispersion function
vsd_predict <- varianceStabilizingTransformation(dds_predict) %>% 
  assay %>% 
  round(., digits=2)

# Table of host genes previously selected by lasso on the full training set (exclude intercept)
lasso_features <- read.csv(paste0(lasso_results_prefix, "coefs.csv"), stringsAsFactors = F)[-1,]

# Fit the random forest on the training samples using only the selected genes
vsd_predict <- t(vsd_predict[lasso_features$gene, ])
rf <- randomForest(vsd_predict[training_samples$sample_name,], as.factor(training_samples$LRTI_adjudication),
                   ntree=n_rf_trees)

# Output rf "votes" of the training samples for use as input features in the integrated classifier 
rf$votes %>%
  as.data.frame() %>% 
  tibble::rownames_to_column("sample_name") %>% 
  select(-`No Evidence`) %>% 
  dplyr::rename("rf_votes_pred"=Definite) %>% 
  write.csv(host_final_def_noevi_rf_votes_csv, row.names = F)

# Return rf predictions on the test samples
pred <- predict(rf, newdata=vsd_predict[prediction_samples$sample_name,], type='prob')[,"Definite"]

prediction_samples %>% 
  mutate("rf_pred"=pred) %>% 
  select(sample_name, LRTI_adjudication, rf_pred) %>% 
  write.csv(host_sus_ind_rf_preds_csv, row.names = F)