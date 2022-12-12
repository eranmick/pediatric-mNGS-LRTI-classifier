## Constants

### CHANGE THIS TO POINT TO THE LOCAL REPO PATH ON YOUR SYSTEM
project_root <- "C:/Users/eranm/Documents/Box Sync/pediatric-mNGS-LRTI-classifier/"

# "raw data" folder
raw_data_dir <- file.path(project_root, "data", "raw")

# "processed data" folder
processed_data_dir <- file.path(project_root, "data", "processed")

# Raw sample metadata
metadata_csv <- file.path(raw_data_dir, "sample_metadata.csv")

# Classifier metadata (with cross-validation folds)
classifier_meta_csv <- file.path(processed_data_dir, "classifier_metadata.csv")

# Raw host gene counts
host_counts_csv <- file.path(raw_data_dir, "host_gene_counts.csv")

# Host gene VST values for cross-validation
cv_vst_csv <- file.path(processed_data_dir, "cv_vst.csv")

# "results" folder
results_dir <- file.path(project_root, "results")

# Filename prefix for logistic-lasso results
lasso_results_prefix <- file.path(results_dir, "lasso_")

# Number of trees to use for random forest
n_rf_trees <- 10000

# Filename prefix for lasso-RF results
lassoRF_results_prefix <- file.path(results_dir, "lasso_RF_")

# Path to raw microbe reports
microbe_reports_csv <- file.path(raw_data_dir, "microbe_reports.csv")

# Minimum number of reads for taxa to be retained
min_reads <- 5

# Path to microbe reports with background filtering stats
microbe_reports_bgfilter_csv <- file.path(processed_data_dir, "microbe_reports_bgfilter.csv")

# Path to list of known respiratory pathogens
known_pathogens_csv <- file.path(raw_data_dir, "known_respiratory_pathogens.csv")

# Path to viral results
viral_hits_csv <- file.path(results_dir, "viral_hits.csv")

# Path to RBM results (bacteria/fungi)
rbm_csv <- file.path(results_dir, "all_rbm_results.csv")

# Path to viral and bacterial/fungal scores
microbe_scores_csv <- file.path(results_dir, "microbial_scores.csv")

# Path to integrated host/microbe cross-validation results
host_microbe_cv_csv <- file.path(results_dir, "host_microbe_cv_preds.csv")

# Filename prefixes for AUC tables and ROC plots of cross-validation results
cv_auc_prefix <- file.path(results_dir, "auc_")
cv_roc_prefix <- file.path(results_dir, "roc_")

# Path to random forest votes for the Definite/No Evidence samples 
# when they are used with the final host classifier as the training 
# set for classifying Suspected/Indeterminate samples
host_final_def_noevi_rf_votes_csv <- file.path(results_dir, "host_final_def_noevi_rf_votes.csv")

# Path to host random forest classifier predictions on Suspected/Indeterminate samples
host_sus_ind_rf_preds_csv <- file.path(results_dir, "host_sus_ind_rf_preds.csv")

# Path to integrated host+microbe classifier predictions on Suspected/Indeterminate samples
integ_sus_ind_preds_csv <- file.path(results_dir, "integ_sus_ind_preds.csv")