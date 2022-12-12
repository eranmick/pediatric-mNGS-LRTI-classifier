library(dplyr)
library(magrittr)
library(DESeq2)
source("constants.R")

# Read classifier metadata table, limiting to Definite/No Evidence samples that are used for CV
classifier_meta <- read.csv(classifier_meta_csv, stringsAsFactors = F) %>% 
  filter(LRTI_adjudication %in% c("Definite", "No Evidence"))

# Read host counts, omitting the first column (gene symbol)
host_counts <- read.csv(host_counts_csv, stringsAsFactors = F, row.names = 1)[,-1]
# Filter host counts to the samples used for CV
host_counts[, classifier_meta$sample_name] ->
  cv_host_counts

# Filter for genes with >10 counts in >20% of samples
keep <- rowSums(cv_host_counts > 10) > 0.2*ncol(cv_host_counts)
cv_host_counts <- cv_host_counts[keep, ]

# Apply variance stabilizing transform
dds <- DESeq2::DESeqDataSetFromMatrix(countData = cv_host_counts,
                                      colData = classifier_meta,
                                      design = ~1)
vsd <- DESeq2::varianceStabilizingTransformation(dds)
# Round transformed values to 2 decimal digits
vsd_mat <- SummarizedExperiment::assay(vsd) %>% 
  round(., digits=2)

# Save VST output
write.csv(vsd_mat, cv_vst_csv)