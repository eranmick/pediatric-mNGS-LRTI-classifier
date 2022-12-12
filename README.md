# pediatric-mNGS-LRTI-classifier
This repo houses data and code for the analyses in the pre-print ["Leveraging the pulmonary immune response and microbiome for improved lower respiratory tract infection diagnosis in critically ill children"](https://doi.org/10.1101/2022.12.01.22282994), focused on development and validation of a classifier for lower respiratory tract infection (LRTI) in critically ill children using host and microbial features from metagenomic next generation sequencing (mNGS) of tracheal aspirate RNA.

## How to get started
1. Download the repo to your local system.
2. Change the `project_root` variable in the script `constants.R` to the path of the main repo folder on your local system.
3. Set the R working directory to the `code` folder in the repo.
4. Run the numbered scripts in the `code` folder in order.

## Raw data

1. `sample_metadata.csv` details all patient (n=261) and water (n=32) samples used in this study including (where applicable) their age, sex, clinical LRTI adjudication, sequencing batch and total nonhost counts from the [CZ-ID](http://czid.org) metagenomic analysis pipeline.
2. `host_gene_counts.csv` includes the raw (unfiltered, unnormalized) host (human) gene counts of all patient samples. The counts were generated by pseudo-alignment with kallisto followed by summarization to the gene level with tximport. Genes are identified by their ENSEMBL ID. See the publication Methods section for full details. These gene counts are also deposited in the NCBI Gene Expression Omnibus (GEO) database under accession [GSE212532](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212532).
3. `microbe_reports.csv` is the nonhost taxa counts table of all the patient and water samples, as generated by the CZ-ID pipeline. 
4. `known_respiratory_pathogens.csv` is a curated list of organisms with known potential to cause LRTI, originally assembled for a [previous publication](https://www.pnas.org/doi/10.1073/pnas.1809700115).

## Scripts

1. `01-generate_cv_folds.R` - Randomly splits the samples with "Definite" or "No Evidence" LRTI status into 5 folds for cross-validation.
2. `02-filter_transform_host_counts_for_cv.R` - Filters the host counts for the cross-validation samples and applies the variance stablizing transformation as implemented in DESeq2.
3. `03-host_lasso_cv.R` - Selects features (genes) for use in a host-based LRTI classifier for each train/test split, as well as for all the Definite/No Evidence samples, using lasso logistic regression.
4. `04-host_lassoRF_cv.R` - Trains a random forest model on the training samples and with the selected host genes of each train/test split, and then generates out-of-fold LRTI probabilities.
5. `05-microbe_background_filtering.R` - Generates background filtering statistics on the microbial taxa using a negative binomial model trained on the water samples, as implemented in `idseqr.R`. [idseqr](https://github.com/czbiohub/idseqr) is a package for working with output from the CZ-ID pipeline in R that is currently in alpha stage.
6. `06-microbe_viral_hits_and_RBM.R` - Outputs i) likely pathogenic viral taxa present in patient samples after background filtering, and ii) microbial/fungal taxa identified by a rules-based model (RBM) as potential pathogens.
7. `07-calc_bac_vir_scores.R` - Calculates summary viral and bacterial/fungal scores for all patient samples.
8. `08-integrated_classifier_cv.R` - Trains a logistic regression model that integrates the i) host LRTI probability, ii) viral score, and iii) bacterial score of the Definite/No Evidence samples in the context of cross-validation, and then generates out-of-fold LRTI probabilities.
9. `09-output_cv_auc_roc` - Generates per-fold AUC tables and ROC curves from the out-of-fold probabilities of the Definite/No Evidence samples based on i) the host random forest classifier and ii) the integrated host+microbial classifier.
10. `10-host_predict_on_sus_ind` - Trains a host random forest model on all the Definite/No Evidence samples, using the genes selected using all those samples, and uses it to classify the "Suspected" and "Indeterminate" LRTI samples.
11. `11-integrated_predict_on_sus_ind` - Trains a integrated logtistic regression model on all the Definite/No Evidence samples using the host probabilities and microbial scores, and uses it to classify the Suspected/Indeterminate samples.
