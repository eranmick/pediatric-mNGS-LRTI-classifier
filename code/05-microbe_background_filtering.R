library(dplyr)
library(magrittr)
source("constants.R")
source("idseqr.R")

# Read all sample metadata
metadata_df <- read.csv(metadata_csv, stringsAsFactors = F)
# Read microbe reports and apriori exclude taxon hits with fewer than 5 NT counts
microbe_reports <- read.csv(microbe_reports_csv, stringsAsFactors = F) %>% 
  filter(nt_count >= min_reads)

# Calculate background filtering stats using the filter_background function from idseqr
controls <- metadata_df %>% filter(is_water) %>% pull(sample_name)
batches <- metadata_df %>% pull(batch)
names(batches) <- metadata_df %>% pull(sample_name)

reports_filterbg <- filter_background(microbe_reports, 
                                      controls=controls,
                                      batches=batches,
                                      normalization=NULL) %>% # the default normalization is based on total non-host counts
        filter(category %in% c("viruses", "bacteria", "eukaryota")) %>% 
        group_by(sample_name) %>%
        mutate("p_adj" = p.adjust(p_val, method = "BH")) %>%  # generate BH adjusted p-values within each sample
        ungroup()

# Output microbe reports with background filtering stats for patient samples
reports_filterbg %>% 
  filter(sample_name %in% (metadata_df %>% filter(!is_water) %>% pull(sample_name))) %>% 
  write.csv(microbe_reports_bgfilter_csv, row.names = F)