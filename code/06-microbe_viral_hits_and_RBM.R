library(dplyr)
library(magrittr)
source("constants.R")

# Read all patient sample metadata
metadata_df <- read.csv(metadata_csv, stringsAsFactors = F) %>% 
  filter(!is_water)
# Read respiratory pathogens list
resp_pathogens <- read.csv(known_pathogens_csv, stringsAsFactors = F)
# Read background filtered microbe reports
microbe_reports_bgf <- read.csv(microbe_reports_bgfilter_csv, stringsAsFactors = F)

### Viral hits
resp_viruses <- resp_pathogens %>% 
  filter(grepl("virus", Pathogen))

# Limit to pathogenic viral hits that passed background filtering
virus_results <- microbe_reports_bgf %>%
  filter(category == "viruses",
         p_adj < 0.05,
         name %in% resp_viruses$Pathogen)

# Write out viral hits
virus_results %>%
  left_join(., y = metadata_df %>% select(sample_name, LRTI_adjudication), by = "sample_name") %>% 
  select(sample_name, LRTI_adjudication, name, tax_id, genus_tax_id, nt_rpm) %>% 
  write.csv(viral_hits_csv, row.names = F)

### Bacteria/fungal hits using rules-based model (RBM)
resp_bf <- resp_pathogens %>% 
  filter(!grepl("virus", Pathogen))

# Limit to bacterial/fungal hits that passed background filtering
bf_reports <- microbe_reports_bgf %>% 
  filter(category == "bacteria" | category == "eukaryota") %>% 
  filter(p_adj < 0.05)

# Function to calculate differences in NT rpM values between taxa in a sample
get_difs <- function(l){
        difs <- list()
        for(i in seq_along(l[1:length(l)-1])){
                difs[[i]] <- l[i] - l[i+1]
        }
        return(unlist(difs))
}

# Function to extract all taxa above the max NT rpM drop-off in a sample
apply_rbm <- function(df){
        
        nt_rpm_values <- df %>%
          select(sample_name, name, tax_id, genus_tax_id, nt_count, nt_rpm, genus_sum_nt_rpm) %>% 
          arrange(desc(nt_rpm))
        
        if(nrow(nt_rpm_values) > 15){
          nt_rpm_values <- nt_rpm_values[seq(1,15),]             # restrict analysis to the top 15 organisms   
        }
        
        x <- get_difs(nt_rpm_values$nt_rpm)       # identify the index of the maximum drop-off in rpM values
        if(is.null(x)){                           # this happens when there is only a single taxon in nt_rpm_values, so differences can't be calculated 
                top_orgs <- nt_rpm_values         # in which case, we return that single taxon
        }else{
                top_orgs <- nt_rpm_values[seq(1, which.max(x)), ]       # o/w, we return all the taxa until the maximum drop-off index
        }
        
        return(top_orgs)
}

# Apply the RBM to all patient samples
samplenames <- metadata_df %>% pull(sample_name)

all_rbm_results <- data.frame()   # data frame that will be populated with the results

# loop through all the samples and apply the RBM
for(sn in samplenames){
  
        bf_reports %>% 
          dplyr::filter(sample_name == sn,
                        genus_tax_id > 0) -> # require a genus taxid (removes "uncultured bacteria")
          report
        
        report %>% 
          group_by(genus_tax_id) %>% 
          mutate(genus_sum_nt_rpm = sum(nt_rpm)) %>% # calculate genus level nt_rpm
          ungroup() ->
          report
        
        report %>%
          group_by(genus_tax_id) %>%          # group by genus
          top_n(1, nt_rpm) %>%                # select the top species within the genus by nt_rpm
          ungroup() ->
          top_species_in_genus_report

        report %>%
          filter(name %in% resp_bf$Pathogen) %>% # filter just for respiratory pathogens
          group_by(genus_tax_id) %>%                    # group by genus
          top_n(1, nt_rpm) %>%                          # select the top (pathogenic) species within the genus by nt_rpm
          ungroup() ->
          top_pathogenic_species_in_genus_report

        top_species_in_genus_report %>%                 # combine selections of top species in genus and top pathogenic species in genus
          bind_rows(top_pathogenic_species_in_genus_report) %>%
          distinct(.keep_all = T) %>% # this prevents duplicated entries for when the top species in the genus was already a pathogen
          as.data.frame() ->
          filtered_report
        
        if(nrow(filtered_report) > 0) {
          result <- apply_rbm(filtered_report)
          all_rbm_results <- rbind(all_rbm_results, result)
        }
        
}

# Write out RBM results
all_rbm_results %>%
  mutate("pathogen" = name %in% resp_bf$Pathogen) %>% 
  dplyr::relocate(pathogen, .after = name) %>% 
  left_join(., y = metadata_df %>% select(sample_name, LRTI_adjudication), by = "sample_name") %>% 
  dplyr::relocate(LRTI_adjudication, .after = sample_name) %>% 
  write.csv(rbm_csv, row.names = F)