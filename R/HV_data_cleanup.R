# This is the script used to cleanup/tidy the data sheets ready for stats

## NOTES ON DATA ##
# MISSING GRASS = MISSING/UNABLE TO LOCATE DUE TO DENSE GRASS/WATER. 
#``               MIGHT BE ALIVE AND HIDING
# MISSING DEAD = MISSING
#               HAD A GOOD LOOK AND SLIM CHANCE OF HIDING SO PRESUMED DEAD. 
#               ALSO USED ON PLANTS THAT WERE FOUND DEAD DURING BASELINE SURVEY
# SURVIVAL STATUS BASELINE = DEAD OR SURVIVAL STATUS RESURVEY
#               DEAD MEANS ACTUALLY FOUND DEAD PLANT

# Clean things up
rm(list=ls())

# Install the necessary packages
#install.packages("tidyverse")

# Load the necessary libraries..
library(tidyverse)

# Read in the data and only include plots that have been resurveyed 
survey_data <- read_csv("Data/minimalist_survey_20220710.csv", 
                        col_types = cols(roots_pot_soil_visible_baseline = 
                                           col_character(),
                                         snapped_height_resurvey_cm = 
                                           col_double())) %>% 
  filter(!is.na(date_resurvey)) %>%
  #Filter out dead/missing/unplanted at initial survey and missing at resurvey
  filter(!plant_id_baseline %in% c("unknown","unplanted","missing") | !plant_id_resurvey 
         %in% c("missing")) %>% 
  filter(!is.na(plant_id_resurvey)) %>%
  #correct species ID from resurvey
  mutate(species = ifelse(plant_id_baseline != plant_id_resurvey & 
                            plant_id_resurvey != "missing", 
                          plant_id_resurvey, plant_id_baseline),
         # replacing ansu and anfl with generic ango
         species = ifelse(species == "ansu" | species == "anfl", "ango", 
                          species)) %>% 
  mutate(survival_status_resurvey = ifelse(survival_status_resurvey == "dead" | 
                                             survival_status_resurvey == "DEAD", 
                                           "dead", survival_status_resurvey )) %>%
  mutate(species = ifelse(species == "unknown", "missing", species)) 
# Add in the rough planting dates.
# Plots which are missing dates: 34, 37-42, 49, 52, 76
# dates could be assumed for 52 and 76 due to plot groupings
# Plots recorded as being planted over 2 different dates: 44, 51, 56
# Replace these all with NA. See orig file from John
planting_dates <- read_csv("Data/edited_greenfleet_planting_dates.csv") 

survey_data <- left_join(survey_data, planting_dates, by = "plot") %>%
  mutate(date_diff_b2r = as.numeric(date_resurvey - date_mapped_baseline)) %>% #days between surveys
  mutate(date_diff_p2b = as.numeric(date_mapped_baseline - date_planted)) %>% #days between planting and baseline
  mutate(date_diff_p2r = as.numeric(date_resurvey - date_planted)) #days between planting and baseline
write.csv(survey_data, "Outputs/HV_survey_data.csv", row.names=F)

# Add in the mixture and treatment information
dry_mixtures <- c("B", "D", "F", "M")
plot_treatments <- read_csv("Data/HV_randomised_treatments_as_planted.csv") %>% 
  filter(is.na (notes) | notes != "insufficient stock so plot removed") %>% 
  mutate(provenance = ifelse(mixture %in% dry_mixtures, "dry", "local")) %>% 
  select(-compost) %>% 
  arrange(plot)

species_mixtures <- read_csv("Data/Species_mixtures_20_10_20.csv") %>% 
  # converts anfl and ansu to ango
  mutate(species = ifelse(species=="ansu" | species=="anfl", "ango", species))

scientific_names <- species_mixtures %>%
  select(scientific_name, species) %>%
  mutate(species = ifelse(species=="eute (dry)", "eute_dry", 
                          species),
         species = ifelse(species=="eucr (dry)", "eucr_dry", 
                          species)) %>%
  mutate(scientific_name = ifelse(scientific_name == "Angophora floribunda", 
                                  "Angophora floribunda/subvelutina", scientific_name),
         scientific_name = ifelse(scientific_name == "Angophora subvelutina", 
                                  "Angophora floribunda/subvelutina", scientific_name)) %>%
  distinct(species, .keep_all=TRUE) 
  

# table of full planting information (species, heights, treatment etc)
planting_details <- left_join(survey_data, plot_treatments, by = "plot") %>%
  mutate(species = ifelse(provenance=="dry" & species=="eute", "eute_dry", 
                          species),
         species = ifelse(provenance=="dry" & species=="eucr", "eucr_dry", 
                          species)) %>%
  rowwise() %>%
  mutate(compost = ifelse(y %in% c(compost_min_row:compost_max_row), 
                          "YES", "NO")) %>%
  # Update unknown certainties from orig survey 
  mutate(id_certainty = 
           ifelse(mapper_baseline == "MMM" & 
                    # Species eute and eucr were most often confused for each 
                    # other and coin and cote were unpredictable mixed up with 
                    # each other and sometimes coci.
                    plant_id_baseline %in% c("eute", "eute_dry", "eucr", 
                                             "eucr_dry", "coin", "cote") & 
                    plant_id_resurvey == "missing" | 
                    plant_id_resurvey == "unknown" , 
                  # Blocks A:D were monocultures so can be confident in IDs
                  ifelse(mixture %in% c("A","B","C","D"), "certain", "uncertain"), "certain"))
# number of days between surveys  
write.csv(planting_details, "Outputs/HV_planting_details.csv", row.names=F)

# elevation data
elevation_data <- read_csv("Data/plot_corner_elevations_LiDAR.csv") 

plot_mean_elevation <- elevation_data%>%
  group_by(plot)%>%
  summarise(mean_elev = mean(elevation),
            cross_fall = max(elevation, na.rm=T) - min(elevation, na.rm=T),
            mean_easting = mean(easting),
            mean_northing = mean(northing),
            n.obs = n())

planting_details <- left_join(planting_details, plot_mean_elevation, 
                              by = "plot") 

avg_trait_vals <- read_csv("Outputs/avg_trait_vals.csv")

planting_details <- left_join(planting_details, avg_trait_vals, 
                              by = "species") 

std_function<-function(x) ((x - mean(x, na.rm=T)) / sd(x, na.rm=T))

nice_data <- planting_details %>%
  filter(id_certainty=="certain")%>%
  filter(block == "A" | block == "B")%>%
  ungroup()%>%
  mutate(std_mean_elev = std_function(mean_elev),
         std_cross_fall = std_function(cross_fall),
         std_sqrt_height = std_function(sqrt(height_baseline_cm)),
         std_date_diff = std_function(date_diff_b2r))


# new data set for survival model - FILTERED FOR IDs CERTAIN
survival_data <- nice_data %>%
  #  filter(survival_status_resurvey == "DEAD") %>%
  mutate(alive_resurvey = ifelse(is.na(height_resurvey), 0, 1)) %>%
  select(block, plot, mean_elev, cross_fall, mixture, innoculation, compost, 
         y, x, mapped_x, species, provenance, trait_pc1, trait_pc2, trait_pc3, 
         la, sla, ldmc, ssd, rtd, log_srl, log_seed_mass, height_baseline_cm, 
         alive_resurvey, date_planted, date_mapped_baseline, date_resurvey, 
         date_diff_b2r, date_diff_p2b, date_diff_p2r, 
         std_mean_elev, std_cross_fall, std_sqrt_height, std_date_diff)

# new data set for growth model. Note: using data set with certain IDs - if 
# plant found at resurvey, then id is certain
growth_data <- nice_data %>%
  # remove data where no height was recorded due to missing at resurvey
  filter(!is.na(height_resurvey) & species != "missing") %>%
  mutate(growth_mm_day = 10*(height_resurvey-height_baseline_cm)/date_diff_b2r) %>%
  select(block, plot, mean_elev, cross_fall, mixture, innoculation, compost, y, 
         x, mapped_x, species, provenance, height_baseline_cm, height_resurvey, 
         growth_mm_day, trait_pc1, trait_pc2, trait_pc3, la, sla, ldmc, ssd, 
         rtd, log_srl, log_seed_mass, date_planted, date_mapped_baseline, 
         date_resurvey, date_diff_b2r, date_diff_p2b, date_diff_p2r, 
         std_mean_elev, std_cross_fall, std_sqrt_height, std_date_diff)

# Filter out negative growth
pos_growth_data <- growth_data %>%
  filter(growth_mm_day >= 0) #Q- what is the effect of inc. 0 growth data vs strictly >0 data?

# Provenance datasets, filtered for eute and eucr, local and dry only
# substring to pull out first 4 characters
provenance_surv_data <- survival_data %>%
  filter(species=="eute" | species=="eute_dry" | species=="eucr" | species=="eucr_dry") %>%
  mutate(species = substring(species, 1, 4), 
         provenance = factor(provenance, levels = c("local", "dry")))

provenance_growth_data <- pos_growth_data %>%
  filter(species=="eute" | species=="eute_dry" | species=="eucr" | species=="eucr_dry") %>%
  mutate(species = substring(species, 1, 4), 
         provenance = factor(provenance, levels = c("local", "dry")))

remove(plot_treatments, species_mixtures, dry_mixtures, plot_mean_elevation, 
       elevation_data, planting_details, 
       survey_data, growth_data)
  