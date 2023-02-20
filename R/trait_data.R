## Variable info:
## species: missing = unable to be located, unplanted = marker or visible tubestock present, unknown = unable to be ID'd
##                    due to death, herbivory etc

## fresh and clean
rm(list=ls())

## install the necessary packages
#install.packages("tidyverse")
#install.packages("lme4")

## load the necessary libraries
library(tidyverse)

## read in the data
numbering <- read_csv("Data/functional_trait_numbering.csv")

leaf_traits <- read_csv("Data/functional_traits_leaves_inc_leaflets.csv")%>%
  filter(is.na (leaflet) | leaflet == "10leaves") %>%  # filter out rachis etc for compound leaves 
  mutate(lamina_area_cm2 = lamina_area_cm2/n_of_leaves, 
         fresh_leaf_mass_g = fresh_leaf_mass_g/n_of_leaves,
         dry_leaf_mass_g = dry_leaf_mass_g/n_of_leaves) %>% # divide values for number of leaves
  select(plant_number, leaf_sample, lamina_area_cm2, fresh_leaf_mass_g, dry_leaf_mass_g)

stem_root_traits <- read_csv("Data/functional_traits_stems_roots.csv")

seeds <- read_csv("Data/functional_trait_seed_mass.csv")

# SSD g cm-3
# LA cm2
# SLA mm^2/mg
# seed mass g
# RTD mg/mm3
# SRL m/g
all_traits <-numbering %>%
  left_join(leaf_traits, by = "plant_number") %>% 
  left_join(stem_root_traits, by = "plant_number") %>%
  left_join(seeds, by = "species") %>%
  relocate(plant_number) %>% 
  mutate(la = lamina_area_cm2,
         sla = 0.1*lamina_area_cm2/dry_leaf_mass_g, 
         ldmc = dry_leaf_mass_g/fresh_leaf_mass_g, 
         ssd = stem_oven_dry_mass_g/stem_displaced_water_volume_g, 
         rtd = 100*dry_root_mass_g/fresh_root_volume_cm3,
         srl = 0.01*fresh_root_length_cm/dry_root_mass_g,
         seed_mass = 0.001*seed_mass_g_per_1000, .after=leaf_sample) %>% 
  select(plant_number, leaf_sample, species, la, sla, ldmc, ssd, rtd, srl, seed_mass)

write.csv(all_traits, "Outputs/all_traits.csv", row.names=F)

avg_trait_vals <- all_traits %>%
  group_by(species) %>%
  summarise(la = mean(la), sla = mean(sla), ldmc = mean(ldmc), 
            ssd = mean(ssd), rtd = mean(rtd), log_srl = log(mean(srl)), log_seed_mass=log(mean(seed_mass))) %>%
  mutate(species_no = row_number(), .after = species) %>% 
  as.data.frame()

# using avg of combined angophora species because we can't distinguish them..
avg_trait_vals <- avg_trait_vals %>%
  mutate(species = ifelse(species == "anfl" | species == "ansu", 
                          "ango", species)) %>%
  group_by(species) %>%
  summarise(la = mean(la), sla = mean(sla), ldmc = mean(ldmc), 
            ssd = mean(ssd), rtd = mean(rtd), log_srl = mean(log_srl), 
            log_seed_mass = mean(log_seed_mass)) %>% 
  as.data.frame()
row.names(avg_trait_vals)<-avg_trait_vals$species

#### PCA based on measurements of functional traits for four replicates of species;...
trait_pca<-princomp(avg_trait_vals[,c(2:8)], cor = TRUE)
#summary(trait_pca) 
#loadings(trait_pca)

avg_trait_vals$trait_pc1<-trait_pca$scores[,1]
avg_trait_vals$trait_pc2<-trait_pca$scores[,2]
avg_trait_vals$trait_pc3<-trait_pca$scores[,3]

write.csv(avg_trait_vals, "Outputs/avg_trait_vals.csv", row.names=F)



