# Libraries ----
library(tidyverse)
library(readxl)
library(vegan)

# Load data ----
benthic_broad <- read_xlsx(path = 'data/primary/CH3-Algal-community-time-0.xlsx', sheet = 1) |>
  select(!Date)
benthic_local <- read_xlsx(path = 'data/primary/CH3-Algal-community-time-0.xlsx', sheet = 2) |>
  select(!Date)
data <- read_xlsx(path = 'data/primary/CH3-recruits-count.xlsx', sheet = 1) |>
  full_join(read_xlsx(path = 'data/primary/CH3-metadata.xlsx', sheet = 2) |> #adding treatments
              select(Tile, Treatment)) 
class <- read_xlsx(path = 'data/primary/Benthic_categories.xlsx', sheet = 1)
str(data)

# Data wrangling ----
# First, declare categorical variables and combine top and bottom recruits as almost no recruits were found on the top of the tiles.
data <- data |>
  filter(!is.na(Unbleached)) |>
  mutate(Tile = factor(Tile),
         Treatment = factor(Treatment, levels = c("No algae", "Only mat", "Only canopy", 
                                                  "All algae"),
                            ordered = TRUE),
         Grazing = factor(Grazing)) |>
  group_by(Tile) |>
  mutate(Total = sum(Unbleached)) |>
  distinct(Tile, .keep_all = TRUE) |>
  select(Tile, Total, Treatment, Grazing)

## Benthic composition ----
benthic_broad <- benthic_broad |>
  mutate(Cover = ifelse(Cover == "+", '0.5', Cover)) |>
  mutate(Cover = as.numeric(Cover)) |> #adjust rare species
  group_by(Tile) |>
  mutate(Freq_broad = Cover / sum(Cover)*100)

benthic_local <- benthic_local |>
  mutate(Perimeter = ifelse(Perimeter == "+", '1', Perimeter)) |>
  mutate(Perimeter = as.numeric(Perimeter)) |> #adjust rare species
  group_by(Tile) |>
  mutate(Freq_local = Perimeter / sum(Perimeter)*100)

## Canopy ----
### Broad ----
canopy_broad <- benthic_broad |>
  filter(Taxa == "Sargassum") |>
  mutate(H_mean_broad = NA,
         H_sd_broad = NA) |>
  rename(Freq_Sarg_broad = Freq_broad)

for (i in 1:nrow(canopy_broad)) {
  canopy_broad[i, 'H_mean_broad'] <- canopy_broad[i, 'Height'] |>
    str_split(', ') |>
    unlist() |>
    as.numeric() |>
    mean()
  
  canopy_broad[i, 'H_sd_broad'] <- canopy_broad[i, 'Height'] |>
    str_split(', ') |>
    unlist() |>
    as.numeric() |>
    sd()
}

### Local ----
canopy_local <- benthic_local |>
  filter(Taxa == "Sargassum") |>
  mutate(H_mean_local = NA,
         H_sd_local = NA,
         Freq_Sarg_local = Freq_local)

for (i in 1:nrow(canopy_local)) {
  canopy_local[i, 'H_mean_local'] <- canopy_local[i, 'Height'] |>
    str_split(', ') |>
    unlist() |>
    as.numeric() |>
    mean()
  
  canopy_local[i, 'H_sd_local'] <- canopy_local[i, 'Height'] |>
    str_split(', ') |>
    unlist() |>
    as.numeric() |>
    sd()
}

## Species Richness and Divesity ----
### Richness ----
R_broad <- benthic_broad |>
  left_join(class) |>
  filter(Classification != 'Other') |>
  group_by(Tile, Classification) |>
  summarise(R_broad = n()) |>
  pivot_wider(names_from = Classification, values_from = R_broad, values_fill = 0) |>
  rename(R_broad_alg = Alga,
         R_broad_coral = Coral)

R_local <- benthic_local |>
  left_join(class) |>
  filter(Classification != 'Other') |>
  group_by(Tile, Classification) |>
  summarise(n = n()) |>
  pivot_wider(names_from = Classification, values_from = n, values_fill = 0) |>
  rename(R_local_alg = Alga,
         R_local_coral = Coral)

Shannon_broad_alg <- benthic_broad |>
  left_join(class) |>
  filter(Classification == 'Alga') |>
  distinct(Tile) |>
  arrange(-Tile) |>
  add_column(benthic_broad |>
               left_join(class) |>
               filter(Classification == 'Alga') |>
               dplyr::select(-Time, -Height, -Density, -Cover, -Classification) |>
               pivot_wider(names_from = Taxa, values_from = Freq_broad, values_fill = 0) |>
               arrange(-Tile) |>
               diversity(index = 'shannon') |>
               as_tibble()) |>
  rename(Shannon_broad_alg = value)

### Diversity ----
Shannon_broad_alg <- benthic_broad |>
  left_join(class) |>
  filter(Classification == 'Alga') |>
  distinct(Tile) |>
  arrange(-Tile) |>
  add_column(benthic_broad |>
               left_join(class) |>
               filter(Classification == 'Alga') |>
               dplyr::select(-Time, -Height, -Density, -Cover, -Classification) |>
               pivot_wider(names_from = Taxa, values_from = Freq_broad, values_fill = 0) |>
               arrange(-Tile) |>
               diversity(index = 'shannon') |>
               as_tibble()) |>
  rename(Shannon_broad_alg = value)

Shannon_local_alg <- benthic_local |>
  left_join(class) |>
  filter(Classification == 'Alga') |>
  distinct(Tile) |>
  arrange(-Tile) |>
  add_column(benthic_local |>
               left_join(class) |>
               filter(Classification == 'Alga') |>
               dplyr::select(-Time, -Height, -Perimeter, -Classification) |>
               pivot_wider(names_from = Taxa, values_from = Freq_local, values_fill = 0) |>
               arrange(-Tile) |>
               diversity(index = 'shannon') |>
               as_tibble()) |>
  rename(Shannon_local_alg = value)

Shannon_broad_cor <- benthic_broad |>
  left_join(class) |>
  filter(Classification == 'Coral') |>
  distinct(Tile) |>
  arrange(-Tile) |>
  add_column(benthic_broad |>
               left_join(class) |>
               filter(Classification == 'Coral') |>
               dplyr::select(-Time, -Height, -Density, -Cover, -Classification) |>
               pivot_wider(names_from = Taxa, values_from = Freq_broad, values_fill = 0) |>
               arrange(-Tile) |>
               diversity(index = 'shannon') |>
               as_tibble()) |>
  rename(Shannon_broad_cor = value)

Shannon_local_cor <- benthic_local |>
  left_join(class) |>
  filter(Classification == 'Coral') |>
  distinct(Tile) |>
  arrange(-Tile) |>
  add_column(benthic_local |>
               left_join(class) |>
               filter(Classification == 'Coral') |>
               dplyr::select(-Time, -Height, -Perimeter, -Classification) |>
               pivot_wider(names_from = Taxa, values_from = Freq_local, values_fill = 0) |>
               arrange(-Tile) |>
               diversity(index = 'shannon') |>
               as_tibble()) |>
  rename(Shannon_local_cor = value)

### Coral cover ----
Tot_broad <- benthic_broad |>
  left_join(class) |>
  filter(Classification != 'Other') |>
  group_by(Tile, Classification) |>
  summarise(n = sum(Freq_broad)) |>
  pivot_wider(names_from = Classification, values_from = n, values_fill = 0) |>
  rename(Tot_broad_alg = Alga,
         Tot_broad_coral = Coral)

Tot_local <- benthic_local |>
  left_join(class) |>
  filter(Classification != 'Other') |>
  group_by(Tile, Classification) |>
  summarise(n = sum(Freq_local)) |>
  pivot_wider(names_from = Classification, values_from = n, values_fill = 0) |>
  rename(Tot_local_alg = Alga,
         Tot_local_coral = Coral)

# Combine datasets ----
recruit <- data |>
  left_join(canopy_broad |>
              mutate(Tile = factor(Tile)) |>
              select(Tile, H_mean_broad, H_sd_broad, Freq_Sarg_broad)) |>
  left_join(canopy_local |>
              mutate(Tile = factor(Tile)) |>
              select(Tile, H_mean_local, H_sd_local, Freq_Sarg_local)) |>
  mutate(Hadj_broad = H_sd_broad/H_mean_broad,
         Hadj_local = H_sd_local/H_mean_local) |>
  left_join(benthic_broad |>
              mutate(Tile = factor(Tile)) |>
              group_by(Tile) |>
              filter(!is.na(Density)) |>
              select(Tile, Density, Freq_broad)) |>
  mutate(D_broad = Density/Freq_broad) |>
  ungroup() |>
  select(-Freq_broad) |>
  left_join(R_broad |>
              mutate(Tile = factor(Tile))) |>
  left_join(R_local |>
              mutate(Tile = factor(Tile))) |>
  left_join(Shannon_broad_alg |>
              mutate(Tile = factor(Tile))) |>
  left_join(Shannon_broad_cor |>
              mutate(Tile = factor(Tile))) |>
  left_join(Shannon_local_alg |>
              mutate(Tile = factor(Tile))) |>
  left_join(Shannon_local_cor |>
              mutate(Tile = factor(Tile))) |>
  left_join(Tot_broad |>
              mutate(Tile = factor(Tile))) |>
  left_join(Tot_local |>
              mutate(Tile = factor(Tile))) |>
  mutate_all(~replace(., is.na(.), 0))

# Save ----
write.csv(recruit, file = '../data/processed/recruit.csv')
