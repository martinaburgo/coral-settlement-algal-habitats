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
tile_benthos <- readxl::read_xlsx(path = "data/primary/CH3-tile-benthos.xlsx", sheet = 1) |>
  mutate(Sediment = factor(Sediment,
                           c('No', 'Low', 'Medium', 'High'), ordered = TRUE),
         Turf_height = factor(Turf_height,
                              c('No', 'Low', 'Medium', 'High'), ordered = TRUE),
         Tile = factor(Tile))
str(data)

# Data wrangling - TILES ----
# First, declare categorical variables and combine top and bottom recruits as almost no recruits were found on the top of the tiles.
data <- data |>
  filter(!is.na(Unbleached)) |>
  mutate(Tile = factor(Tile),
         Treatment = factor(Treatment, levels = c("No algae", "Only mat", "Only canopy", 
                                                  "All algae")),
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
  filter(Category != 'Other') |>
  group_by(Tile, Category) |>
  summarise(R_broad = n()) |>
  pivot_wider(names_from = Category, values_from = R_broad, values_fill = 0) |>
  rename(R_broad_alg = Macroalgae,
         R_broad_coral = Coral)

R_local <- benthic_local |>
  filter(Category != 'Other') |>
  group_by(Tile, Category) |>
  summarise(n = n()) |>
  pivot_wider(names_from = Category, values_from = n, values_fill = 0) |>
  rename(R_local_alg = Macroalgae,
         R_local_coral = Coral)

Shannon_broad_alg <- benthic_broad |>
  filter(Category == 'Macroalgae') |>
  distinct(Tile) |>
  arrange(-Tile) |>
  add_column(benthic_broad |>
               filter(Category == 'Macroalgae') |>
               dplyr::select(-Time, -Height, -Density, -Cover, -Category) |>
               pivot_wider(names_from = Taxa, values_from = Freq_broad, values_fill = 0) |>
               arrange(-Tile) |>
               diversity(index = 'shannon') |>
               as_tibble()) |>
  rename(Shannon_broad_alg = value)

### Diversity ----
Shannon_broad_alg <- benthic_broad |>
  filter(Category == 'Macroalgae') |>
  distinct(Tile) |>
  arrange(-Tile) |>
  add_column(benthic_broad |>
               filter(Category == 'Macroalgae') |>
               dplyr::select(-Time, -Height, -Density, -Cover, -Category) |>
               pivot_wider(names_from = Taxa, values_from = Freq_broad, values_fill = 0) |>
               arrange(-Tile) |>
               diversity(index = 'shannon') |>
               as_tibble()) |>
  rename(Shannon_broad_alg = value)

Shannon_local_alg <- benthic_local |>
  filter(Category == 'Macroalgae') |>
  distinct(Tile) |>
  arrange(-Tile) |>
  add_column(benthic_local |>
               filter(Category == 'Macroalgae') |>
               dplyr::select(-Time, -Height, -Perimeter, -Category) |>
               pivot_wider(names_from = Taxa, values_from = Freq_local, values_fill = 0) |>
               arrange(-Tile) |>
               diversity(index = 'shannon') |>
               as_tibble()) |>
  rename(Shannon_local_alg = value)

Shannon_broad_cor <- benthic_broad |>
  filter(Category == 'Coral') |>
  distinct(Tile) |>
  arrange(-Tile) |>
  add_column(benthic_broad |>
               filter(Category == 'Coral') |>
               dplyr::select(-Time, -Height, -Density, -Cover, -Category) |>
               pivot_wider(names_from = Taxa, values_from = Freq_broad, values_fill = 0) |>
               arrange(-Tile) |>
               diversity(index = 'shannon') |>
               as_tibble()) |>
  rename(Shannon_broad_cor = value)

Shannon_local_cor <- benthic_local |>
  filter(Category == 'Coral') |>
  distinct(Tile) |>
  arrange(-Tile) |>
  add_column(benthic_local |>
               filter(Category == 'Coral') |>
               dplyr::select(-Time, -Height, -Perimeter, -Category) |>
               pivot_wider(names_from = Taxa, values_from = Freq_local, values_fill = 0) |>
               arrange(-Tile) |>
               diversity(index = 'shannon') |>
               as_tibble()) |>
  rename(Shannon_local_cor = value)

### Coral cover ----
Tot_broad <- benthic_broad |>
  filter(Category != 'Other') |>
  group_by(Tile, Category) |>
  summarise(n = sum(Freq_broad)) |>
  pivot_wider(names_from = Category, values_from = n, values_fill = 0) |>
  rename(Tot_broad_alg = Macroalgae,
         Tot_broad_coral = Coral)

Tot_local <- benthic_local |>
  filter(Category != 'Other') |>
  group_by(Tile, Category) |>
  summarise(n = sum(Freq_local)) |>
  pivot_wider(names_from = Category, values_from = n, values_fill = 0) |>
  rename(Tot_local_alg = Macroalgae,
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
  select(-Freq_broad)


recruit <- recruit |>
  left_join(R_broad |>
              mutate(Tile = factor(Tile)) |>
              dplyr::select(Tile, R_broad_coral, R_broad_alg)) |> 
  left_join(R_local |>
              mutate(Tile = factor(Tile)) |>
              dplyr::select(Tile, R_local_coral, R_local_alg)) |> 
  left_join(Shannon_broad_alg |>
              mutate(Tile = factor(Tile))) |> 
  left_join(Shannon_broad_cor |>
              mutate(Tile = factor(Tile))) |>
  left_join(Shannon_local_alg |>
              mutate(Tile = factor(Tile))) |>  
  left_join(Shannon_local_cor |>
              mutate(Tile = factor(Tile))) |> 
  left_join(Tot_broad |>
              mutate(Tile = factor(Tile)) |>
            dplyr::select(Tot_broad_coral, Tot_broad_alg)) |> 
  left_join(Tot_local |>
              mutate(Tile = factor(Tile)) |>
              dplyr::select(Tot_local_coral, Tot_local_alg)) |>
  mutate_all(~replace(., is.na(.), 0)) |>
  left_join(tile_benthos)

# Save ----
write.csv(recruit, file = 'data/processed/recruit.csv')


# Data wrangling - NATURAL RECRUITMENT ----
## Load data ----
corals <- readxl::read_xlsx(path = "data/primary/Maggie-Quads-8.xlsx", sheet = "Florence") |>
  dplyr::select(Habitat, Trip, Period, Transect, Quadrat, Taxa, Interaction, Diameter, Distance, Cover, Height)
canopy <- rbind(readxl::read_xlsx(path = "data/primary/Maggie-Quads-1.xlsx", sheet = "Florence"),
                readxl::read_xlsx(path = "data/primary/Maggie-Quads-2.xlsx", sheet = "Florence"),
                readxl::read_xlsx(path = "data/primary/Maggie-Quads-3.xlsx", sheet = "Florence"),
                readxl::read_xlsx(path = "data/primary/Maggie-Quads-4.xlsx", sheet = "Florence"),
                readxl::read_xlsx(path = "data/primary/Maggie-Quads-5.xlsx", sheet = "Florence"),
                readxl::read_xlsx(path = "data/primary/Maggie-Quads-6.xlsx", sheet = "Florence"),
                readxl::read_xlsx(path = "data/primary/Maggie-Quads-7.xlsx", sheet = "Florence"),
                readxl::read_xlsx(path = "data/primary/Maggie-Quads-8.xlsx", sheet = "Florence")) |>
  dplyr::select(!c(Date, Metre, Taxa, Diameter, Interaction, Distance, Cover)) |>
  dplyr::filter(!is.na(Thalli))

benthos <- rbind(readxl::read_xlsx(path = "data/primary/Maggie-Quads-1.xlsx", sheet = "Florence"),
                 readxl::read_xlsx(path = "data/primary/Maggie-Quads-2.xlsx", sheet = "Florence"),
                 readxl::read_xlsx(path = "data/primary/Maggie-Quads-3.xlsx", sheet = "Florence"),
                 readxl::read_xlsx(path = "data/primary/Maggie-Quads-4.xlsx", sheet = "Florence"),
                 readxl::read_xlsx(path = "data/primary/Maggie-Quads-5.xlsx", sheet = "Florence"),
                 readxl::read_xlsx(path = "data/primary/Maggie-Quads-6.xlsx", sheet = "Florence"),
                 readxl::read_xlsx(path = "data/primary/Maggie-Quads-7.xlsx", sheet = "Florence"),
                 readxl::read_xlsx(path = "data/primary/Maggie-Quads-8.xlsx", sheet = "Florence")) |>
  dplyr::select(!c(Date))

## Final dataset for benthic analysis----

# Shannon's diversity
shannon_div <- benthos |>
  dplyr::filter(Trip == '3' | Trip == '4') |>
  dplyr::filter(Category == 'Alga') |>
  dplyr::select(!c(Island, Site, Metre, Category)) |>
  dplyr::distinct(Habitat, Trip, Transect, Quadrat) |>
  add_column(benthos |>
               dplyr::filter(Trip == '3' | Trip == '4') |>
               dplyr::filter(Category == 'Alga') |>
               pivot_wider(names_from = Taxa, values_from = Cover, values_fill = 0) |>
               dplyr::select(!c(Island, Site, Habitat, Period, Trip, 
                                Transect, Quadrat, Metre, Category)) |>
               diversity(index = 'shannon') |>
               as_tibble()) |>
  rename(Shannon_algae = value) |>
  group_by(Habitat, Transect, Quadrat) |>
  summarise(Shannon_algae = mean(Shannon_algae)) |>
  full_join(benthos |>
              dplyr::filter(Trip == '3' | Trip == '4') |>
              dplyr::filter(Category == 'Coral') |> 
              group_by(Habitat, Trip, Transect, Quadrat, Taxa) |>
              summarise(Cover = mean(Cover)) |>
              ungroup() |>
              dplyr::distinct(Habitat, Trip, Transect, Quadrat) |>
              add_column(benthos |>
                           dplyr::filter(Trip == '3' | Trip == '4') |>
                           dplyr::filter(Category == 'Coral') |> 
                           group_by(Habitat, Trip, Transect, Quadrat, Taxa) |>
                           summarise(Cover = mean(Cover)) |>
                           ungroup() |>
                           pivot_wider(names_from = Taxa, values_from = Cover, values_fill = 0) |>
                           dplyr::select(!c(Habitat, Trip, Transect, Quadrat)) |>
                           diversity(index = 'shannon') |>
                           as_tibble()) |>
              rename(Shannon_coral = value) |>
              group_by(Habitat, Transect, Quadrat) |>
              summarise(Shannon_coral = mean(Shannon_coral)))  |>
  mutate(Shannon_algae = ifelse(Shannon_algae == 0, 0.01, Shannon_algae),
         Shannon_coral = ifelse(Shannon_coral == 0, 0.01, Shannon_coral))


# Height and density of the canopy (trip 3 and 4)
H_D <- benthos |>
  dplyr::filter(Trip == '3' | Trip == '4') |>
  dplyr::filter(Taxa == 'Sargassum') |>
  dplyr::select(Habitat, Trip, Transect, Quadrat, Cover) |>
  full_join(canopy |>
              dplyr::filter(Trip == '3' | Trip == '4')) |>
  mutate(Density = (Thalli/(Cover)*100)) |>
  dplyr::select(!c(Cover, Island, Site, Height, Thalli, Mean, Max)) |>
  group_by(Habitat, Transect, Quadrat) |>
  summarise(Density = mean(Density)) |>
  ungroup() |>
  full_join(benthos |>
              dplyr::filter(Trip == '4') |>
              dplyr::filter(Taxa == 'Sargassum') |>
              dplyr::select(Habitat, Trip, Transect, Quadrat, Cover) |>
              full_join(canopy |>
                          dplyr::filter(Trip == '4')) |>
              mutate(Extra_d = (Thalli/(Cover)*100)) |>
              dplyr::select(!c(Cover, Island, Trip, Period, Site, Height, Thalli, Mean, Max))) |>
  mutate(Density = ifelse(is.na(Density), Extra_d, Density)) |>
  dplyr::select(-Extra_d) |>
  full_join(canopy |>
              dplyr::filter(Trip == '3' | Trip == '4') |>
              dplyr::select(!c(Height, Mean)) |>
              group_by(Habitat, Transect, Quadrat) |>
              summarise(Height = mean(Max)))

data <- H_D |>
  full_join(shannon_div) |>
  right_join(corals |>
               dplyr::filter(!is.na(Diameter)) |>
               dplyr::select(Habitat, Transect, Quadrat, Taxa, Diameter)) |>
  mutate(Shannon_coral = ifelse(is.na(Shannon_coral), 0, Shannon_coral))

# Save dataset ----
write.csv(data, file = 'data/processed/natural_recruitment.csv')

