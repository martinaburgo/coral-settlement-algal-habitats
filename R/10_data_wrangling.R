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

# TIle ----
# First, declare categorical variables and combine top and bottom recruits as almost no recruits were found on the top of the tiles.
data <- data |>
  filter(!is.na(Unbleached)) |>
  mutate(Tile = factor(Tile),
         Treatment = factor(Treatment),
         Grazing = factor(Grazing)) |>
  group_by(Tile) |>
  mutate(Total = sum(Unbleached)) |>
  distinct(Tile, .keep_all = TRUE) |>
  select(Tile, Total, Treatment, Grazing)

## Benthic composition ----
benthic_broad <- benthic_broad |>
  mutate(Cover = ifelse(Cover == "+", '1', Cover)) |>
  mutate(Cover = as.numeric(Cover)) |>
  group_by(Tile) |>
  mutate(Freq_broad = Cover / sum(Cover)*100)

benthic_local <- benthic_local |>
  mutate(Perimeter = ifelse(Perimeter == "+", '1', Perimeter)) |>
  mutate(Perimeter = as.numeric(Perimeter)) |> 
  group_by(Tile) |>
  mutate(Freq_local = Perimeter / sum(Perimeter)*100)

## Canopy ----
### Broad ----
canopy_broad <- benthic_broad |>
  filter(Taxa == "Sargassum") |>
  mutate(H_mean_broad = NA) |>
  rename(Freq_Sarg_broad = Freq_broad)

for (i in 1:nrow(canopy_broad)) {
  canopy_broad[i, 'H_mean_broad'] <- canopy_broad[i, 'Height'] |>
    str_split(', ') |>
    unlist() |>
    as.numeric() |>
    mean()
}

### Local ----
canopy_local <- benthic_local |>
  filter(Taxa == "Sargassum") |>
  mutate(H_mean_local = NA,
         Freq_Sarg_local = Freq_local)

for (i in 1:nrow(canopy_local)) {
  canopy_local[i, 'H_mean_local'] <- canopy_local[i, 'Height'] |>
    str_split(', ') |>
    unlist() |>
    as.numeric() |>
    mean()
}

## Species Richness and Divesity ----
### Richness ----
R_broad <- benthic_broad |>
  filter(Category != 'Other') |>
  group_by(Tile, Category) |>
  summarise(R_broad = n()) |>
  pivot_wider(names_from = Category, values_from = R_broad, values_fill = 0) |>
  rename(R_broad_alg = Macroalgae)

R_local <- benthic_local |>
  filter(Category != 'Other') |>
  group_by(Tile, Category) |>
  summarise(n = n()) |>
  pivot_wider(names_from = Category, values_from = n, values_fill = 0) |>
  rename(R_local_alg = Macroalgae)

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

### Coral cover ----
Tot_broad <- benthic_broad |>
  filter(Category != 'Other') |>
  group_by(Tile, Category) |>
  summarise(n = sum(Freq_broad)) |>
  pivot_wider(names_from = Category, values_from = n, values_fill = 0) |>
  rename(Tot_broad_alg = Macroalgae)

Tot_local <- benthic_local |>
  filter(Category != 'Other') |>
  group_by(Tile, Category) |>
  summarise(n = sum(Freq_local)) |>
  pivot_wider(names_from = Category, values_from = n, values_fill = 0) |>
  rename(Tot_local_alg = Macroalgae)

## Combine datasets ----
recruit <- data |>
  left_join(canopy_broad |>
              mutate(Tile = factor(Tile)) |>
              select(Tile, H_mean_broad, Freq_Sarg_broad)) |>
  left_join(canopy_local |>
              mutate(Tile = factor(Tile)) |>
              select(Tile, H_mean_local, Freq_Sarg_local)) |>
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
              dplyr::select(Tile, R_broad_alg)) |> 
  left_join(R_local |>
              mutate(Tile = factor(Tile)) |>
              dplyr::select(Tile, R_local_alg)) |> 
  left_join(Shannon_broad_alg |>
              mutate(Tile = factor(Tile))) |> 
  left_join(Shannon_local_alg |>
              mutate(Tile = factor(Tile))) |>  
  left_join(Tot_broad |>
              mutate(Tile = factor(Tile)) |>
            dplyr::select(Tot_broad_alg)) |> 
  left_join(Tot_local |>
              mutate(Tile = factor(Tile)) |>
              dplyr::select(Tot_local_alg)) |>
  mutate_all(~replace(., is.na(.), 0)) |>
  left_join(tile_benthos)

### Save ----
#write.csv(recruit, file = 'data/processed/recruit.csv')


# Juveniles density ----
# Load data
corals <- readxl::read_xlsx(path = "data/primary/Maggie-Quads-8.xlsx", sheet = "Florence") |>
  dplyr::select(Habitat, Trip, Period, Transect, Quadrat, Taxa, Interaction, Diameter, Distance, Cover, Height) |>
  dplyr::filter(!is.na(Diameter))
canopy <- rbind(readxl::read_xlsx(path = "data/primary/Maggie-Quads-1.xlsx", sheet = "Florence"),
                readxl::read_xlsx(path = "data/primary/Maggie-Quads-2.xlsx", sheet = "Florence"),
                readxl::read_xlsx(path = "data/primary/Maggie-Quads-3.xlsx", sheet = "Florence"),
                readxl::read_xlsx(path = "data/primary/Maggie-Quads-4.xlsx", sheet = "Florence"),
                readxl::read_xlsx(path = "data/primary/Maggie-Quads-5.xlsx", sheet = "Florence"),
                readxl::read_xlsx(path = "data/primary/Maggie-Quads-6.xlsx", sheet = "Florence"),
                readxl::read_xlsx(path = "data/primary/Maggie-Quads-7.xlsx", sheet = "Florence"),
                readxl::read_xlsx(path = "data/primary/Maggie-Quads-8.xlsx", sheet = "Florence")) |>
  dplyr::select(!c(Date, Metre, Taxa, Diameter, Interaction, Distance)) |>
  dplyr::filter(!is.na(Thalli))
# create coral class sizes
breaks <- c(0, 5, 20,  40, max(corals$Diameter, na.rm = TRUE)) # set up cut-off values 
tags <- c("<5", "6-20", '21-40', '>40') # specify interval/bin labels
group_tags <- cut(corals$Diameter, 
                  breaks = breaks, 
                  include.lowest = TRUE, 
                  right = TRUE, 
                  labels = tags) # bucketing values into bins
summary(group_tags) # inspect bins
class_sizes <- factor(group_tags, 
                      levels = tags,
                      ordered = TRUE)
# Calculate canopy height
for (i in 1:nrow(canopy)) {
  canopy[i, 'Mean'] <- canopy[i, 'Height'] |>
    str_split(', ') |>
    unlist() |>
    as.numeric() |>
    mean()
}

canopy <- canopy |>
  mutate(Mean = as.numeric(Mean))

data <- corals |> 
  mutate(class_sizes = factor(group_tags, 
                              levels = tags,
                              ordered = TRUE)) |>
  dplyr::select(Habitat, Transect, Quadrat, Taxa, Diameter, class_sizes) |>
  left_join(canopy |> 
              dplyr::filter(Trip == '3' | Trip == '2') |>
              dplyr::select(Habitat, Transect, Trip, Quadrat, Mean) |>
              group_by(Habitat, Transect, Quadrat) |>
              summarise(Mean = mean(Mean))) |>
  mutate(across(where(is.numeric), coalesce, 0)) |>
  mutate(Depth = ifelse(Habitat == 'Flat', '0-1 m',
                        ifelse(Habitat == 'Crest', '2-3 m',
                               '4-5 m')))

## Save dataset ----
#write.csv(data, file = 'data/processed/natural_recruitment.csv')

