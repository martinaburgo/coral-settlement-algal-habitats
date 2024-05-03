# Libraries ----
library(tidyverse)   
library(readxl)
library(vegan)
library(stringr)
library(car)
library(brms)
library(patchwork)
library(corrplot)
library(loo)
library(tidybayes)
library(DHARMa)   
library(rstan)    
library(emmeans)
library(bayestestR)
library(ggridges) 

source('R/functions.R')

# Load data ----
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
  dplyr::select(!Date)


# Data wrangling ----
## Benthic data ----
# Adjust cover:
benthos <- benthos |>
  mutate(Cover = as.numeric(ifelse(Cover == "+", '0.5', 
                                   ifelse(Cover == '=', '0.5', Cover)))) |>
  dplyr::filter(!is.na(Cover)) |>
  dplyr::select(!c(Distance, Interaction, Diameter, Height, Thalli)) |>
  group_by(Habitat, Trip, Transect, Quadrat) |>
  mutate(Cover = Cover/sum(Cover)*100) |>
  ungroup()

#Add categories
benthos <- benthos |>
  left_join(read.csv(file = 'data/primary/categories.csv'))


## Canopy height ----

for (i in 1:nrow(canopy)) {
  canopy[i, 'Mean'] <- canopy[i, 'Height'] |>
    str_split(', ') |>
    unlist() |>
    as.numeric() |>
    mean()
  
  canopy[i, 'Max'] <- canopy[i, 'Height'] |>
    str_split(', ') |>
    unlist() |>
    as.numeric() |>
    max()
}

canopy <- canopy |>
  mutate(Mean = as.numeric(Mean),
         Max = as.numeric(Max))


# General benthic composition (Methods) ----
# algae v coral cover trend
benthos |>
  group_by(Habitat, Trip, Transect, Quadrat, Category) |>
  summarise(sum = sum(Cover)) |>
  group_by(Habitat, Trip, Category) |>
  summarise(mean = mean(sum)) |>
  group_by(Habitat, Category) |>
  summarise(mean = mean(mean)) |>
  View()


benthos |>
  dplyr::filter(Taxa == 'Sargassum') |>
  group_by(Habitat, Trip, Transect, Quadrat) |>
  summarise(sum = sum(Cover)) |>
  group_by(Habitat, Trip) |>
  summarise(mean = mean(sum)) |>
  group_by(Habitat) |>
  summarise(mean = mean(mean)) |>
  View()


## Canopy height ----
canopy |>
  dplyr::filter(Habitat == 'Crest') |>
  ggplot(aes(x = Trip, y = Mean)) +
  geom_point() +
  geom_boxplot(aes(group = Trip))

canopy |>
  dplyr::filter(Habitat == 'Crest') |>
  ggplot(aes(x = Trip, y = Max)) +
  geom_point() +
  geom_boxplot(aes(group = Trip))
