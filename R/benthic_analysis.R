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
canopy <- rbind(readxl::read_xlsx(path = "data/primary/Maggie-Quads-3.xlsx", sheet = "Florence"),
                readxl::read_xlsx(path = "data/primary/Maggie-Quads-4.xlsx", sheet = "Florence"),
                readxl::read_xlsx(path = "data/primary/Maggie-Quads-5.xlsx", sheet = "Florence")) |>
  dplyr::select(Habitat, Trip, Transect, Quadrat, Cover, Height, Thalli) |>
  dplyr::filter(!is.na(Thalli))

# Data wrangling ----
canopy |>
  expand(Habitat, Trip, Transect, Quadrat) |>
  full_join(canopy)


for (i in 1:nrow(canopy)) {
  algae_df[i, 'Height'] <- algae_temp[i,] |>
    str_split(', ') |>
    unlist() |>
    as.numeric() |>
    mean()
}

algae_df <- algae_df %>%
  mutate(Height = as.numeric(Height))
plot(x = algae_df$Thalli, y = algae_df$Height)  

summary(lm(formula = 'Height ~ Thalli', data = algae_df))
