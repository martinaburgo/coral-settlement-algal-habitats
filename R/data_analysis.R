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
recruit <- read_csv(file = 'data/processed/recruit.csv', col_select = -1)
tile_benthos <- readxl::read_xlsx(path = "data/primary/CH3-tile-benthos.xlsx", sheet = 1) |>
  mutate(Sediment = factor(Sediment,
                           c('No', 'Low', 'Medium', 'High'), ordered = TRUE),
         Turf_height = factor(Turf_height,
                              c('No', 'Low', 'Medium', 'High'), ordered = TRUE))

recruit <- recruit |>
  mutate(Treatment = factor(Treatment),
         Grazing = factor(Grazing, c('No', 'Light', 'Medium', 'Heavy'), ordered = TRUE)) |>
  full_join(tile_benthos)