# LIBRARIES----
library(clipr)
library(readxl)
library(tidyverse)

# DATA IMPORT ----
tiles <- read_xlsx(path = 'CH3 - Tiles.xlsx', sheet = 1) %>%
  mutate(Tile = factor(Tile))
treat <- read_xlsx(path = 'CH3-Set-up.xlsx', sheet = 2) %>%
  select(Tile, Treatment, T1_Status) %>%
  mutate(Tile = factor(Tile))


data_peri <- read_xlsx(path = '../data/primary/CH3-Algal-community-time-0.xlsx', sheet = 2) %>%
  add_row(., read_xlsx(path = '../data/primary/CH3-Algal-community-time-1.xlsx', sheet = 2)%>%
            mutate(Height = as.character(Height)))

data_peri %>% glimpse() 


read_xlsx(path = '../data/primary/CH3-recruits-count.xlsx', sheet = 1) %>%
  select(!Bleached)