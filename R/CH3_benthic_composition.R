## ---- loadFunctions
source('../R/functions.R')
## ----end

## Benthic composition *****************************************************************

## ---- BenthicComposition readData
data <- read_xlsx(path = '../data/primary/CH3-Algal-community-time-0.xlsx', sheet = 1) %>%
  add_row(., read_xlsx(path = '../data/primary/CH3-Algal-community-time-1.xlsx', sheet = 1)) %>%
  select(!Date)
metadata <- read_xlsx(path = '../data/primary/CH3-metadata.xlsx', sheet = 2) %>%
  select(Tile, Treatment)
## ----end

## ---- BenthicComposition glimpse
data %>% glimpse() 
## ----end

## Process data =================================================================
## ---- BenthicComposition process
algae <- data %>%
  arrange(Tile) %>%
  mutate(Tile = factor(Tile)) %>%
  mutate(Cover = ifelse(Cover == "+", '1', Cover)) %>%
  mutate(Cover = as.numeric(Cover)) %>%
  left_join(., metadata %>%
              mutate(Tile = factor(Tile))) %>%
  select(Time, Tile, Treatment, Taxa, Cover)

canopy <- data %>%
  mutate(Tile = factor(Tile)) %>%
  filter(Taxa == "Sargassum") %>%
  mutate(H_mean = NA,
         H_median = NA,
         H_sd = NA)

for (i in 1:nrow(canopy)) {
  canopy[i, 'H_mean'] <- canopy[i, 'Height'] %>%
    str_split(', ') %>%
    unlist() %>%
    as.numeric() %>%
    mean()
  
  canopy[i, 'H_median'] <- canopy[i, 'Height'] %>%
    str_split(', ') %>%
    unlist() %>%
    as.numeric() %>%
    median()
  
  canopy[i, 'H_sd'] <- canopy[i, 'Height'] %>%
    str_split(', ') %>%
    unlist() %>%
    as.numeric() %>%
    sd()
}


canopy <- algae %>%
  distinct(Time, Tile) %>%
  full_join(., canopy %>% 
              select(!Height)) %>%
  mutate(Taxa = 'Sargassum',
         Cover = ifelse(is.na(Cover), 0, Cover),
         H_mean = ifelse(is.na(H_mean), 0, H_mean),
         H_median = ifelse(is.na(H_median), 0, H_median),
         H_sd = ifelse(is.na(H_sd), 0, H_sd),
         Density = ifelse(is.na(Density), 0, Density)) %>%
  left_join(., metadata %>%
              mutate(Tile = factor(Tile)))
## ----end







data_df <- data %>%
  left_join(., algae_df %>%
              mutate(Perimeter = as.numeric(Perimeter)) %>%
              mutate(Perimeter = Perimeter/40*100) %>%
              select(Tile, Perimeter, Height)) %>%
  mutate(Perimeter = ifelse(!is.na(Perimeter), Perimeter, 0),
         Height = ifelse(!is.na(Height), Height, 0))


summary(glmmTMB(formula = Total ~ Treatment + (1), data = data_df, family = poisson))
summary(glmmTMB(formula = Total ~ Treatment + (1|Perimeter), data = data_df, family = poisson))


dispersion_test <- function(x) 
{
  res <- 1-2 * abs((1 - pchisq((sum((x - mean(x))^2)/mean(x)), length(x) - 1))-0.5)
  
  cat("Dispersion test of count data:\n",
      length(x), " data points.\n",
      "Mean: ",mean(x),"\n",
      "Variance: ",var(x),"\n",
      "Probability of being drawn from Poisson distribution: ", 
      round(res, 3),"\n", sep = "")
  
  invisible(res)
}
dispersion_test(data_df$Total)
