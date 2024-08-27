# Packages ----
library(gbm)         #for gradient boosted models
library(car)
library(dismo) # gradient boosted models and a few additional features
library(pdp) # like a conditional effects plots
library(ggfortify)
library(randomForest) #for running random forests
library(gridExtra)
library(patchwork)
library(easystats)
library(brms)
library(loo)
library(tidybayes)
library(DHARMa)   
library(rstan)    
library(emmeans)
library(bayestestR)
library(ggridges) 
library(tidyverse)
library(readxl)
library(vegan)
library(GGally)
library(corrplot)
library(car)
library(mvabund)
library(scales)
library(ggvegan)
library(ggrepel)
library(glmmTMB)
library(gllvm)
library(EcolUtils)
library(rstan)
library(ggforce)
library(concaveman)

source('R/functions.R')

# Community diversity ----
## Load data ----
data <- read_csv(file = 'data/processed/recruit.csv', col_select = -1) |>
  mutate(Treatment = ifelse(Treatment == 'No algae', 'T1',
                            ifelse(Treatment == 'Only canopy', 'T2',
                                   ifelse(Treatment == 'Only mat', 'T3', 
                                          'T4'))),
         Grazing = factor(Grazing, c('No', 'Light', 'Medium', 'Heavy'), ordered = TRUE)) |>
  mutate(Treatment = factor(Treatment)) |>
  dplyr::select(Tile, Total, Treatment, Grazing, H_mean_broad, D_broad, Shannon_broad_alg) |>
  left_join(read_xlsx(path = 'data/primary/CH3-Algal-community-time-0.xlsx', sheet = 1) |>
              mutate(Taxa = ifelse(Category == 'Coral', 'Coral',
                                   ifelse(Category == 'Other', 'Other', Taxa))) |>
              dplyr::select(Time, Tile, Taxa, Cover) |>
              mutate(Cover = ifelse(Cover == "+", '1', Cover)) |> #adjust rare species 
              mutate(Cover = as.numeric(Cover)) |>
              group_by(Time, Tile, Taxa) |>
              summarise(Cover = sum(Cover)) |>
              ungroup() |>
              group_by(Tile) |>
              mutate(Freq = Cover / sum(Cover)*100) |>
              ungroup() |>
              dplyr::select(-Cover) |>
              tidyr::pivot_wider(names_from = 'Taxa', values_from = 'Freq', values_fill = 0))

data |>
  glimpse()

## nMDS ----
data.mds <- metaMDS(data[,-c(1:5)], k=2,  plot=TRUE) #it always standardise by square root
data.mds #check that the stress is <0.2, the dimensions show to what reduction we drove the data
stressplot(data.mds)

data.mds$stress 
stressplot(data.mds)
data.mds.scores <- data.mds |> 
  fortify() |> 
  mutate(Label = label,
         Score = score) |>
  full_join(data[, 1:5] |> add_rownames(var='Label'))
data.mds.scores.centroids <- data.mds.scores |>
  filter(Score == "sites") |>
  group_by(Treatment) |>
  summarise(across(c(NMDS1, NMDS2), list(c = mean)))
data.mds.scores <- data.mds.scores |>
  full_join(data.mds.scores.centroids)

col_vals <- c("T4" = "#8dd3c7", 
              "T1" = "#fb8072", 
              "T2" = "#d1a95e",
              "T3" = "#a0c58b")

ggplot(data = NULL, aes(y = NMDS2, x = NMDS1)) + 
  ggforce::geom_mark_ellipse(data = data.mds.scores |> 
                               filter(Score=='sites'),
                             aes(y = NMDS2, x = NMDS1, 
                                 fill = Treatment, 
                                 colour = Treatment), expand = 0, alpha = 0.2) + 
  geom_point(data=data.mds.scores |> 
               filter(Score=='sites'),
             aes(color = Treatment, shape = Treatment), size = 2, alpha = 0.6) +
  geom_segment(data = data.mds.scores |> filter(Score == 'species'),
               aes(y = 0, x = 0, yend = NMDS2, xend = NMDS1), alpha = 0.7) +
  geom_text_repel(data = data.mds.scores |> filter(Score == 'species'),
                  aes(y = NMDS2*1.1, x = NMDS1*1.1, label=Label), alpha = 0.7) +
  scale_colour_manual(values = col_vals) + scale_fill_manual(values = col_vals) +
  theme_classic() + theme(legend.position = c(0.85, 0.2))



## brms ----

nMDS.data <- data.mds.scores |> 
  dplyr::filter(score == 'sites') |>
  left_join(data |> 
              dplyr::select(Tile, H_mean_broad, D_broad, Shannon_broad_alg))
nMDS.data |>
  ggplot(mapping = aes(y = NMDS2, x = NMDS1)) +
  geom_point(mapping = aes(size = Total, colour = Shannon_broad_alg)) +
  theme_bw()

nMDS.data |>
  ggplot(mapping = aes(y = NMDS2, x = NMDS1)) +
  geom_point(mapping = aes(size = Total, colour = D_broad)) +
  theme_bw()

# Recruits vs nMDS
data.form <- bf(Total ~ NMDS1 + NMDS2, zi ~ 1, family = zero_inflated_poisson(link = 'log'))
data.form |> get_prior(data = nMDS.data)

priors <- prior(normal(0.5, 3), class = 'Intercept') +
  prior(normal(0, 5), class = 'b')  +
  prior(logistic(0, 1), class = 'Intercept', dpar = 'zi')

nMDS.brm <- brm(data.form, prior = priors, data = nMDS.data, 
                      sample_prior = 'yes', 
                      iter = 5000, 
                      warmup = 1000, 
                      chains = 3, cores = 3, 
                      control = list(adapt_delta = 0.99, 
                                     max_treedepth = 20), 
                      thin = 5, 
                      refresh = 100, 
                      backend = 'rstan') 

nMDS.brm |>
  brms::as_draws_df() |>
  dplyr::select(starts_with('b_')) |>
  mutate(across(everything(), exp)) |>
  summarise_draws(median,
                  HDInterval::hdi,
                  rhat, length, ess_bulk, ess_tail,
                  Pl = ~mean(.x < 1),
                  Pg = ~mean(.x > 1))

# H vs nMDS
data.form <- bf(H_mean_broad ~ NMDS1 + NMDS2, family = gaussian(link = 'log'))
data.form |> get_prior(data = nMDS.data)

priors <- prior(normal(0.5, 3), class = 'Intercept') +
  prior(normal(0, 5), class = 'b')  +
  prior(student_t(3, 0, 1), class = 'sigma')

nMDS.brm2 <- brm(data.form, prior = priors, data = nMDS.data, 
                sample_prior = 'yes', 
                iter = 5000, 
                warmup = 1000, 
                chains = 3, cores = 3, 
                control = list(adapt_delta = 0.99, 
                               max_treedepth = 20), 
                thin = 5, 
                refresh = 100, 
                backend = 'rstan') 

nMDS.brm2 |>
  brms::as_draws_df() |>
  dplyr::select(starts_with('b_')) |>
  mutate(across(everything(), exp)) |>
  summarise_draws(median,
                  HDInterval::hdi,
                  rhat, length, ess_bulk, ess_tail,
                  Pl = ~mean(.x < 1),
                  Pg = ~mean(.x > 1))

# Diversity vs nMDS
data.form <- bf(Shannon_broad_alg ~ NMDS1 + NMDS2, family = gaussian(link = 'log'))
data.form |> get_prior(data = nMDS.data)

priors <- prior(normal(0.5, 3), class = 'Intercept') +
  prior(normal(0, 5), class = 'b')  +
  prior(student_t(3, 0, 1), class = 'sigma')

nMDS.brm3 <- brm(data.form, prior = priors, data = nMDS.data, 
                 sample_prior = 'yes', 
                 iter = 5000, 
                 warmup = 1000, 
                 chains = 3, cores = 3, 
                 control = list(adapt_delta = 0.99, 
                                max_treedepth = 20), 
                 thin = 5, 
                 refresh = 100, 
                 backend = 'rstan') 

nMDS.brm3 |>
  brms::as_draws_df() |>
  dplyr::select(starts_with('b_')) |>
  mutate(across(everything(), exp)) |>
  summarise_draws(median,
                  HDInterval::hdi,
                  rhat, length, ess_bulk, ess_tail,
                  Pl = ~mean(.x < 1),
                  Pg = ~mean(.x > 1))


# Functional diversity ----
## Load data ----
data <- read_csv(file = 'data/processed/recruit.csv', col_select = -1) |>
  mutate(Treatment = ifelse(Treatment == 'No algae', 'T1',
                            ifelse(Treatment == 'Only canopy', 'T2',
                                   ifelse(Treatment == 'Only mat', 'T3', 
                                          'T4'))),
         Grazing = factor(Grazing, c('No', 'Light', 'Medium', 'Heavy'), ordered = TRUE)) |>
  mutate(Treatment = factor(Treatment)) |>
  dplyr::select(Tile, Total, Treatment, Grazing, H_mean_broad, D_broad, Shannon_broad_alg) |>
  left_join(read_xlsx(path = 'data/primary/CH3-Algal-community-time-0.xlsx', sheet = 1) |>
              dplyr::mutate(Taxa = ifelse(Category == 'Coral', 'Coral',
                                   ifelse(Category == 'Other', 'Other', Taxa))) |>
              dplyr::select(Time, Tile, Taxa, Cover) |>
              left_join(read.csv(file = 'data/primary/categories.csv') |>
                          dplyr::select(-Phylum)) |> 
              mutate(Category = ifelse(is.na(Category), Taxa, Category)) |>
              mutate(Taxa = ifelse(Category == 'Alga', Functional.group, Taxa)) |>
              mutate(Cover = ifelse(Cover == "+", '1', Cover)) |>
              mutate(Cover = as.numeric(Cover)) |>
              group_by(Time, Tile, Taxa) |>
              summarise(Cover = sum(Cover)) |>
              ungroup() |>
              group_by(Tile) |>
              mutate(Freq = Cover / sum(Cover)*100) |>
              ungroup() |>
              dplyr::select(-Cover) |>
              tidyr::pivot_wider(names_from = 'Taxa', values_from = 'Freq', values_fill = 0))

data |>
  glimpse()

## nMDS ----
data.mds <- metaMDS(data[,-c(1:8)], k=2,  plot=TRUE) #it always standardise by square root
data.mds #check that the stress is <0.2, the dimensions show to what reduction we drove the data
stressplot(data.mds)

data.mds$stress 
stressplot(data.mds)
data.mds.scores <- data.mds |> 
  fortify() |> 
  mutate(Label = label,
         Score = score) |>
  full_join(data[, 1:8] |> add_rownames(var='Label'))
data.mds.scores.centroids <- data.mds.scores |>
  filter(Score == "sites") |>
  group_by(Treatment) |>
  summarise(across(c(NMDS1, NMDS2), list(c = mean)))
data.mds.scores <- data.mds.scores |>
  full_join(data.mds.scores.centroids)

col_vals <- c("T4" = "#8dd3c7", 
              "T1" = "#fb8072", 
              "T2" = "#d1a95e",
              "T3" = "#a0c58b")

ggplot(data = NULL, aes(y = NMDS2, x = NMDS1)) + 
  ggforce::geom_mark_ellipse(data = data.mds.scores |> 
                               filter(Score=='sites'),
                             aes(y = NMDS2, x = NMDS1, 
                                 fill = Treatment, 
                                 colour = Treatment), expand = 0, alpha = 0.2) + 
  geom_point(data=data.mds.scores |> 
               filter(Score=='sites'),
             aes(color = Treatment, shape = Treatment), size = 2, alpha = 0.6) +
  geom_segment(data = data.mds.scores |> filter(Score == 'species'),
               aes(y = 0, x = 0, yend = NMDS2, xend = NMDS1), alpha = 0.7) +
  geom_text_repel(data = data.mds.scores |> filter(Score == 'species'),
                  aes(y = NMDS2*1.1, x = NMDS1*1.1, label=Label), alpha = 0.7) +
 scale_colour_manual(values = col_vals) + scale_fill_manual(values = col_vals) +
  theme_classic() + theme(legend.position = c(0.85, 0.2))



## brms ----

nMDS.data <- data.mds.scores |> 
  dplyr::filter(score == 'sites') |>
  left_join(data |> 
              dplyr::select(Tile, H_mean_broad, D_broad, Shannon_broad_alg))
nMDS.data |>
  ggplot(mapping = aes(y = NMDS2, x = NMDS1)) +
  geom_point(mapping = aes(size = Total, colour = Shannon_broad_alg)) +
  theme_bw()

nMDS.data |>
  ggplot(mapping = aes(y = NMDS2, x = NMDS1)) +
  geom_point(mapping = aes(size = Total, colour = D_broad)) +
  theme_bw()

nMDS.data |>
  ggplot(mapping = aes(y = NMDS2, x = NMDS1)) +
  geom_point(mapping = aes(size = Total, colour = H_mean_broad)) +
  theme_bw()

# Recruits vs nMDS
data.form <- bf(Total ~ NMDS1 + NMDS2, zi ~ 1, family = zero_inflated_poisson(link = 'log'))
data.form |> get_prior(data = nMDS.data)

priors <- prior(normal(0.5, 3), class = 'Intercept') +
  prior(normal(0, 5), class = 'b')  +
  prior(logistic(0, 1), class = 'Intercept', dpar = 'zi')

nMDS.brm4 <- brm(data.form, prior = priors, data = nMDS.data, 
                sample_prior = 'yes', 
                iter = 5000, 
                warmup = 1000, 
                chains = 3, cores = 3, 
                control = list(adapt_delta = 0.99, 
                               max_treedepth = 20), 
                thin = 5, 
                refresh = 100, 
                backend = 'rstan') 

nMDS.brm4 |>
  brms::as_draws_df() |>
  dplyr::select(starts_with('b_')) |>
  mutate(across(everything(), exp)) |>
  summarise_draws(median,
                  HDInterval::hdi,
                  rhat, length, ess_bulk, ess_tail,
                  Pl = ~mean(.x < 1),
                  Pg = ~mean(.x > 1))

# H vs nMDS
data.form <- bf(H_mean_broad ~ NMDS1 + NMDS2, family = gaussian(link = 'log'))
data.form |> get_prior(data = nMDS.data)

priors <- prior(normal(0.5, 3), class = 'Intercept') +
  prior(normal(0, 5), class = 'b')  +
  prior(student_t(3, 0, 1), class = 'sigma')

nMDS.brm5 <- brm(data.form, prior = priors, data = nMDS.data, 
                 sample_prior = 'yes', 
                 iter = 5000, 
                 warmup = 1000, 
                 chains = 3, cores = 3, 
                 control = list(adapt_delta = 0.99, 
                                max_treedepth = 20), 
                 thin = 5, 
                 refresh = 100, 
                 backend = 'rstan') 

nMDS.brm5 |>
  brms::as_draws_df() |>
  dplyr::select(starts_with('b_')) |>
  mutate(across(everything(), exp)) |>
  summarise_draws(median,
                  HDInterval::hdi,
                  rhat, length, ess_bulk, ess_tail,
                  Pl = ~mean(.x < 1),
                  Pg = ~mean(.x > 1))

# Diversity vs nMDS
data.form <- bf(Shannon_broad_alg ~ NMDS1 + NMDS2, family = gaussian(link = 'log'))
data.form |> get_prior(data = nMDS.data)

priors <- prior(normal(0.5, 3), class = 'Intercept') +
  prior(normal(0, 5), class = 'b')  +
  prior(student_t(3, 0, 1), class = 'sigma')

nMDS.brm6 <- brm(data.form, prior = priors, data = nMDS.data, 
                 sample_prior = 'yes', 
                 iter = 5000, 
                 warmup = 1000, 
                 chains = 3, cores = 3, 
                 control = list(adapt_delta = 0.99, 
                                max_treedepth = 20), 
                 thin = 5, 
                 refresh = 100, 
                 backend = 'rstan') 

 nMDS.brm6 |>
  brms::as_draws_df() |>
  dplyr::select(starts_with('b_')) |>
  mutate(across(everything(), exp)) |>
  summarise_draws(median,
                  HDInterval::hdi,
                  rhat, length, ess_bulk, ess_tail,
                  Pl = ~mean(.x < 1),
                  Pg = ~mean(.x > 1))

