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

# M1: Recruit ~ Treat ----
## Fit model ----
form <- bf(Total ~ Treatment + (1|Sediment), family = poisson(link = 'log')) 
form |> get_prior(data = recruit)

priors <- prior(normal(0.5, 5), class = 'Intercept') +
  prior(normal(0, 10), class = 'b') +
  prior(student_t(3, 0, 3), class = 'sd')

recruit.brm <- brm(form, prior = priors, data = recruit, 
                   sample_prior = 'only', 
                   iter = 5000, 
                   warmup = 1000, 
                   chains = 3, cores = 3, 
                   thin = 5, 
                   refresh = 0, 
                   control = list(adapt_delta = 0.99, 
                                  max_treedepth = 20),
                   backend = 'rstan') 

recruit.brm |> 
  conditional_effects() |> 
  plot(points = TRUE)

recruit.brm2 <- recruit.brm |>
  update(sample_prior = 'yes')

## Diagnostics ----
recruit.brm2 |> 
  SUYR_prior_and_posterior()

## MCMC Sampling diagnostics
(recruit.brm2$fit |> stan_trace()) + (recruit.brm2$fit |> stan_ac()) + (recruit.brm2$fit |> stan_rhat()) + (recruit.brm2$fit |> stan_ess())

## Model validation
### Posterior probability check
recruit.brm2 |> pp_check(type = 'dens_overlay', ndraws = 100)

### Residuals
recruit.resids <- make_brms_dharma_res(recruit.brm2, integerResponse = FALSE)
testUniformity(recruit.resids)
plotResiduals(recruit.resids, form = factor(rep(1, nrow(recruit))))
plotResiduals(recruit.resids, quantreg = TRUE) 
testDispersion(recruit.resids)

## Save model ----
save(recruit.brm2, form, priors, recruit, file = 'data/modelled/M1_base.RData')

## Investigation ----
recruit.brm2 |> 
  conditional_effects() |> 
  plot(points = TRUE)

### Output ----
### ---- M1Output
recruit.brm2 |>
  brms::as_draws_df() |>
  mutate(across(everything(), exp)) |>
  summarise_draws(median,
                  HDInterval::hdi,
                  rhat, length, ess_bulk, ess_tail,
                  Pl = ~mean(.x < 1),
                  Pg = ~mean(.x > 1))
### ----end

recruit.brm2 |>
  brms::as_draws_df() |>
  mutate(across(everything(), exp)) |>
  summarise_draws(median,
                  HDInterval::hdi,
                  rhat, length, ess_bulk, ess_tail,
                  Pl = ~mean(.x < 1),
                  Pg = ~mean(.x > 1)) |>
  mutate(median = round(median, 3),
         lower = round(lower, 3),
         upper = round(upper, 3),
         rhat = round(rhat, 3),
         Pl = round(Pl, 3),
         Pg = round(Pg, 3)) |>
  write.table(file = 'output/tables/M1Output.txt', sep = ",", quote = FALSE, row.names = F)

### Contrast ----
### ---- M1Contrast
recruit.brm2 |> 
  emmeans(~Treatment, type = 'link') |>
  pairs(reverse = TRUE) |>
  gather_emmeans_draws() |> #it's on log scale, so we need to mutate it
  mutate(Response = exp(.value)) |>
  summarise(median_hdci(Response),
            Pl = mean(Response < 1),
            Pg = mean(Response > 1))
### ----end

recruit.brm2 |> 
  emmeans(~Treatment, type = 'link') |>
  pairs(reverse = TRUE) |>
  gather_emmeans_draws() |> #it's on log scale, so we need to mutate it
  mutate(Response = exp(.value)) |>
  summarise(median_hdci(Response),
            Pl = mean(Response < 1),
            Pg = mean(Response > 1)) |>
  write.table(file = 'output/tables/M1Contrast.txt', sep = ",", quote = FALSE, row.names = F)

### ROPE ----
recruit.brm2 |> emmeans(~Treatment, type = 'link') |>
  pairs(reverse = TRUE) |>
  gather_emmeans_draws() |>
  mutate(Response = exp(.value)) |>
  summarise(median_hdci(Response),
            ROPE = rope(Response, range = c(0.9, 1.1))$ROPE_Percentage)


### Figure ----
recruit.em <- recruit.brm2 |> 
  emmeans(~Treatment, type = 'link') |>
  pairs() |>
  gather_emmeans_draws() |>
  mutate(Fit = exp(.value))

M1ContrastFig <- recruit.em |> 
  ggplot(aes(x = Fit)) +
  geom_density_ridges_gradient(aes(y = contrast, fill = stat(x)),
                               alpha = 0.4, color = 'white', 
                               quantile_lines = TRUE,
                               quantiles = c(0.025, 0.975)) + 
  geom_vline(xintercept = 1, linetype = 'dashed') + 
  scale_fill_viridis_c(option = 'C') + 
  scale_x_continuous('', trans = scales::log_trans())  + 
  scale_y_discrete('') +
  theme_classic() +
  theme(text = element_text(colour = 'black'), 
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.5)),
        legend.position = 'none')
M1ContrastFig

ggsave(filename = 'output/figures/M1ContrastFig.png', width = 8, height = 5, dpi = 100)

## Further investigation ----

### ---- PlannedContrast
cmat <- cbind('Mat_Canopy' = c(-1/2, 1/2, -1/2, 1/2),
              'NoAlg_AnyAlg' = c(-1/3, 1, -1/3, -1/3),
              'Mat_AnyAlg' = c(-1/2, 0, -1/2, 1),
              'NoAlg_Canopy' = c(-1/2, 1, -1/2, 0))
recruit.brm2 |>
  emmeans(~Treatment, type = 'link') |>
  contrast(method = list(Treatment = cmat)) |>
  gather_emmeans_draws() |>
  mutate(Response = exp(.value)) |>
  summarise(median_hdci(Response),
            Pl = mean(Response < 1),
            Pg = mean(Response > 1),
            ROPE = rope(Response, range = c(0.9, 1.1))$ROPE_Percentage)
### ----end


recruit.brm2 |>
  emmeans(~Treatment, type = 'link') |>
  contrast(method = list(Treatment = cmat)) |>
  gather_emmeans_draws() |>
  mutate(Response = exp(.value)) |>
  summarise(median_hdci(Response),
            Pl = mean(Response < 1),
            Pg = mean(Response > 1)) |>
  write.table(file = 'output/tables/M1PlannedContrast.txt', sep = ",", quote = FALSE, row.names = F)

# M2: Recruit ~ Height (broad) ----
## Fit model ----
form <- bf(Total ~ H_mean_broad + (1|Sediment), family = poisson(link = 'log')) 
form |> get_prior(data = recruit)

priors <- prior(normal(0.5, 5), class = 'Intercept') +
  prior(normal(6, 5), class = 'b')  +
  prior(student_t(3, 0, 3), class = 'sd')

recruit.brm3 <- brm(form, prior = priors, data = recruit, 
                   sample_prior = 'only', 
                   iter = 5000, 
                   warmup = 1000, 
                   chains = 3, cores = 3, 
                   control = list(adapt_delta = 0.99, 
                                  max_treedepth = 20), 
                   thin = 5, 
                   refresh = 0, 
                   backend = 'rstan') 

recruit.brm3 |> 
  conditional_effects() |> 
  plot(points = TRUE)

recruit.brm4 <- recruit.brm3 |>
  update(sample_prior = 'yes')

## Diagnostics ----
recruit.brm4 |> 
  SUYR_prior_and_posterior()

## MCMC Sampling diagnostics
(recruit.brm4$fit |> stan_trace()) + (recruit.brm4$fit |> stan_ac()) + (recruit.brm4$fit |> stan_rhat()) + (recruit.brm4$fit |> stan_ess())

## Model validation
### Posterior probability check
recruit.brm4 |> pp_check(type = 'dens_overlay', ndraws = 100)

### Residuals
recruit.resids <- make_brms_dharma_res(recruit.brm4, integerResponse = FALSE)
testUniformity(recruit.resids)
plotResiduals(recruit.resids, form = factor(rep(1, nrow(recruit))))
plotResiduals(recruit.resids, quantreg = TRUE) 
testDispersion(recruit.resids)

## Save model ----
save(recruit.brm4, form, priors, recruit, file = 'data/modelled/M2_Height_broad.RData')

## Investigation ----
recruit.brm4 |> 
  conditional_effects() |> 
  plot(points = TRUE)

### Output ----
### ---- M2Output
recruit.brm4 |>
  brms::as_draws_df() |>
  mutate(across(everything(), exp)) |>
  summarise_draws(median,
                  HDInterval::hdi,
                  rhat, length, ess_bulk, ess_tail,
                  Pl = ~mean(.x < 1),
                  Pg = ~mean(.x > 1))
### ----end

recruit.brm4 |>
  brms::as_draws_df() |>
  mutate(across(everything(), exp)) |>
  summarise_draws(median,
                  HDInterval::hdi,
                  rhat, length, ess_bulk, ess_tail,
                  Pl = ~mean(.x < 1),
                  Pg = ~mean(.x > 1)) |>
  mutate(median = round(median, 3),
         lower = round(lower, 3),
         upper = round(upper, 3),
         rhat = round(rhat, 3),
         Pl = round(Pl, 3),
         Pg = round(Pg, 3)) |>
  write.table(file = 'output/tables/M2Output.txt', sep = ",", quote = FALSE, row.names = F)

### Figure ----
newdata <- with(recruit,
                list(H_mean_broad = seq(min(H_mean_broad),
                                        max(H_mean_broad),
                                        len = 50)))
fit <- recruit.brm4 |> 
  emmeans(~H_mean_broad, at = newdata) |>
  gather_emmeans_draws() |>
  mutate(.value = exp(.value)) |>
  group_by(H_mean_broad) |>
  summarise(median_hdci(.value)) |>
  as.data.frame()

M2Figure <- fit |> 
  ggplot(aes(x = H_mean_broad,
             y = y)) + 
  geom_line() + 
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.2) +
  scale_x_continuous(name = expression(paste(italic('Sargassum'), ' height (cm)'))) +
  scale_y_continuous('Total recruits') +
  theme_classic()  +
  theme(text = element_text(colour = 'black'), 
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.5)))
M2Figure

ggsave(filename = 'output/figures/M2Figure.png', width = 8, height = 5, dpi = 100)

# M3: Recruit ~ Height (local) ----
## Fit model ----
form <- bf(Total ~ H_mean_local + (1|Sediment), family = poisson(link = 'log')) 
form |> get_prior(data = recruit)

priors <- prior(normal(0.5, 5), class = 'Intercept') +
  prior(normal(6, 8), class = 'b') + 
  prior(student_t(3, 0, 3), class = 'sd')

recruit.brm5 <- brm(form, prior = priors, data = recruit, 
                    sample_prior = 'only', 
                    iter = 5000, 
                    warmup = 1000, 
                    chains = 3, cores = 3, 
                    thin = 5, 
                    refresh = 0, 
                    control = list(adapt_delta = 0.99, 
                                   max_treedepth = 20), 
                    backend = 'rstan') 

recruit.brm5 |> 
  conditional_effects() |> 
  plot(points = TRUE)

recruit.brm6 <- recruit.brm5 |>
  update(sample_prior = 'yes')

## Diagnostics ----
recruit.brm6 |> 
  SUYR_prior_and_posterior()

## MCMC Sampling diagnostics
(recruit.brm6$fit |> stan_trace()) + (recruit.brm6$fit |> stan_ac()) + (recruit.brm6$fit |> stan_rhat()) + (recruit.brm6$fit |> stan_ess())

## Model validation
### Posterior probability check
recruit.brm6 |> pp_check(type = 'dens_overlay', ndraws = 100)

### Residuals
recruit.resids <- make_brms_dharma_res(recruit.brm6, integerResponse = FALSE)
testUniformity(recruit.resids)
plotResiduals(recruit.resids, form = factor(rep(1, nrow(recruit))))
plotResiduals(recruit.resids, quantreg = TRUE) 
testDispersion(recruit.resids)

## Save model ----
save(recruit.brm6, form, priors, recruit, file = 'data/modelled/M3_Height_local.RData')

## Investigation ----
recruit.brm6 |> 
  conditional_effects() |> 
  plot(points = TRUE)

### Output ----
### ---- M3Output
recruit.brm6 |>
  brms::as_draws_df() |>
  mutate(across(everything(), exp)) |>
  summarise_draws(median,
                  HDInterval::hdi,
                  rhat, length, ess_bulk, ess_tail,
                  Pl = ~mean(.x < 1),
                  Pg = ~mean(.x > 1))
### ----end

recruit.brm6 |>
  brms::as_draws_df() |>
  mutate(across(everything(), exp)) |>
  summarise_draws(median,
                  HDInterval::hdi,
                  rhat, length, ess_bulk, ess_tail,
                  Pl = ~mean(.x < 1),
                  Pg = ~mean(.x > 1)) |>
  mutate(median = round(median, 3),
         lower = round(lower, 3),
         upper = round(upper, 3),
         rhat = round(rhat, 3),
         Pl = round(Pl, 3),
         Pg = round(Pg, 3)) |>
  write.table(file = 'output/tables/M3Output.txt', sep = ",", quote = FALSE, row.names = F)

### Figure ----
newdata <- with(recruit,
                list(H_mean_local = seq(min(H_mean_local),
                                        max(H_mean_local),
                                        len = 50)))
fit <- recruit.brm6 |> 
  emmeans(~H_mean_local, at = newdata) |>
  gather_emmeans_draws() |>
  mutate(.value = exp(.value)) |>
  group_by(H_mean_local) |>
  summarise(median_hdci(.value)) |>
  as.data.frame()

M3Figure <- fit |> 
  ggplot(aes(x = H_mean_local,
             y = y)) + 
  geom_line() + 
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.2) +
  scale_x_continuous(name = expression(paste(italic('Sargassum'), ' height (cm)'))) +
  scale_y_continuous('Total recruits') +
  theme_classic()  +
  theme(text = element_text(colour = 'black'), 
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.5)))
M3Figure

ggsave(filename = 'output/figures/M3Figure.png', width = 8, height = 5, dpi = 100)

# Compare ----
## ---- CompareM2vsM3
loo::loo_compare(loo::loo(recruit.brm4),
                 loo::loo(recruit.brm6))
## ----end


# M4: Recruit ~ Density ----
## Fit model ----
form <- bf(Total ~ D_broad + (1|Sediment), family = poisson(link = 'log')) 
form |> get_prior(data = recruit)

priors <- prior(normal(0.5, 20), class = 'Intercept') +
  prior(normal(1, 3), class = 'b') +
  prior(student_t(3, 0, 3), class = 'sd')

recruit.brm7 <- brm(form, prior = priors, data = recruit, 
                    sample_prior = 'only', 
                    iter = 5000, 
                    warmup = 1000, 
                    chains = 3, cores = 3, 
                    control = list(adapt_delta = 0.99, 
                                   max_treedepth = 20), 
                    thin = 5, 
                    refresh = 0, 
                    backend = 'rstan') 

recruit.brm7 |> 
  conditional_effects() |> 
  plot(points = TRUE)

recruit.brm8 <- recruit.brm7 |>
  update(sample_prior = 'yes')

## Diagnostics ----
recruit.brm8 |> 
  SUYR_prior_and_posterior()

## MCMC Sampling diagnostics
(recruit.brm8$fit |> stan_trace()) + (recruit.brm8$fit |> stan_ac()) + (recruit.brm8$fit |> stan_rhat()) + (recruit.brm8$fit |> stan_ess())

## Model validation
### Posterior probability check
recruit.brm8 |> pp_check(type = 'dens_overlay', ndraws = 100)

### Residuals
recruit.resids <- make_brms_dharma_res(recruit.brm8, integerResponse = FALSE)
testUniformity(recruit.resids)
plotResiduals(recruit.resids, form = factor(rep(1, nrow(recruit))))
plotResiduals(recruit.resids, quantreg = TRUE) 
testDispersion(recruit.resids)

## Save model ----
save(recruit.brm8, form, priors, recruit, file = 'data/modelled/M4_Density.RData')

## Investigation ----
recruit.brm8 |> 
  conditional_effects() |> 
  plot(points = TRUE)

### Output ----
### ---- M4Output
recruit.brm8 |>
  brms::as_draws_df() |>
  mutate(across(everything(), exp)) |>
  summarise_draws(median,
                  HDInterval::hdi,
                  rhat, length, ess_bulk, ess_tail,
                  Pl = ~mean(.x < 1),
                  Pg = ~mean(.x > 1))
### ----end

recruit.brm8 |>
  brms::as_draws_df() |>
  mutate(across(everything(), exp)) |>
  summarise_draws(median,
                  HDInterval::hdi,
                  rhat, length, ess_bulk, ess_tail,
                  Pl = ~mean(.x < 1),
                  Pg = ~mean(.x > 1)) |>
  mutate(median = round(median, 3),
         lower = round(lower, 3),
         upper = round(upper, 3),
         rhat = round(rhat, 3),
         Pl = round(Pl, 3),
         Pg = round(Pg, 3)) |>
  write.table(file = 'output/tables/M4Output.txt', sep = ",", quote = FALSE, row.names = F)

### Figure ----
newdata <- with(recruit,
                list(D_broad = seq(min(D_broad),
                                   max(D_broad),
                                        len = 50)))
fit <- recruit.brm8 |> 
  emmeans(~D_broad, at = newdata) |>
  gather_emmeans_draws() |>
  mutate(.value = exp(.value)) |>
  group_by(D_broad) |>
  summarise(median_hdci(.value)) |>
  as.data.frame()

M4Figure <- fit |> 
  ggplot(aes(x = D_broad,
             y = y)) + 
  geom_line() + 
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.2) +
  scale_x_continuous(name = expression(paste(italic('Sargassum'), ' density'))) +
  scale_y_continuous('Total recruits') +
  theme_classic()  +
  theme(text = element_text(colour = 'black'), 
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.5)))
M4Figure

ggsave(filename = 'output/figures/M4Figure.png', width = 8, height = 5, dpi = 100)

# M5: Recruit ~ Species Richness (broad) ----
## Fit model ----
form <- bf(Total ~ scale(R_broad_alg) + scale(R_broad_coral) + (1|Grazing), family = poisson(link = 'log')) 
form |> get_prior(data = recruit)

priors <- prior(normal(0.5, 10), class = 'Intercept') +
  prior(normal(2, 3), class = 'b') +
  prior(student_t(3, 0, 3), class = 'sd')

recruit.brm9 <- brm(form, prior = priors, data = recruit, 
                    sample_prior = 'only', 
                    iter = 5000, 
                    warmup = 1000, 
                    chains = 3, cores = 3, 
                    control = list(adapt_delta = 0.99, 
                                   max_treedepth = 20), 
                    thin = 5, 
                    refresh = 0, 
                    backend = 'rstan') 

recruit.brm9 |> 
  conditional_effects() |> 
  plot(points = TRUE)

recruit.brm10 <- recruit.brm9 |>
  update(sample_prior = 'yes')

## Diagnostics ----
recruit.brm10 |> 
  SUYR_prior_and_posterior()

## MCMC Sampling diagnostics
(recruit.brm10$fit |> stan_trace()) + (recruit.brm10$fit |> stan_ac()) + (recruit.brm10$fit |> stan_rhat()) + (recruit.brm10$fit |> stan_ess())

## Model validation
### Posterior probability check
recruit.brm10 |> pp_check(type = 'dens_overlay', ndraws = 100)

### Residuals
recruit.resids <- make_brms_dharma_res(recruit.brm10, integerResponse = FALSE)
testUniformity(recruit.resids)
plotResiduals(recruit.resids, form = factor(rep(1, nrow(recruit))))
plotResiduals(recruit.resids, quantreg = TRUE) 
testDispersion(recruit.resids)

## Save model ----
save(recruit.brm10, form, priors, recruit, file = 'data/modelled/M5_Richness_broad.RData')

## Investigation ----
recruit.brm10 |> 
  conditional_effects() |> 
  plot(points = TRUE, ask = FALSE)

### Output ----
### ---- M5Output
recruit.brm10 |>
  brms::as_draws_df() |>
  mutate(across(everything(), exp)) |>
  summarise_draws(median,
                  HDInterval::hdi,
                  rhat, length, ess_bulk, ess_tail,
                  Pl = ~mean(.x < 1),
                  Pg = ~mean(.x > 1))
### ----end

recruit.brm10 |>
  brms::as_draws_df() |>
  mutate(across(everything(), exp)) |>
  summarise_draws(median,
                  HDInterval::hdi,
                  rhat, length, ess_bulk, ess_tail,
                  Pl = ~mean(.x < 1),
                  Pg = ~mean(.x > 1)) |>
  mutate(median = round(median, 3),
         lower = round(lower, 3),
         upper = round(upper, 3),
         rhat = round(rhat, 3),
         Pl = round(Pl, 3),
         Pg = round(Pg, 3)) |>
  write.table(file = 'output/tables/M5Output.txt', sep = ",", quote = FALSE, row.names = F)

# M6: Recruit ~ Species Richness (local) ----
## Fit model ----
form <- bf(Total ~ scale(R_local_alg) + scale(R_local_coral) + (1|Grazing), family = poisson(link = 'log')) 
form |> get_prior(data = recruit)

priors <- prior(normal(0.5, 10), class = 'Intercept') +
  prior(normal(2, 3), class = 'b') +
  prior(student_t(3, 0, 3), class = 'sd')

recruit.brm11 <- brm(form, prior = priors, data = recruit, 
                    sample_prior = 'only', 
                    iter = 5000, 
                    warmup = 1000, 
                    chains = 3, cores = 3, 
                    control = list(adapt_delta = 0.99, 
                                   max_treedepth = 20), 
                    thin = 5, 
                    refresh = 0, 
                    backend = 'rstan') 

recruit.brm11 |> 
  conditional_effects() |> 
  plot(points = TRUE)

recruit.brm12 <- recruit.brm11 |>
  update(sample_prior = 'yes')

## Diagnostics ----
recruit.brm12 |> 
  SUYR_prior_and_posterior()

## MCMC Sampling diagnostics
(recruit.brm12$fit |> stan_trace()) + (recruit.brm12$fit |> stan_ac()) + (recruit.brm12$fit |> stan_rhat()) + (recruit.brm12$fit |> stan_ess())

## Model validation
### Posterior probability check
recruit.brm12 |> pp_check(type = 'dens_overlay', ndraws = 100)

### Residuals
recruit.resids <- make_brms_dharma_res(recruit.brm12, integerResponse = FALSE)
testUniformity(recruit.resids)
plotResiduals(recruit.resids, form = factor(rep(1, nrow(recruit))))
plotResiduals(recruit.resids, quantreg = TRUE) 
testDispersion(recruit.resids)

## Save model ----
save(recruit.brm12, form, priors, recruit, file = 'data/modelled/M6_Richness_local.RData')

## Investigation ----
recruit.brm12 |> 
  conditional_effects() |> 
  plot(points = TRUE, ask = FALSE)

### Output ----
### ---- M6Output
recruit.brm12 |>
  brms::as_draws_df() |>
  mutate(across(everything(), exp)) |>
  summarise_draws(median,
                  HDInterval::hdi,
                  rhat, length, ess_bulk, ess_tail,
                  Pl = ~mean(.x < 1),
                  Pg = ~mean(.x > 1))
### ----end

recruit.brm12 |>
  brms::as_draws_df() |>
  mutate(across(everything(), exp)) |>
  summarise_draws(median,
                  HDInterval::hdi,
                  rhat, length, ess_bulk, ess_tail,
                  Pl = ~mean(.x < 1),
                  Pg = ~mean(.x > 1)) |>
  mutate(median = round(median, 3),
         lower = round(lower, 3),
         upper = round(upper, 3),
         rhat = round(rhat, 3),
         Pl = round(Pl, 3),
         Pg = round(Pg, 3)) |>
  write.table(file = 'output/tables/M6Output.txt', sep = ",", quote = FALSE, row.names = F)

# M7: Recruit ~ Diversity (broad) ----
## Fit model ----
form <- bf(Total ~ scale(Shannon_broad_alg) + scale(Shannon_broad_cor) + (1|Grazing), family = poisson(link = 'log')) 
form |> get_prior(data = recruit)

priors <- prior(normal(0.5, 10), class = 'Intercept') +
  prior(normal(1, 2), class = 'b') + 
  prior(student_t(3, 0, 3), class = 'sd')

recruit.brm13 <- brm(form, prior = priors, data = recruit, 
                     sample_prior = 'only', 
                     iter = 5000, 
                     warmup = 1000, 
                     chains = 3, cores = 3, 
                     control = list(adapt_delta = 0.99, 
                                    max_treedepth = 20), 
                     thin = 5, 
                     refresh = 0, 
                     backend = 'rstan') 

recruit.brm13 |> 
  conditional_effects() |> 
  plot(points = TRUE)

recruit.brm14 <- recruit.brm13 |>
  update(sample_prior = 'yes')

## Diagnostics ----
recruit.brm14 |> 
  SUYR_prior_and_posterior()

## MCMC Sampling diagnostics
(recruit.brm14$fit |> stan_trace()) + (recruit.brm14$fit |> stan_ac()) + (recruit.brm14$fit |> stan_rhat()) + (recruit.brm14$fit |> stan_ess())

## Model validation
### Posterior probability check
recruit.brm14 |> pp_check(type = 'dens_overlay', ndraws = 100)

### Residuals
recruit.resids <- make_brms_dharma_res(recruit.brm14, integerResponse = FALSE)
testUniformity(recruit.resids)
plotResiduals(recruit.resids, form = factor(rep(1, nrow(recruit))))
plotResiduals(recruit.resids, quantreg = TRUE) 
testDispersion(recruit.resids)

## Save model ----
save(recruit.brm14, form, priors, recruit, file = 'data/modelled/M7_Diversity_broad.RData')

## Investigation ----
recruit.brm14 |> 
  conditional_effects() |> 
  plot(points = TRUE, ask = FALSE)

### Output ----
### ---- M7Output
recruit.brm14 |>
  brms::as_draws_df() |>
  mutate(across(everything(), exp)) |>
  summarise_draws(median,
                  HDInterval::hdi,
                  rhat, length, ess_bulk, ess_tail,
                  Pl = ~mean(.x < 1),
                  Pg = ~mean(.x > 1))
### ----end

recruit.brm14 |>
  brms::as_draws_df() |>
  mutate(across(everything(), exp)) |>
  summarise_draws(median,
                  HDInterval::hdi,
                  rhat, length, ess_bulk, ess_tail,
                  Pl = ~mean(.x < 1),
                  Pg = ~mean(.x > 1)) |>
  mutate(median = round(median, 3),
         lower = round(lower, 3),
         upper = round(upper, 3),
         rhat = round(rhat, 3),
         Pl = round(Pl, 3),
         Pg = round(Pg, 3)) |>
  write.table(file = 'output/tables/M7Output.txt', sep = ",", quote = FALSE, row.names = F)

### Figure ----
newdata <- with(recruit,
                list(Shannon_broad_cor = seq(min(Shannon_broad_cor),
                                   max(Shannon_broad_cor),
                                   len = 25)))
fit <- recruit.brm14 |> 
  emmeans(~Shannon_broad_cor, at = newdata) |>
  gather_emmeans_draws() |>
  mutate(.value = exp(.value)) |>
  group_by(Shannon_broad_cor) |>
  summarise(median_hdci(.value)) |>
  as.data.frame()

M7Figure <- fit |> 
  ggplot(aes(x = Shannon_broad_cor,
             y = y)) + 
  geom_line() + 
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.2) +
  scale_x_continuous('Shannon Diversity - corals (broad)') +
  scale_y_continuous('Total recruits') +
  theme_classic()  +
  theme(text = element_text(colour = 'black'), 
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.5)))
M7Figure

ggsave(filename = 'output/figures/M7Figure.png', width = 8, height = 5, dpi = 100)

# M8: Recruit ~ Diversity (local) ----
## Fit model ----
form <- bf(Total ~ scale(Shannon_local_alg) + scale(Shannon_local_cor) + (1|Grazing), family = poisson(link = 'log')) 
form |> get_prior(data = recruit)

priors <- prior(normal(0.5, 10), class = 'Intercept') +
  prior(normal(1, 2), class = 'b') +
  prior(student_t(3, 0, 3), class = 'sd')

recruit.brm15 <- brm(form, prior = priors, data = recruit, 
                     sample_prior = 'only', 
                     iter = 5000, 
                     warmup = 1000, 
                     chains = 3, cores = 3, 
                     control = list(adapt_delta = 0.99, 
                                    max_treedepth = 20), 
                     thin = 5, 
                     refresh = 0, 
                     backend = 'rstan') 

recruit.brm15 |> 
  conditional_effects() |> 
  plot(points = TRUE)

recruit.brm16 <- recruit.brm15 |>
  update(sample_prior = 'yes')

## Diagnostics ----
recruit.brm16 |> 
  SUYR_prior_and_posterior()

## MCMC Sampling diagnostics
(recruit.brm16$fit |> stan_trace()) + (recruit.brm16$fit |> stan_ac()) + (recruit.brm16$fit |> stan_rhat()) + (recruit.brm16$fit |> stan_ess())

## Model validation
### Posterior probability check
recruit.brm16 |> pp_check(type = 'dens_overlay', ndraws = 100)

### Residuals
recruit.resids <- make_brms_dharma_res(recruit.brm16, integerResponse = FALSE)
testUniformity(recruit.resids)
plotResiduals(recruit.resids, form = factor(rep(1, nrow(recruit))))
plotResiduals(recruit.resids, quantreg = TRUE) 
testDispersion(recruit.resids)

## Save model ----
save(recruit.brm16, form, priors, recruit, file = 'data/modelled/M8_Diversity_local.RData')

## Investigation ----
recruit.brm16 |> 
  conditional_effects() |> 
  plot(points = TRUE, ask = FALSE)

### Output ----
### ---- M8Output
recruit.brm16 |>
  brms::as_draws_df() |>
  mutate(across(everything(), exp)) |>
  summarise_draws(median,
                  HDInterval::hdi,
                  rhat, length, ess_bulk, ess_tail,
                  Pl = ~mean(.x < 1),
                  Pg = ~mean(.x > 1))
### ----end

recruit.brm16 |>
  brms::as_draws_df() |>
  mutate(across(everything(), exp)) |>
  summarise_draws(median,
                  HDInterval::hdi,
                  rhat, length, ess_bulk, ess_tail,
                  Pl = ~mean(.x < 1),
                  Pg = ~mean(.x > 1)) |>
  mutate(median = round(median, 3),
         lower = round(lower, 3),
         upper = round(upper, 3),
         rhat = round(rhat, 3),
         Pl = round(Pl, 3),
         Pg = round(Pg, 3)) |>
  write.table(file = 'output/tables/M8Output.txt', sep = ",", quote = FALSE, row.names = F)

# M9: Recruit ~ Total cover (broad) ----
## Fit model ----
form <- bf(Total ~ scale(Tot_broad_alg) + scale(Tot_broad_coral) + (1|Grazing), family = poisson(link = 'log')) 
form |> get_prior(data = recruit)

priors <- prior(normal(0.5, 10), class = 'Intercept') +
  prior(normal(35, 30), class = 'b') +
  prior(student_t(3, 0, 3), class = 'sd')

recruit.brm17 <- brm(form, prior = priors, data = recruit, 
                     sample_prior = 'only', 
                     iter = 5000, 
                     warmup = 1000, 
                     chains = 3, cores = 3, 
                     control = list(adapt_delta = 0.99, 
                                    max_treedepth = 20), 
                     thin = 5, 
                     refresh = 0, 
                     backend = 'rstan') 

recruit.brm17 |> 
  conditional_effects() |> 
  plot(points = TRUE)

recruit.brm18 <- recruit.brm17 |>
  update(sample_prior = 'yes')

## Diagnostics ----
recruit.brm18 |> 
  SUYR_prior_and_posterior()

## MCMC Sampling diagnostics
(recruit.brm18$fit |> stan_trace()) + (recruit.brm18$fit |> stan_ac()) + (recruit.brm18$fit |> stan_rhat()) + (recruit.brm18$fit |> stan_ess())

## Model validation
### Posterior probability check
recruit.brm18 |> pp_check(type = 'dens_overlay', ndraws = 100)

### Residuals
recruit.resids <- make_brms_dharma_res(recruit.brm18, integerResponse = FALSE)
testUniformity(recruit.resids)
plotResiduals(recruit.resids, form = factor(rep(1, nrow(recruit))))
plotResiduals(recruit.resids, quantreg = TRUE) 
testDispersion(recruit.resids)

## Save model ----
save(recruit.brm18, form, priors, recruit, file = 'data/modelled/M9_Cover_broad.RData')

## Investigation ----
recruit.brm18 |> 
  conditional_effects() |> 
  plot(points = TRUE, ask = FALSE)

### Output ----
### ---- M9Output
recruit.brm18 |>
  brms::as_draws_df() |>
  mutate(across(everything(), exp)) |>
  summarise_draws(median,
                  HDInterval::hdi,
                  rhat, length, ess_bulk, ess_tail,
                  Pl = ~mean(.x < 1),
                  Pg = ~mean(.x > 1))
### ----end

recruit.brm18 |>
  brms::as_draws_df() |>
  mutate(across(everything(), exp)) |>
  summarise_draws(median,
                  HDInterval::hdi,
                  rhat, length, ess_bulk, ess_tail,
                  Pl = ~mean(.x < 1),
                  Pg = ~mean(.x > 1)) |>
  mutate(median = round(median, 3),
         lower = round(lower, 3),
         upper = round(upper, 3),
         rhat = round(rhat, 3),
         Pl = round(Pl, 3),
         Pg = round(Pg, 3)) |>
  write.table(file = 'output/tables/M9Output.txt', sep = ",", quote = FALSE, row.names = F)

# M10: Recruit ~ Total cover (local) ----
## Fit model ----
form <- bf(Total ~ scale(Tot_local_alg) + scale(Tot_local_coral) + (1|Grazing), family = poisson(link = 'log')) 
form |> get_prior(data = recruit)

priors <- prior(normal(0.5, 10), class = 'Intercept') +
  prior(normal(35, 30), class = 'b') + 
  prior(student_t(3, 0, 3), class = 'sd')

recruit.brm19 <- brm(form, prior = priors, data = recruit, 
                     sample_prior = 'only', 
                     iter = 5000, 
                     warmup = 1000, 
                     chains = 3, cores = 3, 
                     control = list(adapt_delta = 0.99, 
                                    max_treedepth = 20), 
                     thin = 5, 
                     refresh = 0, 
                     backend = 'rstan') 

recruit.brm19 |> 
  conditional_effects() |> 
  plot(points = TRUE)

recruit.brm20 <- recruit.brm19 |>
  update(sample_prior = 'yes')

## Diagnostics ----
recruit.brm20 |> 
  SUYR_prior_and_posterior()

## MCMC Sampling diagnostics
(recruit.brm20$fit |> stan_trace()) + (recruit.brm20$fit |> stan_ac()) + (recruit.brm20$fit |> stan_rhat()) + (recruit.brm20$fit |> stan_ess())

## Model validation
### Posterior probability check
recruit.brm20 |> pp_check(type = 'dens_overlay', ndraws = 100)

### Residuals
recruit.resids <- make_brms_dharma_res(recruit.brm20, integerResponse = FALSE)
testUniformity(recruit.resids)
plotResiduals(recruit.resids, form = factor(rep(1, nrow(recruit))))
plotResiduals(recruit.resids, quantreg = TRUE) 
testDispersion(recruit.resids)

## Save model ----
save(recruit.brm20, form, priors, recruit, file = 'data/modelled/M10_Cover_local.RData')

## Investigation ----
recruit.brm20 |> 
  conditional_effects() |> 
  plot(points = TRUE, ask = FALSE)

### Output ----
### ---- M9Output
recruit.brm20 |>
  brms::as_draws_df() |>
  mutate(across(everything(), exp)) |>
  summarise_draws(median,
                  HDInterval::hdi,
                  rhat, length, ess_bulk, ess_tail,
                  Pl = ~mean(.x < 1),
                  Pg = ~mean(.x > 1))
### ----end

recruit.brm20 |>
  brms::as_draws_df() |>
  mutate(across(everything(), exp)) |>
  summarise_draws(median,
                  HDInterval::hdi,
                  rhat, length, ess_bulk, ess_tail,
                  Pl = ~mean(.x < 1),
                  Pg = ~mean(.x > 1)) |>
  mutate(median = round(median, 3),
         lower = round(lower, 3),
         upper = round(upper, 3),
         rhat = round(rhat, 3),
         Pl = round(Pl, 3),
         Pg = round(Pg, 3)) |>
  write.table(file = 'output/tables/M10Output.txt', sep = ",", quote = FALSE, row.names = F)

# M11: Recruit ~ Treat + Diversity corals ----
## Fit model ----
form <- bf(Total ~ Treatment + Shannon_broad_cor + (1|Grazing), family = poisson(link = 'log')) 
form |> get_prior(data = recruit)

priors <- prior(normal(0.5, 10), class = 'Intercept') +
  prior(normal(1, 10), class = 'b') +
  prior(student_t(3, 0, 3), class = 'sd')

recruit.brm21 <- brm(form, prior = priors, data = recruit, 
                   sample_prior = 'only', 
                   iter = 5000, 
                   warmup = 1000, 
                   chains = 3, cores = 3, 
                   control = list(adapt_delta = 0.99, 
                                  max_treedepth = 20), 
                   thin = 5, 
                   refresh = 0, 
                   backend = 'rstan') 

recruit.brm21 |> 
  conditional_effects() |> 
  plot(points = TRUE)

recruit.brm22 <- recruit.brm21 |>
  update(sample_prior = 'yes')

## Diagnostics ----
recruit.brm22 |> 
  SUYR_prior_and_posterior()

## MCMC Sampling diagnostics
(recruit.brm22$fit |> stan_trace()) + (recruit.brm22$fit |> stan_ac()) + (recruit.brm22$fit |> stan_rhat()) + (recruit.brm22$fit |> stan_ess())

## Model validation
### Posterior probability check
recruit.brm22 |> pp_check(type = 'dens_overlay', ndraws = 100)

### Residuals
recruit.resids <- make_brms_dharma_res(recruit.brm22, integerResponse = FALSE)
testUniformity(recruit.resids)
plotResiduals(recruit.resids, form = factor(rep(1, nrow(recruit))))
plotResiduals(recruit.resids, quantreg = TRUE) 
testDispersion(recruit.resids)

## Save model ----
save(recruit.brm22, form, priors, recruit, file = 'data/modelled/M11_final.RData')

## Investigation ----
recruit.brm22 |> 
  conditional_effects() |> 
  plot(points = TRUE)

### Output ----
### ---- M11Output
recruit.brm22 |>
  brms::as_draws_df() |>
  mutate(across(everything(), exp)) |>
  summarise_draws(median,
                  HDInterval::hdi,
                  rhat, length, ess_bulk, ess_tail,
                  Pl = ~mean(.x < 1),
                  Pg = ~mean(.x > 1))
### ----end

recruit.brm22 |>
  brms::as_draws_df() |>
  mutate(across(everything(), exp)) |>
  summarise_draws(median,
                  HDInterval::hdi,
                  rhat, length, ess_bulk, ess_tail,
                  Pl = ~mean(.x < 1),
                  Pg = ~mean(.x > 1)) |>
  mutate(median = round(median, 3),
         lower = round(lower, 3),
         upper = round(upper, 3),
         rhat = round(rhat, 3),
         Pl = round(Pl, 3),
         Pg = round(Pg, 3)) |>
  write.table(file = 'output/tables/M11Output.txt', sep = ",", quote = FALSE, row.names = F)

# M12: Recruit ~ Sarg. cover (broad) ----
## Fit model ----
form <- bf(Total ~ Freq_Sarg_broad + (1|Grazing), family = poisson(link = 'log')) 
form |> get_prior(data = recruit)

priors <- prior(normal(0.5, 10), class = 'Intercept') +
  prior(normal(20, 30), class = 'b') + 
  prior(student_t(3, 0, 3), class = 'sd')

recruit.brm23 <- brm(form, prior = priors, data = recruit, 
                    sample_prior = 'only', 
                    iter = 10000, 
                    warmup = 2500, 
                    chains = 3, cores = 3, 
                    control = list(adapt_delta = 0.99, 
                                   max_treedepth = 20), 
                    thin = 5, 
                    refresh = 0, 
                    backend = 'rstan') 

recruit.brm23 |> 
  conditional_effects() |> 
  plot(points = TRUE)

recruit.brm24 <- recruit.brm23 |>
  update(sample_prior = 'yes')

## Diagnostics ----
recruit.brm24 |> 
  SUYR_prior_and_posterior()

## MCMC Sampling diagnostics
(recruit.brm24$fit |> stan_trace()) + (recruit.brm24$fit |> stan_ac()) + (recruit.brm24$fit |> stan_rhat()) + (recruit.brm24$fit |> stan_ess())

## Model validation
### Posterior probability check
recruit.brm24 |> pp_check(type = 'dens_overlay', ndraws = 100)

### Residuals
recruit.resids <- make_brms_dharma_res(recruit.brm24, integerResponse = FALSE)
testUniformity(recruit.resids)
plotResiduals(recruit.resids, form = factor(rep(1, nrow(recruit))))
plotResiduals(recruit.resids, quantreg = TRUE) 
testDispersion(recruit.resids)

## Save model ----
save(recruit.brm24, form, priors, recruit, file = 'data/modelled/M12_Freq_Sarg_broad.RData')

## Investigation ----
recruit.brm24 |> 
  conditional_effects() |> 
  plot(points = TRUE)

### Output ----
### ---- M12Output
recruit.brm24 |>
  brms::as_draws_df() |>
  mutate(across(everything(), exp)) |>
  summarise_draws(median,
                  HDInterval::hdi,
                  rhat, length, ess_bulk, ess_tail,
                  Pl = ~mean(.x < 1),
                  Pg = ~mean(.x > 1))
### ----end

recruit.brm24 |>
  brms::as_draws_df() |>
  mutate(across(everything(), exp)) |>
  summarise_draws(median,
                  HDInterval::hdi,
                  rhat, length, ess_bulk, ess_tail,
                  Pl = ~mean(.x < 1),
                  Pg = ~mean(.x > 1)) |>
  mutate(median = round(median, 3),
         lower = round(lower, 3),
         upper = round(upper, 3),
         rhat = round(rhat, 3),
         Pl = round(Pl, 3),
         Pg = round(Pg, 3)) |>
  write.table(file = 'output/tables/M12Output.txt', sep = ",", quote = FALSE, row.names = F)

### Figure ----
newdata <- with(recruit,
                list(Freq_Sarg_broad = seq(min(Freq_Sarg_broad),
                                        max(Freq_Sarg_broad),
                                        len = 50)))
fit <- recruit.brm24 |> 
  emmeans(~Freq_Sarg_broad, at = newdata) |>
  gather_emmeans_draws() |>
  mutate(.value = exp(.value)) |>
  group_by(Freq_Sarg_broad) |>
  summarise(median_hdci(.value)) |>
  as.data.frame()

M12Figure <- fit |> 
  ggplot(aes(x = Freq_Sarg_broad,
             y = y)) + 
  geom_line() + 
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.2) +
  scale_x_continuous(name = expression(paste(italic('Sargassum'), ' cover'))) +
  scale_y_continuous('Total recruits') +
  theme_classic()  +
  theme(text = element_text(colour = 'black'), 
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.5)))
M12Figure

ggsave(filename = 'output/figures/M12Figure.png', width = 8, height = 5, dpi = 100)

# M13: Recruit ~ Sarg. cover (local) ----
## Fit model ----
form <- bf(Total ~ Freq_Sarg_local + (1|Grazing), family = poisson(link = 'log')) 
form |> get_prior(data = recruit)

priors <- prior(normal(0.5, 3), class = 'Intercept') +
  prior(normal(35, 45), class = 'b') +
  prior(student_t(3, 0, 3), class = 'sd')

recruit.brm25 <- brm(form, prior = priors, data = recruit, 
                     sample_prior = 'only', 
                     iter = 10000, 
                     warmup = 2500, 
                     chains = 3, cores = 3, 
                     control = list(adapt_delta = 0.99, 
                                    max_treedepth = 20), 
                     thin = 5, 
                     refresh = 0, 
                     backend = 'rstan') 

recruit.brm25 |> 
  conditional_effects() |> 
  plot(points = TRUE)

recruit.brm26 <- recruit.brm25 |>
  update(sample_prior = 'yes')

## Diagnostics ----
recruit.brm26 |> 
  SUYR_prior_and_posterior()

## MCMC Sampling diagnostics
(recruit.brm26$fit |> stan_trace()) + (recruit.brm26$fit |> stan_ac()) + (recruit.brm26$fit |> stan_rhat()) + (recruit.brm26$fit |> stan_ess())

## Model validation
### Posterior probability check
recruit.brm26 |> pp_check(type = 'dens_overlay', ndraws = 100)

### Residuals
recruit.resids <- make_brms_dharma_res(recruit.brm26, integerResponse = FALSE)
testUniformity(recruit.resids)
plotResiduals(recruit.resids, form = factor(rep(1, nrow(recruit))))
plotResiduals(recruit.resids, quantreg = TRUE) 
testDispersion(recruit.resids)

## Save model ----
save(recruit.brm26, form, priors, recruit, file = 'data/modelled/M13_Freq_Sarg_local.RData')

## Investigation ----
recruit.brm26 |> 
  conditional_effects() |> 
  plot(points = TRUE)

### Output ----
### ---- M13Output
recruit.brm26 |>
  brms::as_draws_df() |>
  mutate(across(everything(), exp)) |>
  summarise_draws(median,
                  HDInterval::hdi,
                  rhat, length, ess_bulk, ess_tail,
                  Pl = ~mean(.x < 1),
                  Pg = ~mean(.x > 1))
### ----end

recruit.brm26 |>
  brms::as_draws_df() |>
  mutate(across(everything(), exp)) |>
  summarise_draws(median,
                  HDInterval::hdi,
                  rhat, length, ess_bulk, ess_tail,
                  Pl = ~mean(.x < 1),
                  Pg = ~mean(.x > 1)) |>
  mutate(median = round(median, 3),
         lower = round(lower, 3),
         upper = round(upper, 3),
         rhat = round(rhat, 3),
         Pl = round(Pl, 3),
         Pg = round(Pg, 3)) |>
  write.table(file = 'output/tables/M13Output.txt', sep = ",", quote = FALSE, row.names = F)

### Figure ----
newdata <- with(recruit,
                list(Freq_Sarg_local = seq(min(Freq_Sarg_local),
                                           max(Freq_Sarg_local),
                                           len = 50)))
fit <- recruit.brm26 |> 
  emmeans(~Freq_Sarg_local, at = newdata) |>
  gather_emmeans_draws() |>
  mutate(.value = exp(.value)) |>
  group_by(Freq_Sarg_local) |>
  summarise(median_hdci(.value)) |>
  as.data.frame()

M13Figure <- fit |> 
  ggplot(aes(x = Freq_Sarg_local,
             y = y)) + 
  geom_line() + 
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.2) +
  scale_x_continuous(name = expression(paste(italic('Sargassum'), ' cover'))) +
  scale_y_continuous('Total recruits') +
  theme_classic()  +
  theme(text = element_text(colour = 'black'), 
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.5)))
M13Figure

ggsave(filename = 'output/figures/M13Figure.png', width = 8, height = 5, dpi = 100)

# Compare ----
## ---- CompareM12vsM13
loo::loo_compare(loo::loo(recruit.brm24),
                 loo::loo(recruit.brm26))
## ----end


# M14: Recruit ~ Diversities ----
scatterplotMatrix(~Total+R_broad_alg+Shannon_broad_cor, 
                  data = recruit, diagonal = list(method = 'boxplot'))

## Fit model ----
form <- bf(Total ~ scale(R_broad_alg) + scale(Shannon_broad_cor) + (1|Grazing), 
           family = poisson(link = 'log')) 
form |> get_prior(data = recruit)

priors <- prior(normal(0.5, 10), class = 'Intercept') +
  prior(normal(15, 20), class = 'b') +
  prior(student_t(3, 0, 10), class = 'sd')

recruit.brm27 <- brm(form, prior = priors, data = recruit, 
                     sample_prior = 'only', 
                     iter = 10000, 
                     warmup = 2500, 
                     chains = 3, cores = 3, 
                     control = list(adapt_delta = 0.99, 
                                    max_treedepth = 30), 
                     thin = 10, 
                     refresh = 0, 
                     backend = 'rstan') 

recruit.brm27 |> 
  conditional_effects() |> 
  plot(points = TRUE)

recruit.brm28 <- recruit.brm27 |>
  update(sample_prior = 'yes')

## Diagnostics ----
recruit.brm28 |> 
  SUYR_prior_and_posterior()

## MCMC Sampling diagnostics
(recruit.brm28$fit |> stan_trace()) + (recruit.brm28$fit |> stan_ac()) + (recruit.brm28$fit |> stan_rhat()) + (recruit.brm28$fit |> stan_ess())

## Model validation
### Posterior probability check
recruit.brm28 |> pp_check(type = 'dens_overlay', ndraws = 100)

### Residuals
recruit.resids <- make_brms_dharma_res(recruit.brm28, integerResponse = FALSE)
testUniformity(recruit.resids)
plotResiduals(recruit.resids, form = factor(rep(1, nrow(recruit))))
plotResiduals(recruit.resids, quantreg = TRUE) 
testDispersion(recruit.resids)

## Save model ----
save(recruit.brm28, form, priors, recruit, file = 'data/modelled/M14_all_vars.RData')

## Investigation ----
recruit.brm28 |> 
  conditional_effects() |> 
  plot(points = TRUE)

### Output ----
### ---- M14Output
recruit.brm28 |>
  brms::as_draws_df() |>
  mutate(across(everything(), exp)) |>
  summarise_draws(median,
                  HDInterval::hdi,
                  rhat, length, ess_bulk, ess_tail,
                  Pl = ~mean(.x < 1),
                  Pg = ~mean(.x > 1))
### ----end

recruit.brm28 |>
  brms::as_draws_df() |>
  mutate(across(everything(), exp)) |>
  summarise_draws(median,
                  HDInterval::hdi,
                  rhat, length, ess_bulk, ess_tail,
                  Pl = ~mean(.x < 1),
                  Pg = ~mean(.x > 1)) |>
  mutate(median = round(median, 3),
         lower = round(lower, 3),
         upper = round(upper, 3),
         rhat = round(rhat, 3),
         Pl = round(Pl, 3),
         Pg = round(Pg, 3)) |>
  write.table(file = 'output/tables/M14Output.txt', sep = ",", quote = FALSE, row.names = F)

### Contrast ----
### ---- M14Contrast
recruit.brm28 |> 
  emmeans(~Treatment, type = 'link') |>
  pairs(reverse = TRUE) |>
  gather_emmeans_draws() |> #it's on log scale, so we need to mutate it
  mutate(Response = exp(.value)) |>
  summarise(median_hdci(Response),
            Pl = mean(Response < 1),
            Pg = mean(Response > 1))
### ----end



# Comparisons ----
loo::loo_compare(loo::loo(recruit.brm2),
                 loo::loo(recruit.brm22),
                 loo::loo(recruit.brm28))

# Thresholds ----
recruit.brm4 |> 
  emmeans(~H_mean_broad, at = with(recruit,
                                   list(H_mean_broad = seq(min(H_mean_broad),
                                                           max(H_mean_broad),
                                                           len = 50)))) |>
  gather_emmeans_draws() |>
  mutate(.value = exp(.value)) |>
  filter(.value > 0) |>
  ungroup() |>
  summarise(median_hdci(H_mean_broad),
            Pl = mean(H_mean_broad < 1),
            Pg = mean(H_mean_broad > 1))

recruit.brm8 |> 
  emmeans(~D_broad, at = with(recruit,
                                   list(D_broad = seq(min(D_broad),
                                                           max(D_broad),
                                                           len = 50)))) |>
  gather_emmeans_draws() |>
  mutate(.value = exp(.value)) |>
  filter(.value > 0) |>
  ungroup() |>
  summarise(median_hdci(D_broad),
            Pl = mean(D_broad < 1),
            Pg = mean(D_broad > 1))

recruit.brm24 |> 
  emmeans(~Freq_Sarg_broad, at = with(recruit,
                              list(Freq_Sarg_broad = seq(min(Freq_Sarg_broad),
                                                 max(Freq_Sarg_broad),
                                                 len = 50)))) |>
  gather_emmeans_draws() |>
  mutate(.value = exp(.value)) |>
  filter(.value > 0) |>
  ungroup() |>
  summarise(median_hdci(Freq_Sarg_broad),
            Pl = mean(Freq_Sarg_broad < 1),
            Pg = mean(Freq_Sarg_broad > 1))


recruit.brm14 |> 
  emmeans(~Shannon_broad_cor, at = with(recruit,
                                      list(Shannon_broad_cor = seq(min(Shannon_broad_cor),
                                                                 max(Shannon_broad_cor),
                                                                 len = 50)))) |>
  gather_emmeans_draws() |>
  mutate(.value = exp(.value)) |>
  filter(.value > 0) |>
  ungroup() |>
  summarise(median_hdci(Shannon_broad_cor),
            Pl = mean(Shannon_broad_cor < 1),
            Pg = mean(Shannon_broad_cor > 1))


# Supplementary checks ----
form <- bf(Shannon_broad_cor ~ Treatment, family = gaussian()) 
form |> get_prior(data = recruit)

priors <- prior(normal(0.2, 1), class = 'Intercept') +
  prior(normal(0, 5), class = 'b')

diversity.brm <- brm(form, prior = priors, data = recruit, 
                   sample_prior = 'only', 
                   iter = 5000, 
                   warmup = 1000, 
                   chains = 3, cores = 3, 
                   thin = 5, 
                   refresh = 0, 
                   control = list(adapt_delta = 0.99, 
                                  max_treedepth = 20),
                   backend = 'rstan') 

diversity.brm |> 
  conditional_effects() |> 
  plot(points = TRUE)

diversity.brm2 <- diversity.brm |>
  update(sample_prior = 'yes')

## Diagnostics ----
diversity.brm2 |> 
  SUYR_prior_and_posterior()

## MCMC Sampling diagnostics
(diversity.brm2$fit |> stan_trace()) + (diversity.brm2$fit |> stan_ac()) + (diversity.brm2$fit |> stan_rhat()) + (diversity.brm2$fit |> stan_ess())

## Model validation
### Posterior probability check
diversity.brm2 |> pp_check(type = 'dens_overlay', ndraws = 100)

### Residuals
diversity.resids <- make_brms_dharma_res(diversity.brm2, integerResponse = FALSE)
testUniformity(diversity.resids)
plotResiduals(diversity.resids, form = factor(rep(1, nrow(recruit))))
plotResiduals(diversity.resids, quantreg = TRUE) 
testDispersion(diversity.resids)

## Investigation ----
diversity.brm2 |> 
  conditional_effects() |> 
  plot(points = TRUE)

### Output ----
diversity.brm2 |>
  summarise_draws(median,
                  HDInterval::hdi,
                  rhat, length, ess_bulk, ess_tail,
                  Pl = ~mean(.x < 0),
                  Pg = ~mean(.x > 0)) |>
  mutate(median = round(median, 3),
         lower = round(lower, 3),
         upper = round(upper, 3),
         rhat = round(rhat, 3),
         Pl = round(Pl, 3),
         Pg = round(Pg, 3))


# M15: Recruit ~ Tile benthos ----
scatterplotMatrix(~Total+Sediment+Turf_height+Turf_cover+CCA_cover, 
                  data = recruit, diagonal = list(method = 'boxplot'))

## Fit model ----
form <- bf(Total ~ Sediment*Turf_height, 
           family = poisson(link = 'log')) 
form |> get_prior(data = recruit)

priors <- prior(normal(0.5, 20), class = 'Intercept') +
  prior(normal(15, 100), class = 'b')

recruit.brm29 <- brm(form, prior = priors, data = recruit, 
                     sample_prior = 'only', 
                     iter = 10000, 
                     warmup = 2500, 
                     chains = 3, cores = 3, 
                     control = list(adapt_delta = 0.99, 
                                    max_treedepth = 30), 
                     thin = 10, 
                     refresh = 0, 
                     backend = 'rstan') 

recruit.brm29 |> 
  conditional_effects() |> 
  plot(points = TRUE)

recruit.brm30 <- recruit.brm29 |>
  update(sample_prior = 'yes')

## Diagnostics ----
recruit.brm30 |> 
  SUYR_prior_and_posterior()

## MCMC Sampling diagnostics
(recruit.brm30$fit |> stan_trace()) + (recruit.brm30$fit |> stan_ac()) + (recruit.brm30$fit |> stan_rhat()) + (recruit.brm30$fit |> stan_ess())

## Model validation
### Posterior probability check
recruit.brm30 |> pp_check(type = 'dens_overlay', ndraws = 100)

### Residuals
recruit.resids <- make_brms_dharma_res(recruit.brm30, integerResponse = FALSE)
testUniformity(recruit.resids)
plotResiduals(recruit.resids, form = factor(rep(1, nrow(recruit))))
plotResiduals(recruit.resids, quantreg = TRUE) 
testDispersion(recruit.resids)

## Save model ----
save(recruit.brm30, form, priors, recruit, file = 'data/modelled/M15_all_vars.RData')

## Investigation ----
recruit.brm30 |> 
  conditional_effects() |> 
  plot(points = TRUE)

### Output ----
### ---- M15Output
recruit.brm30 |>
  brms::as_draws_df() |>
  mutate(across(everything(), exp)) |>
  summarise_draws(median,
                  HDInterval::hdi,
                  rhat, length, ess_bulk, ess_tail,
                  Pl = ~mean(.x < 1),
                  Pg = ~mean(.x > 1))
### ----end