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
recruit <- read_csv(file = 'data/processed/recruit.csv', col_select = -1) |>
  mutate(Sediment = factor(Sediment,
                           c('No', 'Low', 'Medium', 'High'), ordered = TRUE),
         Turf_height = factor(Turf_height,
                              c('No', 'Low', 'Medium', 'High'), ordered = TRUE),
         Treatment = ifelse(Treatment == 'No algae', 'T1',
                            ifelse(Treatment == 'Only canopy', 'T2',
                                   ifelse(Treatment == 'Only mat', 'T3', 
                                          'T4'))),
         Grazing = factor(Grazing, c('No', 'Light', 'Medium', 'Heavy'), ordered = TRUE)) |>
  mutate(Treatment = factor(Treatment))

recruit |>
  str()

# Exploratory analysis ----
recruit |> 
  ggplot(aes(y = Total, x = Treatment)) +
  geom_boxplot() +
  geom_point(color = 'black') + 
  theme_classic() +
  theme(text = element_text(colour = 'black'), 
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.5))) +
  scale_y_continuous(name = 'Total recruits', n.breaks = 3)

# scatterplot
scatterplotMatrix(~Total + H_mean_broad + H_mean_local + D_broad +
                    Shannon_broad_alg + Shannon_local_alg, 
                  data = recruit, diagonal = list(method = 'boxplot'))

# MZI1 - Treatments ----

## Fit ----

## ---- MZI1Fit
recruit.form <- bf(Total ~ Treatment + Turf_height + (1|Grazing), 
                zi ~ 1, 
                family = zero_inflated_poisson(link = 'log'))
## ----end

priors <- prior(normal(0.5, 3), class = 'Intercept') +
  prior(normal(0, 5), class = 'b') +
  prior(student_t(3, 0, 2), class = 'sd') +
  prior(logistic(0, 1), class = 'Intercept', dpar = 'zi')

recruitZI.brm <- brm(recruit.form, prior = priors, data = recruit, 
                 sample_prior = 'yes', 
                 iter = 5000, 
                 warmup = 1000, 
                 chains = 3, cores = 3, 
                 thin = 5, 
                 control = list(adapt_delta = 0.99, max_treedepth = 20),
                 refresh = 100, 
                 backend = 'rstan') 

## Diagnostics ----

recruitZI.brm |>
  conditional_effects() |>
  plot(points = TRUE)

recruitZI.brm |> 
  SUYR_prior_and_posterior()   +
  theme_classic() + theme(text = element_text(colour = 'black'), 
                          axis.text = element_text(size = rel(1.2)),
                          axis.title = element_text(size = rel(1.5)),
                          legend.text = element_text(size = rel(1.2)),
                          legend.title = element_text(size = rel(1.5)),
                          legend.position = 'bottom')

ggsave(file = paste0(FIGS_PATH, "/MZI1PriorsPosteriors.png"), 
       width = 10,
       height = 8, 
       dpi = 300)

(recruitZI.brm$fit |> stan_trace()) + (recruitZI.brm$fit |> stan_ac()) + 
(recruitZI.brm$fit |> stan_rhat()) + (recruitZI.brm$fit |> stan_ess())

ggsave(file = paste0(FIGS_PATH, "/M1ZIMCMC.png"), 
       width = 10,
       height = 8, 
       dpi = 300)

recruitZI.brm |> 
  pp_check(type = 'dens_overlay', 
           ndraws = 100)

ggsave(file = paste0(FIGS_PATH, "/MZI1PPCheck.png"), 
       width = 10,
       height = 8, 
       dpi = 300)

### Residuals
recruit.resids <- make_brms_dharma_res(recruitZI.brm, 
                                       integerResponse = TRUE)
testUniformity(recruit.resids)
plotResiduals(recruit.resids, form = factor(rep(1, nrow(recruit))))
plotResiduals(recruit.resids, quantreg = TRUE) 
testDispersion(recruit.resids)

ggsave(filename = paste0(FIGS_PATH, '/MZI1DHARMa.png'),
       wrap_elements(~testUniformity(recruit.resids)) +
         wrap_elements(~plotResiduals(recruit.resids, form = factor(rep(1, nrow(recruit))))) +
         wrap_elements(~testDispersion(recruit.resids)),
       width = 12,
       height = 4,
       dpi = 300)

## Save model ----
#save(recruitZI.brm, recruit.form, priors, recruit, file = 'data/modelled/MZI1_base.RData')

## Investigation ----
recruitZI.brm |> 
  conditional_effects() |> 
  plot(points = TRUE)

### Output ----
### ---- MZI1Output
recruitZI.brm |>
  brms::as_draws_df() |>
  dplyr::select(starts_with('b_')) |>
  mutate(across(everything(), exp)) |>
  summarise_draws(median,
                  HDInterval::hdi,
                  rhat, length, ess_bulk, ess_tail,
                  Pl = ~mean(.x < 1),
                  Pg = ~mean(.x > 1))
### ----end

recruitZI.brm |>
  brms::as_draws_df() |>
  dplyr::select(starts_with('b_')) |>
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
  write.table(file = 'output/tables/MZI1Output.txt', sep = ",", quote = FALSE, row.names = F)

### Contrast ----
### ---- MZI1Contrast
recruitZI.brm |> 
  emmeans(~Treatment, type = 'link') |>
  pairs(reverse = TRUE) |>
  gather_emmeans_draws() |> 
  mutate(Response = exp(.value)) |> #it's on log scale, so we need to mutate it
  summarise(median_hdci(Response),
            Pl = mean(Response < 1),
            Pg = mean(Response > 1))
### ----end

recruitZI.brm |> 
  emmeans(~Treatment, type = 'link') |>
  pairs(reverse = TRUE) |>
  gather_emmeans_draws() |>
  mutate(Response = exp(.value)) |>
  summarise(median_hdci(Response),
            Pl = mean(Response < 1),
            Pg = mean(Response > 1)) |>
  write.table(file = 'output/tables/MZI1Contrast.txt', sep = ",", quote = FALSE, row.names = F)

### ROPE 
recruitZI.brm |> 
  emmeans(~Treatment, type = 'link') |>
  pairs(reverse = TRUE) |>
  gather_emmeans_draws() |>
  mutate(Response = exp(.value)) |>
  summarise(median_hdci(Response),
            ROPE = rope(Response, range = c(0.9, 1.1))$ROPE_Percentage)


## Further investigation ----

### ---- ZIPlannedContrast
cmat <- cbind('Mat_Canopy' = c(1/2, -1/2, 1/2, -1/2))
recruitZI.brm |>
  emmeans(~Treatment, type = 'link') |>
  contrast(method = list(Treatment = cmat)) |>
  gather_emmeans_draws() |>
  mutate(Response = exp(.value)) |>
  summarise(median_hdci(Response),
            Pl = mean(Response < 1),
            Pg = mean(Response > 1),
            ROPE = rope(Response, range = c(0.9, 1.1))$ROPE_Percentage)
### ----end

recruitZI.brm |>
  emmeans(~Treatment, type = 'link') |>
  contrast(method = list(Treatment = cmat)) |>
  gather_emmeans_draws() |>
  mutate(Response = exp(.value)) |>
  summarise(median_hdci(Response),
            Pl = mean(Response < 1),
            Pg = mean(Response > 1)) |>
  write.table(file = 'output/tables/MZI1PlannedContrast.txt', sep = ",", quote = FALSE, row.names = F)

# in the actual number of recruits:
recruitZI.brm |> 
  emmeans(~Treatment) |>
  regrid() |> #turns values onto the response scale, so that pairs() which is for contrasts, is on the response scale
  pairs(reverse = TRUE) |>
  gather_emmeans_draws() |>
  summarise(median_hdci(.value))

## Summary figures ----

# Estimated number of recruits and their CI
recruitZI.brm |> 
  emmeans(~Treatment, type = 'response') |>
  as.tibble() |>
  ggplot(aes(x = Treatment, y = rate)) +
  geom_pointrange(aes(ymin = lower.HPD, ymax = upper.HPD)) +
  theme_classic() +
  ylab('Coral recruitment (recruits/tile)') +
  geom_point(data = recruit, aes(x = Treatment, y = Total), alpha = 0.4,
             position = position_jitter(width = 0.1, height = 0.1))

ggsave(file = paste0(FIGS_PATH, "/MZI1_Treat_fig.png"), 
       width = 160, 
       height = 160/1.6, 
       units = "mm", 
       dpi = 300)

# Contrasts figure
recruitZI.brm |> 
  emmeans(~Treatment, type = 'link') |>
  pairs() |>
  gather_emmeans_draws() |>
  mutate(Fit = exp(.value)) |> 
  ggplot() +
  stat_slab(aes(x = Fit, y = contrast,
                fill = stat(ggdist::cut_cdf_qi(cdf,
                                               .width = c(0.5, 0.8, 0.95),
                                               labels = scales::percent_format()))), 
            color = 'black', linewidth = 0.3) +
  geom_vline(xintercept = 1, linetype = 'dashed') +
  scale_fill_brewer('Interval', 
                    direction  =  -1, 
                    na.translate = FALSE) +
  scale_x_continuous('Effect', trans = scales::log2_trans(), labels = scales::comma) +
  theme_classic() +
  theme(legend.position = 'bottom')

ggsave(file = paste0(FIGS_PATH, "/MZI1_Contrast_fig.png"), 
       width = 160, 
       height = 160/1.6, 
       units = "mm", 
       dpi = 300)


# MZI2: Recruit ~ H (broad) ----
recruit |> dplyr::filter(H_mean_broad != 0) |> dplyr::select(H_mean_broad) |> min()/2

## Fit model ----

## ---- MZI2Fit
recruit.form <- bf(Total ~ log(H_mean_broad + 2.2) + (1|Grazing), 
           zi ~ 1, 
           family = zero_inflated_poisson(link = 'log'))
## ----end


recruit.form |> get_prior(data = recruit)

priors <- prior(normal(0.5, 3), class = 'Intercept') +
  prior(normal(0, 5), class = 'b')  +
  prior(student_t(3, 0, 2), class = 'sd') +
  prior(logistic(0, 1), class = 'Intercept', dpar = 'zi')

recruitZI.brm2 <- brm(recruit.form, prior = priors, data = recruit, 
                    sample_prior = 'yes', 
                    iter = 5000, 
                    warmup = 1000, 
                    chains = 3, cores = 3, 
                    control = list(adapt_delta = 0.99, 
                                   max_treedepth = 20), 
                    thin = 5, 
                    refresh = 100, 
                    backend = 'rstan') 

## Diagnostics ----
recruitZI.brm2 |> 
  SUYR_prior_and_posterior() +
  theme_classic() + theme(text = element_text(colour = 'black'), 
                          axis.text = element_text(size = rel(1.2)),
                          axis.title = element_text(size = rel(1.5)),
                          legend.text = element_text(size = rel(1.2)),
                          legend.title = element_text(size = rel(1.5)),
                          legend.position = 'bottom')

ggsave(file = paste0(FIGS_PATH, "/MZI2PriorsPosteriors.png"), 
       width = 10,
       height = 8, 
       dpi = 300)

## MCMC Sampling diagnostics
(recruitZI.brm2$fit |> stan_trace()) + (recruitZI.brm2$fit |> stan_ac()) + 
(recruitZI.brm2$fit |> stan_rhat()) + (recruitZI.brm2$fit |> stan_ess())

ggsave(file = paste0(FIGS_PATH, "/M2ZIMCMC.png"), 
       width = 10,
       height = 8, 
       dpi = 300)

## Model validation
### Posterior probability check
recruitZI.brm2 |> 
  pp_check(type = 'dens_overlay', ndraws = 100)

ggsave(file = paste0(FIGS_PATH, "/MZI2PPCheck.png"), 
       width = 10,
       height = 8, 
       dpi = 300)

### Residuals
recruit.resids <- make_brms_dharma_res(recruitZI.brm2, 
                                       integerResponse = TRUE)
testUniformity(recruit.resids)
plotResiduals(recruit.resids, form = factor(rep(1, nrow(recruit))))
plotResiduals(recruit.resids, quantreg = TRUE) 
testDispersion(recruit.resids)

ggsave(filename = paste0(FIGS_PATH, '/MZI2DHARMa.png'),
       wrap_elements(~testUniformity(recruit.resids)) +
         wrap_elements(~plotResiduals(recruit.resids, form = factor(rep(1, nrow(recruit))))) +
         wrap_elements(~testDispersion(recruit.resids)),
       width = 12,
       height = 4,
       dpi = 300)

## Save model ----
#save(recruitZI.brm2, recruit.form, priors, recruit, file = 'data/modelled/MZI2_Height_broad.RData')

## Investigation ----
recruitZI.brm2 |> 
  conditional_effects() |> 
  plot(points = TRUE)

### Output ----
### ---- MZI2Output
recruitZI.brm2 |>
  brms::as_draws_df() |>
  dplyr::select(starts_with('b_')) |>
  mutate(across(everything(), exp)) |>
  summarise_draws(median,
                  HDInterval::hdi,
                  rhat, length, ess_bulk, ess_tail,
                  Pl = ~mean(.x < 1),
                  Pg = ~mean(.x > 1))
### ----end

recruitZI.brm2 |>
  brms::as_draws_df() |>
  dplyr::select(starts_with('b_')) |>
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
  write.table(file = 'output/tables/MZI2Output.txt', sep = ",", quote = FALSE, row.names = F)

### Figure ----
newdata <- with(recruit,
                list(H_mean_broad = seq(min(H_mean_broad),
                                        max(H_mean_broad),
                                        len = 50)))
fit <- recruitZI.brm2 |> 
  emmeans(~H_mean_broad, at = newdata) |>
  gather_emmeans_draws() |>
  mutate(.value = exp(.value)) |>
  group_by(H_mean_broad) |>
  summarise(median_hdci(.value)) |>
  as.data.frame()

MZI2Figure <- fit |> 
  ggplot(aes(x = H_mean_broad,
             y = y)) + 
  geom_line() + 
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.2) +
  geom_point(data = recruit, aes(x = H_mean_broad, y = Total), alpha = 0.4, size = 2, 
             position = position_jitter(width = 0.1, height = 0.1)) +
  scale_x_continuous(name = expression(paste(italic('Sargassum'), ' height (cm)'))) +
  scale_y_continuous('Coral recruitment (recruits/tile)') +
  theme_classic()  +
  theme(text = element_text(colour = 'black'), 
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.5)))
MZI2Figure


# MZI3: Recruit ~ H (local) ----
recruit |> dplyr::filter(H_mean_local != 0) |> dplyr::select(H_mean_local) |> min()/2

## Fit model ----
## ---- MZI3Fit
recruit.form <- bf(Total ~ log(H_mean_local + 2) + (1|Grazing), 
           zi ~ 1, 
           family = zero_inflated_poisson(link = 'log'))
## ----end

recruit.form |> get_prior(data = recruit)

priors <- prior(normal(0.5, 3), class = 'Intercept') +
  prior(normal(0, 5), class = 'b') + 
  prior(student_t(3, 0, 2), class = 'sd') +
  prior(logistic(0, 1), class = 'Intercept', dpar = 'zi')

recruitZI.brm3 <- brm(recruit.form, prior = priors, data = recruit, 
                      sample_prior = 'yes', 
                      iter = 5000, 
                      warmup = 1000, 
                      chains = 3, cores = 3, 
                      control = list(adapt_delta = 0.99, 
                                     max_treedepth = 20), 
                      thin = 5, 
                      refresh = 0, 
                      backend = 'rstan') 

## Diagnostics ----
recruitZI.brm3 |> 
  SUYR_prior_and_posterior()  +
  theme_classic() + theme(text = element_text(colour = 'black'), 
                          axis.text = element_text(size = rel(1.2)),
                          axis.title = element_text(size = rel(1.5)),
                          legend.text = element_text(size = rel(1.2)),
                          legend.title = element_text(size = rel(1.5)),
                          legend.position = 'bottom')

ggsave(file = paste0(FIGS_PATH, "/MZI3PriorsPosteriors.png"), 
       width = 10,
       height = 8, 
       dpi = 300)

## MCMC Sampling diagnostics
(recruitZI.brm3$fit |> stan_trace()) + (recruitZI.brm3$fit |> stan_ac()) + 
  (recruitZI.brm3$fit |> stan_rhat()) + (recruitZI.brm3$fit |> stan_ess())

ggsave(file = paste0(FIGS_PATH, "/M3ZIMCMC.png"), 
       width = 10,
       height = 8, 
       dpi = 300)

## Model validation
### Posterior probability check
recruitZI.brm3 |> 
  pp_check(type = 'dens_overlay', ndraws = 100)

ggsave(file = paste0(FIGS_PATH, "/MZI3PPCheck.png"), 
       width = 10,
       height = 8, 
       dpi = 300)

### Residuals
recruit.resids <- make_brms_dharma_res(recruitZI.brm3, integerResponse = FALSE)
testUniformity(recruit.resids)
plotResiduals(recruit.resids, form = factor(rep(1, nrow(recruit))))
plotResiduals(recruit.resids, quantreg = TRUE) 
testDispersion(recruit.resids)

ggsave(filename = paste0(FIGS_PATH, '/MZI3DHARMa.png'),
       wrap_elements(~testUniformity(recruit.resids)) +
         wrap_elements(~plotResiduals(recruit.resids, form = factor(rep(1, nrow(recruit))))) +
         wrap_elements(~testDispersion(recruit.resids)),
       width = 12,
       height = 4,
       dpi = 300)

## Save model ----
#save(recruitZI.brm3, recruit.form, priors, recruit, file = 'data/modelled/MZI3_Height_local.RData')

## Investigation ----
recruitZI.brm3 |> 
  conditional_effects() |> 
  plot(points = TRUE)

### Output ----
### ---- MZI3Output
recruitZI.brm3 |>
  brms::as_draws_df() |>
  dplyr::select(starts_with('b_')) |>
  mutate(across(everything(), exp)) |>
  summarise_draws(median,
                  HDInterval::hdi,
                  rhat, length, ess_bulk, ess_tail,
                  Pl = ~mean(.x < 1),
                  Pg = ~mean(.x > 1))
### ----end

recruitZI.brm3 |>
  brms::as_draws_df() |>
  dplyr::select(starts_with('b_')) |>
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
  write.table(file = 'output/tables/MZI3Output.txt', sep = ",", quote = FALSE, row.names = F)

### Figure ----
newdata <- with(recruit,
                list(H_mean_local = seq(min(H_mean_local),
                                        max(H_mean_local),
                                        len = 50)))
fit <- recruitZI.brm3 |> 
  emmeans(~H_mean_local, at = newdata) |>
  gather_emmeans_draws() |>
  mutate(.value = exp(.value)) |>
  group_by(H_mean_local) |>
  summarise(median_hdci(.value)) |>
  as.data.frame()

MZI3Figure <- fit |> 
  ggplot(aes(x = H_mean_local,
             y = y)) + 
  geom_line() + 
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.2) +
  geom_point(data = recruit, aes(x = H_mean_local, y = Total), alpha = 0.4, size = 2, 
             position = position_jitter(width = 0.1, height = 0.1))  +
  scale_x_continuous(name = expression(paste(italic('Sargassum'), ' height (cm)'))) +
  scale_y_continuous('Coral recruitment (recruits/tile)') +
  labs(tag = 'a') +
  theme_classic()  +
  theme(text = element_text(colour = 'black'), 
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.5)))
MZI3Figure

ggsave(file = paste0(FIGS_PATH, "/MZI3_H_local_fig.png"), 
       width = 160, 
       height = 160/1.6, 
       units = "mm", 
       dpi = 300)

# Compare ----
## ---- CompareMZI2vsMZI3
loo::loo_compare(loo::loo(recruitZI.brm2),
                 loo::loo(recruitZI.brm3))
## ----end


# MZI4: Recruit ~ Density ----
recruit |> dplyr::filter(D_broad != 0) |> dplyr::select(D_broad) |> min()/2

## Fit model ----
## ---- MZI4Fit
recruit.form <- bf(Total ~ log(D_broad + 0.3) + (1|Grazing), 
                   zi ~ 1, 
                   family = zero_inflated_poisson(link = 'log'))
## ----end
recruit.form |> get_prior(data = recruit)

priors <- prior(normal(0.5, 3), class = 'Intercept') +
  prior(normal(0, 5), class = 'b') +
  prior(student_t(3, 0, 2), class = 'sd')  +
  prior(logistic(0, 1), class = 'Intercept', dpar = 'zi')

recruitZI.brm4 <- brm(recruit.form, prior = priors, data = recruit, 
                    sample_prior = 'yes', 
                    iter = 5000, 
                    warmup = 1000, 
                    chains = 3, cores = 3, 
                    control = list(adapt_delta = 0.99, 
                                   max_treedepth = 20), 
                    thin = 5, 
                    refresh = 0, 
                    backend = 'rstan') 

## Diagnostics ----
recruitZI.brm4 |> 
  SUYR_prior_and_posterior()   +
  theme_classic() + theme(text = element_text(colour = 'black'), 
                          axis.text = element_text(size = rel(1.2)),
                          axis.title = element_text(size = rel(1.5)),
                          legend.text = element_text(size = rel(1.2)),
                          legend.title = element_text(size = rel(1.5)),
                          legend.position = 'bottom')

ggsave(file = paste0(FIGS_PATH, "/MZI4PriorsPosteriors.png"), 
       width = 10,
       height = 8, 
       dpi = 300)

## MCMC Sampling diagnostics
(recruitZI.brm4$fit |> stan_trace()) + (recruitZI.brm4$fit |> stan_ac()) + 
  (recruitZI.brm4$fit |> stan_rhat()) + (recruitZI.brm4$fit |> stan_ess())

ggsave(file = paste0(FIGS_PATH, "/M4ZIMCMC.png"), 
       width = 10,
       height = 8, 
       dpi = 300)

## Model validation
### Posterior probability check
recruitZI.brm4 |> 
  pp_check(type = 'dens_overlay', ndraws = 100)

ggsave(file = paste0(FIGS_PATH, "/MZI4PPCheck.png"), 
       width = 10,
       height = 8, 
       dpi = 300)

### Residuals
recruit.resids <- make_brms_dharma_res(recruitZI.brm4, 
                                       integerResponse = TRUE)
testUniformity(recruit.resids)
plotResiduals(recruit.resids, form = factor(rep(1, nrow(recruit))))
plotResiduals(recruit.resids, quantreg = TRUE) 
testDispersion(recruit.resids)

ggsave(filename = paste0(FIGS_PATH, '/MZI4DHARMa.png'),
       wrap_elements(~testUniformity(recruit.resids)) +
         wrap_elements(~plotResiduals(recruit.resids, form = factor(rep(1, nrow(recruit))))) +
         wrap_elements(~testDispersion(recruit.resids)),
       width = 12,
       height = 4,
       dpi = 300)

## Save model ----
#save(recruitZI.brm4, recruit.form, priors, recruit, file = 'data/modelled/MZI4_Density.RData')

## Investigation ----
recruitZI.brm4 |> 
  conditional_effects() |> 
  plot(points = TRUE)

### Output ----
### ---- MZI4Output
recruitZI.brm4 |>
  brms::as_draws_df() |>
  dplyr::select(starts_with('b_')) |>
  mutate(across(everything(), exp)) |>
  summarise_draws(median,
                  HDInterval::hdi,
                  rhat, length, ess_bulk, ess_tail,
                  Pl = ~mean(.x < 1),
                  Pg = ~mean(.x > 1))
### ----end

recruitZI.brm4 |>
  brms::as_draws_df() |>
  dplyr::select(starts_with('b_')) |>
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
  write.table(file = 'output/tables/MZI4Output.txt', sep = ",", quote = FALSE, row.names = F)

### Figure ----
newdata <- with(recruit,
                list(D_broad = seq(min(D_broad),
                                   max(D_broad),
                                   len = 50)))
fit <- recruitZI.brm4 |> 
  emmeans(~D_broad, at = newdata) |>
  gather_emmeans_draws() |>
  mutate(.value = exp(.value)) |>
  group_by(D_broad) |>
  summarise(median_hdci(.value)) |>
  as.data.frame()

MZI4Figure <- fit |> 
  ggplot(aes(x = D_broad*100,
             y = y)) + 
  geom_line() + 
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.2) +
  geom_point(data = recruit, aes(x = D_broad*100, y = Total), alpha = 0.4, size = 2,
             position = position_jitter(width = 0.1, height = 0.1))   +
  scale_x_continuous(name = expression(paste(italic('Sargassum'), ' density (thalli/', m^{2}, ')'))) +
  scale_y_continuous('Coral recruitment (recruits/tile)') +
  labs(tag = 'b') +
  theme_classic()  +
  theme(text = element_text(colour = 'black'), 
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.5)))
MZI4Figure

ggsave(file = paste0(FIGS_PATH, "/MZI4_Density_fig.png"), 
       width = 160, 
       height = 160/1.6, 
       units = "mm", 
       dpi = 300)

# Summary fig (H + D) ----
MZI3Figure / MZI4Figure +
  plot_layout(axis_titles = 'collect')

ggsave(file = paste0(FIGS_PATH, "/MZI_H_D_fig.png"), 
       width = 200/1.2, 
       height = 200, 
       units = "mm", 
       dpi = 300)

# Compare H vs D ----
## ---- CompareMZI3vsMZI4
loo::loo_compare(loo::loo(recruitZI.brm2),
                 loo::loo(recruitZI.brm4))
## ----end

# Thresholds ----
recruitZI.brm2 |> 
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

recruitZI.brm3 |> 
  emmeans(~H_mean_local, at = with(recruit,
                                   list(H_mean_local = seq(min(H_mean_local),
                                                           max(H_mean_local),
                                                           len = 50)))) |>
  gather_emmeans_draws() |>
  mutate(.value = exp(.value)) |>
  filter(.value > 0) |>
  ungroup() |>
  summarise(median_hdci(H_mean_local),
            Pl = mean(H_mean_local < 1),
            Pg = mean(H_mean_local > 1))

recruitZI.brm4 |> 
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

# MZI5: Recruit ~ Diversity (local) ----
recruit |> dplyr::filter(Shannon_local_alg != 0) |> dplyr::select(Shannon_local_alg) |> min()/2

## Fit model ----

## ---- MZI5Fit
recruit.form <- bf(Total ~ log(Shannon_local_alg + 0.01) + 
                     (1|Grazing), 
           zi ~ 1, 
           family = zero_inflated_poisson(link = 'log'))
## ----end
recruit.form |> get_prior(data = recruit)

priors <- prior(normal(0.5, 3), class = 'Intercept') +
  prior(normal(0, 5), class = 'b') +
  prior(student_t(3, 0, 2), class = 'sd') +
  prior(logistic(0, 1), class = 'Intercept', dpar = 'zi')

recruitZI.brm5 <- brm(recruit.form, prior = priors, data = recruit, 
                    sample_prior = 'yes', 
                    iter = 5000, 
                    warmup = 1000, 
                    chains = 3, cores = 3, 
                    control = list(adapt_delta = 0.99, 
                                   max_treedepth = 20), 
                    thin = 5, 
                    refresh = 0, 
                    backend = 'rstan') 

## Diagnostics ----
recruitZI.brm5 |> 
  SUYR_prior_and_posterior()   +
  theme_classic() + theme(text = element_text(colour = 'black'), 
                          axis.text = element_text(size = rel(1.2)),
                          axis.title = element_text(size = rel(1.5)),
                          legend.text = element_text(size = rel(1.2)),
                          legend.title = element_text(size = rel(1.5)),
                          legend.position = 'bottom')

ggsave(file = paste0(FIGS_PATH, "/MZI5PriorsPosteriors.png"), 
       width = 10,
       height = 8, 
       dpi = 300)

## MCMC Sampling diagnostics
(recruitZI.brm5$fit |> stan_trace()) + (recruitZI.brm5$fit |> stan_ac()) + 
  (recruitZI.brm5$fit |> stan_rhat()) + (recruitZI.brm5$fit |> stan_ess())

ggsave(file = paste0(FIGS_PATH, "/M5ZIMCMC.png"), 
       width = 10,
       height = 8, 
       dpi = 300)

## Model validation
### Posterior probability check
recruitZI.brm5 |> 
  pp_check(type = 'dens_overlay', ndraws = 100)

ggsave(file = paste0(FIGS_PATH, "/MZI5PPCheck.png"), 
       width = 10,
       height = 8, 
       dpi = 300)

### Residuals
recruit.resids <- make_brms_dharma_res(recruitZI.brm5, 
                                       integerResponse = TRUE)
testUniformity(recruit.resids)
plotResiduals(recruit.resids, form = factor(rep(1, nrow(recruit))))
plotResiduals(recruit.resids, quantreg = TRUE) 
testDispersion(recruit.resids)

ggsave(filename = paste0(FIGS_PATH, '/MZI5DHARMa.png'),
       wrap_elements(~testUniformity(recruit.resids)) +
         wrap_elements(~plotResiduals(recruit.resids, form = factor(rep(1, nrow(recruit))))) +
         wrap_elements(~testDispersion(recruit.resids)),
       width = 12,
       height = 4,
       dpi = 300)

## Save model ----
#save(recruitZI.brm5, recruit.form, priors, recruit, file = 'data/modelled/MZI5_Local_diversity.RData')

## Investigation ----
recruitZI.brm5 |> 
  conditional_effects() |> 
  plot(points = TRUE)

### Output ----
### ---- MZI5Output
recruitZI.brm5 |>
  brms::as_draws_df() |>
  dplyr::select(starts_with('b_')) |>
  mutate(across(everything(), exp)) |>
  summarise_draws(median,
                  HDInterval::hdi,
                  rhat, length, ess_bulk, ess_tail,
                  Pl = ~mean(.x < 1),
                  Pg = ~mean(.x > 1))
### ----end

recruitZI.brm5 |>
  brms::as_draws_df() |>
  dplyr::select(starts_with('b_')) |>
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
  write.table(file = 'output/tables/MZI5Output.txt', sep = ",", quote = FALSE, row.names = F)

# MZI6: Recruit ~ Diversity (broad) ----
recruit |> dplyr::filter(Shannon_broad_alg != 0) |> dplyr::select(Shannon_broad_alg) |> min()/2

## Fit model ----
## ---- MZI16Fit
recruit.form <- bf(Total ~ log(Shannon_broad_alg + 0.03) + 
                     (1|Grazing), 
           zi ~ 1, 
           family = zero_inflated_poisson(link = 'log'))
## ----end
recruit.form |> get_prior(data = recruit)

priors <- prior(normal(0.5, 3), class = 'Intercept') +
  prior(normal(0, 5), class = 'b') +
  prior(student_t(3, 0, 2), class = 'sd') +
  prior(logistic(0, 1), class = 'Intercept', dpar = 'zi')

recruitZI.brm6 <- brm(recruit.form, prior = priors, data = recruit, 
                     sample_prior = 'yes', 
                     iter = 5000, 
                     warmup = 1000, 
                     chains = 3, cores = 3, 
                     control = list(adapt_delta = 0.99, 
                                    max_treedepth = 20), 
                     thin = 5, 
                     refresh = 0, 
                     backend = 'rstan') 

## Diagnostics ----
recruitZI.brm6 |> 
  SUYR_prior_and_posterior()   +
  theme_classic() + theme(text = element_text(colour = 'black'), 
                          axis.text = element_text(size = rel(1.2)),
                          axis.title = element_text(size = rel(1.5)),
                          legend.text = element_text(size = rel(1.2)),
                          legend.title = element_text(size = rel(1.5)),
                          legend.position = 'bottom')

ggsave(file = paste0(FIGS_PATH, "/MZI6PriorsPosteriors.png"), 
       width = 10,
       height = 8, 
       dpi = 300)

## MCMC Sampling diagnostics
(recruitZI.brm6$fit |> stan_trace()) + (recruitZI.brm6$fit |> stan_ac()) + 
  (recruitZI.brm6$fit |> stan_rhat()) + (recruitZI.brm6$fit |> stan_ess())

ggsave(file = paste0(FIGS_PATH, "/M6ZIMCMC.png"), 
       width = 10,
       height = 8, 
       dpi = 300)

## Model validation
### Posterior probability check
recruitZI.brm6 |> 
  pp_check(type = 'dens_overlay', ndraws = 100)

ggsave(file = paste0(FIGS_PATH, "/MZI6PPCheck.png"), 
       width = 10,
       height = 8, 
       dpi = 300)

### Residuals
recruit.resids <- make_brms_dharma_res(recruitZI.brm6, 
                                       integerResponse = TRUE)
testUniformity(recruit.resids)
plotResiduals(recruit.resids, form = factor(rep(1, nrow(recruit))))
plotResiduals(recruit.resids, quantreg = TRUE) 
testDispersion(recruit.resids)

ggsave(filename = paste0(FIGS_PATH, '/MZI6DHARMa.png'),
       wrap_elements(~testUniformity(recruit.resids)) +
         wrap_elements(~plotResiduals(recruit.resids, form = factor(rep(1, nrow(recruit))))) +
         wrap_elements(~testDispersion(recruit.resids)),
       width = 12,
       height = 4,
       dpi = 300)

## Save model ----
#save(recruitZI.brm6, recruit.form, priors, recruit, file = 'data/modelled/MZI6_Broad_diversity.RData')

## Investigation ----
recruitZI.brm6 |> 
  conditional_effects() |> 
  plot(points = TRUE)

### Output ----
### ---- MZI6Output
recruitZI.brm6 |>
  brms::as_draws_df() |>
  dplyr::select(starts_with('b_')) |>
  mutate(across(everything(), exp)) |>
  summarise_draws(median,
                  HDInterval::hdi,
                  rhat, length, ess_bulk, ess_tail,
                  Pl = ~mean(.x < 1),
                  Pg = ~mean(.x > 1))
### ----end

recruitZI.brm6 |>
  brms::as_draws_df() |>
  dplyr::select(starts_with('b_')) |>
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
  write.table(file = 'output/tables/MZI6Output.txt', sep = ",", quote = FALSE, row.names = F)

# Compare ----
## ---- CompareMZI5vsMZI6
loo::loo_compare(loo::loo(recruitZI.brm5),
                 loo::loo(recruitZI.brm6))
## ----end


loo::loo_compare(loo::loo(recruitZI.brm2),
                 loo::loo(recruitZI.brm3),
                 loo::loo(recruitZI.brm4),
                 loo::loo(recruitZI.brm5),
                 loo::loo(recruitZI.brm6))
