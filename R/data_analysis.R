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
  mutate(Treatment = factor(Treatment)) |>
  dplyr::filter(!is.na(Turf_height))

recruit |>
  str()

# Exploratory analysis ----
recruit |> 
  ggplot(aes(y = Total, x = Treatment)) +
  geom_violin() +
  geom_point(color = 'black') + 
  theme_classic() +
  theme(text = element_text(colour = 'black'), 
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.5))) +
  scale_y_continuous(name = 'Total recruits', n.breaks = 3)

# correlation plot
recruit |> 
  dplyr::select(H_mean_broad, Freq_Sarg_broad, H_mean_local, Freq_Sarg_local, D_broad, 
                R_broad_alg, R_broad_coral, R_local_alg, R_local_coral, 
                Shannon_broad_alg, Shannon_broad_cor, Shannon_local_alg, Shannon_local_cor, 
                Tot_broad_alg, Tot_broad_coral, Tot_local_alg, Tot_local_coral) |> 
  cor() |> 
  corrplot()

# M1: Recruit ~ Treatment ----
## Fit model ----
form <- bf(Total ~ Treatment + (1|Turf_height), family = poisson(link = 'log')) 
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
  SUYR_prior_and_posterior() +
  theme_classic() + theme(text = element_text(colour = 'black'), 
                          axis.text = element_text(size = rel(1.2)),
                          axis.title = element_text(size = rel(1.5)),
                          legend.text = element_text(size = rel(1.2)),
                          legend.title = element_text(size = rel(1.5)),
                          legend.position = 'bottom')

ggsave(file = paste0(FIGS_PATH, "/M1PriorsPosteriors.png"), 
       width = 10,
       height = 8, 
       dpi = 300)

## MCMC Sampling diagnostics
(recruit.brm2$fit |> stan_trace()) + (recruit.brm2$fit |> stan_ac()) + (recruit.brm2$fit |> stan_rhat()) + (recruit.brm2$fit |> stan_ess())

ggsave(file = paste0(FIGS_PATH, "/M1MCMC.png"), 
       width = 10,
       height = 8, 
       dpi = 300)

## Model validation
### Posterior probability check
recruit.brm2 |> pp_check(type = 'dens_overlay', ndraws = 100)

ggsave(file = paste0(FIGS_PATH, "/M1PPCheck.png"), 
       width = 10,
       height = 8, 
       dpi = 300)

### Residuals
recruit.resids <- make_brms_dharma_res(recruit.brm2, integerResponse = FALSE)
testUniformity(recruit.resids)
plotResiduals(recruit.resids, form = factor(rep(1, nrow(recruit))))
plotResiduals(recruit.resids, quantreg = TRUE) 
testDispersion(recruit.resids)

ggsave(filename = paste0(FIGS_PATH, '/M1DHARMa.png'),
       wrap_elements(~testUniformity(recruit.resids)) +
         wrap_elements(~plotResiduals(recruit.resids, form = factor(rep(1, nrow(recruit))))) +
         wrap_elements(~testDispersion(recruit.resids)),
       width = 12,
       height = 4,
       dpi = 300)

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
recruit.brm2 |> 
  emmeans(~Treatment, type = 'link') |>
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

## Further investigation ----

### ---- PlannedContrast
cmat <- cbind('Mat_Canopy' = c(1/2, -1/2, 1/2, -1/2),
              'NoAlg_AnyAlg' = c(1, -1/3, -1/3, -1/3),
              'Mat_AnyAlg' = c(0, -1/2, 1, -1/2),
              'NoAlg_Canopy' = c(1, -1/2, 0, -1/2))
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

# in the actual number of recruits:
recruit.brm2 |> 
  emmeans(~Treatment) |>
  regrid() |> #turns values onto the response scale, so that pairs() which is for contrasts, is on the response scale
  pairs(reverse = TRUE) |>
  gather_emmeans_draws() |>
  summarise(median_hdci(.value))

## Summary figures ----

# Estimated number of recruits and their CI
recruit.brm2 |> 
  emmeans(~Treatment, type = 'response') |>
  as.tibble() |>
  ggplot(aes(x = Treatment, y = rate)) +
  geom_pointrange(aes(ymin = lower.HPD, ymax = upper.HPD)) +
  theme_classic() +
  ylab('Number of coral recruits') +
  geom_point(data = recruit, aes(x = Treatment, y = Total), alpha = 0.4,
             position = position_jitter(width = 0.1, height = 0.1))

ggsave(file = paste0(FIGS_PATH, "/M1_Treat_fig.png"), 
       width = 160, 
       height = 160/1.6, 
       units = "mm", 
       dpi = 300)
  
# Planned contrasts figure
recruit.brm2 |>
  emmeans(~Treatment, type = 'link') |>
  contrast(method = list(Treatment = cmat)) |>
  gather_emmeans_draws() |>
  mutate(Response = exp(.value)) |> 
  ggplot(aes(x = Response)) +
  geom_density_ridges_gradient(aes(y = contrast, fill = stat(x)),
                               alpha = 0.4, color = 'white', 
                               quantile_lines = TRUE,
                               quantiles = c(0.025, 0.975)) + 
  geom_vline(xintercept = 1, linetype = 'dashed') + 
  scale_fill_viridis_c(option = 'C') + 
  scale_x_continuous('', trans = scales::log2_trans())  + 
  scale_y_discrete('') +
  theme_classic() +
  theme(text = element_text(colour = 'black'), 
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.5)),
        legend.position = 'none')

ggsave(file = paste0(FIGS_PATH, "/M1_Planned_contrast_fig.png"), 
       width = 160, 
       height = 160/1.6, 
       units = "mm", 
       dpi = 300)


# M2: Recruit ~ H (broad) ----
## Fit model ----
form <- bf(Total ~ H_mean_broad + (1|Turf_height), 
           family = poisson(link = 'log')) 
form |> get_prior(data = recruit)

priors <- prior(normal(0.5, 5), class = 'Intercept') +
  prior(normal(5, 6), class = 'b')  +
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
  SUYR_prior_and_posterior() +
  theme_classic() + theme(text = element_text(colour = 'black'), 
                          axis.text = element_text(size = rel(1.2)),
                          axis.title = element_text(size = rel(1.5)),
                          legend.text = element_text(size = rel(1.2)),
                          legend.title = element_text(size = rel(1.5)),
                          legend.position = 'bottom')

ggsave(file = paste0(FIGS_PATH, "/M2PriorsPosteriors.png"), 
       width = 10,
       height = 8, 
       dpi = 300)

## MCMC Sampling diagnostics
(recruit.brm4$fit |> stan_trace()) + (recruit.brm4$fit |> stan_ac()) + (recruit.brm4$fit |> stan_rhat()) + (recruit.brm4$fit |> stan_ess())
ggsave(file = paste0(FIGS_PATH, "/M2MCMC.png"), 
       width = 10,
       height = 8, 
       dpi = 300)

## Model validation
### Posterior probability check
recruit.brm4 |> 
  pp_check(type = 'dens_overlay', ndraws = 100)

ggsave(file = paste0(FIGS_PATH, "/M2PPCheck.png"), 
       width = 10,
       height = 8, 
       dpi = 300)

### Residuals
recruit.resids <- make_brms_dharma_res(recruit.brm4, integerResponse = FALSE)
testUniformity(recruit.resids)
plotResiduals(recruit.resids, form = factor(rep(1, nrow(recruit))))
plotResiduals(recruit.resids, quantreg = TRUE) 
testDispersion(recruit.resids)

ggsave(filename = paste0(FIGS_PATH, '/M2DHARMa.png'),
       wrap_elements(~testUniformity(recruit.resids)) +
         wrap_elements(~plotResiduals(recruit.resids, form = factor(rep(1, nrow(recruit))))) +
         wrap_elements(~testDispersion(recruit.resids)),
       width = 12,
       height = 4,
       dpi = 300)

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
  geom_point(data = recruit, aes(x = H_mean_broad, y = Total), alpha = 0.4,
             position = position_jitter(width = 0.1, height = 0.1)) +
  scale_x_continuous(name = expression(paste(italic('Sargassum'), ' height (cm)'))) +
  scale_y_continuous('Number of coral recruits') +
  theme_classic()  +
  theme(text = element_text(colour = 'black'), 
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.5)))
M2Figure

# M3: Recruit ~ H (local) ----
## Fit model ----
form <- bf(Total ~ H_mean_local + (1|Turf_height), 
           family = poisson(link = 'log'))  
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
  SUYR_prior_and_posterior()  +
  theme_classic() + theme(text = element_text(colour = 'black'), 
                          axis.text = element_text(size = rel(1.2)),
                          axis.title = element_text(size = rel(1.5)),
                          legend.text = element_text(size = rel(1.2)),
                          legend.title = element_text(size = rel(1.5)),
                          legend.position = 'bottom')

ggsave(file = paste0(FIGS_PATH, "/M3PriorsPosteriors.png"), 
       width = 10,
       height = 8, 
       dpi = 300)

## MCMC Sampling diagnostics
(recruit.brm6$fit |> stan_trace()) + (recruit.brm6$fit |> stan_ac()) + (recruit.brm6$fit |> stan_rhat()) + (recruit.brm6$fit |> stan_ess())
ggsave(file = paste0(FIGS_PATH, "/M3MCMC.png"), 
       width = 10,
       height = 8, 
       dpi = 300)


## Model validation
### Posterior probability check
recruit.brm6 |> 
  pp_check(type = 'dens_overlay', ndraws = 100)

ggsave(file = paste0(FIGS_PATH, "/M3PPCheck.png"), 
       width = 10,
       height = 8, 
       dpi = 300)

### Residuals
recruit.resids <- make_brms_dharma_res(recruit.brm6, integerResponse = FALSE)
testUniformity(recruit.resids)
plotResiduals(recruit.resids, form = factor(rep(1, nrow(recruit))))
plotResiduals(recruit.resids, quantreg = TRUE) 
testDispersion(recruit.resids)

ggsave(filename = paste0(FIGS_PATH, '/M3DHARMa.png'),
       wrap_elements(~testUniformity(recruit.resids)) +
         wrap_elements(~plotResiduals(recruit.resids, form = factor(rep(1, nrow(recruit))))) +
         wrap_elements(~testDispersion(recruit.resids)),
       width = 12,
       height = 4,
       dpi = 300)

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
  geom_point(data = recruit, aes(x = H_mean_local, y = Total), alpha = 0.4,
             position = position_jitter(width = 0.1, height = 0.1))  +
  scale_x_continuous(name = expression(paste(italic('Sargassum'), ' height (cm)'))) +
  scale_y_continuous('Number of coral recruits') +
  labs(tag = 'a') +
  theme_classic()  +
  theme(text = element_text(colour = 'black'), 
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.5)))
M3Figure

ggsave(file = paste0(FIGS_PATH, "/M3_H_local_fig.png"), 
       width = 160, 
       height = 160/1.6, 
       units = "mm", 
       dpi = 300)

# Compare ----
## ---- CompareM2vsM3
loo::loo_compare(loo::loo(recruit.brm4),
                 loo::loo(recruit.brm6))
## ----end


# M4: Recruit ~ Density ----
## Fit model ----
form <- bf(Total ~ D_broad + (1|Turf_height), family = poisson(link = 'log')) 
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
  SUYR_prior_and_posterior()   +
  theme_classic() + theme(text = element_text(colour = 'black'), 
                          axis.text = element_text(size = rel(1.2)),
                          axis.title = element_text(size = rel(1.5)),
                          legend.text = element_text(size = rel(1.2)),
                          legend.title = element_text(size = rel(1.5)),
                          legend.position = 'bottom')

ggsave(file = paste0(FIGS_PATH, "/M4PriorsPosteriors.png"), 
       width = 10,
       height = 8, 
       dpi = 300)

## MCMC Sampling diagnostics
(recruit.brm8$fit |> stan_trace()) + (recruit.brm8$fit |> stan_ac()) + (recruit.brm8$fit |> stan_rhat()) + (recruit.brm8$fit |> stan_ess())
ggsave(file = paste0(FIGS_PATH, "/M4MCMC.png"), 
       width = 10,
       height = 8, 
       dpi = 300)

## Model validation
### Posterior probability check
recruit.brm8 |> 
  pp_check(type = 'dens_overlay', ndraws = 100)

ggsave(file = paste0(FIGS_PATH, "/M4PPCheck.png"), 
       width = 10,
       height = 8, 
       dpi = 300)

### Residuals
recruit.resids <- make_brms_dharma_res(recruit.brm8, integerResponse = FALSE)
testUniformity(recruit.resids)
plotResiduals(recruit.resids, form = factor(rep(1, nrow(recruit))))
plotResiduals(recruit.resids, quantreg = TRUE) 
testDispersion(recruit.resids)

ggsave(filename = paste0(FIGS_PATH, '/M4DHARMa.png'),
       wrap_elements(~testUniformity(recruit.resids)) +
         wrap_elements(~plotResiduals(recruit.resids, form = factor(rep(1, nrow(recruit))))) +
         wrap_elements(~testDispersion(recruit.resids)),
       width = 12,
       height = 4,
       dpi = 300)

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
  ggplot(aes(x = D_broad*100,
             y = y)) + 
  geom_line() + 
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.2) +
  geom_point(data = recruit, aes(x = D_broad*100, y = Total), alpha = 0.4,
             position = position_jitter(width = 0.1, height = 0.1))   +
  scale_x_continuous(name = expression(paste(italic('Sargassum'), ' density (thalli/', m^{2}, ')'))) +
  scale_y_continuous('Number of coral recruits') +
  labs(tag = 'b') +
  theme_classic()  +
  theme(text = element_text(colour = 'black'), 
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.5)))
M4Figure

ggsave(file = paste0(FIGS_PATH, "/M4_Density_fig.png"), 
       width = 160, 
       height = 160/1.6, 
       units = "mm", 
       dpi = 300)

# Summary fig (H + D) ----
M3Figure / M4Figure

ggsave(file = paste0(FIGS_PATH, "/M_H_D_fig.png"), 
       width = 200/1.2, 
       height = 200, 
       units = "mm", 
       dpi = 300)

# Compare ALL ----
loo::loo_compare(loo::loo(recruit.brm2),
                 loo::loo(recruit.brm4),
                 loo::loo(recruit.brm6),
                 loo::loo(recruit.brm8))

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

recruit.brm6 |> 
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



# M5: Recruit ~ Diversity (local) ----
## Fit model ----
form <- bf(Total ~ scale(Shannon_local_cor) + scale(Shannon_local_alg) + (1|Turf_height), 
           family = poisson(link = 'log')) 
form |> get_prior(data = recruit)

priors <- prior(normal(0.5, 20), class = 'Intercept') +
  prior(normal(1, 5), class = 'b') +
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
  SUYR_prior_and_posterior()   +
  theme_classic() + theme(text = element_text(colour = 'black'), 
                          axis.text = element_text(size = rel(1.2)),
                          axis.title = element_text(size = rel(1.5)),
                          legend.text = element_text(size = rel(1.2)),
                          legend.title = element_text(size = rel(1.5)),
                          legend.position = 'bottom')

ggsave(file = paste0(FIGS_PATH, "/M5PriorsPosteriors.png"), 
       width = 10,
       height = 8, 
       dpi = 300)

## MCMC Sampling diagnostics
(recruit.brm10$fit |> stan_trace()) + (recruit.brm10$fit |> stan_ac()) + (recruit.brm10$fit |> stan_rhat()) + (recruit.brm10$fit |> stan_ess())
ggsave(file = paste0(FIGS_PATH, "/M5MCMC.png"), 
       width = 10,
       height = 8, 
       dpi = 300)

## Model validation
### Posterior probability check
recruit.brm10 |> 
  pp_check(type = 'dens_overlay', ndraws = 100)

ggsave(file = paste0(FIGS_PATH, "/M5PPCheck.png"), 
       width = 10,
       height = 8, 
       dpi = 300)

### Residuals
recruit.resids <- make_brms_dharma_res(recruit.brm10, integerResponse = FALSE)
testUniformity(recruit.resids)
plotResiduals(recruit.resids, form = factor(rep(1, nrow(recruit))))
plotResiduals(recruit.resids, quantreg = TRUE) 
testDispersion(recruit.resids)

ggsave(filename = paste0(FIGS_PATH, '/M5DHARMa.png'),
       wrap_elements(~testUniformity(recruit.resids)) +
         wrap_elements(~plotResiduals(recruit.resids, form = factor(rep(1, nrow(recruit))))) +
         wrap_elements(~testDispersion(recruit.resids)),
       width = 12,
       height = 4,
       dpi = 300)

## Save model ----
save(recruit.brm10, form, priors, recruit, file = 'data/modelled/M5_Local_diversity.RData')

## Investigation ----
recruit.brm10 |> 
  conditional_effects() |> 
  plot(points = TRUE)

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

# M6: Recruit ~ Diversity (broad) ----
## Fit model ----
form <- bf(Total ~ scale(Shannon_broad_cor) + scale(Shannon_broad_alg) + (1|Turf_height), 
           family = poisson(link = 'log')) 
form |> get_prior(data = recruit)

priors <- prior(normal(0.5, 15), class = 'Intercept') +
  prior(normal(1, 10), class = 'b') +
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
  SUYR_prior_and_posterior()   +
  theme_classic() + theme(text = element_text(colour = 'black'), 
                          axis.text = element_text(size = rel(1.2)),
                          axis.title = element_text(size = rel(1.5)),
                          legend.text = element_text(size = rel(1.2)),
                          legend.title = element_text(size = rel(1.5)),
                          legend.position = 'bottom')

ggsave(file = paste0(FIGS_PATH, "/M6PriorsPosteriors.png"), 
       width = 10,
       height = 8, 
       dpi = 300)

## MCMC Sampling diagnostics
(recruit.brm12$fit |> stan_trace()) + (recruit.brm12$fit |> stan_ac()) + (recruit.brm12$fit |> stan_rhat()) + (recruit.brm12$fit |> stan_ess())
ggsave(file = paste0(FIGS_PATH, "/M6MCMC.png"), 
       width = 10,
       height = 8, 
       dpi = 300)

## Model validation
### Posterior probability check
recruit.brm12 |> 
  pp_check(type = 'dens_overlay', ndraws = 100)

ggsave(file = paste0(FIGS_PATH, "/M6PPCheck.png"), 
       width = 10,
       height = 8, 
       dpi = 300)

### Residuals
recruit.resids <- make_brms_dharma_res(recruit.brm12, integerResponse = FALSE)
testUniformity(recruit.resids)
plotResiduals(recruit.resids, form = factor(rep(1, nrow(recruit))))
plotResiduals(recruit.resids, quantreg = TRUE) 
testDispersion(recruit.resids)

ggsave(filename = paste0(FIGS_PATH, '/M6DHARMa.png'),
       wrap_elements(~testUniformity(recruit.resids)) +
         wrap_elements(~plotResiduals(recruit.resids, form = factor(rep(1, nrow(recruit))))) +
         wrap_elements(~testDispersion(recruit.resids)),
       width = 12,
       height = 4,
       dpi = 300)

## Save model ----
save(recruit.brm12, form, priors, recruit, file = 'data/modelled/M6_Broad_diversity.RData')

## Investigation ----
recruit.brm12 |> 
  conditional_effects() |> 
  plot(points = TRUE)

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

# Compare ----
## ---- CompareM5vsM6
loo::loo_compare(loo::loo(recruit.brm10),
                 loo::loo(recruit.brm12))
## ----end


loo::loo_compare(loo::loo(recruit.brm4),
                 loo::loo(recruit.brm6),
                 loo::loo(recruit.brm8),
                 loo::loo(recruit.brm10),
                 loo::loo(recruit.brm12))
