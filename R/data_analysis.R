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

recruit <- recruit |>
  mutate(Treatment = factor(Treatment),
         Grazing = factor(Grazing, c('No', 'Light', 'Medium', 'Heavy'), ordered = TRUE))

# M1: Recruit ~ Treat ----
## Fit model ----
form <- bf(Total ~ Treatment, family = poisson(link = 'log')) 
form |> get_prior(data = recruit)

priors <- prior(normal(0.5, 1), class = 'Intercept') +
  prior(normal(0, 3), class = 'b')

recruit.brm <- brm(form, prior = priors, data = recruit, 
                   sample_prior = 'only', 
                   iter = 5000, 
                   warmup = 1000, 
                   chains = 3, cores = 3, 
                   thin = 5, 
                   refresh = 0, 
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
              'Mat_AnyAlg' = c(-1/2, 0, -1/2, 1))
recruit.brm2 |>
  emmeans(~Treatment, type = 'link') |>
  contrast(method = list(Treatment = cmat)) |>
  gather_emmeans_draws() |>
  mutate(Response = exp(.value)) |>
  summarise(median_hdci(Response),
            Pl = mean(Response < 1),
            Pg = mean(Response > 1))
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

# M2: ----

