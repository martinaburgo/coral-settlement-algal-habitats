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
data <- read.csv(file = 'data/processed/natural_recruitment.csv') |>
  dplyr::select(-X)

data |> glimpse()

# Data prep ----
data <- data |>
  mutate(Habitat = factor(Habitat, levels = c('Flat', 'Crest', 'Slope'), 
                        ordered = TRUE),
         class_sizes = factor(class_sizes, levels = c('<5', '6-20', '21-40', '>40'), 
                              ordered = TRUE),
         Depth = factor(factor(Depth, levels = c('0-1 m', '2-3 m', '4-5 m'), 
                               ordered = TRUE)))

# Data analysis ----
## EDA ----
data |>
  ggplot(aes(y = Diameter, x = Mean)) +
  geom_point(aes(colour = class_sizes), size = 2, position = position_jitterdodge()) +
  scale_color_viridis_d(option = 'D', name = 'Coral class sizes') +
  scale_x_continuous(expression('Average '*italic(Sargassum)*' height (cm)')) +
  ylab('Coral class size (cm)') +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal')

ggsave(file = paste0(FIGS_PATH, "/EDA_natural_juvs.png"), 
       width = 160, 
       height = 160/1.6, 
       units = "mm", 
       dpi = 300)

data |>
  ggplot(aes(y = Diameter, x = Density)) +
  geom_point(aes(colour = class_sizes), size = 2, position = position_jitterdodge()) +
  scale_color_viridis_d(option = 'D', name = 'Coral class sizes') +
  scale_x_continuous(expression('Average '*italic(Sargassum)*' density (thalli/m'^'2'*')')) +
  ylab('Coral class size (cm)') +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal')

data |>
  ggplot(aes(y = Diameter, x = Diversity)) +
  geom_point(aes(colour = class_sizes), size = 2, position = position_jitterdodge()) +
  scale_color_viridis_d(option = 'D', name = 'Coral class sizes') +
  scale_x_continuous('Average macroalgal diversity') +
  ylab('Coral class size (cm)') +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal')

data |>
  ggplot(aes(y = Diameter, x = Sargassum)) +
  geom_point(aes(colour = class_sizes), size = 2, position = position_jitterdodge()) +
  scale_color_viridis_d(option = 'D', name = 'Coral class sizes') +
  scale_x_continuous(expression('Average '*italic(Sargassum)*' cover (%)')) +
  ylab('Coral class size (cm)') +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal')

data |>
  ggplot(aes(y = Diameter, x = Lobophora)) +
  geom_point(aes(colour = class_sizes), size = 2, position = position_jitterdodge()) +
  scale_color_viridis_d(option = 'D', name = 'Coral class sizes') +
  scale_x_continuous(expression('Average '*italic(Lobophora)*' cover (%)')) +
  ylab('Coral class size (cm)') +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal')


# Data analysis ----
cor(data[, -c(1:6, 9)]) |> corrplot(method = 'number')

## MN2 - Height + Density ----
### Fit ----
## ----MN2Fit
size.form <- bf(class_sizes ~ scale(Mean) + scale(Density) + (1|Taxa),
                family = cumulative(link = 'logit', 
                                      threshold = 'flexible'))
## ----end

priors <- prior(normal(0, 1), class = 'Intercept') +
  prior(normal(0, 1), class = 'b') +
  prior(student_t(3, 0, 1), class = 'sd')


size.brm <- brm(size.form,
                  data = data,
                  prior = priors,
                  sample_prior = 'only',
                  iter = 5000,
                  warmup = 1000,
                  chains = 3, cores = 3,
                  thin = 5,
                  backend = 'rstan')

#check priors
size.brm |>
  conditional_effects(categorical = TRUE) |>
  plot(poits = TRUE, ask = FALSE)

#add data 
size.brm2 <- update(size.brm,
                    sample_prior = 'yes',
                    control = list(adapt_delta = 0.99),
                    refresh = 100)

#save model
#save(size.brm2, size.form, priors, data, file = 'data/modelled/MN2_HandD.RData')

size.brm2 |>
  conditional_effects('Mean',
                      categorical = TRUE)
size.brm2 |>
  conditional_effects('Density',
                      categorical = TRUE)

### Diagnostics ----

# MCMC Sampling diagnostics
pars <- size.brm2 |> get_variables() |>
  str_subset("^b_.*|^sd_.*")

(size.brm2$fit |> stan_trace(pars = pars)) + (size.brm2$fit |> stan_ac(pars= pars)) + 
  (size.brm2$fit |> stan_rhat()) + (size.brm2$fit |> stan_ess())

ggsave(file = paste0(FIGS_PATH, "/MN2MCMC.png"), 
       width = 10,
       height = 8, 
       dpi = 300)

## Model validation
### Posterior probability check
size.brm2 |> 
  pp_check(type = 'dens_overlay', ndraws = 100)

ggsave(file = paste0(FIGS_PATH, "/MN2PPCheck.png"), 
       width = 10,
       height = 8, 
       dpi = 300)

# Residuals
recruit.resids <- make_brms_dharma_res(size.brm2, 
                                       integerResponse = TRUE)
testUniformity(recruit.resids)
plotResiduals(recruit.resids, quantreg = FALSE) 
testDispersion(recruit.resids)

ggsave(filename = paste0(FIGS_PATH, '/MN2DHARMa.png'),
       wrap_elements(~testUniformity(recruit.resids)) +
         wrap_elements(~plotResiduals(recruit.resids)) +
         wrap_elements(~testDispersion(recruit.resids)),
       width = 12,
       height = 4,
       dpi = 300)

### Output ----
size.brm2 |>
  as_draws_df() |>
  dplyr::select(matches("^b_.*|^sd_.*")) |> #they're on logit scale
  mutate(across(matches("^b_.*"), exp)) |>
  summarise_draws(median,
                  HDInterval::hdi,
                  Pl = ~ mean(.x < 1),
                  Pg = ~ mean(.x > 1),
                  rhat,
                  ess_bulk,
                  ess_tail) |>
  mutate(across(where(is.numeric), ~ round(., 2))) |>
  #write.table(file = 'output/tables/MN2Output.txt', sep = ",", quote = FALSE, row.names = F)
  
size.brm2 |>
  emmeans(~Mean, type = 'response', at = list(Mean = c(0, 10, 20, 30, 40, 50, 60, 70, 80))) |>
  pairs(rev = TRUE) |>
  tidy_draws() |>
 # exp() |> 
  summarise_draws(median,
                  HDInterval::hdi,
                  Pl = ~ mean(.x < 1),
                  Pg = ~ mean(.x > 1))

## Figure ----
preds <- add_epred_draws(size.brm2, 
                         newdata = expand.grid(Mean = seq(min(data$Mean), max(data$Mean), len = 100),
                                               Density = mean(data$Density)),
                re_formula = NA) |>
  group_by(Mean, .category) |>
  summarise(median_hdci(.epred))

coeff <- 150

ggplot() +
  geom_ribbon(data = preds, 
              aes(x = Mean, 
                  ymin = ymin, 
                  ymax = ymax, 
                  fill = .category), 
              alpha = 0.2) +
  geom_line(data = preds, 
            aes(y = y,
                x = Mean, 
                colour = .category), 
            size = 0.8) +
  geom_point(data = data, 
             mapping = aes(y = Diameter/coeff, 
                           x = Mean),
             size = 1.5,
             alpha = 0.5,
             position = position_jitter(width = 0.1, height = -0.1)) +
  labs(col = 'Coral class size (cm)', 
       fill = 'Coral class size (cm)', 
       x = expression(paste(italic('Sargassum'), ' height (cm)')))  +
  scale_y_continuous(name = "Probability",
                     sec.axis = sec_axis(~.*coeff, 
                                         name="Coral diameter (cm)")) +
  scale_fill_viridis_d(option = 'D') +
  scale_colour_viridis_d(option = 'D') +
  theme_classic() +
  theme(text = element_text(colour = 'black'), 
        axis.text = element_text(size = rel(1.1)),
        axis.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.1)),
        legend.title = element_text(size = rel(1.2)),
        legend.position = 'bottom')

ggsave(file = paste0(FIGS_PATH, "/MN2_height.png"), 
       width = 160, 
       height = 160/1.6, 
       units = "mm", 
       dpi = 300)

## V2
MN2Figure <- ggplot() +
  geom_ribbon(data = preds, 
              aes(x = Mean, 
                  ymin = ymin, 
                  ymax = ymax, 
                  fill = .category), 
              alpha = 0.2) +
  geom_line(data = preds, 
            aes(y = y,
                x = Mean, 
                colour = .category), 
            size = 0.8) +
  labs(col = 'Coral class size (cm)', 
       fill = 'Coral class size (cm)', 
       x = expression(paste('Canopy height (cm)')),
       tag = 'a')  +
  scale_y_continuous(name = "Probability") +
  scale_fill_viridis_d(option = 'D') +
  scale_colour_viridis_d(option = 'D') +
  theme_classic() +
  theme(text = element_text(colour = 'black'), 
        axis.text = element_text(size = rel(1.1)),
        axis.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.1)),
        legend.title = element_text(size = rel(1.2)),
        legend.position = 'bottom')

ggsave(file = paste0(FIGS_PATH, "/MN2_HandD.png"), 
       width = 160, 
       height = 160/1.4, 
       units = "mm", 
       dpi = 300)


# Thresholds ----
add_epred_draws(size.brm2, 
                newdata = expand.grid(Mean = seq(min(data$Mean),
                                                 max(data$Mean),
                                                 len = 100),
                                      Density = mean(data$Density)), 
                re_formula = NA) |>
  group_by(Mean, .category) |>
  filter(.epred < 0.1 & .category == '<5') |>
  ungroup() |>
  summarise(median_hdci(Mean),
            Pl = mean(Mean < 1),
            Pg = mean(Mean > 1))

add_epred_draws(size.brm2, 
                newdata = expand.grid(Mean = seq(min(data$Mean),
                                                 max(data$Mean),
                                                 len = 100),
                                      Density = mean(data$Density)), 
                re_formula = NA) |>
  group_by(Mean, .category) |>
  filter(.epred > 0.85 & .category == '>40' | .category == '21-40') |>
  ungroup() |>
  summarise(median_hdci(Mean),
            Pl = mean(Mean < 1),
            Pg = mean(Mean > 1))

## MN3 - Sargassum + Lobophora ----
### Fit ----
## ----MN3Fit
size.form <- bf(class_sizes ~ scale(Sargassum) + scale(Lobophora) + (1|Taxa),
                family = cumulative(link = 'logit', 
                                    threshold = 'flexible'))
## ----end

priors <- prior(normal(0, 1), class = 'Intercept') +
  prior(normal(0, 1), class = 'b') +
  prior(student_t(3, 0, 1), class = 'sd')


size.brm3 <- brm(size.form,
                data = data,
                prior = priors,
                sample_prior = 'only',
                iter = 5000,
                warmup = 1000,
                chains = 3, cores = 3,
                thin = 5,
                backend = 'rstan')

#check priors
size.brm3 |>
  conditional_effects(categorical = TRUE) |>
  plot(poits = TRUE, ask = FALSE)

#add data 
size.brm3 <- update(size.brm3,
                    sample_prior = 'yes',
                    control = list(adapt_delta = 0.99),
                    refresh = 100)

#save model
#save(size.brm3, size.form, priors, data, file = 'data/modelled/MN3_density.RData')

size.brm3 |>
  conditional_effects('Sargassum',
                      categorical = TRUE)
size.brm3 |>
  conditional_effects('Lobophora',
                      categorical = TRUE)

### Diagnostics ----

# MCMC Sampling diagnostics
pars <- size.brm3 |> get_variables() |>
  str_subset("^b_.*|^sd_.*")

(size.brm3$fit |> stan_trace(pars = pars)) + (size.brm3$fit |> stan_ac(pars= pars)) + 
  (size.brm3$fit |> stan_rhat()) + (size.brm3$fit |> stan_ess())

ggsave(file = paste0(FIGS_PATH, "/MN3MCMC.png"), 
       width = 10,
       height = 8, 
       dpi = 300)

## Model validation
### Posterior probability check
size.brm3 |> 
  pp_check(type = 'dens_overlay', ndraws = 100)

ggsave(file = paste0(FIGS_PATH, "/MN3PPCheck.png"), 
       width = 10,
       height = 8, 
       dpi = 300)

# Residuals
recruit.resids <- make_brms_dharma_res(size.brm3, 
                                       integerResponse = TRUE)
testUniformity(recruit.resids)
plotResiduals(recruit.resids, quantreg = FALSE) 
testDispersion(recruit.resids)

ggsave(filename = paste0(FIGS_PATH, '/MN3DHARMa.png'),
       wrap_elements(~testUniformity(recruit.resids)) +
         wrap_elements(~plotResiduals(recruit.resids)) +
         wrap_elements(~testDispersion(recruit.resids)),
       width = 12,
       height = 4,
       dpi = 300)

### Output ----
size.brm3 |>
  as_draws_df() |>
  dplyr::select(matches("^b_.*|^sd_.*")) |> #they're on logit scale
  mutate(across(matches("^b_.*"), exp)) |>
  summarise_draws(median,
                  HDInterval::hdi,
                  Pl = ~ mean(.x < 1),
                  Pg = ~ mean(.x > 1),
                  rhat,
                  ess_bulk,
                  ess_tail) |>
  mutate(across(where(is.numeric), ~ round(., 2))) |>
  #write.table(file = 'output/tables/MN3Output.txt', sep = ",", quote = FALSE, row.names = F)
  
  size.brm3 |>
  emmeans(~Sargassum, type = 'response', at = list(Sargassum = c(8, 10, 12, 14, 16, 18, 20))) |>
  pairs(rev = TRUE) |>
  tidy_draws() |>
  # exp() |> 
  summarise_draws(median,
                  HDInterval::hdi,
                  Pl = ~ mean(.x < 1),
                  Pg = ~ mean(.x > 1))

## Figure ----
preds <- add_epred_draws(size.brm3, 
                         newdata = expand.grid(Sargassum = seq(min(data$Sargassum),
                                                          max(data$Sargassum),
                                                          len = 100),
                                               Lobophora = mean(data$Lobophora)), 
                         re_formula = NA) |>
  group_by(Sargassum, .category) |>
  summarise(median_hdci(.epred))

MN3Figure <- ggplot() +
  geom_ribbon(data = preds, 
              aes(x = Sargassum, 
                  ymin = ymin, 
                  ymax = ymax, 
                  fill = .category), 
              alpha = 0.2) +
  geom_line(data = preds, 
            aes(y = y,
                x = Sargassum, 
                colour = .category), 
            size = 0.8) +
  labs(col = 'Coral class size (cm)', 
       fill = 'Coral class size (cm)', 
       x = expression(italic(Sargassum)*' cover (%)'),
       tag = 'b')  +
  scale_y_continuous(name = "Probability") +
  scale_fill_viridis_d(option = 'D') +
  scale_colour_viridis_d(option = 'D') +
  theme_classic() +
  theme(text = element_text(colour = 'black'), 
        axis.text = element_text(size = rel(1.1)),
        axis.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.1)),
        legend.title = element_text(size = rel(1.2)),
        legend.position = 'bottom')

ggsave(file = paste0(FIGS_PATH, "/MN3_sargassum.png"), 
       width = 160, 
       height = 160/1.4, 
       units = "mm", 
       dpi = 300)

# Combined figure ----
MN2Figure / MN3Figure +
  plot_layout(axis_titles = 'collect', guides = 'collect') &
  theme(legend.position = 'bottom')

ggsave(file = paste0(FIGS_PATH, "/MN_combined.png"), 
       width = 220/1.4, 
       height = 220, 
       units = "mm", 
       dpi = 300)

## MN4 - Diversity ----
### Fit ----
## ----MN4Fit
size.form <- bf(class_sizes ~ Diversity + (1|Taxa),
                family = cumulative(link = 'logit', 
                                    threshold = 'flexible'))
## ----end

priors <- prior(normal(0, 1), class = 'Intercept') +
  prior(normal(0, 1), class = 'b') +
  prior(student_t(3, 0, 1), class = 'sd')


size.brm4 <- brm(size.form,
                 data = data,
                 prior = priors,
                 sample_prior = 'only',
                 iter = 5000,
                 warmup = 1000,
                 chains = 3, cores = 3,
                 thin = 5,
                 backend = 'rstan')

#check priors
size.brm4 |>
  conditional_effects(categorical = TRUE) |>
  plot(poits = TRUE, ask = FALSE)

#add data 
size.brm4 <- update(size.brm4,
                    sample_prior = 'yes',
                    control = list(adapt_delta = 0.99),
                    refresh = 100)

#save model
#save(size.brm4, size.form, priors, data, file = 'data/modelled/MN4_diversity.RData')

size.brm4 |>
  conditional_effects('Diversity',
                      categorical = TRUE)

### Diagnostics ----

# MCMC Sampling diagnostics
pars <- size.brm4 |> get_variables() |>
  str_subset("^b_.*|^sd_.*")

(size.brm4$fit |> stan_trace(pars = pars)) + (size.brm4$fit |> stan_ac(pars= pars)) + 
  (size.brm4$fit |> stan_rhat()) + (size.brm4$fit |> stan_ess())

ggsave(file = paste0(FIGS_PATH, "/MN4MCMC.png"), 
       width = 10,
       height = 8, 
       dpi = 300)

## Model validation
### Posterior probability check
size.brm4 |> 
  pp_check(type = 'dens_overlay', ndraws = 100)

ggsave(file = paste0(FIGS_PATH, "/MN4PPCheck.png"), 
       width = 10,
       height = 8, 
       dpi = 300)

# Residuals
recruit.resids <- make_brms_dharma_res(size.brm4, 
                                       integerResponse = TRUE)
testUniformity(recruit.resids)
plotResiduals(recruit.resids, quantreg = FALSE) 
testDispersion(recruit.resids)

ggsave(filename = paste0(FIGS_PATH, '/MN4DHARMa.png'),
       wrap_elements(~testUniformity(recruit.resids)) +
         wrap_elements(~plotResiduals(recruit.resids)) +
         wrap_elements(~testDispersion(recruit.resids)),
       width = 12,
       height = 4,
       dpi = 300)

### Output ----
size.brm4 |>
  as_draws_df() |>
  dplyr::select(matches("^b_.*|^sd_.*")) |> #they're on logit scale
  mutate(across(matches("^b_.*"), exp)) |>
  summarise_draws(median,
                  HDInterval::hdi,
                  Pl = ~ mean(.x < 1),
                  Pg = ~ mean(.x > 1),
                  rhat,
                  ess_bulk,
                  ess_tail) |>
  mutate(across(where(is.numeric), ~ round(., 2))) |>
  #write.table(file = 'output/tables/MN3Output.txt', sep = ",", quote = FALSE, row.names = F)
  
  size.brm4 |>
  emmeans(~Diversity, type = 'response', at = list(Diversity = c(0, 0.2, 0.4, 0.6, 0.8, 1))) |>
  pairs(rev = TRUE) |>
  tidy_draws() |>
  # exp() |> 
  summarise_draws(median,
                  HDInterval::hdi,
                  Pl = ~ mean(.x < 1),
                  Pg = ~ mean(.x > 1))

## Figure ----
preds <- add_epred_draws(size.brm4, 
                         newdata = expand.grid(Diversity = seq(min(data$Diversity),
                                                             max(data$Diversity),
                                                             len = 100)), 
                         re_formula = NA) |>
  group_by(Diversity, .category) |>
  summarise(median_hdci(.epred))

ggplot() +
  geom_ribbon(data = preds, 
              aes(x = Diversity, 
                  ymin = ymin, 
                  ymax = ymax, 
                  fill = .category), 
              alpha = 0.2) +
  geom_line(data = preds, 
            aes(y = y,
                x = Diversity, 
                colour = .category), 
            size = 0.8) +
  labs(col = 'Coral class size (cm)', 
       fill = 'Coral class size (cm)', 
       x = 'Macroalgal diversity')  +
  scale_y_continuous(name = "Probability") +
  scale_fill_viridis_d(option = 'D') +
  scale_colour_viridis_d(option = 'D') +
  theme_classic() +
  theme(text = element_text(colour = 'black'), 
        axis.text = element_text(size = rel(1.1)),
        axis.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.1)),
        legend.title = element_text(size = rel(1.2)),
        legend.position = 'bottom')

ggsave(file = paste0(FIGS_PATH, "/MN4_diversity.png"), 
       width = 160, 
       height = 160/1.4, 
       units = "mm", 
       dpi = 300)



