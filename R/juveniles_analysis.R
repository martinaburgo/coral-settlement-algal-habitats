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
  ggplot(aes(y = class_sizes, x = Mean)) +
  geom_point(aes(colour = Depth), size = , position = position_jitterdodge()) +
  scale_color_viridis_d(option = 'G') +
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


# Data analysis ----

## MN2 - Class sizes ----
### Fit ----
## ----MN2Fit
size.form <- bf(class_sizes ~ Mean + (1|Taxa),
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
save(size.brm2, size.form, priors, dat2, file = 'data/modelled/MN2_height.RData')

size.brm2 |>
  conditional_effects('Mean',
                      categorical = TRUE)


## MCMC Sampling diagnostics
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

### Residuals
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
  mutate(across(where(is.numeric), ~ round(., 2)))

#Figure
newdata <- with(dat2, 
                expand.grid(Mean = seq(min(dat2$Mean),
                                       max(dat2$Mean),
                                       len = 100)))

add_epred_draws(size.brm2, 
                newdata = newdata, 
                re_formula = NA) |>
  group_by(Mean, .category) |>
  summarise(median_hdci(.epred)) |> 
  ggplot(aes(y = y, 
             x = Mean, 
             colour = .category, 
             fill = .category)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.2) +
  geom_line(size = 0.8) +
  labs(col = 'Coral class size (cm)', 
       fill = 'Coral class size (cm)', 
       y = 'Probability', 
       x = expression(paste(italic('Sargassum'), ' height (cm)')))  +
  theme_classic() +
  theme(legend.position = 'bottom')

ggsave(file = paste0(FIGS_PATH, "/MN2_height.png"), 
       width = 250, 
       height = 250/1.6, 
       units = "mm", 
       dpi = 300)

# Old analyses ----

## MN1 - Height ----

### EDA ----
dat2 |>
  ggplot(aes(y = class_sizes, x = Mean)) +
  geom_point(position = position_dodge2(width = 0.3)) +
  theme_classic() +
  scale_x_continuous(expression('Average '*italic(Sargassum)*' height (cm)')) +
  ylab('Coral class size (cm)')


# Fit model
priors <- prior(normal(3, 2), class = 'Intercept') +
  prior(normal(0, 1), class = 'b') +
  prior(student_t(3, 0, 1), class = 'sd')  +
  prior(student_t(3, 0, 1), class = 'sigma')

form <- bf(Diameter ~ log(Height) + (1|Taxa), family = gaussian(link = 'log')) 
form |>
  get_prior(data = data)

nat.brm <- brm(form, prior = priors, data = data, 
               sample_prior = 'only', 
               iter = 5000, 
               warmup = 1000, 
               chains = 3, cores = 3, 
               thin = 5, 
               refresh = 0, 
               control = list(adapt_delta = 0.99, 
                              max_treedepth = 20),
               backend = 'rstan') 
nat.brm |> 
  conditional_effects() |> 
  plot(points = TRUE)

nat.brm2 <- nat.brm |>
  update(sample_prior = 'yes')

### Diagnostics ----
nat.brm2 |> 
  SUYR_prior_and_posterior() +
  theme_classic() + theme(text = element_text(colour = 'black'), 
                          axis.text = element_text(size = rel(1.2)),
                          axis.title = element_text(size = rel(1.5)),
                          legend.text = element_text(size = rel(1.2)),
                          legend.title = element_text(size = rel(1.5)),
                          legend.position = 'bottom')

ggsave(file = paste0(FIGS_PATH, "/MN1PriorsPosteriors.png"), 
       width = 10,
       height = 8, 
       dpi = 300)

## MCMC Sampling diagnostics
(nat.brm2$fit |> stan_trace()) + (nat.brm2$fit |> stan_ac()) + (nat.brm2$fit |> stan_rhat()) + (nat.brm2$fit |> stan_ess())

ggsave(file = paste0(FIGS_PATH, "/MN1MCMC.png"), 
       width = 10,
       height = 8, 
       dpi = 300)

## Model validation
### Posterior probability check
nat.brm2 |> 
  pp_check(type = 'dens_overlay', ndraws = 100)

ggsave(file = paste0(FIGS_PATH, "/MN1PPCheck.png"), 
       width = 10,
       height = 8, 
       dpi = 300)

### Residuals
nat.resids <- make_brms_dharma_res(nat.brm2, integerResponse = FALSE)
testUniformity(nat.resids)
plotResiduals(nat.resids, form = factor(rep(1, nrow(data))))
plotResiduals(nat.resids, quantreg = TRUE) 
testDispersion(nat.resids)

ggsave(filename = paste0(FIGS_PATH, '/MN1DHARMa.png'),
       wrap_elements(~testUniformity(nat.resids)) +
         wrap_elements(~plotResiduals(nat.resids, form = factor(rep(1, nrow(data))))) +
         wrap_elements(~testDispersion(nat.resids)),
       width = 12,
       height = 4,
       dpi = 300)

### Save model ----
save(nat.brm2, form, priors, data, file = 'data/modelled/MN1_Height.RData')

### Investigation ----
nat.brm2 |> 
  conditional_effects() |> 
  plot(points = TRUE)

#### Output ----
### ---- MN1Output
nat.brm2 |>
  summarise_draws(median,
                  ~ HDInterval::hdi(.x),
                  rhat, length, ess_bulk, ess_tail,
                  Pl = ~mean(.x < 1),
                  Pg = ~mean(.x > 1))
### ----end

nat.brm2 |>
  summarise_draws(median,
                  ~ HDInterval::hdi(.x),
                  rhat, length, ess_bulk, ess_tail,
                  Pl = ~mean(.x < 1),
                  Pg = ~mean(.x > 1)) |>
  mutate(median = round(median, 3),
         lower = round(lower, 3),
         upper = round(upper, 3),
         rhat = round(rhat, 3),
         Pl = round(Pl, 3),
         Pg = round(Pg, 3)) |>
  write.table(file = 'output/tables/MN1Output.txt', sep = ",", quote = FALSE, row.names = F)

### Depth ----

d.clmm2 <- ordinal::clmm(class_sizes ~ Depth + (1|Taxa),
                         data = dat2) #Fit model
summary(d.clmm2)
emmeans(d.clmm2, ~Depth, mode = 'mean.class') |>
  pairs()

emmeans(d.clmm2, ~class_sizes|Depth, mode = 'prob') |>
  as.data.frame() |>
  ggplot(aes(y = prob, x = Depth,
             colour = class_sizes)) +
  geom_pointrange(aes(ymin = abs(asymp.LCL), ymax = asymp.UCL)) +
  geom_point() +
  geom_line(aes(group = class_sizes)) +
  theme_classic() +
  theme(legend.position = 'bottom',
        text = element_text(size = 20)) + labs(colour = 'Coral class size') +
  ylab('Probability') +
  scale_color_viridis_d(labels = c('1' = '> 5 cm', 
                                   '2' = '6-20 cm',
                                   '3' = '21-40 cm',
                                   '4' = '> 40 cm'), 
                        option = 'viridis')

ggsave(file = paste0(FIGS_PATH, "/CLMM_depth.png"), 
       width = 200, 
       height = 200/1.6, 
       units = "mm", 
       dpi = 300)

#sarg height

d.clmm <- ordinal::clmm(class_sizes ~ Mean + (1|Taxa),
                        data = dat2) #Fit model
summary(d.clmm)


sjPlot::plot_model(d.clmm, type = 'pred')$data |>
  as.data.frame() |>
  ggplot(aes(x = x, y = predicted, colour = response.level, fill = response.level)) +
  geom_line() +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha = 0.2) +
  theme_classic() +
  theme(legend.position = 'bottom',
        text = element_text(size = 20)) + labs(colour = 'Coral class size (cm)') +
  ylab('Probability') +
  scale_x_continuous(expression('Average '*italic(Sargassum)*' height (cm)')) +
  scale_color_viridis_d(labels = c('1' = '< 5', 
                                   '2' = '6-20',
                                   '3' = '21-40',
                                   '4' = '> 40'), 
                        option = 'viridis')  +
  scale_fill_viridis_d(alpha = 0) +
  guides(fill = FALSE)

ggsave(file = paste0(FIGS_PATH, "/CLMM_height.png"), 
       width = 200, 
       height = 200/1.6, 
       units = "mm", 
       dpi = 300)

