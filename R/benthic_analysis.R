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
corals <- readxl::read_xlsx(path = "data/primary/Maggie-Quads-8.xlsx", sheet = "Florence") |>
  dplyr::select(Habitat, Trip, Period, Transect, Quadrat, Taxa, Interaction, Diameter, Distance, Cover, Height)
canopy <- rbind(readxl::read_xlsx(path = "data/primary/Maggie-Quads-1.xlsx", sheet = "Florence"),
                readxl::read_xlsx(path = "data/primary/Maggie-Quads-2.xlsx", sheet = "Florence"),
                readxl::read_xlsx(path = "data/primary/Maggie-Quads-3.xlsx", sheet = "Florence"),
                readxl::read_xlsx(path = "data/primary/Maggie-Quads-4.xlsx", sheet = "Florence"),
                readxl::read_xlsx(path = "data/primary/Maggie-Quads-5.xlsx", sheet = "Florence"),
                readxl::read_xlsx(path = "data/primary/Maggie-Quads-6.xlsx", sheet = "Florence"),
                readxl::read_xlsx(path = "data/primary/Maggie-Quads-7.xlsx", sheet = "Florence"),
                readxl::read_xlsx(path = "data/primary/Maggie-Quads-8.xlsx", sheet = "Florence")) |>
  dplyr::select(!c(Date, Metre, Taxa, Diameter, Interaction, Distance, Cover)) |>
  dplyr::filter(!is.na(Thalli))

benthos <- rbind(readxl::read_xlsx(path = "data/primary/Maggie-Quads-1.xlsx", sheet = "Florence"),
                 readxl::read_xlsx(path = "data/primary/Maggie-Quads-2.xlsx", sheet = "Florence"),
                 readxl::read_xlsx(path = "data/primary/Maggie-Quads-3.xlsx", sheet = "Florence"),
                 readxl::read_xlsx(path = "data/primary/Maggie-Quads-4.xlsx", sheet = "Florence"),
                 readxl::read_xlsx(path = "data/primary/Maggie-Quads-5.xlsx", sheet = "Florence"),
                 readxl::read_xlsx(path = "data/primary/Maggie-Quads-6.xlsx", sheet = "Florence"),
                 readxl::read_xlsx(path = "data/primary/Maggie-Quads-7.xlsx", sheet = "Florence"),
                 readxl::read_xlsx(path = "data/primary/Maggie-Quads-8.xlsx", sheet = "Florence")) |>
  dplyr::select(!c(Date))

data <- read.csv(file = 'data/processed/natural_recruitment.csv') |>
  dplyr::select(-X)

corals <- corals |>
  dplyr::filter(!is.na(Diameter))

breaks <- c(0, 5, 20,  40, max(corals$Diameter, na.rm = TRUE)) # set up cut-off values 
tags <- c("<5", "6-20", '21-40', '>40') # specify interval/bin labels
group_tags <- cut(corals$Diameter, 
                  breaks = breaks, 
                  include.lowest = TRUE, 
                  right = FALSE, 
                  labels = tags) # bucketing values into bins
summary(group_tags) # inspect bins
class_sizes <- factor(group_tags, 
                      levels = tags,
                      ordered = TRUE)


# Explore benthic data ----
## Benthic data ----
# Adjust cover:
benthos <- benthos |>
  mutate(Cover = as.numeric(ifelse(Cover == "+", '1', 
                                   ifelse(Cover == '=', '1', Cover)))) |>
  dplyr::filter(!is.na(Cover)) |>
  dplyr::select(!c(Distance, Interaction, Diameter, Height, Thalli)) |>
  group_by(Habitat, Trip, Transect, Quadrat) |>
  mutate(Cover = Cover/sum(Cover)*100) |>
  ungroup()

#Add categories
benthos <- benthos |>
  left_join(read.csv(file = 'data/primary/categories.csv'))


## Canopy height ----

for (i in 1:nrow(canopy)) {
  canopy[i, 'Mean'] <- canopy[i, 'Height'] |>
    str_split(', ') |>
    unlist() |>
    as.numeric() |>
    mean()
  
  canopy[i, 'Max'] <- canopy[i, 'Height'] |>
    str_split(', ') |>
    unlist() |>
    as.numeric() |>
    max()
}

canopy <- canopy |>
  mutate(Mean = as.numeric(Mean),
         Max = as.numeric(Max))


## General benthic composition (Methods) ----
# algae v coral cover trend
benthos |>
  group_by(Habitat, Trip, Transect, Quadrat, Category) |>
  summarise(sum = sum(Cover)) |>
  group_by(Habitat, Trip, Category) |>
  summarise(mean = mean(sum)) |>
  group_by(Habitat, Category) |>
  summarise(mean = mean(mean)) |>
  View()


benthos |>
  dplyr::filter(Taxa == 'Sargassum') |>
  group_by(Habitat, Trip, Transect, Quadrat) |>
  summarise(sum = sum(Cover)) |>
  group_by(Habitat, Trip) |>
  summarise(mean = mean(sum)) |>
  group_by(Habitat) |>
  summarise(mean = mean(mean)) |>
  View()

benthos |>
  dplyr::filter(Category == 'Alga') |>
  dplyr::filter(Habitat == 'Crest') |>
  group_by(Trip, Transect, Taxa) |>
  summarise(sum = mean(Cover)) |>
  group_by(Trip, Taxa) |>
  summarise(mean = mean(sum)) |>
  ungroup() |>
  dplyr::filter(Taxa == 'Halymenia' | Taxa == 'Galaxaura' | Taxa == 'Neomeris') |>
  ggplot(aes(y = mean, x = factor(Trip), fill = Taxa)) +
  geom_histogram(stat = 'identity')


### Canopy height ----
canopy |>
  mutate(Habitat = factor(Habitat, levels = c('Flat', 'Crest', 'Slope'), ordered = TRUE)) |>
  ggplot(aes(x = Trip, y = Mean)) +
  geom_violin(aes(group = Trip)) +
  facet_wrap(~Habitat)

#### Supp. Fig. ----

canopy |>
  mutate(Habitat = factor(Habitat, levels = c('Flat', 'Crest', 'Slope'), ordered = TRUE)) |>
  ggplot(aes(y = Mean, x = factor(Trip), color = Habitat)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.02, dodge.width = 0.9),
             alpha = 0.2) +
  geom_boxplot(fill = NA) +
  theme_classic() +
  theme(legend.position = c(0.9, 0.85),
        text = element_text(size = 20)) +
  xlab('Period') +
  scale_y_continuous(expression(italic(Sargassum)*' height (cm)'))  +
  scale_x_discrete(breaks = c("1","2","3", "4", "5", "6", "7", "8"),
                   labels = c("Jul/Aug","Sep/Oct","Nov/Dec", "Jan", "Feb/Mar", "Apr/May", "Jun", "Jul/Aug")) +
  scale_color_viridis_d(labels = c('Flat' = '0-1 m',
                                  'Crest' = '2-3 m',
                                  'Slope' = '4-5 m'),
                        name = 'Depth')

ggsave(file = paste0(FIGS_PATH, "/SF1_canopy_height.png"), 
       width = 250, 
       height = 250/1.6, 
       units = "mm", 
       dpi = 300)


canopy |> 
  mutate(Habitat = factor(Habitat, levels = c('Flat', 'Crest', 'Slope'), ordered = TRUE)) |>
  dplyr::filter(Trip == '3') |>
  ggplot(aes(x = Habitat, y = Mean)) +
  geom_violin() +
  geom_point() +
  theme_classic()


### Coral sizes ----
corals |>
  mutate(Habitat = factor(Habitat, levels = c('Flat', 'Crest', 'Slope'), ordered = TRUE)) |>
  dplyr::filter(!is.na(Diameter)) |>
  ggplot(aes(x = Habitat, y = Diameter)) +
  geom_violin() +
  geom_point(position = position_jitter(width = 0.2, height = 0.1)) +
  theme_classic()


corals |> 
  mutate(Habitat = factor(Habitat, levels = c('Flat', 'Crest', 'Slope'), ordered = TRUE)) |>
  mutate(class_sizes = factor(group_tags, 
                                  levels = tags,
                                  ordered = TRUE)) |>
  count(Habitat, class_sizes) |>
  ggplot() + geom_bar(mapping = aes(x = class_sizes, y = n,
                                    fill = Habitat), 
                      stat = 'identity', color = 'black') +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     text = element_text(size = 20),
                     legend.direction = 'vertical',
                     legend.position = c(0.85, 0.9))  +
  scale_fill_manual(values = c("white", "grey65", 'grey20')) +
  ylab('Number of corals') + xlab('Coral class size (cm)')


corals.sum <- corals |> 
  mutate(Habitat = factor(Habitat, levels = c('Flat', 'Crest', 'Slope'), ordered = TRUE)) |>
  mutate(class_sizes = factor(group_tags, 
                              levels = tags,
                              ordered = TRUE)) |>
  count(Habitat, class_sizes) |>
  group_by(Habitat) |>
  mutate(prop = prop.table(n))

ggplot(data = corals.sum, 
       aes(y = prop, x = Habitat)) +
  geom_bar(stat = 'Identity', #'Identity' means don't do anything on it, multiply by 1
           aes(fill = class_sizes), 
           color = 'black') + 
  theme_classic() +
  theme(panel.spacing.y = unit(10, 'pt'),
        legend.position = 'bottom') + labs(fill = 'Coral class sizes')  + 
  ylab('Proportion')


# Data analysis ----
data <- data |>
  mutate(Density = (ifelse(Density > 150, 150, Density)/100),
         Diameter = Diameter/100,
         Height = Height/100)

## Exploratory ----
ggplot(data = data) +
  geom_point(aes(y = Diameter, x = Height, colour = Habitat))

ggplot(data = data) +
  geom_point(aes(y = Diameter, x = Density, colour = Habitat))

ggplot(data = data) +
  geom_point(aes(y = Diameter, x = Shannon_algae, colour = Habitat))

ggplot(data = data) +
  geom_point(aes(y = Diameter, x = Shannon_coral, colour = Habitat))


data |> 
  dplyr::select(Density, Height, Shannon_algae, Shannon_coral) |> 
  cor() |> 
  corrplot()

ggplot(data = data, aes(Diameter)) +
  geom_density()


# Analysis

## MN1 - Height ----
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

## MN2 - Class sizes ----
dat2 <- corals |> 
  mutate(Habitat = factor(Habitat, levels = c('Flat', 'Crest', 'Slope'), ordered = TRUE)) |>
  mutate(class_sizes = factor(group_tags, 
                              levels = tags,
                              ordered = TRUE)) |>
  dplyr::select(Habitat, Transect, Quadrat, Taxa, class_sizes) |>
  left_join(canopy |> 
              mutate(Habitat = factor(Habitat, levels = c('Flat', 'Crest', 'Slope'), ordered = TRUE)) |>
              dplyr::filter(Trip == '3' | Trip == '2') |>
              dplyr::select(Habitat, Transect, Trip, Quadrat, Mean) |>
              group_by(Habitat, Transect, Quadrat) |>
              summarise(Mean = mean(Mean))) |>
  mutate(across(where(is.numeric), coalesce, 0)) |>
  mutate(Depth = ifelse(Habitat == 'Flat', '0-1 m',
                        ifelse(Habitat == 'Crest', '2-3 m',
                               '4-5 m')))

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

### Sargassum height ----

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

