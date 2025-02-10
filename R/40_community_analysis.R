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
              dplyr::select(Tile, Taxa, Cover) |>
              mutate(Cover = ifelse(Cover == "+", '1', Cover)) |> #adjust rare species 
              mutate(Cover = as.numeric(Cover)) |>
              group_by(Tile, Taxa) |>
              summarise(Cover = sum(Cover)) |>
              ungroup() |>
              group_by(Tile) |>
              mutate(Freq = Cover / sum(Cover)*100) |>
              ungroup() |>
              dplyr::select(-Cover) |>
              tidyr::pivot_wider(names_from = 'Taxa', values_from = 'Freq', values_fill = 0))

data |>
  glimpse()

data[,-c(1:8)]^0.25 |> #square root
  ggpairs(lower = list(continuous = "smooth"),
          diag = list(continuous = "density"),
          axisLabels = "show")

data.std <- (data[,-c(1:8)]^0.25) |> 
  wisconsin()
data.std

## PCA ----
data.rda <- rda(data.std, 
                scale = FALSE) 

summary(data.rda, 
        display = NULL)

data.rda$tot.chi/data.rda$CA$v |> nrow() #threshold for what to keep

autoplot(data.rda)
autoplot(data.rda) + theme_bw()
autoplot(data.rda, geom = 'text') + theme_bw()

data.rda.scores <- data.rda |> 
  fortify() |> #prepare for ggplot
  full_join(data[, 1:3] |> add_column(data.rda |> 
                                        fortify() |> 
                                        as.data.frame() |> 
                                        filter(score == 'sites') |> 
                                        dplyr::select(label)), 
            by = 'label')

ggplot(data = NULL, aes(y = PC2, x = PC1)) + 
  geom_point(data = data.rda.scores |> 
               filter(score=='sites'),
             mapping = aes(size = Total), alpha = 0.7) +
  geom_segment(data = data.rda.scores |> filter(score == 'species') |> filter(label != 'Other'),
               aes(y = 0, x = 0, yend = PC2, xend = PC1), alpha = 0.6) +
  geom_text_repel(data = data.rda.scores |> filter(score == 'species') |> filter(label != 'Other'),
                  aes(y = PC2*1.1, x = PC1*1.1, label=label), alpha = 0.7) +
  #scale_x_continuous(limits = c(-0.7, 0.7)) +
  #scale_y_continuous(limits = c(-0.7, 0.7)) +
  scale_size('Coral recruitment \n (recruits/tile)') +
  geom_hline(aes(yintercept = 0), linetype = 2, alpha = 0.3) +
  geom_vline(aes(xintercept = 0), linetype = 2, alpha = 0.3) +
  theme_classic() + theme(legend.position = c(0.9, 0.85), 
                          legend.title = element_text(hjust = 1))

## CCA ----
data.cca <- cca(data.std ~ H_mean_broad  + Shannon_broad_alg, data=data[, c(1:7)], scale=FALSE)

summary(data.cca, display=NULL)
anova(data.cca)

autoplot(data.cca)
vif.cca(data.cca)
#overall test
anova(data.cca)
anova(data.cca, by='axis')
anova(data.cca, by='margin')

coef(data.cca)

RsquareAdj(data.cca)

screeplot(data.cca)

data.cca.scores <- data.cca |> 
  fortify() |> #prepare for ggplot
  full_join(data[, 1:7] |> add_column(data.cca |> 
                                             fortify() |> 
                                             as.data.frame() |> 
                                             filter(score == 'sites') |> 
                                             dplyr::select(label)), 
            by = 'label')

ggplot(data = NULL, aes(y = CCA2, x = CCA1)) + 
  geom_point(data = data.cca.scores |> 
               filter(score=='sites'),
             mapping = aes(size = Total), alpha = 0.7) +
  geom_segment(data = data.cca.scores |> filter(score == 'species') |> filter(label != 'Other'),
               aes(y = 0, x = 0, yend = CCA2, xend = CCA1), alpha = 0.6) +
  geom_text_repel(data = data.cca.scores |> filter(score == 'species') |> filter(label != 'Other'),
                  aes(y = CCA2*1.1, x = CCA1*1.1, label=label), alpha = 0.7) +
  scale_size('Coral recruitment \n (recruits/tile)') +
  geom_hline(aes(yintercept = 0), linetype = 2, alpha = 0.3) +
  geom_vline(aes(xintercept = 0), linetype = 2, alpha = 0.3) +
  theme_classic() + theme(legend.position = c(0.9, 0.85), 
                          legend.title = element_text(hjust = 1),
                          legend.background = element_blank()) +
  scale_y_continuous(limits = c(-2.2, 2.2)) +
  scale_x_continuous(limits = c(-2.2, 2.2)) + 
  geom_segment(data = data.cca.scores |> 
                 filter(score=='centroids'),
               aes(y=0, x=0, yend=CCA2, xend=CCA1),
               arrow=arrow(length=unit(0.3,'lines')), color='blue') +
  geom_text(data=data.cca.scores |> 
              filter(score=='centroids'),
            aes(y=CCA2*1.1, x=CCA1*1.1, label=label), color='blue')

data[, c(2, 5:7)] |>
  cor() |> 
  corrplot(type = 'upper',
           order = 'FPC',
           diag = FALSE)

data[, c(2, 5:7)] |> 
  ggpairs(lower = list(continuous = "smooth"),
          diag = list(continuous = "density"),
          axisLabels = "show")

data.envfit <- envfit(data.cca, env = data[, c(2:3,5:7)])
data.envfit 

data.env.scores <- data.envfit |> 
  fortify()

cca1 <- data.cca.scores |> 
  filter(score=='sites') |> 
  pull(CCA1)
cca2 <- data.cca.scores |> 
  filter(score=='sites') |> 
  pull(CCA2)


lm(1:nrow(data[, c(2, 5:7)]) ~ H_mean_broad  + D_broad  + Shannon_broad_alg , 
   data = data[, c(2, 5:7)]) |>
  vif() #to calculate which predictors cannot be modeled together
#remove some:
lm(1:nrow(data[, c(2, 5:7)]) ~ H_mean_broad  + Shannon_broad_alg, 
   data = data[, c(2, 5:7)]) |>
  vif()

#you could run it Bayesian!
lm(cca1 ~ H_mean_broad  + Shannon_broad_alg, 
   data = data[, c(2, 5:7)]) |> 
  summary()

lm(cca2 ~ H_mean_broad + Shannon_broad_alg , 
   data = data[, c(2, 5:7)]) |> 
  summary()

lm(cca1 ~ Shannon_broad_alg, 
   data = data[, c(2, 5:7)]) |> 
  summary()

lm(cca2 ~ Shannon_broad_alg, 
   data = data[, c(2, 5:7)]) |> 
  summary()


#Genus analysis ----
scatterplotMatrix(~Total+Sargassum+Colpomenia+Lobophora+Hypnea+Padina+Dictyota+Halimeda+Galaxaura+Halymenia+Dictyopteris+Gelidiopsis+Asparagopsis, 
                  data = data,
                  diagonal = list(method='boxplot'))

## Model ----

vif(lm(Total ~ Sargassum + Lobophora + Hypnea + Colpomenia + Turf, data = data))

data |> dplyr::filter(Sargassum != 0) |> dplyr::select(Sargassum) |> min()/2
data |> dplyr::filter(Lobophora != 0) |> dplyr::select(Lobophora) |> min()/2
data |> dplyr::filter(Hypnea != 0) |> dplyr::select(Hypnea) |> min()/2
data |> dplyr::filter(Colpomenia != 0) |> dplyr::select(Colpomenia) |> min()/2
data |> dplyr::filter(Galaxaura != 0) |> dplyr::select(Galaxaura) |> min()/2
data |> dplyr::filter(Dictyota != 0) |> dplyr::select(Dictyota) |> min()/2

### Fit model ----

## ---- MZI7Fit
data.form <- bf(Total ~ scale(log(Sargassum + 5)) + scale(log(Lobophora + 0.49)) +
                  scale(log(Hypnea + 0.49)) + scale(log(Colpomenia + 0.47)) +
                  (1|Grazing), 
                zi ~ 1, 
                family = zero_inflated_poisson(link = 'log'))
## ----end


data.form |> get_prior(data = data)

priors <- prior(normal(0.5, 3), class = 'Intercept') +
  prior(normal(0, 5), class = 'b')  +
  prior(student_t(3, 0, 2), class = 'sd') +
  prior(logistic(0, 1), class = 'Intercept', dpar = 'zi')

dataZI.brm <- brm(data.form, prior = priors, data = data, 
                  sample_prior = 'yes', 
                  iter = 5000, 
                  warmup = 1000, 
                  chains = 3, cores = 3, 
                  control = list(adapt_delta = 0.99, 
                                 max_treedepth = 20), 
                  thin = 5, 
                  refresh = 100, 
                  backend = 'rstan') 

### Diagnostics ----
dataZI.brm |> 
  SUYR_prior_and_posterior() +
  theme_classic() + theme(text = element_text(colour = 'black'), 
                          axis.text = element_text(size = rel(1.2)),
                          axis.title = element_text(size = rel(1.5)),
                          legend.text = element_text(size = rel(1.2)),
                          legend.title = element_text(size = rel(1.5)),
                          legend.position = 'bottom')


(dataZI.brm$fit |> stan_trace()) + (dataZI.brm$fit |> stan_ac()) + 
  (dataZI.brm$fit |> stan_rhat()) + (dataZI.brm$fit |> stan_ess())


## Model validation
### Posterior probability check
dataZI.brm |> 
  pp_check(type = 'dens_overlay', ndraws = 100)


### Residuals
data.resids <- make_brms_dharma_res(dataZI.brm, 
                                    integerResponse = TRUE)
testUniformity(data.resids)
plotResiduals(data.resids, form = factor(rep(1, nrow(data))))
plotResiduals(data.resids, quantreg = TRUE) 
testDispersion(data.resids)

### Save model ----
#save(dataZI.brm, data.form, priors, data, file = 'data/modelled/MZI7_algae.RData')

### Investigation ----
dataZI.brm |> 
  conditional_effects() |> 
  plot(points = TRUE, ask = FALSE)

### ---- MZI7Output
dataZI.brm |>
  brms::as_draws_df() |>
  dplyr::select(starts_with('b_')) |>
  mutate(across(everything(), exp)) |>
  summarise_draws(median,
                  HDInterval::hdi,
                  rhat, length, ess_bulk, ess_tail,
                  Pl = ~mean(.x < 1),
                  Pg = ~mean(.x > 1))
### ----end

### Figure ----

Sarg <- dataZI.brm |> 
  emmeans(~Sargassum, at = with(data,
                                   list(Sargassum = seq(min(Sargassum),
                                                           max(Sargassum),
                                                           len = 50)))) |>
  gather_emmeans_draws() |>
  mutate(.value = exp(.value)) |>
  group_by(Sargassum) |>
  summarise(median_hdci(.value)) |>
  as.data.frame() |> 
  ggplot(aes(x = Sargassum,
             y = y)) + 
  geom_line() + 
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.2) +
  geom_point(data = data, aes(x = Sargassum, y = Total), alpha = 0.4, size = 2, 
             position = position_jitter(width = 0.1, height = 0.1)) +
  scale_x_continuous(name = expression(paste(italic('Sargassum'), ' cover (%)'))) +
  scale_y_continuous('Coral recruitment (recruits/tile)') +
  labs(tag = 'a')  +
  theme_classic()  +
  theme(text = element_text(colour = 'black'), 
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.5)))

Lobo <- dataZI.brm |> 
  emmeans(~Lobophora, at = with(data,
                                list(Lobophora = seq(min(Lobophora),
                                                     max(Lobophora),
                                                     len = 50)))) |>
  gather_emmeans_draws() |>
  mutate(.value = exp(.value)) |>
  group_by(Lobophora) |>
  summarise(median_hdci(.value)) |>
  as.data.frame() |> 
  ggplot(aes(x = Lobophora,
             y = y)) + 
  geom_line() + 
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.2) +
  geom_point(data = data, aes(x = Lobophora, y = Total), alpha = 0.4, size = 2, 
             position = position_jitter(width = 0.1, height = 0.1)) +
  scale_x_continuous(name = expression(paste(italic('Lobophora'), ' cover (%)'))) +
  scale_y_continuous('Coral recruitment (recruits/tile)') +
  labs(tag = 'b')  +
  theme_classic()  +
  theme(text = element_text(colour = 'black'), 
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.5)))

Hyp <- dataZI.brm |> 
  emmeans(~Hypnea, at = with(data,
                                list(Hypnea = seq(min(Hypnea),
                                                     max(Hypnea),
                                                     len = 50)))) |>
  gather_emmeans_draws() |>
  mutate(.value = exp(.value)) |>
  group_by(Hypnea) |>
  summarise(median_hdci(.value)) |>
  as.data.frame() |> 
  ggplot(aes(x = Hypnea,
             y = y)) + 
  geom_line() + 
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.2) +
  geom_point(data = data, aes(x = Hypnea, y = Total), alpha = 0.4, size = 2, 
             position = position_jitter(width = 0.1, height = 0.1)) +
  scale_x_continuous(name = expression(paste(italic('Hypnea'), ' cover (%)'))) +
  scale_y_continuous('Coral recruitment (recruits/tile)') +
  labs(tag = 'c')  +
  theme_classic()  +
  theme(text = element_text(colour = 'black'), 
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.5)))

Colp <- dataZI.brm |> 
  emmeans(~Colpomenia, at = with(data,
                             list(Colpomenia = seq(min(Colpomenia),
                                               max(Colpomenia),
                                               len = 50)))) |>
  gather_emmeans_draws() |>
  mutate(.value = exp(.value)) |>
  group_by(Colpomenia) |>
  summarise(median_hdci(.value)) |>
  as.data.frame() |> 
  ggplot(aes(x = Colpomenia,
             y = y)) + 
  geom_line() + 
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.2) +
  geom_point(data = data, aes(x = Colpomenia, y = Total), alpha = 0.4, size = 2, 
             position = position_jitter(width = 0.1, height = 0.1)) +
  scale_x_continuous(name = expression(paste(italic('Colpomenia'), ' cover (%)'))) +
  scale_y_continuous('Coral recruitment (recruits/tile)') +
  labs(tag = 'd')  +
  theme_classic()  +
  theme(text = element_text(colour = 'black'), 
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.5)))

(Sarg + Lobo + Hyp + Colp) +
  plot_layout(axis_titles = 'collect', nrow = 2, )

ggsave(file = paste0(FIGS_PATH, "/MZI_algae_fig.png"), 
       width = 200, 
       height = 200/1.2, 
       units = "mm", 
       dpi = 300)

# Functional diversity ----
## Load data ----
func_data <- read_csv(file = 'data/processed/recruit.csv', col_select = -1) |>
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
              dplyr::select(Tile, Taxa, Cover) |>
              mutate(Cover = ifelse(Cover == "+", '1', Cover)) |> #adjust rare species 
              mutate(Cover = as.numeric(Cover)) |>
              left_join(read.csv(file = 'data/primary/categories.csv')) |> 
              mutate(spp = ifelse(!is.na(Functional.group), Functional.group, Taxa)) |> 
              group_by(Tile, spp) |>
              summarise(Cover = sum(Cover)) |>
              ungroup() |>
              group_by(Tile) |>
              mutate(Freq = Cover / sum(Cover)*100) |>
              ungroup() |>
              dplyr::select(-Cover) |>
              tidyr::pivot_wider(names_from = 'spp', values_from = 'Freq', values_fill = 0))

func_data |>
  glimpse()

func_data[,-c(1:7)]^0.25 |> #square root
  ggpairs(lower = list(continuous = "smooth"),
          diag = list(continuous = "density"),
          axisLabels = "show")

data.std <- (func_data[,-c(1:7)]^0.25) |> 
  wisconsin()
data.std

## PCA ----
data.rda <- rda(data.std, 
                scale = FALSE) 

summary(data.rda, 
        display = NULL)

data.rda$tot.chi/data.rda$CA$v |> nrow() #threshold for what to keep

autoplot(data.rda)
autoplot(data.rda) + theme_bw()
autoplot(data.rda, geom = 'text') + theme_bw()

data.rda.scores <- data.rda |> 
  fortify() |> #prepare for ggplot
  full_join(func_data[, 1:3] |> add_column(data.rda |> 
                                        fortify() |> 
                                        as.data.frame() |> 
                                        filter(score == 'sites') |> 
                                        dplyr::select(label)), 
            by = 'label')

ggplot(data = NULL, aes(y = PC2, x = PC1)) + 
  geom_point(data = data.rda.scores |> 
               filter(score=='sites'),
             mapping = aes(size = Total), alpha = 0.7) +
  geom_segment(data = data.rda.scores |> filter(score == 'species') |> filter(label != 'Other'),
               aes(y = 0, x = 0, yend = PC2, xend = PC1), alpha = 0.6) +
  geom_text_repel(data = data.rda.scores |> filter(score == 'species') |> filter(label != 'Other'),
                  aes(y = PC2*1.1, x = PC1*1.1, label=label), alpha = 0.7) +
  scale_size('Coral recruitment \n (recruits/tile)') +
  geom_hline(aes(yintercept = 0), linetype = 2, alpha = 0.3) +
  geom_vline(aes(xintercept = 0), linetype = 2, alpha = 0.3) +
  theme_classic() + theme(legend.position = c(0.9, 0.85), 
                          legend.title = element_text(hjust = 1)) +
  scale_y_continuous(paste(names(data.rda$CA$eig[2]), sprintf('(%0.1f%% explained var.)',
                                                  100 * data.rda$CA$eig[2]/data.rda$CA$tot.chi)),
                     limits = c(-0.5, 0.5))+
  scale_x_continuous(paste(names(data.rda$CA$eig[1]), sprintf('(%0.1f%% explained var.)',
                                                  100 * data.rda$CA$eig[1]/data.rda$CA$tot.chi)),
                     limits = c(-0.5, 0.5))

func_data[, c(5:7)] |>
  cor() |> 
  corrplot(type = 'upper',
           order = 'FPC',
           diag = FALSE)

func_data[, c(5:7)] |> 
  ggpairs(lower = list(continuous = "smooth"),
          diag = list(continuous = "density"),
          axisLabels = "show")

data.envfit <- envfit(data.rda, env = func_data[, c(5:7)])
data.envfit 

data.env.scores <- data.envfit |> 
  fortify()

ggplot(data = NULL, aes(y = PC2, x = PC1)) + 
  geom_point(data = data.rda.scores |> 
               filter(score=='sites'),
             mapping = aes(size = Total), alpha = 0.7) +
  geom_segment(data = data.rda.scores |> filter(score == 'species') |> filter(label != 'Other'),
               aes(y = 0, x = 0, yend = PC2, xend = PC1), alpha = 0.6) +
  geom_text_repel(data = data.rda.scores |> filter(score == 'species') |> filter(label != 'Other'),
                  aes(y = PC2*1.1, x = PC1*1.1, label=label), alpha = 0.7) +
  #scale_x_continuous(limits = c(-0.7, 0.7)) +
  #scale_y_continuous(limits = c(-0.7, 0.7)) +
  scale_size('Coral recruitment \n (recruits/tile)') +
  geom_hline(aes(yintercept = 0), linetype = 2, alpha = 0.3) +
  geom_vline(aes(xintercept = 0), linetype = 2, alpha = 0.3) +
  theme_classic() + theme(legend.position = c(0.9, 0.85), 
                          legend.title = element_text(hjust = 1)) +
  scale_y_continuous(paste(names(data.rda$CA$eig[2]), sprintf('(%0.1f%% explained var.)',
                                                              100 * data.rda$CA$eig[2]/data.rda$CA$tot.chi)))+
  scale_x_continuous(paste(names(data.rda$CA$eig[1]), sprintf('(%0.1f%% explained var.)',
                                                              100 * data.rda$CA$eig[1]/data.rda$CA$tot.chi))) + 
  geom_segment(data = data.env.scores,
               aes(y = 0, x = 0, yend = PC2, xend = PC1),
               arrow = arrow(length = unit(0.3,'lines')), color = 'blue') +
  geom_text(data = data.env.scores,
            aes(y = PC2*1.1, x = PC1*1.1, label = label), color = 'blue')

pc1 <- data.rda.scores |> 
  filter(score=='sites') |> 
  pull(PC1)
pc2 <- data.rda.scores |> 
  filter(score=='sites') |> 
  pull(PC2)


lm(1:nrow(func_data[, c(5:7)]) ~ H_mean_broad  + D_broad  + Shannon_broad_alg , 
   data = func_data[, c(5:7)]) |>
  vif() #to calculate which predictors cannot be modeled together
#remove some:
lm(1:nrow(func_data[, c(5:7)]) ~ D_broad  + Shannon_broad_alg , 
   data = func_data[, c(5:7)]) |>
  vif()

lm(pc1 ~ H_mean_broad  + D_broad  + Shannon_broad_alg, 
   data = func_data[, c(5:7)]) |> 
  summary()

lm(pc2 ~ H_mean_broad  + D_broad  + Shannon_broad_alg, 
   data = func_data[, c(5:7)]) |> 
  summary()

## CCA ----
data.cca <- cca(data.std ~ H_mean_broad + Shannon_broad_alg, data=func_data[, c(1:7)], scale=FALSE)

summary(data.cca, display=NULL)
anova(data.cca)

autoplot(data.cca)
vif.cca(data.cca)
#overall test
anova(data.cca)
anova(data.cca, by='axis')
anova(data.cca, by='margin')

coef(data.cca)

RsquareAdj(data.cca)

screeplot(data.cca)

data.cca.scores <- data.cca |> 
  fortify() |> #prepare for ggplot
  full_join(func_data[, 1:7] |> add_column(data.cca |> 
                                             fortify() |> 
                                             as.data.frame() |> 
                                             filter(score == 'sites') |> 
                                             dplyr::select(label)), 
            by = 'label')

ggplot(data = NULL, aes(y = CCA2, x = CCA1)) + 
  geom_point(data = data.cca.scores |> 
               filter(score=='sites'),
             mapping = aes(size = Total), alpha = 0.7) +
  geom_segment(data = data.cca.scores |> filter(score == 'species') |> filter(label != 'Other'),
               aes(y = 0, x = 0, yend = CCA2, xend = CCA1), alpha = 0.6) +
  geom_text_repel(data = data.cca.scores |> filter(score == 'species') |> filter(label != 'Other'),
                  aes(y = CCA2*1.1, x = CCA1*1.1, label=label), alpha = 0.7) +
  scale_size('Coral recruitment \n (recruits/tile)') +
  geom_hline(aes(yintercept = 0), linetype = 2, alpha = 0.3) +
  geom_vline(aes(xintercept = 0), linetype = 2, alpha = 0.3) +
  theme_classic() + theme(legend.position = c(0.9, 0.85), 
                          legend.title = element_text(hjust = 1),
                          legend.background = element_blank()) +
  scale_y_continuous(limits = c(-2.2, 2.2))+
  scale_x_continuous(limits = c(-2.2, 2.2))

func_data[, c(2, 5:7)] |>
  cor() |> 
  corrplot(type = 'upper',
           order = 'FPC',
           diag = FALSE)

func_data[, c(2, 5:7)] |> 
  ggpairs(lower = list(continuous = "smooth"),
          diag = list(continuous = "density"),
          axisLabels = "show")

data.envfit <- envfit(data.cca, env = func_data[, c(2, 5:7)])
data.envfit 

data.env.scores <- data.envfit |> 
  fortify()

cca1 <- data.cca.scores |> 
  filter(score=='sites') |> 
  pull(CCA1)
cca2 <- data.cca.scores |> 
  filter(score=='sites') |> 
  pull(CCA2)


lm(1:nrow(data[, c(2, 5:7)]) ~ H_mean_broad  + D_broad  + Shannon_broad_alg , 
   data = data[, c(2, 5:7)]) |>
  vif() #to calculate which predictors cannot be modeled together
#remove some:
lm(1:nrow(data[, c(2, 5:7)]) ~ H_mean_broad  + Shannon_broad_alg, 
   data = data[, c(2, 5:7)]) |>
  vif()

#you could run it Bayesian!
lm(cca1 ~ H_mean_broad  + Shannon_broad_alg, 
   data = data[, c(2, 5:7)]) |> 
  summary()

lm(cca2 ~ H_mean_broad + Shannon_broad_alg , 
   data = data[, c(2, 5:7)]) |> 
  summary()

lm(cca1 ~ Shannon_broad_alg, 
   data = data[, c(2, 5:7)]) |> 
  summary()

lm(cca2 ~ Shannon_broad_alg, 
   data = data[, c(2, 5:7)]) |> 
  summary()

