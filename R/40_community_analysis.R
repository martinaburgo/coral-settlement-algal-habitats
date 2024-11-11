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
data.cca <- cca(data.std~H_mean_broad + D_broad + Shannon_broad_alg, data=data[, c(1:7)], scale=FALSE)

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

data.envfit <- envfit(data.cca, env = func_data[, c(2,5:7)])
data.envfit 

data.env.scores <- data.envfit |> 
  fortify()

cca1 <- data.cca.scores |> 
  filter(score=='sites') |> 
  pull(CCA1)
cca2 <- data.cca.scores |> 
  filter(score=='sites') |> 
  pull(CCA2)


lm(1:nrow(func_data[, c(2, 5:7)]) ~ H_mean_broad  + D_broad  + Shannon_broad_alg , 
   data = func_data[, c(2, 5:7)]) |>
  vif() #to calculate which predictors cannot be modeled together
#remove some:
lm(1:nrow(func_data[, c(2, 5:7)]) ~ H_mean_broad  + D_broad, 
   data = func_data[, c(2, 5:7)]) |>
  vif()

#you could run it Bayesian!
lm(cca1 ~ H_mean_broad  +D_broad, 
   data = func_data[, c(2, 5:7)]) |> 
  summary()

lm(cca2 ~ H_mean_broad  +D_broad , 
   data = func_data[, c(2, 5:7)]) |> 
  summary()

lm(cca1 ~ Shannon_broad_alg, 
   data = func_data[, c(2, 5:7)]) |> 
  summary()

lm(cca2 ~ Shannon_broad_alg, 
   data = func_data[, c(2, 5:7)]) |> 
  summary()