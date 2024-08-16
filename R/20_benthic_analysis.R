# PACKAGES ----

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

source(file = 'R/functions.R')

# READ DATA ----
data <- read_xlsx(path = 'data/primary/CH3-recruits-count.xlsx', sheet = 1) |>
  full_join(read_xlsx(path = 'data/primary/CH3-metadata.xlsx', sheet = 2) |>
              dplyr::select(Tile, Treatment, T1_survey) |>
              dplyr::filter(T1_survey == 'Done')) |>
  full_join(read_xlsx(path = 'data/primary/CH3-tile-benthos.xlsx', sheet = 1)) |>
  group_by(Tile) |>
  mutate(Total = sum(Unbleached)) |>
  dplyr::distinct(Tile, .keep_all = TRUE) |>
  ungroup() |>
  dplyr::select(Tile, Grazing, Treatment, T1_survey,
                Sediment, Turf_height, Turf_cover, CCA_cover, Algae_cover, Total) |>
  left_join(read_xlsx(path = 'data/primary/CH3-Algal-community-time-0.xlsx', sheet = 1) |>
              dplyr::filter(Category != "Coral") |>
              dplyr::filter(Taxa != 'Zooanthids') |>
              dplyr::select(Time, Tile, Taxa, Cover) |>
              mutate(Cover = ifelse(Cover == "+", '0.5', Cover)) |> #adjust rare species 
              mutate(Cover = as.numeric(Cover)) |>
              group_by(Tile) |>
              mutate(Freq = Cover / sum(Cover)*100) |>
              ungroup() |>
              dplyr::select(-Cover) |>
              tidyr::pivot_wider(names_from = 'Taxa', values_from = 'Freq', values_fill = 0)) |>
  dplyr::filter(T1_survey == 'Done') |>
  dplyr::select(!T1_survey) |>
  mutate(Treatment = ifelse(Treatment == 'No algae', 'T1',
                            ifelse(Treatment == 'Only canopy', 'T2',
                                   ifelse(Treatment == 'Only mat', 'T3', 
                                          'T4'))))
#write.csv(data, file = 'data/processed/multivariate_data.csv')

# for summary:
benthos <- rbind(readxl::read_xlsx(path = "data/primary/Maggie-Quads-1.xlsx", sheet = "Florence"),
                 readxl::read_xlsx(path = "data/primary/Maggie-Quads-2.xlsx", sheet = "Florence"),
                 readxl::read_xlsx(path = "data/primary/Maggie-Quads-3.xlsx", sheet = "Florence"),
                 readxl::read_xlsx(path = "data/primary/Maggie-Quads-4.xlsx", sheet = "Florence"),
                 readxl::read_xlsx(path = "data/primary/Maggie-Quads-5.xlsx", sheet = "Florence"),
                 readxl::read_xlsx(path = "data/primary/Maggie-Quads-6.xlsx", sheet = "Florence"),
                 readxl::read_xlsx(path = "data/primary/Maggie-Quads-7.xlsx", sheet = "Florence"),
                 readxl::read_xlsx(path = "data/primary/Maggie-Quads-8.xlsx", sheet = "Florence")) |>
  dplyr::select(!c(Date))

canopy <- read.csv('data/primary/canopy_data.csv')

# MDS ----
## MDS Plot ----
## ---- MDS fit
#non-metric multidimensional scaling
data.mds <- metaMDS(data[,-c(1:10)], k=2,  plot=TRUE) #it always standardise by square root
data.mds #check that the stress is <0.2, the dimensions show to what reduction we drove the data
stressplot(data.mds) #N.B. is usually a non linear trend
## ----end

data.mds.scores <- data.mds |> 
  fortify() |> 
  mutate(Label = label,
         Score = score) |>
  full_join(data[, 1:10] |> add_rownames(var='Label'))
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

### Supplementary Figure 1 ----
nMDS_plot <- ggplot(data = NULL, aes(y = NMDS2, x = NMDS1)) + 
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
nMDS_plot

# save plot
ggsave(file = paste0(FIGS_PATH, "/supp_fig_1.png"), 
       plot = nMDS_plot, 
       width = 160, 
       height = 160/1.6, 
       units = "mm", 
       dpi = 300)

# plot showing just the taxa:
ggplot() + 
  geom_segment(data = data.mds.scores |> filter(Score == 'species'),
               aes(y = 0, x = 0, yend = NMDS2, xend = NMDS1)) +
  geom_text(data = data.mds.scores |> filter(Score == 'species'),
            aes(y = NMDS2*1.1, x = NMDS1*1.1, label=Label))  +
  theme_classic() +
  ylab("NMDS2") + xlab("NMDS1")

## ---- ANOVA
data.dist <- vegdist(wisconsin(data[,-c(1:10)]^0.25),"bray")
data.adonis<-adonis2(data.dist~Treatment,  data=data)
data.adonis
## ----end

Xmat <- model.matrix(~-1+Sediment+Turf_height+Algae_cover+CCA_cover, data=data) #need to dummy code because it's categorical variables
colnames(Xmat) <-gsub("Treatment","",colnames(Xmat))
envfit <- envfit(data.mds, env=Xmat)
envfit #shows how different each variable is from the centroid, so a significant values is whether that variable is moved from the 'general' community


# METHODS ----
data |>
  dplyr::select(!c(Grazing, Sediment, Turf_height, 
                   Turf_cover, CCA_cover, Algae_cover, Total, Time)) |>
  pivot_longer(!c(Tile, Treatment)) |>
  group_by(Treatment, name) |>
  summarise(mean = mean(value)) |>
  dplyr::filter(mean != 0) |>
  arrange(desc(mean), .by_group = TRUE) |>
  View()

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


# algal genera in each treatment:
data |>
  pivot_longer(names_to = 'Taxa', values_to = 'Cover', cols = 11:26) |>
  dplyr::select(Tile, Treatment, Taxa, Cover) |>
  group_by(Treatment, Taxa) |>
  summarise(avg = round(mean(Cover), 1)) |>
  arrange(Treatment, desc(avg)) |>
  View()

### Canopy height ----

for (i in 1:nrow(canopy)) {
  canopy[i, 'Mean'] <- canopy[i, 'Height'] |>
    str_split(', ') |>
    unlist() |>
    as.numeric() |>
    mean()
}

canopy <- canopy |>
  mutate(Mean = as.numeric(Mean))


canopy |>
  mutate(Habitat = factor(Habitat, levels = c('Flat', 'Crest', 'Slope'), ordered = TRUE)) |>
  dplyr::filter(Site == 'Florence') |>
  droplevels() |>
  ggplot(aes(x = Trip, y = Mean)) +
  geom_violin(aes(group = Trip)) +
  facet_wrap(~Habitat)

#### Supp. Fig. ----

canopy |>
  mutate(Habitat = factor(Habitat, levels = c('Flat', 'Crest', 'Slope'), ordered = TRUE)) |>
  dplyr::filter(Site == 'Florence') |>
  droplevels() |>
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
  dplyr::filter(Site == 'Florence') |>
  droplevels() |>
  dplyr::filter(Trip == '3') |>
  ggplot(aes(x = Habitat, y = Mean)) +
  geom_violin() +
  geom_point() +
  scale_y_continuous(expression(italic(Sargassum)*' height (cm)'))  +
  theme_classic()

ggsave(file = paste0(FIGS_PATH, "/SF2_canopy_height_spawning.png"), 
       width = 250, 
       height = 250/1.6, 
       units = "mm", 
       dpi = 300)

