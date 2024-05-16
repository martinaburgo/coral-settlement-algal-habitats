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
  mutate(Treatment = ifelse(Treatment == 'No algae', 'Control',
                            ifelse(Treatment == 'Only canopy', 'Canopy-forming',
                                   ifelse(Treatment == 'Only mat', 'Mat-forming', Treatment))))

# MDS ----
## MDS Plot ----
## ---- MDS fit
#non-metric multidimensional scaling
data.mds <- metaMDS(data[,-c(1:10)], k=2,  plot=TRUE) #it always standardise by square root
data.mds #check that the stress is <0.2, the dimensions show to what reduction we drove the data
stressplot(data.mds) #N.B. is usually a non linear trend
## ----end

data.env.scores <- envfit |>
  fortify() |> 
  mutate(Label = label)

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

col_vals <- c("All algae" = "#8dd3c7", 
              "Control" = "#fb8072", 
              "Canopy-forming" = "#d1a95e",
              "Mat-forming" = "#a0c58b")

### Supplementary Figure 1 ----
nMDS_plot <- ggplot(data = NULL, aes(y = NMDS2, x = NMDS1)) + 
    geom_segment(data = data.mds.scores |> filter(Score == 'species'),
                 aes(y = 0, x = 0, yend = NMDS2, xend = NMDS1), alpha = 0.7) +
    geom_text_repel(data = data.mds.scores |> filter(Score == 'species'),
              aes(y = NMDS2*1.1, x = NMDS1*1.1, label=Label), alpha = 0.7) +
  geom_point(data=data.mds.scores |> 
               filter(Score=='sites'),
             aes(color = Treatment, shape = Treatment), size = 2, alpha = 0.6) +
    ggforce::geom_mark_ellipse(data = data.mds.scores |> 
                                 filter(Score=='sites'),
                               aes(y = NMDS2, x = NMDS1, 
                                   fill = Treatment, 
                                   colour = Treatment), expand = 0, alpha = 0.2) + 
    scale_colour_manual(values = col_vals) + scale_fill_manual(values = col_vals) +
    theme_classic() + theme(legend.position = c(0.9, 0.2))
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

## ANOVA ----

data.dist <- vegdist(wisconsin(data[,-c(1:10)]^0.25),"bray")
data.adonis<-adonis2(data.dist~Treatment,  data=data)
data.adonis


Xmat <- model.matrix(~-1+Sediment+Turf_height+Algae_cover+CCA_cover, data=data) #need to dummy code because it's categorical variables
colnames(Xmat) <-gsub("Treatment","",colnames(Xmat))
envfit <- envfit(data.mds, env=Xmat)
envfit #shows how different each variable is from the centroid, so a significant values is whether that variable is moved from the 'general' community
