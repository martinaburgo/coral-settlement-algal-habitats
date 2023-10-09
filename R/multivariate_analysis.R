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


data <- read_xlsx(path = 'data/primary/CH3-recruits-count.xlsx', sheet = 1) |>
  full_join(read_xlsx(path = 'data/primary/CH3-metadata.xlsx', sheet = 2) |> #adding treatments
              dplyr::select(Tile, Treatment)) |>
  group_by(Tile) |>
  mutate(Total = sum(Unbleached)) |>
  dplyr::distinct(Tile, .keep_all = TRUE) |>
  ungroup() |>
  dplyr::select(Tile, Grazing, Treatment, Total) |> 
  full_join(read_xlsx(path = 'data/primary/CH3-Algal-community-time-0.xlsx', sheet = 1) |>
  dplyr::select(Tile, Taxa, Cover) |>
  mutate(Cover = ifelse(Cover == "+", '0.5', Cover)) |>
  mutate(Cover = as.numeric(Cover)) |> #adjust rare species
  group_by(Tile) |>
  mutate(Freq = Cover / sum(Cover)*100) |>
  ungroup() |>
  dplyr::select(-Cover) |>
  pivot_wider(names_from = Taxa, values_from = Freq, values_fill = 0)) |>
  dplyr::filter(!is.na(Total))

head(data)

## ---- MDS fit
#non-metric multidimensional scaling
data.mds <- metaMDS(data[,-c(1:4)], k=2,  plot=TRUE) #it always standardise by square root
data.mds #check that the stress is <0.2, the dimensions show to what reduction we drove the data
stressplot(data.mds) #N.B. is usually a non linear trend
## ----end


# ---- MDS plot part 1
data.mds.scores <- data.mds |> 
  fortify() |> 
  mutate(Label = label,
         Score = score) |>
  full_join(data |> add_rownames(var='Label'))
data.mds.scores.centroids <- data.mds.scores |>
  filter(Score == "sites") |>
  group_by(Treatment) |>
  summarise(across(c(NMDS1, NMDS2), list(c = mean)))
data.mds.scores <- data.mds.scores |>
  full_join(data.mds.scores.centroids)

g <-
  ggplot(data = NULL, aes(y=NMDS2, x=NMDS1)) +
  geom_hline(yintercept=0, linetype='dotted') +
  geom_vline(xintercept=0, linetype='dotted') +
  geom_point(data=data.mds.scores %>% filter(Score=='sites'),
             aes(color=Treatment)) +
  geom_text(data=data.mds.scores %>% filter(Score=='sites'),
            aes(label=Label, color=Treatment), hjust=-0.2) +
  geom_segment(data=data.mds.scores %>% filter(Score=='species'),
               aes(y=0, x=0, yend=NMDS2, xend=NMDS1),
               arrow=arrow(length=unit(0.3,'lines')), color='red',
               alpha =  0.2) +
  geom_text(data=data.mds.scores %>% filter(Score=='species'),
            aes(y=NMDS2*1.1, x=NMDS1*1.1, label=Label), color='red',
            alpha =  0.2) + ggforce::geom_mark_ellipse(data=data.mds.scores %>% filter(Score=='sites'),
                                                       aes(y=NMDS2, x=NMDS1, fill=Treatment), expand=0) + 
  geom_segment(data = data.mds.scores,
               aes(x = NMDS1_c, xend = NMDS1, y = NMDS2_c, yend = NMDS2, colour = Treatment))
g
## ----end

Xmat <- model.matrix(~-1+Treatment, data=data) #need to dummy code because it's categorical variables
colnames(Xmat) <-gsub("Treatment","",colnames(Xmat))
envfit <- envfit(data.mds, env=Xmat)
envfit #shows how different each variable is from the centroid, so a significsnt values is whethet that variable is moved from the 'general' community

data.env.scores <- envfit |>
  fortify() |> 
  mutate(Label = label)

g <- ggplot(data = NULL, aes(y=NMDS2, x=NMDS1)) +
  geom_hline(yintercept=0, linetype='dotted') +
  geom_vline(xintercept=0, linetype='dotted') +
  geom_point(data=data.mds.scores %>% filter(Score=='sites'),
             aes(color=Treatment)) +
  geom_text(data=data.mds.scores %>% filter(Score=='sites'),
            aes(label=Label, color=Treatment), hjust=-0.2) +
  geom_segment(data=data.mds.scores %>% filter(Score=='species'),
               aes(y=0, x=0, yend=NMDS2, xend=NMDS1),
               arrow=arrow(length=unit(0.3,'lines')), color='red',
               alpha =  0.2) +
  geom_text(data=data.mds.scores %>% filter(Score=='species'),
            aes(y=NMDS2*1.1, x=NMDS1*1.1, label=Label), color='red',
            alpha =  0.2) + ggforce::geom_mark_ellipse(data=data.mds.scores %>% filter(Score=='sites'),
                                                       aes(y=NMDS2, x=NMDS1, fill=Treatment), expand=0) + 
  geom_segment(data = data.mds.scores,
               aes(x = NMDS1_c, xend = NMDS1, y = NMDS2_c, yend = NMDS2, colour = Treatment)) + 
  geom_segment(data=data.env.scores,
               aes(y=0, x=0, yend=NMDS2, xend=NMDS1),
               arrow=arrow(length=unit(0.3,'lines')), color='blue') +
  geom_text(data=data.env.scores,
            aes(y=NMDS2*1.1, x=NMDS1*1.1, label=Label), color='blue')
g