library(sf)
library(dplyr)
library(ozmaps)
library(ggplot2)
library(ggrepel)
library(scales)
library(OpenStreetMap)


labs <- data.frame(
  long = c(145.754120, 146.816956, 153.021072),
  lat = c(-16.925491, -19.258965, -27.470125),
  cities = c("Cairns", "Townsville", "Brisbane"),
  stringsAsFactors = FALSE) 

sf_oz <- ozmap("states")
subset(sf_oz, NAME == "Queensland") |>
  ggplot() + geom_sf() + 
  geom_point(data = labs, aes(long, lat), size = 2) +
  geom_text_repel(data = labs, aes(long, lat, label = cities), size = 5) +
  ylab("") + xlab("") +
  theme_classic() +
  theme(text = element_text(colour = 'black'), 
        axis.text = element_text(size = rel(1.5)))


labs <- data.frame(
  long = c(146.816956, 146.8809),
  lat = c(-19.258965, -19.1228),
  cities = c("Townsville", "Florence Bay"),
  stringsAsFactors = FALSE) 

lga_sf <- ozmap_data("abs_lga")
subset(lga_sf, NAME == "Townsville (C)") |>
 ggplot() + geom_sf() +
  coord_sf(
    xlim = c(146.7, 147.05),
    ylim = c(-19.32, -19.05)
  ) + 
  geom_point(data = labs, aes(long, lat), size = 2) +
  geom_text_repel(data = labs, aes(long, lat, label = cities), size = 5) +
  ylab("") + xlab("") +
  theme_classic() +
  theme(text = element_text(colour = 'black'), 
        axis.text = element_text(size = rel(1.5)))
