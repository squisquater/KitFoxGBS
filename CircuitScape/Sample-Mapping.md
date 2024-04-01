# Mapping Sampling Units for Landscape Genetics

This script takes a file of genetic sampels and groups them by Sampling Unit (i.e. Bakersfield, Panoche, etc) and then generates a single sampling coordinate for each unit that is a centroid of all datapoints associated with that unit. 
```
library(dplyr)
library(sf)
library(rnaturalearth)
library(ggplot2)
library(readr)
library(tidyr)

##set working directory
setwd("/Users/sophiepq/Desktop/KitFox/Manuscript-GBS")

##input samples with location data
data <- read.csv("StructureResults_K2-7_20xn230.csv", header = TRUE)

##remove missing data
data <- data %>%
  filter(!is.na(Lat))

##make sure lat/long data is numeric
data$Long <- as.numeric(data$Long)
data$Lat <- as.numeric(data$Lat)

# Convert data to an sf object
data_sf <- st_as_sf(data, coords = c("Long", "Lat"), crs = 4326)

# Calculate centroids for sampling regions
centroids <- data_sf %>%
  group_by(MiniPip) %>%
  summarize(geometry = st_centroid(st_combine(geometry))) %>%
  mutate(longitude = st_coordinates(geometry)[, 1],
         latitude = st_coordinates(geometry)[, 2])
         
## generate the basemap
NAmap <- ne_states(c("united states of america"), returnclass = "sf")
NAmap <- st_transform(NAmap, st_crs(4326))

## constrain the map
bbox <- st_bbox(c(xmin = -122, xmax = -118, ymin = 34.5, ymax = 37), crs = st_crs(NAmap))
filtered_map <- st_crop(NAmap, bbox)

## plot it
map <- ggplot() +
  geom_sf(data = NAmap, fill = "white", color = "black") +
  geom_sf(data = centroids, color = "red", size = 1) +
  ggtitle("Kit Fox Sampling Units") +
  theme_minimal() +
  coord_sf(xlim = c(-122, -118), ylim = c(34.5, 37), expand = FALSE)
```
