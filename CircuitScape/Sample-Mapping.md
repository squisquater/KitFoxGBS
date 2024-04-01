# Sample Mapping

library(dplyr)
library(sf)
library(rnaturalearth)
library(ggplot2)
library(readr)
library(tidyr)

##set working directory
setwd("/Users/sophiepq/Desktop/")

##input samples with location data
data <- read.csv("LynxGeneticDatabase-20231214.csv", header = TRUE)

##remove poor DNA samples
data <- data %>%
  filter(DNA.Individual.ID != "poor DNA")

##make sure lat/long data is numeric
data$Longitude <- as.numeric(data$Longitude)
data$Latitude <- as.numeric(data$Latitude)

# Convert data to an sf object
data_sf <- st_as_sf(data, coords = c("Longitude", "Latitude"), crs = 4326)

# Calculate centroids for replicate samples
centroids <- data_sf %>%
  group_by(DNA.Individual.ID) %>%
  summarize(geometry = st_centroid(st_combine(geometry))) %>%
  mutate(longitude = st_coordinates(geometry)[, 1],
         latitude = st_coordinates(geometry)[, 2])
