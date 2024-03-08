Running all this directly on farm for computational power

Create and activate a conda environment
```
micromamba create --name LandGen

micromamba activate LandGen
```

Install devtools using conda
```
micromamba install conda-forge::r-devtools
micromamba install conda-forge::r-gdistance

```
Activate R
```
R
```
Install & load required packages
```
devtools::install_github("hhwagner1/LandGenCourse") # note that I opted not to update any existing packages when prompted to do so.

library(LandGenCourse)
library(terra)
library(gdistance)
```

Load and plot the raster file
```
myRaster <- rast("/group/ctbrowngrp2/sophiepq/KitFoxGBS/Omniscape/KitFox-ESARPmodel-Raster30x30-NEWmod01.tif")
png("myRasterPlot.png", width = 800, height = 600)
terra::plot(myRaster)
dev.off()
```

let's make some mock data points that will be "centroids" of the subpops. I'll replace these later.
```
df <- data.frame(lon = c(-119.063, -119.462, -119.607, -119.579, -120.043, -120.301, -120.263, -120.749, -120.878), lat = c(35.366, 35.139, 35.375, 35.691, 35.372, 35.850, 36.187, 36.569, 36.644))
```
Convert to spatial points
```
points <- vect(df, geom=c("lon", "lat"), crs="+proj=longlat +datum=NAD83")
```

Project them to match crs of the raster
```
points <- project(points, crs(myRaster))
```

plot it
```
png("myRasterPlot_data.png", width = 800, height = 600)
terra::plot(myRaster[[1]])
terra::points(points, pch=21, col="black", bg="white", cex=2)
dev.off()
```
Create a cost surface
```
kf.cost <- (100-myRaster[[1]])/100

png("KF.CostSurface.png", width = 800, height = 600)
terra::plot(kf.cost[[1]])
terra::points(points, pch=21, col="black", bg="white", cex=2)
dev.off()
```
Convert KF conductance into effective distance
```
kf.cost1 <- gdistance::transition(raster::raster(kf.cost), transitionFunction=mean, directions=8) 
```

Correct for geometric distortion
```
kf.cost1 <- gdistance::geoCorrection(kf.cost1,type = "c",multpl=FALSE)

png("KF.EffectiveDistance.png", width = 800, height = 600)
raster::plot(raster::raster(kf.cost1))
dev.off()
```
