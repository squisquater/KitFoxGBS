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
library(sf)
library(spdep)

Load and plot the raster files
```
#myRaster <- rast("/group/ctbrowngrp2/sophiepq/KitFoxGBS/Omniscape/KitFox-ESARPmodel-Raster30x30-NEWmod01.tif")
#KFRaster <- rast("/group/ctbrowngrp2/sophiepq/KitFoxGBS/LandGenR/KitFox-ESARPmodel-Raster500x500.tif")
kfsuit <- rast("/group/ctbrowngrp2/sophiepq/KitFoxGBS/LandGenR/ESRP-kfsuit-continuous-modified-reprojected-1000x1000.tif")
png("kfsuit.png", width = 800, height = 600)
terra::plot(kfsuit)
dev.off()

roads <- rast("/group/ctbrowngrp2/sophiepq/KitFoxGBS/LandGenR/InterstateHwy5_rasterized_rescaled_road10.tif")
png("roads.png", width = 800, height = 600)
terra::plot(roads)
dev.off()
```

merge the two raster layers together
```
merged_raster <- merge(kfsuit, roads)
png("merged_raster.png", width = 800, height = 600)
terra::plot(merged_raster)
dev.off()
```

let's make some mock data points that will be "centroids" of the subpops. I'll replace these later.
```
df <- data.frame(lon = c(-119.063, -118.769, -120.301, -119.836, -120.263, -119.607, -120.878, -120.749, -119.579, -119.462, -120.043), lat = c(35.366, 35.354, 35.85, 35.174, 36.187, 35.375, 36.644, 36.569, 35.691, 35.139, 35.372))
```
Convert to spatial points
```
points <- vect(df, geom=c("lon", "lat"), crs="+proj=longlat +datum=NAD83")
```
Project them to match crs of the raster
```
points <- project(points, crs(merged_raster))
```

plot it
```
png("myRasterPlot_20240703.png", width = 800, height = 600)
terra::plot(merged_raster[[1]])
terra::points(points, pch=21, col="black", bg="white", cex=2)
dev.off()
```
Create a cost surface
```
kfsuit.cost <- (100-kfsuit[[1]])*0.5
roads.cost <- roads[[1]]*5

png("kfsuit.CostSurface.png", width = 800, height = 600)
terra::plot(kfsuit.cost[[1]])
terra::points(points, pch=21, col="black", bg="white", cex=2)
dev.off()

png("roads.CostSurface.png", width = 800, height = 600)
terra::plot(roads.cost[[1]])
terra::points(points, pch=21, col="black", bg="white", cex=2)
dev.off()
```
Create a single landscape conductance raster and plot it
```
cost1 <- (kfsuit.cost + roads.cost)

png("KF.CostSurface-20240704.png", width = 800, height = 600)
terra::plot(cost1[[1]])
terra::points(points, pch=21, col="black", bg="white", cex=2)
dev.off()
```

tr.cost1 <- gdistance::transition(raster::raster(cost1), transitionFunction=mean, directions=8) 
tr.cost1

png("KF.TransissionLayer-20240704.png", width = 800, height = 600)
par(mar=c(2,2,1,1))
raster::plot(raster::raster(tr.cost1))
dev.off()

Correct for geometric distortion
tr.cost1 <- gdistance::geoCorrection(tr.cost1,type = "c",multpl=FALSE)

Plot shortest paths in space

Convert SpatVector to sf object
```
points_sf <- st_as_sf(points)
```
Convert sf object to SpatialPointsDataFrame
```
sites.sp <- sf::as_Spatial(points_sf)
```

Plot the shortest path between two points
```
png("kf-cost-shortestpath-20240703.png", width = 800, height = 600)
par(mar=c(2,2,1,2))
AtoB <- gdistance::shortestPath(tr.cost1, origin=sites.sp[1,], 
                                goal=sites.sp[2,], output="SpatialLines")
raster::plot(raster::raster(tr.cost1), xlab="x coordinate (m)", 
             ylab="y coordinate (m)",legend.lab="Conductance")
lines(AtoB, col="red", lwd=2)
points(sites.sp[1:2,])
dev.off()
```
What about for all of them!
```
png("kf-cost-shortestpathAll-20240703.png", width = 800, height = 600)
par(mar=c(2,2,1,2))
raster::plot(raster::raster(tr.cost1), xlab="x coordinate (m)", 
             ylab="y coordinate (m)", legend.lab="Conductance")
points(sites.sp)

Neighbours <- spdep::tri2nb(sites.sp@coords, row.names = sites.sp$SiteName)

plot(Neighbours, sites.sp@coords, col="darkgrey", add=TRUE)
for(i in 1:length(Neighbours))
{
  for(j in Neighbours[[i]][Neighbours[[i]] > i])
  {
    AtoB <- gdistance::shortestPath(tr.cost1, origin=sites.sp[i,], 
                                goal=sites.sp[j,], output="SpatialLines")
    lines(AtoB, col="red", lwd=1.5)
  }
}
dev.off()
```





###### OLD

#Convert SpatRaster to RasterLayer
cost1_raster <- raster::raster(cost1)

Convert KF conductance into effective distance
```
tr.cost1 <- gdistance::transition(cost1_raster, transitionFunction=mean, directions=8) 
```
Visually inspect the raster

tr.cost1_raster <- raster::raster(tr.cost1)

raster::plot(raster::raster(tr.cost1))
```
png("KF.EffectiveDistance-20240704.png", width = 800, height = 600)
par(mar=c(2,2,1,1))
raster::plot(raster::raster(tr.cost1))
dev.off()
```
Correct for geometric distortion
```
tr.cost1 <- gdistance::geoCorrection(tr.cost1,type = "c",multpl=FALSE)
```
png("KF.EffectiveDistance-20240703.png", width = 800, height = 600)
raster::plot(raster::raster(tr.cost1))
dev.off()

Plot shortest paths in space

Convert SpatVector to sf object
```
points_sf <- st_as_sf(points)
```
Convert sf object to SpatialPointsDataFrame
```
sites.sp <- sf::as_Spatial(points_sf)
```

Plot the shortest path between two points
```
png("kf-cost-shortestpath-20240703.png", width = 800, height = 600)
par(mar=c(2,2,1,2))
AtoB <- gdistance::shortestPath(tr.cost1, origin=sites.sp[1,], 
                                goal=sites.sp[2,], output="SpatialLines")
raster::plot(raster::raster(tr.cost1), xlab="x coordinate (m)", 
             ylab="y coordinate (m)",legend.lab="Conductance")
lines(AtoB, col="red", lwd=2)
points(sites.sp[1:2,])
dev.off()
```
What about for all of them!
```
png("kf-cost-shortestpathAll-20240703.png", width = 800, height = 600)
par(mar=c(2,2,1,2))
raster::plot(raster::raster(tr.cost1), xlab="x coordinate (m)", 
             ylab="y coordinate (m)", legend.lab="Conductance")
points(sites.sp)

Neighbours <- spdep::tri2nb(sites.sp@coords, row.names = sites.sp$SiteName)

plot(Neighbours, sites.sp@coords, col="darkgrey", add=TRUE)
for(i in 1:length(Neighbours))
{
  for(j in Neighbours[[i]][Neighbours[[i]] > i])
  {
    AtoB <- gdistance::shortestPath(tr.cost1, origin=sites.sp[i,], 
                                goal=sites.sp[j,], output="SpatialLines")
    lines(AtoB, col="red", lwd=1.5)
  }
}
dev.off()
```


