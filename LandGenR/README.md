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
```
Load and plot the raster files
```
#myRaster <- rast("/group/ctbrowngrp2/sophiepq/KitFoxGBS/Omniscape/KitFox-ESARPmodel-Raster30x30-NEWmod01.tif")
#kfsuit <- rast("/group/ctbrowngrp2/sophiepq/KitFoxGBS/LandGenR/Radish/KitFox-ESARPmodel-Raster1000x1000-FillAllCells.tif")
#kfsuit <- rast("/group/ctbrowngrp2/sophiepq/KitFoxGBS/LandGenR/ESRP-kfsuit-continuous-modified-reprojected-1000x1000.tif")
kfsuit <- rast("/group/ctbrowngrp2/sophiepq/KitFoxGBS/LandGenR/ESRP-kfsuit-continuous-modified-gdal-reprojected-1000x1000.tif")
png("kfsuit-20240708.png", width = 800, height = 600)
terra::plot(kfsuit)
dev.off()

#roads <- rast("/group/ctbrowngrp2/sophiepq/KitFoxGBS/LandGenR/InterstateHwy5_rasterized_rescaled_road10.tif")
roads <- rast("/group/ctbrowngrp2/sophiepq/KitFoxGBS/LandGenR/InterstateHwy5_rasterized.tif")
png("roads-20240708.png", width = 800, height = 600)
terra::plot(roads)
dev.off()
```

merge the two raster layers together
```
merged_raster <- merge(kfsuit, roads)
png("merged_raster-20240708.png", width = 800, height = 600)
terra::plot(merged_raster)
dev.off()
```

let's make some mock data points that will be "centroids" of the subpops. I'll replace these later.
```
df <- data.frame(lon = c(-119.063, -118.769, -120.301, -119.836, -120.263, -119.607, -120.878, -120.749, -119.579, -119.462, -120.043), lat = c(35.366, 35.354, 35.85, 35.174, 36.187, 35.375, 36.644, 36.569, 35.691, 35.139, 35.372), SiteName = c("Bakersfield", "Bena", "CalFlats", "Carrizo", "Coalinga", "Lokern", "Panoche-North", "Panoche-South", "Semitropic", "Taft", "Topaz"))
```
Convert to spatial points
```
points <- vect(df, geom=c("lon", "lat"), crs="+proj=longlat +datum=NAD83")
```
Project them to match crs of the raster
```
points <- project(points, crs(merged_raster))
```
Convert to SpatialPoints
```
coords <- crds(points)
sites.sp <- SpatialPoints(coords, proj4string=CRS("+proj=aea +lat_0=0 +lon_0=-120 +lat_1=34 +lat_2=40.5 +x_0=0 +y_0=-4000000 +datum=NAD83 +units=m +no_defs"))
```
Convert to SpatialPointsDataFrame
```
sites.spdf <- SpatialPointsDataFrame(sites.sp, data=df)
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
roads.cost <- roads[[1]]*50

png("kfsuit.CostSurface-20240708.png", width = 800, height = 600)
terra::plot(kfsuit.cost[[1]])
terra::points(points, pch=21, col="black", bg="white", cex=2)
dev.off()

png("roads.CostSurface-20240708.png", width = 800, height = 600)
terra::plot(roads.cost[[1]])
terra::points(points, pch=21, col="black", bg="white", cex=2)
dev.off()
```
Create a conductance surface
```
kfsuit.conductance <- (kfsuit[[1]])*0.5
roads.conductance <- (1-roads[[1]])*50

png("kfsuit.ConductanceSurface-20240708.png", width = 800, height = 600)
terra::plot(kfsuit.conductance[[1]])
terra::points(points, pch=21, col="black", bg="white", cex=2)
dev.off()

png("roads.ConductanceSurface-20240708.png", width = 800, height = 600)
terra::plot(roads.conductance[[1]])
terra::points(points, pch=21, col="black", bg="white", cex=2)
dev.off()
```
Create a single landscape conductance raster and plot it
```
conductance1 <- (kfsuit.conductance + roads.conductance)
conductance1 <- (kfsuit.conductance)

png("KF.ConductanceSurface-20240708.png", width = 800, height = 600)
terra::plot(conductance1[[1]])
terra::points(points, pch=21, col="black", bg="white", cex=2)
dev.off()
```

## Read in Genetic Distance
```
genDist <- read.table("/group/ctbrowngrp2/sophiepq/KitFoxGBS/LandGenR/Radish/KF_NeisD_matrix.txt", sep = "\t", header = TRUE, row.names = 1)
genDist <- as.dist(genDist)
```
## Create a geographic distance matrix
```
coords <- df[, 1:2]
geoDist <- pointDistance(coords, longlat = TRUE)
```
*How well does geographic distance predict genetic distance?*
```
cor(genDist, geoDist)
[1] 0.2498877
```
## Let's see if least-cost distance or resistance distance are better predictors of genetic distance

### Create a transition layer
```
tr.cost1 <- gdistance::transition(raster::raster(conductance1), transitionFunction=mean, directions=8) 
tr.cost1

png("KF.TransissionLayer-20240708.png", width = 800, height = 600)
par(mar=c(2,2,1,1))
raster::plot(raster::raster(tr.cost1))
dev.off()
```
### Correct for geometric distortion
```
trC.cost1 <- gdistance::geoCorrection(tr.cost1,type = "c",multpl=FALSE)
trR.cost1 <- gdistance::geoCorrection(tr.cost1,type = "r",multpl=FALSE)
```
### Plot shortest paths in space

Start by plotting the shortest path between two points
```
png("kf-cost-shortestpath-20240708.png", width = 800, height = 600)
par(mar=c(2,2,1,2))
AtoB <- gdistance::shortestPath(trC.cost1, origin=sites.spdf[1,], 
                                goal=sites.spdf[2,], output="SpatialLines")
raster::plot(raster::raster(trC.cost1), xlab="x coordinate (m)", 
             ylab="y coordinate (m)",legend.lab="Conductance")
lines(AtoB, col="red", lwd=2)
points(sites.spdf[1:2,])
dev.off()
```
What about for all of them!
```
png("kf-cost-shortestpathAll-20240708.png", width = 800, height = 600)
par(mar=c(2,2,1,2))
raster::plot(raster::raster(trC.cost1), xlab="x coordinate (m)", 
             ylab="y coordinate (m)", legend.lab="Conductance")
points(sites.spdf)

Neighbours <- spdep::tri2nb(sites.spdf@coords, row.names = sites.spdf$SiteName)

plot(Neighbours, sites.spdf@coords, col="darkgrey", add=TRUE)
for(i in 1:length(Neighbours))
{
  for(j in Neighbours[[i]][Neighbours[[i]] > i])
  {
    AtoB <- gdistance::shortestPath(trC.cost1, origin=sites.spdf[i,], 
                                goal=sites.spdf[j,], output="SpatialLines")
    lines(AtoB, col="red", lwd=1.5)
  }
}
dev.off()
```

### Least Cost Distance
```
cost1.dist <- gdistance::costDistance(trC.cost1,sites.sp)
```
### Cost-distance matrix based on random paths (similar to Circuitscape)
I might need to use a different geocorrection (R) for this one -- that wasn't the issue. I think there's something with NA values in my raster files.... I'll explore more. For now I'm jsut working with an older kfsuit layer until I problem solve the existing layers.
```
comm1.dist <- gdistance::commuteDistance(x = trR.cost1, coords = sites.sp)
```
Compare cost distances
```
dist_df <- data.frame("cost1.dist"=as.numeric(cost1.dist),
                      "comm1.dist"=as.numeric(comm1.dist))
```
Look at the corelation between the two
```
corr.LCD.comm <- cor(dist_df$cost1.dist, dist_df$comm1.dist, method = "spearman")
corr.LCD.comm
[1] 0.9154401

# Convert the matrices to vectors for plotting

cost1.dist_vector <- as.numeric(cost1.dist)
comm1.dist_vector <- as.numeric(comm1.dist)

# Create a data frame for comparison
dist_df <- data.frame("cost1.dist" = cost1.dist_vector,
                      "comm1.dist" = comm1.dist_vector)

png("CorrPlot.png", width = 800, height = 600)
plot(dist_df$comm1.dist, dist_df$cost1.dist, xlab = "Commute Distance", ylab = "Cost Distance", main = "Correlation Plot")
dev.off()
```

### How well do these predict genetic distance?
cor(cost1.dist, geoDist)
[1] 0.9990757

cor(comm1.dist, geoDist)
[1] 0.9726797

## Mantel Tests

