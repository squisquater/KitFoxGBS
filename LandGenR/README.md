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
kfsuit <- rast("/group/ctbrowngrp2/sophiepq/KitFoxGBS/LandGenR/ESRP-kfsuit-continuous-modified-gdal-5-reprojected-1000x1000-nodata.tif")
png("kfsuit-20240723.png", width = 800, height = 600)
terra::plot(kfsuit)
dev.off()

roads <- rast("/group/ctbrowngrp2/sophiepq/KitFoxGBS/LandGenR/InterstateHwy5_rasterized.tif")
png("roads-20240723.png", width = 800, height = 600)
terra::plot(roads)
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

Create a cost surface
```
kfsuit.cost <- (100-kfsuit[[1]])*0.5
roads.cost <- roads[[1]]*50

png("kfsuit.CostSurface-20240722.png", width = 800, height = 600)
terra::plot(kfsuit.cost[[1]])
terra::points(points, pch=21, col="black", bg="white", cex=2)
dev.off()

png("roads.CostSurface-20240722.png", width = 800, height = 600)
terra::plot(roads.cost[[1]])
terra::points(points, pch=21, col="black", bg="white", cex=2)
dev.off()
```
Create a conductance surface
```
kfsuit.conductance <- (kfsuit[[1]])*0.5
roads.conductance <- (1-roads[[1]])*50

png("kfsuit.ConductanceSurface-20240722.png", width = 800, height = 600)
terra::plot(kfsuit.conductance[[1]])
terra::points(points, pch=21, col="black", bg="white", cex=2)
dev.off()

png("roads.ConductanceSurface-20240722.png", width = 800, height = 600)
terra::plot(roads.conductance[[1]])
terra::points(points, pch=21, col="black", bg="white", cex=2)
dev.off()
```
Create a single landscape conductance raster and plot it
```
conductance1 <- (kfsuit.conductance + roads.conductance)
#conductance1 <- (kfsuit.conductance)

png("KF.ConductanceSurface-FullModel-20240722.png", width = 800, height = 600)
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
geoDist <- as.dist(geoDist)
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

png("KF.TransissionLayer-20240722.png", width = 800, height = 600)
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
```
png("kf-cost-shortestpathAll-20240722.png", width = 800, height = 600)
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
[1] 0.9627706 - full model (50:50)
[1] 0.9332612 - kfsuit model
[1] 0.9930014 - roads model

# Convert the matrices to vectors for plotting

cost1.dist_vector <- as.numeric(cost1.dist)
comm1.dist_vector <- as.numeric(comm1.dist)

# Create a data frame for comparison
dist_df <- data.frame("cost1.dist" = cost1.dist_vector,
                      "comm1.dist" = comm1.dist_vector)

png("CorrPlot-CostComm-20240722.png", width = 800, height = 600)
plot(dist_df$comm1.dist, dist_df$cost1.dist, xlab = "Commute Distance", ylab = "Cost Distance", main = "Correlation Plot")
dev.off()
```

### How well do these predict genetic distance?
To reiterate geographic vs genetic
cor(genDist, geoDist)
[1] 0.2498877

cor(genDist,cost1.dist)
[1] 0.2546895 - full model
[1] 0.2533487 - kfsuit model
[1] 0.2598093 - roads model

cor(genDist,comm1.dist)
[1] 0.2370583 - full model
[1] 0.1791814 - kfsuit model
[1] 0.3315328 - roads model *

## Perform Mantel tests
install.packages("vegan")
library(vegan)

### Mantel test for geographic distance vs genetic distance
```
mantel_geo_gen <- mantel(geoDist, genDist, method="pearson", permutations=9999)
print(mantel_geo_gen)
> mantel_geo_gen

Mantel statistic based on Pearson's product-moment correlation 

Call:
mantel(xdis = geoDist, ydis = genDist, method = "pearson", permutations = 9999) 

Mantel statistic r: 0.2499 
      Significance: 0.0869 

Upper quantiles of permutations (null model):
  90%   95% 97.5%   99% 
0.235 0.297 0.351 0.397 
Permutation: free
Number of permutations: 9999
```
### Mantel test for least-cost distance vs genetic distance
```
mantel_cost_gen <- mantel(as.dist(cost1.dist), genDist, method="pearson", permutations=9999)
print(mantel_cost_gen)
> mantel_cost_gen (full model)
Mantel statistic based on Pearson's product-moment correlation 

Call:
mantel(xdis = as.dist(cost1.dist), ydis = genDist, method = "pearson",      permutations = 9999) 

Mantel statistic r: 0.2547 
      Significance: 0.0942 

Upper quantiles of permutations (null model):
  90%   95% 97.5%   99% 
0.249 0.311 0.361 0.414 
Permutation: free
Number of permutations: 9999

> mantel_cost_gen (kfsuit model)

Mantel statistic based on Pearson's product-moment correlation 

Call:
mantel(xdis = as.dist(cost1.dist), ydis = genDist, method = "pearson",      permutations = 9999) 

Mantel statistic r: 0.2533 
      Significance: 0.1005 

Upper quantiles of permutations (null model):
  90%   95% 97.5%   99% 
0.254 0.320 0.371 0.419 
Permutation: free
Number of permutations: 9999

> mantel_cost_gen (roads model)

Mantel statistic based on Pearson's product-moment correlation 

Call:
mantel(xdis = as.dist(cost1.dist), ydis = genDist, method = "pearson",      permutations = 9999) 

Mantel statistic r: 0.2598 
      Significance: 0.0794 

Upper quantiles of permutations (null model):
  90%   95% 97.5%   99% 
0.237 0.299 0.345 0.396 
Permutation: free
Number of permutations: 9999
```
### Mantel test for commute distance vs genetic distance
```
mantel_comm_gen <- mantel(as.dist(comm1.dist), genDist, method="pearson", permutations=9999)
print(mantel_comm_gen)
> mantel_comm_gen (full model)

Mantel statistic based on Pearson's product-moment correlation 

Call:
mantel(xdis = as.dist(comm1.dist), ydis = genDist, method = "pearson",      permutations = 9999) 

Mantel statistic r: 0.2371 
      Significance: 0.118 

Upper quantiles of permutations (null model):
  90%   95% 97.5%   99% 
0.255 0.319 0.371 0.419 
Permutation: free
Number of permutations: 9999

> mantel_comm_gen (kfsuit model)

Mantel statistic based on Pearson's product-moment correlation 

Call:
mantel(xdis = as.dist(comm1.dist), ydis = genDist, method = "pearson",      permutations = 9999) 

Mantel statistic r: 0.1792 
      Significance: 0.2127 

Upper quantiles of permutations (null model):
  90%   95% 97.5%   99% 
0.284 0.359 0.412 0.472 
Permutation: free
Number of permutations: 9999

> mantel_comm_gen (roads model)

Mantel statistic based on Pearson's product-moment correlation 

Call:
mantel(xdis = as.dist(comm1.dist), ydis = genDist, method = "pearson",      permutations = 9999) 

Mantel statistic r: 0.3315 
      Significance: 0.0151 

Upper quantiles of permutations (null model):
  90%   95% 97.5%   99% 
0.209 0.263 0.308 0.351 
Permutation: free
Number of permutations: 9999
```

## Compare Models
```
# Perform Mantel tests for least-cost distance and commute distance vs genetic distance
mantel_geo_gen <- mantel(as.dist(geoDist), genDist, method="pearson", permutations=999)
print(mantel_geo_gen)

mantel_cost_gen <- mantel(as.dist(cost1.dist), genDist, method="pearson", permutations=9999)
print(mantel_cost_gen)

mantel_comm_gen <- mantel(as.dist(comm1.dist), genDist, method="pearson", permutations=9999)
print(mantel_comm_gen)

# Statistical comparison of Mantel statistics
set.seed(123)
permutations <- 999

mantel_diff <- function(matrix1, model1, model2, permutations) {
  matrix1 <- as.matrix(matrix1) # Convert distance objects to matrices
  model1 <- as.matrix(model1)
  model2 <- as.matrix(model2)
  diff_stats <- numeric(permutations)
  for (i in 1:permutations) {
    permuted <- sample(nrow(matrix1))
    mantel1 <- mantel(as.dist(matrix1[permuted, permuted]), as.dist(model1), method = "pearson", permutations = 0)$statistic
    mantel2 <- mantel(as.dist(matrix1[permuted, permuted]), as.dist(model2), method = "pearson", permutations = 0)$statistic
    diff_stats[i] <- mantel1 - mantel2
  }
  return(diff_stats)
}

### Compute the difference in Mantel statistics with bootstrapping
diff_stats_geo_cost <- mantel_diff(genDist, as.dist(geoDist), as.dist(cost1.dist), permutations)
diff_stats_geo_comm <- mantel_diff(genDist, as.dist(geoDist), as.dist(comm1.dist), permutations)
diff_stats_cost_comm <- mantel_diff(genDist, as.dist(cost1.dist), as.dist(comm1.dist), permutations)

### Observed difference in Mantel statistics
observed_diff_geo_cost <- mantel_geo_gen$statistic - mantel_cost_gen$statistic
observed_diff_geo_comm <- mantel_geo_gen$statistic - mantel_comm_gen$statistic
observed_diff_cost_comm <- mantel_cost_gen$statistic - mantel_comm_gen$statistic

### Calculate p-value
p_value_geo_cost <- sum(abs(diff_stats_geo_cost) >= abs(observed_diff_geo_cost)) / permutations
p_value_geo_comm <- sum(abs(diff_stats_geo_comm) >= abs(observed_diff_geo_comm)) / permutations
p_value_cost_comm <- sum(abs(diff_stats_cost_comm) >= abs(observed_diff_cost_comm)) / permutations

### Print the observed difference and p-value
cat("P-value (geo vs cost):", p_value_geo_cost, "\n")
# P-value (geo vs cost): 0.8218218 
cat("Observed difference in Mantel statistics (geo vs commute):", observed_diff_geo_comm, "\n")
# Observed difference in Mantel statistics (geo vs commute): 0.01282938 
cat("P-value (geo vs commute):", p_value_geo_comm, "\n")
# P-value (geo vs commute): 0.9049049 
cat("Observed difference in Mantel statistics (cost vs commute):", observed_diff_cost_comm, "\n")
# Observed difference in Mantel statistics (cost vs commute): 0.01763121 
cat("P-value (cost vs commute):", p_value_cost_comm, "\n")
# P-value (cost vs commute): 0.8598599 
```
**Summary:**  For the Full Model (kfsuit + roads) neither cost distance nor commute distance predicts genetic distance better than the null model. 


# What if I remove baskersfield and bena???
```
df <- data.frame(lon = c(-120.301, -119.836, -120.263, -119.607, -120.878, -120.749, -119.579, -119.462, -120.043), lat = c(35.85, 35.174, 36.187, 35.375, 36.644, 36.569, 35.691, 35.139, 35.372), SiteName = c("CalFlats", "Carrizo", "Coalinga", "Lokern", "Panoche-North", "Panoche-South", "Semitropic", "Taft", "Topaz"))
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

## Read in Genetic Distance (without Bakersfield and Bena)
```
genDist <- read.table("/group/ctbrowngrp2/sophiepq/KitFoxGBS/LandGenR/Radish/KF_NeisD_matrix.txt", sep = "\t", header = TRUE, row.names = 1)

genDist <- genDist[-c(1, 2), -c(1, 2)]

genDist <- as.dist(genDist)
```
## Create a geographic distance matrix
```
coords <- df[, 1:2]
geoDist <- pointDistance(coords, longlat = TRUE)
geoDist <- as.dist(geoDist)
```
*How well does geographic distance predict genetic distance?*
```
cor(genDist, geoDist)
[1] 0.2498877 (all sites)
[1] 0.2433318 (no bakersfield/bena)
```

```
kfsuit.conductance <- (kfsuit[[1]])*0.5
roads.conductance <- (1-roads[[1]])*50
conductance1 <- (kfsuit.conductance + roads.conductance)
tr.cost1 <- gdistance::transition(raster::raster(conductance1), transitionFunction=mean, directions=8) 
trC.cost1 <- gdistance::geoCorrection(tr.cost1,type = "c",multpl=FALSE)
trR.cost1 <- gdistance::geoCorrection(tr.cost1,type = "r",multpl=FALSE)
cost1.dist <- gdistance::costDistance(trC.cost1,sites.sp)
comm1.dist <- gdistance::commuteDistance(x = trR.cost1, coords = sites.sp)
```

### Testing correlations and comparing to model with all geographic locations
To reiterate geographic vs genetic
```
cor(genDist, geoDist)
[1] 0.2498877 (all sites)
[1] 0.2433318 (no bakersfield/bena)

cor(genDist,cost1.dist)
[1] 0.2546895 - full model (all sites)
[1] 0.2584517 - full model (no Bakers/Bena)
[1] 0.2533487 - kfsuit model (all sites)
[1] 0.2838693 - kfsuit model (no Bakers/Bena)
[1] 0.2598093 - roads model (all sites)
[1] 0.2447101 - roads model (no Bakers/Bena)

cor(genDist,comm1.dist)
[1] 0.2370583 - full model (all sites)
[1] 0.3002291 - full model (no Bakers/Bena) *
[1] 0.1791814 - kfsuit model (all sites)
[1] 0.2198528 - kfsuit model (no Bakers/Bena)
[1] 0.3315328 - roads model (all sites) *
[1] 0.3404909 - roads model (no Bakers/Bena) **
```

### Mantel test for geographic distance vs genetic distance
```
mantel_geo_gen <- mantel(geoDist, genDist, method="pearson", permutations=9999)
print(mantel_geo_gen)
> mantel_geo_gen

Mantel statistic based on Pearson's product-moment correlation 

Call:
mantel(xdis = geoDist, ydis = genDist, method = "pearson", permutations = 9999) 

Mantel statistic r: 0.2433 
      Significance: 0.1003 

Upper quantiles of permutations (null model):
  90%   95% 97.5%   99% 
0.243 0.305 0.353 0.406 
Permutation: free
Number of permutations: 9999
```
### Mantel test for least-cost distance vs genetic distance
```
mantel_cost_gen <- mantel(as.dist(cost1.dist), genDist, method="pearson", permutations=9999)
print(mantel_cost_gen)
> mantel_cost_gen (full model)

Mantel statistic based on Pearson's product-moment correlation 

Call:
mantel(xdis = as.dist(cost1.dist), ydis = genDist, method = "pearson",      permutations = 9999) 

Mantel statistic r: 0.2585 
      Significance: 0.0959 

Upper quantiles of permutations (null model):
  90%   95% 97.5%   99% 
0.254 0.320 0.367 0.418 
Permutation: free
Number of permutations: 9999

> mantel_cost_gen (kfsuit model)

Mantel statistic based on Pearson's product-moment correlation 

Call:
mantel(xdis = as.dist(cost1.dist), ydis = genDist, method = "pearson",      permutations = 9999) 

Mantel statistic r: 0.2839 
      Significance: 0.07 

Upper quantiles of permutations (null model):
  90%   95% 97.5%   99% 
0.247 0.313 0.366 0.412 
Permutation: free
Number of permutations: 9999

> mantel_cost_gen (roads model)

Mantel statistic based on Pearson's product-moment correlation 

Call:
mantel(xdis = as.dist(cost1.dist), ydis = genDist, method = "pearson",      permutations = 9999) 

Mantel statistic r: 0.2447 
      Significance: 0.0987 

Upper quantiles of permutations (null model):
  90%   95% 97.5%   99% 
0.243 0.303 0.348 0.402 
Permutation: free
Number of permutations: 9999
```

### Mantel test for commute distance vs genetic distance
```
mantel_comm_gen <- mantel(as.dist(comm1.dist), genDist, method="pearson", permutations=9999)
print(mantel_comm_gen)
> mantel_comm_gen (full model)

Mantel statistic based on Pearson's product-moment correlation 

Call:
mantel(xdis = as.dist(comm1.dist), ydis = genDist, method = "pearson",      permutations = 9999) 

Mantel statistic r: 0.3002 
      Significance: 0.0476 

Upper quantiles of permutations (null model):
  90%   95% 97.5%   99% 
0.239 0.296 0.342 0.382 
Permutation: free
Number of permutations: 9999

> mantel_comm_gen (kfsuit model)

Mantel statistic based on Pearson's product-moment correlation 

Call:
mantel(xdis = as.dist(comm1.dist), ydis = genDist, method = "pearson",      permutations = 9999) 

Mantel statistic r: 0.2199 
      Significance: 0.1695 

Upper quantiles of permutations (null model):
  90%   95% 97.5%   99% 
0.296 0.369 0.421 0.470 
Permutation: free
Number of permutations: 9999

> mantel_comm_gen (roads model)

Mantel statistic based on Pearson's product-moment correlation 

Call:
mantel(xdis = as.dist(comm1.dist), ydis = genDist, method = "pearson",      permutations = 9999) 

Mantel statistic r: 0.3405 
      Significance: 0.005 

Upper quantiles of permutations (null model):
  90%   95% 97.5%   99% 
0.199 0.250 0.286 0.319 
Permutation: free
Number of permutations: 9999
```

Plot the best model to observe fit
```
png("GenDist_vs_CommDist_roads_noBakersBena_20240711.png", width = 800, height = 600)
plot(comm1.dist, genDist, xlab = "Commute Distance", ylab = "Genetic Distance", main = "Genetic Distance vs. Commute Distance")
dev.off()

png("GenDist_vs_CommDist_fullmodel_noBakersBena_20240711.png", width = 800, height = 600)
plot(comm1.dist, genDist, xlab = "Commute Distance", ylab = "Genetic Distance", main = "Genetic Distance vs. Commute Distance")
dev.off()
```

## Compare Models Using Partial Mantel Test
```
result <- mantel.partial(genDist, geoDist, cost1.dist, method = "pearson",  permutations = 999)

```
