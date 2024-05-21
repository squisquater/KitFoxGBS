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

Load and plot the raster files
```
#myRaster <- rast("/group/ctbrowngrp2/sophiepq/KitFoxGBS/Omniscape/KitFox-ESARPmodel-Raster30x30-NEWmod01.tif")
KFRaster <- rast("/group/ctbrowngrp2/sophiepq/KitFoxGBS/LandGenR/KitFox-ESARPmodel-Raster500x500.tif")
png("myRasterPlot500.png", width = 800, height = 600)
terra::plot(KFRaster)
dev.off()

roads <- rast("/group/ctbrowngrp2/sophiepq/KitFoxGBS/LandGenR/CaliforniaRoads-Reprojected-500x500.tif")
png("Roads_raster500.png", width = 800, height = 600)
terra::plot(roads)
dev.off()

KFRaster <- rast("/group/ctbrowngrp2/sophiepq/KitFoxGBS/LandGenR/KitFox-ESARPmodel-Raster1000x1000.tif")
png("myRasterPlot1000.png", width = 800, height = 600)
terra::plot(KFRaster)
dev.off()

roads <- rast("/group/ctbrowngrp2/sophiepq/KitFoxGBS/LandGenR/CaliforniaRoads-Reprojected-1000x1000.tif")
png("Roads_raster1000.png", width = 800, height = 600)
terra::plot(roads)
dev.off()
```

merge the two raster layers together
```
merged_raster <- merge(KFRaster, roads)
png("merged_raster1000.png", width = 800, height = 600)
terra::plot(merged_raster)
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
points <- project(points, crs(merged_raster))
```

plot it
```
png("myRasterPlot_data1000.png", width = 800, height = 600)
terra::plot(merged_raster[[1]])
terra::points(points, pch=21, col="black", bg="white", cex=2)
dev.off()
```
Create a cost surface
```
kf.cost <- (100-KFRaster[[1]])/100

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



GENERATE A CHORD (DISTANCE) MATRIX

I can't find a good tool for data conversion so I'm tweaking the format the plink file myself.
Note: For some weird reason I have logical values (TRUE only) in my plink file. I converted all of these to "T".
```
# Load the necessary data if not already loaded
plink.ped <- read.table("/group/ctbrowngrp2/sophiepq/KitFoxGBS/stacks/SJKF_metapop_LandGen_n223_20x/populations.plink_noheader.ped", header=FALSE, stringsAsFactors = FALSE)
plink.map <- read.table("/group/ctbrowngrp2/sophiepq/KitFoxGBS/stacks/SJKF_metapop_LandGen_n223_20x/populations.plink_noheader.map", header=FALSE, stringsAsFactors = FALSE)

# Extract FamilyID, IndividualID, and SNP columns
individual_info <- plink.ped[, 1:2]
num_snps <- (ncol(plink.ped) - 6) / 2
snp_data <- plink.ped[, 7:ncol(plink.ped)]

# Convert logical columns to character to ensure no logical values
snp_data <- apply(snp_data, 2, function(col) {
  if (is.logical(col)) {
    col <- as.character(col)
  }
  return(col)
})

# Replace "TRUE" with "T" in the snp_data
snp_data[snp_data == "TRUE"] <- "T"

# Check for any logical values in the cleaned snp_data
logical_values_after <- unique(as.vector(snp_data[sapply(snp_data, is.logical)]))
if (length(logical_values_after) > 0) {
  print("Logical values found in snp_data after cleaning:")
  print(logical_values_after)
} else {
  print("No logical values found in snp_data after cleaning.")
}

# Check unique values in snp_data before merging to identify problematic entries
unique_snp_values_before <- unique(as.vector(snp_data))
print("Unique values in snp_data before merging:")
print(unique_snp_values_before)

# Function to merge SNP alleles into a single genotype and convert "00" to NA
merge_snps <- function(df) {
  merged_snps <- apply(df, 1, function(row) {
    snps <- as.character(row[seq(1, length(row), by=2)])
    alleles <- as.character(row[seq(2, length(row), by=2)])
    genotypes <- paste(snps, alleles, sep="")
    genotypes[genotypes == "00"] <- NA
    return(genotypes)
  })
  return(t(merged_snps))
}

# Merge the SNP columns
merged_snps <- merge_snps(snp_data)

# Check unique values in merged_snps to verify format
unique_values <- unique(as.vector(merged_snps))
print("Unique values in merged_snps:")
print(unique_values)

# Assign column names to merged SNPs using SNP IDs from .map file
snp_ids <- plink.map$V2
colnames(merged_snps) <- snp_ids

# Combine the individual info with merged SNP data
final_data <- cbind(individual_info, merged_snps)

# Print a subset of final_data to check the combined data format
print(head(final_data, 10))  # Print a larger sample of rows

# Import the data frame into a genind object
data.genind <- adegenet::df2genind(X=final_data[,-c(1:2)], sep="", ncode=1,   
                          ind.names=final_data$V2, loc.names=NULL, 
                          pop=final_data$V1, NA.char="NA", ploidy=2, 
                          type="codom", strata=NULL, hierarchy=NULL)
data.genind

genpop_obj <- genind2genpop(data.genind)

chord_dist <- dist.genpop(genpop_obj, method = 2, diag = TRUE, upper = TRUE)



