```
library(adegenet)
library(radish)
library(raster)

GENERATE A CHORD (DISTANCE) MATRIX

I can't find a good tool for data conversion so I'm tweaking the format the plink file myself.
Note: For some weird reason I have logical values (TRUE only) in my plink file. I converted all of these to "T".

# Load the necessary data if not already loaded - make sure if this is coming from stacks there is no header
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

# convert to matrix
chord_dist_matrix <- as.matrix(chord_dist)

#################################################################################

# Run Radish
myRaster <- raster("/group/ctbrowngrp2/sophiepq/KitFoxGBS/LandGenR/KitFox-ESARPmodel-Raster1000x1000.tif")

# scaling spatial covariates helps avoid numeric overflow
covariates <- raster::stack(list(kfsuit = raster::scale(myRaster)))

#Load in spatial points dataframe
df <- data.frame(lon = c(-119.063, -118.769, -120.301, -119.836, -120.263, -119.607, -120.878, -120.749, -119.579, -119.462, -120.043), lat = c(35.366, 35.354, 35.85, 35.174, 36.187, 35.375, 36.644, 36.569, 35.691, 35.139, 35.372))
#Convert to spatial points dataframe
coordinates(df) <- ~ lon + lat
proj4string(df) <- CRS("+proj=longlat +datum=WGS84 +no_defs")

# Extract the CRS from the existing raster
raster_crs <- crs(covariates)

# Reproject the spatial points to match the raster CRS
df_projected <- spTransform(df, CRSobj = raster_crs)

#Extract to a spatial points object
df_points <- SpatialPoints(coordinates(df_projected), proj4string = CRS(proj4string(df_projected)))

png("KFsuit.png", width = 800, height = 600)
plot(covariates[["kfsuit"]])
points(df_points, pch = 19)
dev.off()

surface <- conductance_surface(covariates, df_points, directions = 8)

fit_nnls <- radish(melip.Fst ~ forestcover + altitude, surface, 
                   radish::loglinear_conductance, radish::leastsquares)
summary(fit_nnls)

# refit with with a different measurement model that models
# dependence among pairwise measurements (radish::mlpe)
fit_nnls <- radish(chord_dist ~ kfsuit, surface, 
                   radish::loglinear_conductance, radish::leastsquares)
summary(fit_mlpe)

# visualisation:
png("Optimized Resistance Distance.png", width = 800, height = 600)
plot(fitted(fit_mlpe, "distance"), melip.Fst, pch = 19,
     xlab = "Optimized resistance distance", ylab = "Fst")
dev.off()

# visualise estimated conductance surface and asymptotic confidence intervals
fitted_conductance <- conductance(surface, fit_mlpe, quantile = 0.95)

png("FittedConductance.png", width = 800, height = 600)
plot(fitted_conductance[["est"]], 
     main = "Fitted conductance surface\n(forestcover + altitude)")
plot(fitted_conductance[["lower95"]], 
     main = "Fitted conductance surface\n(lower 95% CI)")
plot(fitted_conductance[["upper95"]], main = 
     "Fitted conductance surface\n(upper 95% CI)")
dev.off()

# visualise likelihood surface across grid (takes awhile)
theta <- as.matrix(expand.grid(forestcover=seq(-1,1,length.out=21), 
                               altitude=seq(-1,1,length.out=21)))
grid <- radish_grid(theta, melip.Fst ~ forestcover + altitude, surface,
                    radish::loglinear_conductance, radish::mlpe)

library(ggplot2)
png("LikelihoodSurface.png", width = 800, height = 600)
ggplot(data.frame(loglik=grid$loglik, grid$theta), 
       aes(x=forestcover, y=altitude)) + 
  geom_tile(aes(fill=loglik)) + 
  geom_contour(aes(z=loglik), color="black") +
  annotate(geom = "point", colour = "red",
           x = coef(fit_mlpe)["forestcover"], 
           y = coef(fit_mlpe)["altitude"]) +
  theme_bw() +
  xlab(expression(theta[altitude])) +
  ylab(expression(theta[forestcover]))
dev.off()

# calculate resistance distances across grid
distances <- radish_distance(theta, ~forestcover + altitude, 
                             surface, radish::loglinear_conductance)

ibd <- which(theta[,1] == 0 & theta[,2] == 0)

png("Null Resistance Distance.png", width = 800, height = 600)
plot(distances$distance[,,ibd], melip.Fst, pch = 19, 
     xlab = "Null resistance distance (IBD)", ylab = "Fst")
dev.off()
# model selection:
# fit a reduced model without "forestcover" covariate, and compare to 
# full model via a likelihood ratio test
fit_mlpe_reduced <- radish(melip.Fst ~ altitude, surface, 
                           radish::loglinear_conductance, radish::mlpe)
anova(fit_mlpe, fit_mlpe_reduced)

# test for an interaction
fit_mlpe_interaction <- radish(melip.Fst ~ forestcover * altitude, surface, 
                               radish::loglinear_conductance, radish::mlpe)
anova(fit_mlpe, fit_mlpe_interaction)

# test against null model of IBD
fit_mlpe_ibd <- radish(melip.Fst ~ 1, surface, 
                       radish::loglinear_conductance, radish::mlpe)
anova(fit_mlpe, fit_mlpe_ibd)
```
