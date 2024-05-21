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

> summary(fit_nnls)
Conductance surface with 69553 vertices (11 focal) estimated by maximum likelihood
Call:   radish(formula = chord_dist_matrix ~ kfsuit, data = surface, 
    conductance_model = radish::loglinear_conductance, measurement_model = radish::leastsquares)

Loglikelihood: 154.8743 (4 degrees freedom)
AIC: -301.7485 

Number of function calls: 32 
Number of Newton-Raphson steps: 7 
Norm of gradient at MLE: 2.415845e-13 

Nuisance parameters:
 alpha    beta     tau  
0.1329  0.1571  6.6318  

Coefficients:
       Estimate Std. Error z value Pr(>|z|)   
kfsuit  -0.6854     0.2205  -3.108  0.00188 **
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# refit with with a different measurement model that models
# dependence among pairwise measurements (radish::mlpe)
fit_mlpe <- radish(chord_dist_matrix ~ kfsuit, surface, 
                   radish::loglinear_conductance, radish::mlpe)
summary(fit_mlpe)

> summary(fit_mlpe)
Conductance surface with 69553 vertices (11 focal) estimated by maximum likelihood
Call:   radish(formula = chord_dist_matrix ~ kfsuit, data = surface, 
    conductance_model = radish::loglinear_conductance, measurement_model = radish::mlpe)

Loglikelihood: 200.6914 (5 degrees freedom)
AIC: -391.3828 

Number of function calls: 13 
Number of Newton-Raphson steps: 4 
Norm of gradient at MLE: 4.950834e-07 

Nuisance parameters:
 alpha    beta     tau     rho  
0.1580  0.3843  5.9903  3.2324  

Coefficients:
       Estimate Std. Error z value Pr(>|z|)    
kfsuit   1.0167     0.3046   3.338 0.000844 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# visualisation:
png("KF Optimized Resistance Distance.png", width = 800, height = 600)
plot(fitted(fit_mlpe, "distance"), chord_dist_matrix, pch = 19,
     xlab = "Optimized resistance distance", ylab = "chord distance")
dev.off()

# visualise estimated conductance surface and asymptotic confidence intervals
fitted_conductance <- conductance(surface, fit_mlpe, quantile = 0.95)

png("KF FittedConductance.png", width = 800, height = 600)
plot(fitted_conductance[["est"]], 
     main = "Fitted conductance surface\n(kfsuit)")
dev.off()

png("KF FittedConductance Lower95.png", width = 800, height = 600)
plot(fitted_conductance[["lower95"]], 
     main = "Fitted conductance surface\n(lower 95% CI)")
dev.off()

png("KF FittedConductance Upper95.png", width = 800, height = 600)
plot(fitted_conductance[["upper95"]], main = 
     "Fitted conductance surface\n(upper 95% CI)")
dev.off()

# visualise likelihood surface across grid (takes awhile)
theta <- as.matrix(expand.grid(kfsuit=seq(-1,1,length.out=21)))
grid <- radish_grid(theta, chord_dist_matrix ~ kfsuit, surface,
                    radish::loglinear_conductance, radish::mlpe)

library(ggplot2)
png("KF LikelihoodSurface.png", width = 800, height = 600)
ggplot(data.frame(loglik=grid$loglik, grid$theta), 
       aes(x=kfsuit, y=kfsuit)) + 
  geom_tile(aes(fill=loglik)) + 
  geom_contour(aes(z=loglik), color="black") +
  annotate(geom = "point", colour = "red",
           x = coef(fit_mlpe)["kfsuit"], 
           y = coef(fit_mlpe)["kfsuit"]) +
  theme_bw() +
  xlab(expression(theta[kfsuit])) +
  ylab(expression(theta[kfsuit]))
dev.off()

# calculate resistance distances across grid
distances <- radish_distance(theta, ~ kfsuit, 
                             surface, radish::loglinear_conductance)

#ibd <- which(theta[,1] == 0 & theta[,2] == 0) ## this would be if I had two layers
ibd <- which(theta[,1] == 0)

png("KF Null Resistance Distance.png", width = 800, height = 600)
plot(distances$distance[,,ibd], chord_dist_matrix, pch = 19, 
     xlab = "Null resistance distance (IBD)", ylab = "Chord Distance")
dev.off()

#### Again I'll do this when I have multiple layers
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

#####
# test against null model of IBD
fit_mlpe_ibd <- radish(chord_dist_matrix ~ 1, surface, 
                       radish::loglinear_conductance, radish::mlpe)

> summary(fit_mlpe_ibd)
Conductance surface with 69553 vertices (11 focal) estimated by maximum likelihood
Call:   radish(formula = chord_dist_matrix ~ 1, data = surface, conductance_model = radish::loglinear_conductance, 
    measurement_model = radish::mlpe)

Loglikelihood: 197.4194 (4 degrees freedom)
AIC: -386.8387 

Number of function calls: 1 
Number of Newton-Raphson steps: 0 
Norm of gradient at MLE: NA 

Nuisance parameters:
 alpha    beta     tau     rho  
0.1288  0.2663  6.3520  2.5914  

No coefficients

> anova(fit_mlpe, fit_mlpe_ibd)
Likelihood ratio test
Null: ~ 1
Alt: ~ kfsuit
     logLik Df  ChiSq Df(ChiSq) Pr(>Chi)  
Null 197.42  4                            
Alt  200.69  5 6.5441         1  0.01052 *
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
