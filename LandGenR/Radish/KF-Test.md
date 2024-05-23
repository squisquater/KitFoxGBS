Some comments on generating raster layers. I did this in QGIS. made sure they were the same extent and that cells were the same size (1km). I also needed to add in background values (0) for the roads layer and did this using Processing Toolbox > Raster Tools > Fill NoData cells. 

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

genpop_obj

chord_dist <- dist.genpop(genpop_obj, method = 2, diag = TRUE, upper = TRUE)

# convert to matrix
chord_dist_matrix <- as.matrix(chord_dist)

# write to an external file for future use.
write.table(chord_dist_matrix, file = "KF_chord_dist_matrix.txt", sep = "\t", row.names = TRUE, quote = FALSE)

#################################################################################

# Run Radish
myRaster <- raster("/group/ctbrowngrp2/sophiepq/KitFoxGBS/LandGenR/Radish/KitFox-ESARPmodel-Raster1000x1000-FillAllCells.tif")

roadRaster <- raster("/group/ctbrowngrp2/sophiepq/KitFoxGBS/LandGenR/Radish/CaliforniaRoads-Reprojected-1000x1000-FillAllCells.tif")
roadRaster <- raster("/group/ctbrowngrp2/sophiepq/KitFoxGBS/LandGenR/Radish/CaliforniaRoads-Reprojected-1000x1000-MajorRoadsOnly.tif")



# scaling spatial covariates helps avoid numeric overflow
covariates <- raster::stack(list(kfsuit = raster::scale(myRaster),
                                roads = raster::scale(roadRaster)))


## Trying just the Major roads because I was getting a weird inverse relationship compared to what I expected.
covariates <- raster::stack(list(roads = raster::scale(roadRaster)))
                                

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

png("Majorroads.png", width = 800, height = 600)
plot(covariates[["roads"]])
points(df_points, pch = 19)
dev.off()

surface <- conductance_surface(covariates, df_points, directions = 8)

fit_nnls <- radish(chord_dist_matrix ~ kfsuit + roads, surface, 
                   radish::loglinear_conductance, radish::leastsquares)

summary(fit_nnls)


### This is the model summary for just the major roads alone (no kf suit or minor roads)

### This is the model summary for the kfsuit and roads (Major + Minor) combined ###
> summary(fit_nnls)
Conductance surface with 146970 vertices (11 focal) estimated by maximum likelihood
Call:   radish(formula = chord_dist_matrix ~ kfsuit + roads, data = surface, 
    conductance_model = radish::loglinear_conductance, measurement_model = radish::leastsquares)

Loglikelihood: 163.9698 (5 degrees freedom)
AIC: -317.9395 

Number of function calls: 43 
Number of Newton-Raphson steps: 10 
Norm of gradient at MLE: 1.947029e-07 

Nuisance parameters:
  alpha     beta      tau  
 0.1726  13.7574   6.9625  

Coefficients:
       Estimate Std. Error z value Pr(>|z|)   
kfsuit   2.0990     1.1358   1.848  0.06460 . 
roads    2.6079     0.8918   2.924  0.00345 **
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Correlation of Coefficients:
         kfsuit
roads 0.9366792

# refit with with a different measurement model that models
# dependence among pairwise measurements (radish::mlpe)
fit_mlpe <- radish(chord_dist_matrix ~ kfsuit + roads, surface, 
                   radish::loglinear_conductance, radish::mlpe)
summary(fit_mlpe)

### This is the model summary for just the major roads alone (no kf suit or minor roads)
> summary(fit_mlpe)
Conductance surface with 146970 vertices (11 focal) estimated by maximum likelihood
Call:   radish(formula = chord_dist_matrix ~ roads, data = surface, conductance_model = radish::loglinear_conductance, 
    measurement_model = radish::mlpe)

Loglikelihood: 199.5061 (5 degrees freedom)
AIC: -389.0122 

Number of function calls: 5 
Number of Newton-Raphson steps: 2 
Norm of gradient at MLE: 3.680808e-11 

Nuisance parameters:
  alpha     beta      tau      rho  
0.09507  0.62003  6.39888  2.63087  

Coefficients:
      Estimate Std. Error z value Pr(>|z|)
roads   -4.622        NaN     NaN      NaN
Warning messages:
1: In sqrt(diag(solve(x$fit$hessian))) : NaNs produced
2: In summary.radish(fit_mlpe) :
  Hessian matrix has negative eigenvalues: possibly a saddle point
3: In sqrt(1/diag(V)) : NaNs produced
4: In cov2cor(vcov) :
  diag(.) had 0 or NA entries; non-finite result is doubtful

### This is the model summary for the kfsuit and Major roads (no minor roads)###
> summary(fit_mlpe)
Conductance surface with 146970 vertices (11 focal) estimated by maximum likelihood
Call:   radish(formula = chord_dist_matrix ~ kfsuit + roads, data = surface, 
    conductance_model = radish::loglinear_conductance, measurement_model = radish::mlpe)

Loglikelihood: 201.0082 (6 degrees freedom)
AIC: -390.0164 

Number of function calls: 25 
Number of Newton-Raphson steps: 7 
Norm of gradient at MLE: 1.417013e-08 

Nuisance parameters:
 alpha    beta     tau     rho  
0.1531  0.8361  5.9783  3.2631  

Coefficients:
       Estimate Std. Error z value Pr(>|z|)   
kfsuit   0.9316     0.3159   2.949  0.00318 **
roads   -0.1823     1.2182  -0.150  0.88102   
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Correlation of Coefficients:
         kfsuit
roads 0.2689689

### This is the model summary for the kfsuit and roads (Major + Minor) combined ###
> summary(fit_mlpe)
Conductance surface with 146970 vertices (11 focal) estimated by maximum likelihood
Call:   radish(formula = chord_dist_matrix ~ kfsuit + roads, data = surface, 
    conductance_model = radish::loglinear_conductance, measurement_model = radish::mlpe)

Loglikelihood: 203.0221 (6 degrees freedom)
AIC: -394.0442 

Number of function calls: 33 
Number of Newton-Raphson steps: 9 
Norm of gradient at MLE: 2.579128e-10 

Nuisance parameters:
  alpha     beta      tau      rho  
0.08522  0.64523  6.58993  2.54468  

Coefficients:
       Estimate Std. Error z value Pr(>|z|)    
kfsuit   0.2086     0.2817   0.741    0.459    
roads    0.7371     0.1307   5.639 1.71e-08 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Correlation of Coefficients:
         kfsuit
roads 0.3542544

# visualisation:
png("KF Optimized Resistance Distance MajorRoads + KFSuit Model.png", width = 800, height = 600)
plot(fitted(fit_mlpe, "distance"), chord_dist_matrix, pch = 19,
     xlab = "Optimized resistance distance", ylab = "chord distance")
dev.off()

# visualise estimated conductance surface and asymptotic confidence intervals
fitted_conductance <- conductance(surface, fit_mlpe, quantile = 0.95)

png("KF FittedConductance  MajorRoads + KFSuit Model.png", width = 800, height = 600)
plot(fitted_conductance[["est"]], 
     main = "Fitted conductance surface\n(MajorRoads + KFSuit)")
points(df_points, pch = 19)
dev.off()

png("KF FittedConductance Lower95 MajorRoads Model.png", width = 800, height = 600)
plot(fitted_conductance[["lower95"]], 
     main = "Fitted conductance surface\n(roads)")
points(df_points, pch = 19)
dev.off()

png("KF FittedConductance Upper95 MajorRoads Model.png", width = 800, height = 600)
plot(fitted_conductance[["upper95"]], main = 
     "Fitted conductance surface\n(roads)")
points(df_points, pch = 19)
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
distances <- radish_distance(theta, ~ kfsuit + roads, 
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
