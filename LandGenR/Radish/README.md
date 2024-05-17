```
library(radish)
library(raster)

data(melip)

# scaling spatial covariates helps avoid numeric overflow
covariates <- raster::stack(list(altitude = raster::scale(melip.altitude), 
                                 forestcover = raster::scale(melip.forestcover)))
png("altitude.png", width = 800, height = 600)
plot(covariates[["altitude"]])
points(melip.coords, pch = 19)
dev.off()

png("forestcover.png", width = 800, height = 600)
plot(covariates[["forestcover"]])
points(melip.coords, pch = 19)
dev.off()

surface <- conductance_surface(covariates, melip.coords, directions = 8)

fit_nnls <- radish(melip.Fst ~ forestcover + altitude, surface, 
                   radish::loglinear_conductance, radish::leastsquares)
summary(fit_nnls)

# refit with with a different measurement model that models
# dependence among pairwise measurements (radish::mlpe)
fit_mlpe <- radish(melip.Fst ~ forestcover + altitude, surface, 
                   radish::loglinear_conductance, radish::mlpe)
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
