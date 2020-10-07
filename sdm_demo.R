# This script demonstrates some methods of using R in SDM, and attempts to teach 
# SDM concepts through demonstration this is more or less the first draft of the 
# RMarkdown File in this repo

# For downloading occurrence data from multiple sources
library(spocc)
# For data manipulation and presentation
library(tidyr)
library(dplyr)
library(ggplot2)
# Provides various functions relating to spatial analysis
library(maptools)
library(sp)
# Reads and manipulates vector data
library(rgdal)
# Reads and manipulates raster data
library(raster)
# SDM packages
library(dismo)
library(biomod2)


# Let's ready some political boundaries in for visualization
usa.shp <- readOGR("USA")
CA_OR.shp <- usa.shp[usa.shp@data$NAME %in% c("California", "Oregon"),]

# Source occurrence data from iNaturalist
dc <- occ(query = "Darlingtonia californica", from = 'inat')
# dc <- occ(query = "Darlingtonia californica", from = 'gbif')

# Horrible nasty function that makes species occurrence data into lon, lat dataframe
dc.df <-  as.data.frame(dc$inat$data$Darlingtonia_californica%>%
  dplyr::filter(identifications_most_agree == TRUE) %>%
  dplyr::select(location)) %>% 
  filter(!is.na(location)) %>% 
  separate(location, c("lat", "lon"), ",") %>% 
  select(lon, lat)

# Lat and Lon columns must be specified as numeric values
dc.df$lon <- as.numeric(dc.df$lon)
dc.df$lat <- as.numeric(dc.df$lat)

# Let's take a look at our raw occurrence points
sp::plot(CA_OR.shp, main = "Darlingtonia californica (iNaturalist)")
points(dc.df$lon, dc.df$lat)

# Download bioclimatic variables from worldclim
# run ??getData to see other options for climate data
bioclim_global.stack <- raster::getData(name = "worldclim",
                        var = 'bio', 
                        res = 2.5)

# This replaces bioX name with more descriptive acronym
# We choose 4 climatic variables here more for examplary purposes than biologically
# sound reason
names(bioclim_global.stack)[names(bioclim_global.stack) %in% c("bio1", "bio10", "bio11", "bio12")] <- 
  c("MAT", "MTWQ", "MTCQ", "MAP")

# This crops the extent to Oregon and California
bioclim.stack <- crop(bioclim_global.stack, extent(CA_OR.shp))

# Let's only look at these 4 variables
predictors <- raster::subset(bioclim.stack, c("MAT", "MTWQ", "MTCQ", "MAP"))
# bio1 = Mean Annual Temperature (MAT)
# bio10 = Mean Temperature of Warmest Quarter (MTWQ)
# bio11 = Mean Temperature of Coldest Quarter (MTCQ)
# bio12 = Annual Precipitation (MAP)
# etc.

# Let's see what Mean Annual temperature raster looks like
sp::plot(predictors$MAT, main = "Mean Annual Temperature")


# Organizing species presence and environmental predictors into data.frames

# These are the predictor values at locations of species presence
presvals <- raster::extract(predictors, dc.df)

# These are 500 random locations, used as in place of absence values as 
# 'pseudoabsences' (the species probably doesn't occur at any random point)
backgr <- randomPoints(predictors, 500)

# predictor values at random locations
absvals <- raster::extract(predictors, backgr)

# We know that expected habitat suitability (Ey) is 1 for areas where the species
# was found, and we assume it's 0 for the random background points
Ey <- c(rep(1, nrow(presvals)), rep(0, nrow(absvals)))

# Now here we have a dataframe with the response variable (Ey) and corresponding
# predictor values
sdmdata <- data.frame(cbind(Ey, rbind(presvals, absvals)))
View(sdmdata)

# Look for collinearity using a scatterplot matrix
pairs(sdmdata[,2:length(sdmdata)], cex = 0.1, fig=TRUE)

# MAT and MTWQ are collinear across this area!
# So let's remove bio1 and try again
predictors <- raster::subset(bioclim.stack, c("MTWQ", "MTCQ", "MAP"))
presvals <- raster::extract(predictors, dc.df)
backgr <- randomPoints(predictors, 500)
absvals <- raster::extract(predictors, backgr)
Ey <- c(rep(1, nrow(presvals)), rep(0, nrow(absvals)))
sdmdata <- data.frame(cbind(Ey, rbind(presvals, absvals)))

pairs(sdmdata[,2:length(sdmdata)], cex = 0.1, fig=TRUE)
# Better! A popular metric of collinearity analysis is Variance Inlation Factor
# (see https://www.rdocumentation.org/packages/car/versions/3.0-9/topics/vif)


# Let's start Modeling! Here's the most basic SDM

## Linear Model SDM
# E(y) = mx + b + e
# E(y) is expected habitat suitability
# b = intercept
# m = coefficient
# x = predictor
# e = error

# Assumes: Normal distributed error
#          mean = 0
#          Variance(y) is constant

sdm_lm <- lm(Ey ~ MTWQ, data = sdmdata)
plot(sdmdata$MTWQ, sdmdata$Ey, 
     main = "D. californica Linear SDM",
     xlab = "Mean Temp of Warmest Quarter",
     ylab = "Habitat Suitability"
     )
abline(sdm_lm)

# The simplest SDM ever! Do this for each predictor and average habitat 
# suitability in each cell (I think?)


## Generalized Linear Model SDM
# g(E(y)) = LP = mx + b + e
# g() is link function
# LP = Linear Predictor

# Does not have to have normal distribution of error
# does not have to have homogenous variance. Options of different link functions.
# Different link functions have their own assumptions

plot(sdmdata$MTWQ, sdmdata$Ey, 
     main = "D. californica GLM SDM",
     xlab = "Mean Temp of Warmest Quarter",
     ylab = "Habitat Suitability"
)

# A GLM with a Gaussian Link Function is a Linear Model
sdm_lm <- glm(Ey ~ MTWQ, data = sdmdata, family = gaussian)
abline(sdm_lm, lwd = 3, col = "blue")
legend("topright", legend = c("Linear Model"), col = c("blue"), lty=1, lwd = 3)

# In SDM we use the logit link function, to turn the binary presence/absence
# into a continuous response
# GLM: g(E(y)) = LP = log10(Ey / (1 - Ey)) = mx + b + e
# Solve for Ey:
# Ey = e^(mx + b + e) / (1 + e^(mx + b + e))
# Go to slide >>>

sdm_glm <- glm(Ey ~ MTWQ, data = sdmdata, family = binomial)

# This is the equation of line 142, written as an R function
Ey_calc <- function(x){
  y_int <- as.numeric(sdm_glm$coefficients["(Intercept)"])
  m <- as.numeric(sdm_glm$coefficients["MTWQ"])
  LP = m*x + y_int
  return(exp(LP) / (1 + exp(LP)))
}

# Need to sort in order for plot to add line properly
MTWQ_sorted <- sort(sdmdata$MTWQ)
Ey <- Ey_calc(MTWQ_sorted)

lines(MTWQ_sorted, Ey, lwd = 3, col = "red")
legend("topright", legend = c("Linear Model", "GLM Model"), col = c("blue", "red"), lty=1, lwd = 3)

# In practice, GLM SDMs are made using machine learning, which makes 
# many versions with slight tweaks to LP and equations that relate the dif. 
# predictor responses. Here's a tweak to the LP:
Ey_calc2 <- function(x){
  y_int <- as.numeric(sdm_glm$coefficients["(Intercept)"])
  m <- as.numeric(sdm_glm$coefficients["MTWQ"])
  LP = m*x + x^.05 + y_int
  return(exp(LP) / (1 + exp(LP)))
}
Ey <- Ey_calc2(MTWQ_sorted)

lines(MTWQ_sorted, Ey, lwd = 3, col = "green")
legend("topright", legend = c("Linear Model", "GLM Model", "GLM Model 2"), 
       col = c("blue", "red", "green"), lty=1, lwd = 3)



# Here's a function we'll use to plot SDM projections
project.sdm <- function(prediction, plotName){
  sp::plot(prediction, main = plotName)
  sp::plot(CA_OR.shp, add = T)
  points(dc.df, pch = 16, cex = .2)
  #legend("bottomright", legend = "obs. occurrences", pch = 16, cex=.2)
}

# Here's a map of projection of the SDM model
prediction_glm <- raster::predict(bioclim.stack, sdm_glm)
# we need to log transform it to give us the expected presence
# We have g(E(x)), where E(x) = e^(LP) / (1 + e^LP)
# where E(x) is expected presence or absence of species
prediction_glm.Ey <- exp(prediction_glm) / (1+exp(prediction_glm))

par(mfrow=c(1,2))
plot(MTWQ_sorted, Ey, lwd = 3, col = "green", type = 'l',
     main = "GLM",
     xlab = "Mean Temp of Warmest Quarter",
     ylab = "Habitat Suitability")

project.sdm(prediction_glm.Ey, "GLM SDM (D. californica), MTWQ only")

# Let's make  GLM with all 3 variables
par(mfrow=c(1,1))
sdm_glm3 <- glm(Ey ~ MTWQ + MTCQ + MAP, data = sdmdata, family = binomial)
prediction_glm3 <- raster::predict(bioclim.stack, sdm_glm3)
prediction_glm3.Ey <- exp(prediction_glm3) / (1+exp(prediction_glm3))
project.sdm(prediction_glm3.Ey, "GLM SDM (D. californica)")



# Now let's use real software based SDMs
# Here's Machine Learning GLM
# Let's try sdm package
library(sdm)
# sdm packages can do GLM, GAM, GBM, RF, TREE, MARS, SVM

# Prep data format for the sdm package
sdm.pkg.df_pres <- cbind(dc.df, presvals)
sdm.pkg.df_pres$Ey <- 1
names(sdm.pkg.df_pres)[1:2] <- c("x", "y")
sdm.pkg.df_abs <- data.frame(cbind(backgr, absvals))
sdm.pkg.df_abs$Ey <- 0
sdmdf_sdmpkg <- rbind(sdm.pkg.df_pres, sdm.pkg.df_abs)
sdmdata_sdmpkg <- sdmData(Ey ~ MTWQ + MTCQ + MAP, train = sdmdf_sdmpkg)

# Run the model and project
sdm_ml.glm <- sdm::sdm(Ey ~ MTWQ + MTCQ + MAP, data = sdmdata_sdmpkg, methods=c("glm"))
prediction_ml.glm <- raster::predict(sdm_ml.glm, bioclim.stack)
project.sdm(prediction_ml.glm, "ML GLM SDM (D. californica)")

# Use this function to determine which variables were more important
getVarImp(sdm_ml.glm)


# Let's try GAM
# Run the model and project
sdm_ml.gam <- sdm::sdm(Ey ~ MTWQ + MTCQ + MAP, data = sdmdata_sdmpkg, methods=c("gam"))
prediction_ml.gam <- raster::predict(sdm_ml.gam, bioclim.stack)
project.sdm(prediction_ml.gam, "ML GAM SDM (D. californica)")

getVarImp(sdm_ml.gam)


# Now Let's try more strictly Machine Learning: Random Forest
# Run the model and project
sdm_rf <- sdm::sdm(Ey ~ MTWQ + MTCQ + MAP, data = sdmdata_sdmpkg, methods=c("rf"))
prediction_rf <- raster::predict(sdm_rf, bioclim.stack)
project.sdm(prediction_rf, "Random Forest SDM (D. californica)")

getVarImp(sdm_rf)


# ENSEMBLE let's make an ensemble model of all of them
# a number of methods can be used, see documentation, but they include 
# weighted mean, unweighted mean, median, entropy, etc
sdm_glm.gam.rf <- sdm::sdm(Ey ~ MTWQ + MTCQ + MAP, data = sdmdata_sdmpkg, methods=c("glm","gam","rf"))
sdm_ensemble <- sdm::ensemble(sdm_glm.gam.rf, bioclim.stack, 
                              setting=list(method="weighted", stat="TSS"))
project.sdm(sdm_ensemble, "Ensemble (D. californica)")

getVarImp(sdm_glm.gam.rf)



# Let's see what dismo has
# BIOCLIM, DOMAIN, MaxEnt

# but powerful model evaluation software, and software to test the importance of 
# particular variables, producing bioclim variables for monthly climate data (CMIP5 is great)
# niche equivalency, niche overlap, a tool that finds the best threshold based on 
# particular accuracy metrics

# Here's BIOCLIM, one of the oldest models, and EE
# Puts a convex hull around all the points in n-dimensional niche space
# Not very powerful at extrapolation

sdm_bioclim <- bioclim(presvals)
response(sdm_bioclim)
prediction_bioclim <- dismo::predict(sdm_bioclim, bioclim.stack)
project.sdm(prediction_bioclim, "BIOCLIM SDM (D. californica)")


# Maxent, need to install maxent (https://biodiversityinformatics.amnh.org/open_source/maxent/)
# and place it here:
system.file("java", package="dismo")

sdm_maxent <- maxent(predictors, dc.df)
prediction_maxent <- dismo::predict(sdm_maxent, bioclim.stack)
project.sdm(prediction_maxent, "MaxEnt SDM (D. californica)")

# Look at response for each predictor
response(sdm_maxent)

# Look at variable contribution for maxent
plot(sdm_maxent)



# Make an ensemble model!!
ensemble_bioclim.maxent <- mean(prediction_bioclim, prediction_maxent)
# Seems less sophisticated than the 'sdm' package!
plot(ensemble_bioclim.maxent, main = "MaxEnt BIOCLIM ensemble")
points(dc.df, pch = 16, cex = .2)
legend("bottomright", legend = "obs. occurrences", pch = 16)


# Amalgamation of all models:

par(mfrow=c(3,3))
project.sdm(prediction_glm.Ey, "GLM SDM (D. californica), MTWQ only")
project.sdm(prediction_glm3.Ey, "GLM SDM (D. californica)")
project.sdm(prediction_ml.glm, "ML GLM SDM (D. californica)")
project.sdm(prediction_ml.gam, "ML GAM SDM (D. californica)")
project.sdm(prediction_rf, "Random Forest SDM (D. californica)")
project.sdm(sdm_ensemble, "Ensemble (D. californica)")
project.sdm(prediction_bioclim, "BIOCLIM SDM (D. californica)")
project.sdm(prediction_maxent, "MaxEnt SDM (D. californica)")
plot(ensemble_bioclim.maxent, main = "MaxEnt BIOCLIM ensemble")
points(dc.df, pch = 16, cex = .2)
par(mfrow=c(1,1))
