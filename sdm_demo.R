# This script demonstrates some methods of using R in SDM, and attempts to teach
# SDM concepts through demonstration this is more or less the first draft of the
# RMarkdown File in this repo
#
# NOTE: Updated to use 'terra' for raster and vector data, and 'sf' for vector data, replacing 'raster' and 'rgdal'.
#
# Minimal changes were made to preserve the educational structure.

# For downloading occurrence data from multiple sources
library(spocc)
# For data manipulation and presentation
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
# Provides various functions relating to spatial analysis
# library(maptools)
# library(sp)
# Reads and manipulates vector data
# library(rgdal) # replaced by sf/terra
# Reads and manipulates raster data
# library(raster) # replaced by terra
library(terra)
library(tidyterra)
# library(sf)
# SDM packages
library(dismo)
library(sdm)
library(geodata)
library(patchwork)


# Let's ready some political boundaries in for visualization
# usa.shp <- readOGR("USA")
usa_borders <- sf::st_read("data/usa_borders.gpkg")
CA_OR_borders <- usa_borders[usa_borders$NAME %in% c("California", "Oregon"), ]

# Source occurrence data from iNaturalist
dc <- spocc::occ(query = "Darlingtonia californica", from = "inat")
# Could just as easily source from gbif:
# dc <- occ(query = "Darlingtonia californica", from = 'gbif')

# Horrible nasty function that makes species occurrence data into lon, lat dataframe
dc.df <- as.data.frame(dc$inat$data$Darlingtonia_californica %>%
  dplyr::filter(identifications_most_agree == TRUE) %>%
  dplyr::select(location)) %>%
  filter(!is.na(location)) %>%
  separate(location, c("lat", "lon"), ",") %>%
  dplyr::select(lon, lat)

# Lat and Lon columns must be specified as numeric values
dc.df$lon <- as.numeric(dc.df$lon)
dc.df$lat <- as.numeric(dc.df$lat)

# Let's take a look at our raw occurrence points
ggplot() +
  geom_sf(data = CA_OR_borders) +
  geom_point(data = dc.df, aes(x = lon, y = lat)) +
  labs(title = "Darlingtonia californica (iNaturalist)") +
  theme_minimal()
# Wow there's a D. californica waayyy to the east, which must be erroneous.

# Let's filter out only areas in Oregon and California
dc_target.df <- dc.df |>
  # Convert to sf object for spatial filtering
  sf::st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) %>%
  sf::st_intersection(CA_OR_borders) |>
  # Back to data.frame
  sf::st_drop_geometry() |>
  as_tibble()

ggplot() +
  geom_sf(data = CA_OR_borders) +
  geom_point(data = dc.df, aes(x = lon, y = lat)) +
  labs(title = "Darlingtonia californica (iNaturalist)") +
  theme_minimal()
# Great.

# Now let's look at the bioclimatic variables

# Download bioclimatic variables from worldclim
bioclim_world.stack <- geodata::worldclim_global(var = "bio", res = 10, path = "data")
# This replaces bioX name with more descriptive acronym
# We choose 4 climatic variables here more for examplary purposes than biologically
# sound reason
names(bioclim_world.stack)[str_which(names(bioclim_world.stack), pattern = "bio_(1|10|11|12)$")] <-
  c("MAT", "MTWQ", "MTCQ", "MAP")

# This crops the extent to Oregon and California
# Need to mask the raster to the actual shape of CA/OR to avoid NaN values outside the states
bioclim.stack <- terra::crop(bioclim_world.stack, CA_OR.shp, mask = TRUE)

# Let's only look at these 4 variables
predictors <- terra::subset(bioclim.stack, c("MAT", "MTWQ", "MTCQ", "MAP"))
# bio1 = Mean Annual Temperature (MAT)
# bio10 = Mean Temperature of Warmest Quarter (MTWQ)
# bio11 = Mean Temperature of Coldest Quarter (MTCQ)
# bio12 = Annual Precipitation (MAP)
# etc.

# Let's see what Mean Annual temperature raster looks like

ggplot() +
  tidyterra::geom_spatraster(data = predictors[["MAT"]]) +
  scale_fill_viridis_c() +
  labs(title = "Mean Annual Temperature") +
  theme_minimal()


# Organizing species presence and environmental predictors into data.frames

# These are the predictor values at locations of species presence
presvals <- terra::extract(predictors, dc_target.df[, c("lon", "lat")])[, -1]

# These are 500 random locations, used as in place of absence values as
# 'pseudoabsences' (the species probably doesn't occur at any random point)
set.seed(1)
backgr <- as.data.frame(terra::spatSample(predictors, 500, method = "random", xy = TRUE)[, c("x", "y")])
colnames(backgr) <- c("lon", "lat")

# predictor values at random locations
absvals <- terra::extract(predictors, backgr)[, -1]

# We know that expected habitat suitability (Ey) is 1 for areas where the species
# was found, and we assume it's 0 for the random background points
Ey <- c(rep(1, nrow(presvals)), rep(0, nrow(absvals)))

# Now here we have a dataframe with the response variable (Ey) and corresponding
# predictor values
sdmdata <- data.frame(cbind(Ey, rbind(presvals, absvals)))
View(sdmdata)

# Look for collinearity using a scatterplot matrix
pairs(sdmdata[, 2:length(sdmdata)], cex = 0.1, fig = TRUE)

# MAT and MTWQ are collinear across this area!
# So let's remove bio1 and try again
predictors <- terra::subset(bioclim.stack, c("MTWQ", "MTCQ", "MAP"))
presvals <- terra::extract(predictors, dc_target.df[, c("lon", "lat")])[, -1]
set.seed(1)
backgr <- as.data.frame(terra::spatSample(predictors, 500, method = "random", xy = TRUE)[, c("x", "y")])
colnames(backgr) <- c("lon", "lat")
absvals <- terra::extract(predictors, backgr)[, -1]
Ey <- c(rep(1, nrow(presvals)), rep(0, nrow(absvals)))
sdmdata <- data.frame(cbind(Ey, rbind(presvals, absvals)))

pairs(sdmdata[, 2:length(sdmdata)], cex = 0.1, fig = TRUE)
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

# Fit linear model
sdm_lm <- lm(Ey ~ MTWQ, data = sdmdata)

# Now get Ey for every point in the dataset
sdmdata$Ey_pred <- predict(sdm_lm, sdmdata)

plot_lm <- ggplot(sdmdata) +
  geom_point(aes(x = MTWQ, y = Ey)) +
  geom_line(aes(x = MTWQ, y = Ey_pred), color = "red") +
  labs(
    title = "D. californica Linear SDM",
    x = "Mean Temp of Warmest Quarter",
    y = "Habitat Suitability"
  ) +
  theme_minimal()

plot_lm


## Generalized Linear Model SDM
# g(E(y)) = LP = mx + b + e
# g() is link function
# LP = Linear Predictor

# Does not have to have normal distribution of error
# does not have to have homogenous variance. Options of different link functions.
# Different link functions have their own assumptions

# A GLM with a Gaussian Link Function is a Linear Model
# Solve for Ey:
# Ey = mx + b + e
sdm_glm_gaussian <- glm(Ey ~ MTWQ, data = sdmdata, family = gaussian)

Ey_calc_gaussian <- function(x) {
  y_int <- as.numeric(sdm_glm_gaussian$coefficients["(Intercept)"])
  m <- as.numeric(sdm_glm_gaussian$coefficients["MTWQ"])
  LP <- (m * x) + y_int
  return(LP)
}

sdmdata$Ey_pred_glm_gauss <- Ey_calc_gaussian(sdmdata$MTWQ)

plot_glm_gaussian <- ggplot(sdmdata) +
  geom_point(aes(x = MTWQ, y = Ey)) +
  geom_line(
    aes(x = MTWQ, y = Ey_pred_glm_gauss, color = "Gaussian"),
    linewidth = 1
  ) +
  labs(
    title = "D. californica GLM SDM",
    x = "Mean Temp of Warmest Quarter",
    y = "Habitat Suitability"
  ) +
  scale_color_manual(name = "Link Function", values = c("Gaussian" = "blue")) +
  theme_minimal()

plot_glm_gaussian

# In SDM we can use the logit link function, to turn the binary presence/absence
# into a continuous response
# GLM: g(E(y)) = LP = log10(Ey / (1 - Ey)) = mx + b + e
# Solve for Ey:
# Ey = e^(mx + b + e) / (1 + e^(mx + b + e))
# Go to slide >>>

sdm_glm_logit <- glm(Ey ~ MTWQ, data = sdmdata, family = binomial)

# This is the equation of line 142, written as an R function
Ey_calc_logit <- function(x) {
  y_int <- as.numeric(sdm_glm_logit$coefficients["(Intercept)"])
  m <- as.numeric(sdm_glm_logit$coefficients["MTWQ"])
  LP <- m * x + y_int
  return(exp(LP) / (1 + exp(LP)))
}

# Need to sort in order for plot to add line properly
sdmdata$Ey_pred_glm_logit <- Ey_calc_logit(sdmdata$MTWQ)

plot_glm_logit <- plot_glm_gaussian +
  geom_line(
    data = sdmdata,
    aes(x = MTWQ, y = Ey_pred_glm_logit, color = "Logit"), linewidth = 1
  ) +
  scale_color_manual(name = "Link Function", values = c("Gaussian" = "blue", "Logit" = "red")) +
  labs(title = "D. californica GLM SDM") +
  theme_minimal()

plot_glm_logit

# Here's a function we'll use to plot SDM projections
project.sdm <- function(prediction, plotName) {
  ggplot() +
    geom_spatraster(data = prediction) +
    scale_fill_viridis_c() +
    geom_sf(data = sf::st_as_sf(CA_OR.shp), fill = NA) +
    geom_point(data = dc_target.df, aes(x = lon, y = lat), size = 0.2) +
    labs(title = plotName) +
    theme_minimal()
}

# Here's a map of projection of the SDM model
prediction_glm <- terra::predict(bioclim.stack, sdm_glm_logit)
# we need to log transform it to give us the expected presence
# We have g(E(x)), where E(x) = e^(LP) / (1 + e^LP)
# where E(x) is expected presence or absence of species
prediction_glm.Ey <- exp(prediction_glm) / (1 + exp(prediction_glm))

par(mfrow = c(1, 2))
plot(MTWQ_sorted, Ey,
  lwd = 3, col = "green", type = "l",
  main = "GLM",
  xlab = "Mean Temp of Warmest Quarter",
  ylab = "Habitat Suitability"
)


# Put these two plots together with patchwork package
project.sdm(prediction_glm.Ey, "GLM SDM (D. californica), MTWQ only") +
  ggplot(sdmdata) +
  geom_line(aes(x = MTWQ, y = Ey_pred_glm_logit), linewidth = 1) + theme_minimal()

# Let's make  GLM with all 3 variables
sdm_glm3 <- glm(Ey ~ MTWQ + MTCQ + MAP, data = sdmdata, family = binomial)
prediction_glm3 <- terra::predict(bioclim.stack, sdm_glm3)
prediction_glm3.Ey <- exp(prediction_glm3) / (1 + exp(prediction_glm3))

project.sdm(prediction_glm3.Ey, "GLM SDM (D. californica)")



# Now let's use real software based SDMs
# Here's Machine Learning GLM
# Let's try sdm package
# sdm packages can do GLM, GAM, GBM, RF, TREE, MARS, SVM

# Prep data format for the sdm package
sdm.pkg.df_pres <- cbind(dc_target.df |> dplyr::select(lon, lat), presvals)
sdm.pkg.df_pres$Ey <- 1
sdm.pkg.df_abs <- data.frame(cbind(backgr, absvals))
sdm.pkg.df_abs$Ey <- 0
sdmdf_sdmpkg <- rbind(sdm.pkg.df_pres, sdm.pkg.df_abs)
sdmdata_sdmpkg <- sdmData(Ey ~ MTWQ + MTCQ + MAP, train = sdmdf_sdmpkg)

# Run the model and project
sdm_ml.glm <- sdm::sdm(Ey ~ MTWQ + MTCQ + MAP, data = sdmdata_sdmpkg, methods = c("glm"))
prediction_ml.glm <- predict(sdm_ml.glm, bioclim.stack)
project.sdm(prediction_ml.glm, "ML GLM SDM (D. californica)")

# Use this function to determine which variables were more important
getVarImp(sdm_ml.glm)


# Let's try GAM
# Run the model and project
sdm_ml.gam <- sdm::sdm(Ey ~ MTWQ + MTCQ + MAP, data = sdmdata_sdmpkg, methods = c("gam"))
prediction_ml.gam <- predict(sdm_ml.gam, bioclim.stack)
project.sdm(prediction_ml.gam, "ML GAM SDM (D. californica)")

getVarImp(sdm_ml.gam)


# Now Let's try more strictly Machine Learning: Random Forest
# Run the model and project
sdm_rf <- sdm::sdm(Ey ~ MTWQ + MTCQ + MAP, data = sdmdata_sdmpkg, methods = c("rf"))
prediction_rf <- predict(sdm_rf, bioclim.stack)
project.sdm(prediction_rf, "Random Forest SDM (D. californica)")

getVarImp(sdm_rf)


# ENSEMBLE let's make an ensemble model of all of them
# a number of methods can be used, see documentation, but they include
# weighted mean, unweighted mean, median, entropy, etc
sdm_glm.gam.rf <- sdm::sdm(Ey ~ MTWQ + MTCQ + MAP, data = sdmdata_sdmpkg, methods = c("glm", "gam", "rf"))
sdm_ensemble <- sdm::ensemble(sdm_glm.gam.rf, bioclim.stack,
  setting = list(method = "weighted", stat = "TSS")
)
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
# prediction_bioclim <- dismo::predict(sdm_bioclim, bioclim.stack)
prediction_bioclim <- terra::predict(bioclim.stack, sdm_bioclim)
project.sdm(prediction_bioclim, "BIOCLIM SDM (D. californica)")


# Maxent, need to install maxent (https://biodiversityinformatics.amnh.org/open_source/maxent/)
# and place it here:
# Dismo maxent doesn't use terra but uses the older raster package, so i'm disabling this for now
# system.file("java", package = "dismo")

# sdm_maxent <- dismo::maxent(predictors, dc.df)
# prediction_maxent <- terra::predict(bioclim.stack, sdm_maxent)
# project.sdm(prediction_maxent, "MaxEnt SDM (D. californica)")

# # Look at response for each predictor
# response(sdm_maxent)

# # Look at variable contribution for maxent
# # Convert maxent plot to ggplot
# maxent_contrib <- data.frame(
#   variable = names(sdm_maxent@results),
#   contribution = sdm_maxent@results
# )
# ggplot(maxent_contrib, aes(x = reorder(variable, contribution), y = contribution)) +
#   geom_bar(stat = "identity") +
#   coord_flip() +
#   labs(title = "MaxEnt Variable Contributions") +
#   theme_minimal()



# # Make an ensemble model!!
# ensemble_bioclim.maxent <- (prediction_bioclim + prediction_maxent) / 2
# # Seems less sophisticated than the 'sdm' package!
# plot(ensemble_bioclim.maxent, main = "MaxEnt BIOCLIM ensemble")
# points(dc_target.df, pch = 16, cex = .2)
# legend("bottomright", legend = "obs. occurrences", pch = 16)


# Amalgamation of all models:
# Create a grid of plots using patchwork

p1 <- project.sdm(prediction_glm.Ey, "GLM SDM (D. californica), MTWQ only")
p2 <- project.sdm(prediction_glm3.Ey, "GLM SDM (D. californica)")
p3 <- project.sdm(prediction_ml.glm, "ML GLM SDM (D. californica)")
p4 <- project.sdm(prediction_ml.gam, "ML GAM SDM (D. californica)")
p5 <- project.sdm(prediction_rf, "Random Forest SDM (D. californica)")
p6 <- project.sdm(sdm_ensemble, "Ensemble (D. californica)")
p7 <- project.sdm(prediction_bioclim, "BIOCLIM SDM (D. californica)")
# p8 <- project.sdm(prediction_maxent, "MaxEnt SDM (D. californica)")
# p9 <- ggplot() +
#   geom_spatraster(data = ensemble_bioclim.maxent) +
#   scale_fill_viridis_c() +
#   geom_point(data = dc_target.df, aes(x = lon, y = lat), size = 0.2) +
#   labs(title = "MaxEnt BIOCLIM ensemble") +
#   theme_minimal()

# Arrange plots in a 3x3 grid
(p1 + p2 + p3) / (p4 + p5 + p6) / (p7) # + p8 + p9)
