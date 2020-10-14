An overview of Species Distribution Modeling in R
=================================================

### This demo was largely based on [Hijmans & Elith, 2013](https://cran.r-project.org/web/packages/dismo/vignettes/sdm.pdf), [Franklin 2010](https://doi.org/10.1017/CBO9780511810602) and the [sdm](https://cran.r-project.org/web/packages/sdm/sdm.pdf) and [dismo](https://cran.r-project.org/web/packages/dismo/dismo.pdf) documentation pages. This serves to provide a walkthrough of SDM and popular implementations in R as of Fall 2020. I’ve included the extent of my understanding of how the different methods work.

### If you already have an understanding of the motivation of SDM, continue here to learn more about implementation. If you’d like a brief refresher on the motivation and general application of SDM take a look at the powerpoint in this git repository.

### Feel free to download the sdm_demo.R script, install the necessary packages, and run and manipulate this code as you see fit.

Let’s start by reading in the necessary packages

``` r
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
library(sdm)
library(usdm) # For vifstep function
```

Data Preparation
================

Species occurrence data
-----------------------

We’ll use the west coast carnivorous plant, [*Darlingtonia
californica*](https://en.wikipedia.org/wiki/Darlingtonia_californica),
as our species of interest. ![*Darlingtonia
Californica*](/includes/dar.jpg)

``` r
# Source occurrence data from iNaturalist
dc <- occ(query = "Darlingtonia californica", from = 'inat')
```

    ## 2 
    ## 3

``` r
# Could just as easily source from gbif:
# dc <- occ(query = "Darlingtonia californica", from = 'gbif')

# Horrible nasty function that makes species occurrence data into lon, lat dataframe
# This is specific to the format that iNaturalist puts data in
dc.df <-  as.data.frame(dc$inat$data$Darlingtonia_californica%>%
  dplyr::filter(identifications_most_agree == TRUE) %>%
  dplyr::select(location)) %>% 
  filter(!is.na(location)) %>% 
  separate(location, c("lat", "lon"), ",") %>% 
  select(lon, lat)

# Lat and Lon columns must be specified as numeric values
dc.df$lon <- as.numeric(dc.df$lon)
dc.df$lat <- as.numeric(dc.df$lat)
head(dc.df)
```

    ##         lon      lat
    ## 1 -123.9390 41.67295
    ## 2 -124.0979 44.04678
    ## 3 -124.0980 44.04678
    ## 4 -120.9603 39.81911
    ## 5 -120.8786 40.05013
    ## 6 -120.8630 40.07212

``` r
# Let's read some political boundaries in for visualization
usa.shp <- readOGR("USA", verbose=F)
CA_OR.shp <- usa.shp[usa.shp@data$NAME %in% c("California", "Oregon"),]

# Let's take a look at our raw occurrence points
sp::plot(CA_OR.shp, main = "Darlingtonia californica (iNaturalist)")
points(dc.df$lon, dc.df$lat)
```

![](README_files/figure-markdown_github/occurrence%20data-1.png)

Environmental predictor data
----------------------------

``` r
# Download bioclimatic variables from worldclim
# run ??getData to see other options for climate data
bioclim_global.stack <- raster::getData(name = "worldclim",
                        var = 'bio', 
                        res = 2.5)

# This replaces bioX name with more descriptive acronym
# We choose 4 climatic variables here more for exemplary purposes than biologically
# sound reasoning, but all predictors chosen should be motivated/informed
# by biological knowledge
names(bioclim_global.stack)[names(bioclim_global.stack) %in% c("bio1", "bio10", "bio11", "bio12")] <- 
  c("MAT", "MTWQ", "MTCQ", "MAP")
# bio1 = Mean Annual Temperature (MAT)
# bio10 = Mean Temperature of Warmest Quarter (MTWQ)
# bio11 = Mean Temperature of Coldest Quarter (MTCQ)
# bio12 = Annual Precipitation (MAP)
# etc.

# This crops the extent to Oregon and California
bioclim.stack <- crop(bioclim_global.stack, extent(CA_OR.shp))

# Let's only look at these 4 variables
predictors <- raster::subset(bioclim.stack, c("MAT", "MTWQ", "MTCQ", "MAP"))

# Let's see what Mean Annual temperature raster looks like
sp::plot(predictors$MAT, main = "Mean Annual Temperature")
```

![](README_files/figure-markdown_github/predictor%20data-1.png)

SDM-ready dataframe of predictor values and species presence/(pseudo)absence
----------------------------------------------------------------------------

``` r
# These are the predictor values at locations of species presence
presvals <- raster::extract(predictors, dc.df)

# These are 500 random locations, used as in place of absence values as 
# 'pseudoabsences' (the species probably doesn't occur at any random point)
# True 'absences' are better, but folks don't often record all the places where 
# they don't find species
backgr <- randomPoints(predictors, 500)

# predictor values at random locations
absvals <- raster::extract(predictors, backgr)

# We know that probability of presence is 1 for areas where the species
# was found, and we assume it's 0 for the random background points
Y <- c(rep(1, nrow(presvals)), rep(0, nrow(absvals)))

# Now here we have a dataframe with the response variable (Y) and corresponding
# predictor values
sdmdata <- data.frame(cbind(Y, rbind(presvals, absvals)))
head(sdmdata)
```

    ##   Y MAT MTWQ MTCQ  MAP
    ## 1 1 102  158   52 1833
    ## 2 1 111  156   70 2003
    ## 3 1 111  156   70 2003
    ## 4 1  88  170   16 1213
    ## 5 1  75  160   -2  956
    ## 6 1  86  172    6  915

``` r
# Collinearity can cause problems in Species Distribution Models. My understanding
# is that if two or more predictors are collinear across the environmental space, then 
# it's difficult to determine which predictor is actually influencing the distribution
# of the species (if not both). If the predictors are not collinear in the areas
# to which the SDM is projected then the model won't know how to assign probability
# of habitat suitability to the independent predictors
# We can check for collinearity by visually inspecting a scatterplot matrix
pairs(sdmdata[,2:length(sdmdata)], cex = 0.1, fig=TRUE)
```

![](README_files/figure-markdown_github/SDM%20prep-1.png)

``` r
# MAT and MTWQ appear to be collinear across this area!
# There are more quantitative assessments of collinearity (like using the Variance Inlation Factor)
# vifstep() from the usdm package uses VIF to identify the most collinear predictors
# From the vifstep() documentation:
# "vifstep calculate VIF for all variables, exclude one with highest VIF (greater than threshold), 
#repeat the procedure until no variables with VIF greater than th remains."
vifstep(select(sdmdata, -Y), th=10)
```

    ## 1 variables from the 4 input variables have collinearity problem: 
    ##  
    ## MAT 
    ## 
    ## After excluding the collinear variables, the linear correlation coefficients ranges between: 
    ## min correlation ( MAP ~ MTCQ ):  0.2449545 
    ## max correlation ( MTCQ ~ MTWQ ):  0.6056576 
    ## 
    ## ---------- VIFs of the remained variables -------- 
    ##   Variables      VIF
    ## 1      MTWQ 6.146709
    ## 2      MTCQ 4.793498
    ## 3       MAP 4.140401

``` r
# It appears that when MAT is removed. There aren't significant collinearity 
# problems among the predictors

# So let's remove MAT and try again
predictors <- raster::subset(bioclim.stack, c("MTWQ", "MTCQ", "MAP"))
presvals <- raster::extract(predictors, dc.df)
backgr <- randomPoints(predictors, 500)
absvals <- raster::extract(predictors, backgr)
Y <- c(rep(1, nrow(presvals)), rep(0, nrow(absvals)))
sdmdata <- data.frame(cbind(Y, rbind(presvals, absvals)))

pairs(sdmdata[,2:length(sdmdata)], cex = 0.1, fig=TRUE)
```

![](README_files/figure-markdown_github/SDM%20prep-2.png)

``` r
# Better! 
```

Manual (base R) SDM
===================

Let’s do some SDMs using linear models and Generalized Linear Models–
these are very rudimentary and not used in published papers, but I think
these examples are important for understanding how some of the simplest
SDMs (under the family of regression methods) generally work

Let’s start Modeling! Here’s the most basic regression SDM, a linear
model

Linear model
------------

*E*<sub>*Y*</sub> = *m* \* *x* + *b* + *e*

*E*<sub>*Y*</sub> is expected habitat suitability

*b* = intercept

*m* = coefficient

*x* = predictor

*e* = error

Assumes: Normal distributed error mean = 0 Variance(*E*<sub>*Y*</sub>)
is constant

``` r
sdm_lm <- lm(Y ~ MTWQ, data = sdmdata)
plot(sdmdata$MTWQ, sdmdata$Y, 
     main = "D. californica Linear SDM",
     xlab = "Mean Temp of Warmest Quarter",
     ylab = "Presence/Absence"
     )
abline(sdm_lm)
```

![](README_files/figure-markdown_github/Linear%20Model-1.png) The
simplest SDM ever! As Mean Temperature of the Warmest Quarter increases,
the more absent the species is. This relationship can be mapped out in
geographic space, but would not be too informative because the model fit
is poor and we know that other climatic variables impact species
distributions.

Generalized Linear Model SDM
----------------------------

*g*(*E*<sub>*Y*</sub>) = *L**P* = *m**x* + *b* + *e*

*g*() is link function

*L**P* = Linear Predictor

GLMs Do not have to have: - normal distribution of error - homogenous
variance Different link functions can be used, which each have their own
assumptions. A GLM with a Gaussian link function is simply a Linear
Model– it carries the same assumptions as a linear model. We can make
the same linear model as above using the glm() function

``` r
plot(sdmdata$MTWQ, sdmdata$Y, 
     main = "D. californica Gaussian GLM SDM",
     xlab = "Mean Temp of Warmest Quarter",
     ylab = "Presence/Absence"
)

sdm_lm <- glm(Y ~ MTWQ, data = sdmdata, family = gaussian)
abline(sdm_lm, lwd = 3, col = "blue")
```

![](README_files/figure-markdown_github/Gaussian%20GLM-1.png)

The Gaussian link function doesn’t make too much sense for SDM because
SDM has a binary response variable: presence or absence. When trying to
interpret the above model we run into a problem– At MTWQ = 20.0°C, what
does .4 Presence/Absence mean? How can we have a plant .4 present at a
location? This does not make sense, and we must transform our model to
provide the *probability* of presence (arguably better interpreted as
“habitat suitability”) as the response variable. This would give meaning
to values between 1 and 0. We can do this with the logit link function.

*g*(*E*<sub>*Y*</sub>) = *L**P* = *m**x* + *b* + *e*

*L**P* = *l**o**g*10(*E*<sub>*Y*</sub>/(1 − *E*<sub>*Y*</sub>))

solve for *E*<sub>*Y*</sub>

Expanded Logit GLM:
*E*<sub>*Y*</sub> = *e*<sup>*m**x* + *b* + *e*</sup>/(1 + *e*<sup>*m**x* + *b* + *e*</sup>)

``` r
sdm_glm <- glm(Y ~ MTWQ, data = sdmdata, family = binomial)
# This is the expanded logit GLM equation, written as an R function
Ey_calc <- function(x){
  y_int <- as.numeric(sdm_glm$coefficients["(Intercept)"])
  m <- as.numeric(sdm_glm$coefficients["MTWQ"])
  LP = m*x + y_int
  return(exp(LP) / (1 + exp(LP)))
}

MTWQ_sorted <- sort(sdmdata$MTWQ)
Ey <- Ey_calc(MTWQ_sorted)

plot(sdmdata$MTWQ, sdmdata$Y, 
     main = "D. californica Logit GLM SDM",
     xlab = "Mean Temp of Warmest Quarter",
     ylab = "Probability of Presence"
)

lines(MTWQ_sorted, Ey, lwd = 3, col = "red")
```

![](README_files/figure-markdown_github/Logit%20GLM-1.png)

Now let’s project this onto geographic space, taking known values of the
environment from raster files and using our GLM to predict the
‘probability of presence’ in each raster cell

``` r
# Here's a function we'll use to plot SDM projections
project.sdm <- function(prediction, plotName){
  sp::plot(prediction, main = plotName)
  sp::plot(CA_OR.shp, add = T)
  points(dc.df, pch = 16, cex = .2)
  legend("bottomright", legend = "D. californica occ.", pch = 16, cex=.4)
}

# Here's a map of projection of the Logit GLM SDM model. We can use the raster 
# 'predict()' function to produce a raster layer of the predicted values
prediction_glm <- raster::predict(bioclim.stack, sdm_glm)

project.sdm(prediction_glm, "Logit GLM SDM (D. californica), MTWQ only")
```

![](README_files/figure-markdown_github/projection1-1.png) Notice how
the raster legend shows values between -6 and 4. These aren’t
probabilities! This happened because we’re plotting the logit
transformed response variable
*g*(*E*<sub>*Y*</sub>) = *m**x* + *b* + *e* where
*g*(*E*<sub>*Y*</sub>) = *l**o**g*10(*E*<sub>*Y*</sub>/(1 − *E*<sub>*Y*</sub>))
Solve for *E*<sub>*Y*</sub> to get the probability value between 0 and
1:
*E*<sub>*Y*</sub> = *e*<sup>(</sup>*g*(*E*<sub>*Y*</sub>)/(1 + *e*<sup>(</sup>*g*(*E*<sub>*Y*</sub>)))

``` r
prediction_glm.Ey <- exp(prediction_glm) / (1 + exp(prediction_glm))
project.sdm(prediction_glm.Ey, "Logit GLM SDM (D. californica), MTWQ only")
```

![](README_files/figure-markdown_github/projection2-1.png)

Looks like a real SDM! Although you can see that the model predicts high
habitat suitability in areas where *D. californica* wasn’t actually
observed. Maybe if we throw the other climatic predictors in it will
help.

``` r
sdm_glm2 <- glm(Y ~ MTWQ + MTCQ + MAP, data = sdmdata, family = binomial)
prediction_glm2 <- raster::predict(bioclim.stack, sdm_glm2)
prediction_glm2.Ey <- exp(prediction_glm2) / (1+exp(prediction_glm2))
project.sdm(prediction_glm2.Ey, "GLM SDM (D. californica)")
```

![](README_files/figure-markdown_github/logitGLM2-1.png)

Notice how the area where the model predicts high habitat suitability
where there weren’t actual occurrences is much lower. We can still do
better! Let’s begin looking at some of the explicit SDM packages and
functions that are popular in R.

Popular SDM methods
===================

Rather than organize these methods by the 3 general categories of SDM

1.  Profiling Methods: distance and envelope based methods
