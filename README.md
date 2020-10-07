An overview of Species Distribution Modeling in R
=================================================

### This demo was largely based on [Hijmans & Elith, 2013](https://cran.r-project.org/web/packages/dismo/vignettes/sdm.pdf) and the [sdm](https://cran.r-project.org/web/packages/sdm/sdm.pdf) and [dismo](https://cran.r-project.org/web/packages/dismo/dismo.pdf) documentation pages. This serves to provide a walkthrough of SDM and popular implementations in R as of Fall 2020. I’ve included the extent of my understanding of how the different methods work.

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
```

Prepare species occurrence data
===============================

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
    ## 1 -123.6607 42.23318
    ## 2 -124.1054 43.78552
    ## 3 -122.9737 40.97181
    ## 4 -124.0976 44.04676
    ## 5 -123.9310 41.78430
    ## 6 -111.6500 40.24918

``` r
# Let's ready some political boundaries in for visualization
usa.shp <- readOGR("USA", verbose=F)
CA_OR.shp <- usa.shp[usa.shp@data$NAME %in% c("California", "Oregon"),]

# Let's take a look at our raw occurrence points
sp::plot(CA_OR.shp, main = "Darlingtonia californica (iNaturalist)")
points(dc.df$lon, dc.df$lat)
```

![](README_files/figure-markdown_github/occurrence%20data-1.png)

Prepare environmental predictor data
====================================

``` r
# Download bioclimatic variables from worldclim
# run ??getData to see other options for climate data
bioclim_global.stack <- raster::getData(name = "worldclim",
                        var = 'bio', 
                        res = 2.5)

# This replaces bioX name with more descriptive acronym
# We choose 4 climatic variables here more for examplary purposes than biologically
# sound reasoning
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

Build SDM-ready dataframe of predictor values and species presence/(pseudo)absence
==================================================================================

``` r
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
head(sdmdata)
```

    ##   Ey MAT MTWQ MTCQ  MAP
    ## 1  1 111  185   44 1378
    ## 2  1 111  156   70 1819
    ## 3  1  52  127   -6 1180
    ## 4  1 111  156   70 2003
    ## 5  1  99  158   48 1793
    ## 6  1  NA   NA   NA   NA

``` r
# Collinearity can cause problems in Species Distribution Models. My understanding
# is that if two or more predictors are collinear across the environmental space, then 
# it's difficult to determine which predictor is actually influencing the distribution
# of the species (if not both). If the predictors are not collinear in the areas
# to which the SDM is projected then the model won't know how to assign probability
# of habitat suitability to the independent predictors
# We can do this by visually inspecting collinearity with a scatterplot matrix
pairs(sdmdata[,2:length(sdmdata)], cex = 0.1, fig=TRUE)
```

![](README_files/figure-markdown_github/SDM%20prep-1.png)

``` r
# MAT and MTWQ appear to be collinear across this area!
# There are more quantitative assessments of collinearity (like using the Variance Inlation Factor, see https://www.rdocumentation.org/packages/car/versions/3.0-9/topics/vif) that I won't do here

# So let's remove MAT and try again
predictors <- raster::subset(bioclim.stack, c("MTWQ", "MTCQ", "MAP"))
presvals <- raster::extract(predictors, dc.df)
backgr <- randomPoints(predictors, 500)
absvals <- raster::extract(predictors, backgr)
Ey <- c(rep(1, nrow(presvals)), rep(0, nrow(absvals)))
sdmdata <- data.frame(cbind(Ey, rbind(presvals, absvals)))

pairs(sdmdata[,2:length(sdmdata)], cex = 0.1, fig=TRUE)
```

![](README_files/figure-markdown_github/SDM%20prep-2.png)

``` r
# Better! 
```
