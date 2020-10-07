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

Prepare occurrence data
=======================

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
usa.shp <- readOGR("USA")
```

    ## OGR data source with driver: ESRI Shapefile 
    ## Source: "/home/avery/SDM_Demo/USA", layer: "USA"
    ## with 400 features
    ## It has 5 fields

``` r
CA_OR.shp <- usa.shp[usa.shp@data$NAME %in% c("California", "Oregon"),]

# Let's take a look at our raw occurrence points
sp::plot(CA_OR.shp, main = "Darlingtonia californica (iNaturalist)")
points(dc.df$lon, dc.df$lat)
```

![](README_files/figure-markdown_github/occurrence%20data-1.png)
