
library(dlstats)
library(pkgsearch)

# Which are the most popular packages?
sdm_packages <- pkg_search("Species Distribution Model")
sdm_package_names <- sdm_packages$package
x <- cran_stats(sdm_package_names, use_cache = F)

x_recent <- x %>% 
  dplyr::filter(end > as.Date("2019-01-01"))

ggplot(x_recent, aes(end, downloads, group=package, color=package)) +
  geom_line() + geom_point(aes(shape=package))

