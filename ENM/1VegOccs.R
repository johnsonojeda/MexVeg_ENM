##############################################################################
###           Create Occurrence points from land cover polygons           ####
##############################################################################

##' This script is to generate occurrence records for vegetation type based on the INEGI land cover polygons 

#Loading the necessary packages
library(raster)
library(rgdal)
library(rgeos)
library(dplyr)

##########################################
#STEP 1: loading all the data 
##########################################

#Specifying the directories where the vegetation shapefiles and Chelsa rasters are located. 
veg.dir <- veg.dir <- "data/Veg/Veg_CarsoHuasteco"
env.dir <- "data/Env/Env_CarsoHuasteco"

#creating lists of files for
veg.list <- file.path(veg.dir, dir(veg.dir, pattern = ".shp$"))
env.list <- file.path(env.dir, dir(env.dir, pattern = ".tif$"))

#reading in vegetation shapefiles
veg.shp <- lapply(veg.list, function(x){
  a<- raster::shapefile(x)
  return(a)
})

#loading the environmental data 
env <- raster::stack(env.list)

##############################################
#STEP 2: Creating occ points for each veg type
##############################################

# Rasterizing the veg polygons using Chelsa bioclim layers as a template
veg.ras <- raster::stack(unlist(lapply(veg.shp, 
                                       raster::rasterize, env[[1]])))
names(veg.ras) <- c("CF", "PO", "SS")

# Now converting the veg rasters to spatial points for my occurrences
veg.occ <- raster::rasterToPoints(veg.ras, spatial = T)

# Converting the points to a spatial points dataframe with lat/long columns
veg.occ.xy <- tbl_df(cbind(veg.occ@coords, veg.occ@data)) %>% 
  rename(scientific_name = name, longitude = x, latitude = y) 

#Exporting occ data for each veg type individually
CF.occ <- veg.occ.xy %>% filter(!is.na(CF)) %>% 
  mutate(scientific_name = rep("Cloud_forest", nrow(.))) %>%
  select(scientific_name, longitude, latitude)

PO.occ <- veg.occ.xy %>% filter(!is.na(PO)) %>% 
  mutate(scientific_name = rep("PineOak_forest", nrow(.))) %>%
  select(scientific_name, longitude, latitude)


SS.occ <- veg.occ.xy %>% filter(!is.na(SS)) %>% 
  mutate(scientific_name = rep("Submontane_scrub", nrow(.))) %>%
  select(scientific_name, longitude, latitude)


write.csv(CF.occ, "data/OccData/FullOccs/CarsoHuasteco/CF_Occ_CH.csv")
write.csv(PO.occ, "data/OccData/FullOccs/CarsoHuasteco/PO_Occ_CH.csv")
write.csv(SS.occ, "data/OccData/FullOccs/CarsoHuasteco/SS_Occ_CH.csv")
