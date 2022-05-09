##############################################################################
###          RASTERIZING INEGI SOIL LAYERS                                ####
##############################################################################

##' Converting the INEGI Edaphic vector files to rasters matching CHELSA resolution for ENMs


#loading the necessary R packages
library(rgdal)
library(rgeos)
library(raster)

###############################################
# LOADING SOIL AND CLIMATE LAYERS
###############################################
#loading shapefile for Mexico
MEX <- readOGR("D:/PhD_Research/Mex_VegModel/Analyses/CountryData/Mex_shp/gadm36_MEX_0.shp", layer = "gadm36_MEX_0")

#Loading Chelsa layer to use as template for rasterization
chelsa <- raster("data/Old/CHELSA_bio/CHELSA_bio10_05.tif")

#cropping Chelsa layer to Mexico
chelsa <- crop(chelsa, MEX)
chelsa <- mask(chelsa, MEX)

#Loading soil mosaic data
soil <- readOGR("D:/PhD_Research/Mex_VegModel/Analyses/CountryData/MEX_Edaf2007/SoilMos_MEXutm.shp", layer = "SoilMos_MEXutm")

###############################################
# REPROJECTING AND RASTERIZING SOIL DATA
###############################################
#1. Defining the original and desired CRS 
proj1 <- CRS("+proj=utm +zone=14 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs") #Original CRS for soil data ITRF 92 Zone 14N

proj2 <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") #No projection WGS84

#Assiging the original projection to shapefile
crs(soil) <- proj1 

#Reprojecting to WGS 84 to match the chelsa data
soil <- spTransform(soil, proj2)

#2. Rasterizing the soil shapefile using chelsa as a template

#Assigning a number value to each soil type to be able to rasterize  
SoilType <- unique(soil$GRUPO1)
SoilType.df <- data.frame(ID = 1:length(SoilType), SoilType = SoilType)

#Adding dictionary to soil spatial polygons dataframe
soil$ID <- SoilType.df$ID[match(soil$GRUPO1, SoilType.df$SoilType)]

#Rasterizing the soil polygon by soil type using the Chelsa layer as a template
soil_ras <- rasterize(soil, chelsa, field = "ID")

#Reclassifying NA values (now NA = 8) as NA
soil_ras[soil_ras == 8] <- NA

# cropping to remove areas far from the data 
e <- c(--103.0034, -95.58052, 15.34395, 27.34527)
soil_ras <- crop(soil_ras, e)

#Exporting the soil raster for eastern Mexico 
writeRaster(soil_ras, "data/Soil_MEX.tif", format = "GTiff", overwrite = T)

###################################
# CLIPPING TO THE CARSO HUASTECO 
###################################

#importing shapefiles for the Carso Huasteco 
CH <- readOGR("data/CarsoHuasteco_WGS.shp", layer = "CarsoHuasteco_WGS")

#Cropping the soil raster to the study region
CH.soil <- crop(soil_ras, CH)
CH.soil <- mask(CH.soil, CH)

#Exporting the cropped soil raster
writeRaster(CH.soil, "data/Env_CarsoHuasteco/Soil.tif", format = "GTiff", overwrite = T)

