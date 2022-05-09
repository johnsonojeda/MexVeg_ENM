##############################################################################
###             Occurrence point forest cover thresholding                 ###
##############################################################################

##' Here we're looking to remove occurrence points falling in areas below/above a forest cover threshold prior to building our veg ENMs

#loading the necessary R packages
library(raster)
library(rgdal)
library(rgeos)

##################################
## 1. UPLOADING ALL THE DATA
##################################
#loading veg type occurrence points drawn from the INEGI land cover polygons
occ.dir <- "data/OccData/FullOccs/CarsoHuasteco"

#Carso Huasteco
CH.occ <- lapply(file.path(occ.dir, dir(occ.dir,pattern = ".csv")), read.csv)

#Loading the Carso Huasteco shapefile
CH <- shapefile("data/CarsoHuasteco.shp")
  
#loading the averaged 2016-2018 MODIS VCF raster 
VCF <- raster("/home/erica/PhD/Mex_VegModel/Analyses/ForestMask/data/MODIS_VCF/MODISVCF_2016-2018Avg.tif")

#Reprojecting the Carso Huasteco shp to match MODIS VCF 
CH<- spTransform(CH, CRS(projection(VCF)))

#Cropping the VCF layer to the extent of the Carso Huasteco
VCF <- crop(VCF, CH)
VCF<- mask(VCF, CH)

#########################################
## 2. EXTRACT RASTER VALUES AT OCC POINTS
#########################################
#' Now we will get  the raster values at each occurrence point per veg type and 
#' remove points falling below specified forest cover threshold for each veg type

#lapply allows to use extract function  to every element in the list

Occs.fc <- lapply(CH.occ,function(x){ 
  #Converting occurrences to spatial points dataframe
  occs.sp <- SpatialPointsDataFrame(x[2:3], x, proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))
  ForestCov<- extract(VCF, occs.sp[2:3]) #cols 2 & 3 of each df represent long/lat
  fcov <- cbind(x, ForestCov) #appending forest cover values to the occs df
  fcov <- if(unique(fcov$scientific_name) == "Cloud_forest"){ #removing points by cover threshold
    subset(fcov, fcov$ForestCov >= 60) #Cloud forest must have >=60% tree cover
  } else if (unique(fcov$scientific_name) == "PineOak_forest"){
    subset(fcov, fcov$ForestCov >= 40)#Pine-oak forest must have >40% tree cover
  } else {
    subset(fcov, fcov$ForestCov < 20) #Scrublands must have <20% tree cover
  }
  fcov <- fcov[1:3]
  return(fcov)
})

#Exporting each of the filtered occurrence points as a unique csv file
out.dir<- "data/OccData/VCF_thresholded/"
write.csv(Occs.fc[[1]], paste0(out.dir, "CarsoHuasteco/CF_Occ_VCF.csv"))
write.csv(Occs.fc[[2]], paste0(out.dir, "CarsoHuasteco/PO_Occ_VCF.csv"))
write.csv(Occs.fc[[3]], paste0(out.dir, "CarsoHuasteco/SS_Occ_VCF.csv"))
  