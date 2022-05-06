##########################################################
###               PERCENT TREE COVER MASKING          ####
##########################################################
##' Removing areas falling below % tree cover threshold from vegetation ENMs

library(raster)
library(rgdal)
library(rgeos)

##################################
## 1. UPLOADING ENMs and VCF DATA
##################################
# Loading unprocessed and post-processed ENMs

#Continuous ENMs
CH.enm.dir <- '/home/erica/PhD/Mex_VegModel/Analyses/ENM/Outputs/RasPreds/CarsoHuasteco'
CH.enm.path <- file.path(CH.enm.dir, dir(CH.enm.dir, pattern = "25_mtp.tif$"))
CH.enm <- raster::stack(CH.enm.path)

#WBC masked ENMs
CH.wbc.dir <- '/home/erica/PhD/Mex_VegModel/Analyses/SVM/outputs/CarsoHuasteco/WBC'
CH.wbc.path <- file.path(CH.wbc.dir, dir(CH.wbc.dir, pattern = "WBC_mtp.tif$"))
CH.wbc <- raster::stack(CH.wbc.path)

#SVMsp masked ENMs
CH.SVMsp.dir <- '/home/erica/PhD/Mex_VegModel/Analyses/SVM/outputs/CarsoHuasteco/SVMsp'
CH.SVMsp.path <- file.path(CH.SVMsp.dir, dir(CH.SVMsp.dir, pattern = "SVMsp_mtp.tif$"))
CH.SVMsp <- raster::stack(CH.SVMsp.path)

#SVMhyb masked ENMs
CH.SVMhyb.dir <- '/home/erica/PhD/Mex_VegModel/Analyses/SVM/outputs/CarsoHuasteco/SVMhyb'
CH.SVMhyb.path <- file.path(CH.SVMhyb.dir, dir(CH.SVMhyb.dir, pattern = "SVMhyb_mtp.tif$"))
CH.SVMhyb <- raster::stack(CH.SVMhyb.path)

#Importing the 3-year MODIS VCF average raster for the Carso Huasteco
VCF <- raster("data/MODIS_VCF/MODISVCF_2016-2018Avg.tif")

#Loading the Carso Huasteco shapefile matching it to the the VCF raster's spatial reference system
CH <- rgdal::readOGR("/home/erica/PhD/Mex_VegModel/Analyses/ENM/data/CarsoHuasteco.shp", 
                     layer = "CarsoHuasteco")

CHt <- spTransform(CH, CRS(projection(VCF)))

#Cropping the VCF layer to the extent of the Carso Huasteco
VCF <- crop(VCF, CHt)
VCF<- mask(VCF, CHt)

############################################
# 2. PERCENT TREE COVER ENM THRESHOLDING 
############################################

#' @param enm a raster object of enm predictions
#' @param mask raster to mask the enm predictions. Must be same extent and projection.

#In WGS 84 (no projection) for visualization
FC.thresh <- function(enm, mask, out.dir){
  CF <- subset(enm, which(grepl("CF_", names(enm)) == TRUE))
  PO <- subset(enm, which(grepl("PO_", names(enm)) == TRUE))
  
  #subsetting the mask by forest cover threshold
  CFm <- mask >= 60
  CFm[CFm == 0]<- NA
  
  # Reprojecting to match extent & resolution of ENMs
  CFm <- projectRaster(CFm, CF)
  CFm <- crop(CFm, CF)
  
  POm <- mask >= 40
  POm[POm == 0]<- NA
  POm <- projectRaster(POm, PO)
  POm <- crop(POm, PO)
  
  # Masking the vegetation ENMs
  CFmsk <- mask(CF, CFm)
  POmsk <- mask(PO, POm)
  
  #Creating a single raster stack of masked ENMs
  enm.msk <- raster::stack(CFmsk, POmsk)
  
  #Exporting masked ENM raster files 
  writeRaster(enm.msk, paste(out.dir, paste0(names(enm), "_FC"),
                             sep = "/"), bylayer = T, format = "GTiff", overwrite = T) 
  return(enm.msk)
}

#In sinusoidal for areal calculations
FC.thresh2 <- function(enm, mask, out.dir){
  CF <- subset(enm, which(grepl("CF_", names(enm)) == TRUE))
  PO <- subset(enm, which(grepl("PO_", names(enm)) == TRUE))
  
  #Subsetting the mask by percent tree cover threshold
  CFm <- mask >= 60
  CFm[CFm == 0]<- NA
  
  # Reprojecting to match extent & resolution of ENM
  CF <- projectRaster(CF, CFm)
  
  POm <- mask >= 40 
  POm[POm == 0]<- NA
  PO <- projectRaster(PO, POm)
  
  #Masking the vegetation ENMs
  CFmsk <- mask(CF, CFm)
  POmsk <- mask(PO, POm)
  
  #Creating a single raster stack of masked ENMs
  enm.msk <- raster::stack(CFmsk, POmsk)
  
  #Exporting masked ENM raster files   
  writeRaster(enm.msk, paste(out.dir, 
                             paste0(names(enm), "_FC"), sep = "/"), 
              bylayer = T, format = "GTiff", overwrite = T) 
  
  return(enm.msk)
}

############################################
# 3. MASKING ENMs WITH TREE COVER TRHESHOLDS
############################################

#creating a list with all predictions
ENM.list <- list(CH.enm, CH.wbc, CH.SVMsp, CH.SVMhyb)

#applying the forest cover threshold to all pre and post-processed predictions
ENM.th <- lapply(ENM.list, FC.thresh, mask = VCF, out.dir = "outputs/WGS84")

#applying the forest cover threshold to all predictions (in sinusoidal proj)
ENM.th2 <- lapply(ENM.list, FC.thresh2, mask = VCF, out.dir = "outputs/Sinusoidal")
