##########################################################
###               PERCENT TREE COVER MASKING          ####
##########################################################
##' Removing areas falling below % tree cover threshold from vegetation ENMs

library(terra)

########################################################
## 1. UPLOADING ENMs, OTHER LAND COVER, & VCF DATA
########################################################
# Loading unprocessed and post-processed ENMs

#Continuous ENMs
CH.enm.dir <- '/home/erica/PhD/Mex_VegModel/Analyses/ENM/Outputs/RasPreds/CarsoHuasteco'
CH.enm.path <- file.path(CH.enm.dir, dir(CH.enm.dir, pattern = "25_mtp.tif$"))
CH.enm <- terra::rast(CH.enm.path)

#WBC masked ENMs
CH.wbc.dir <- '/home/erica/PhD/Mex_VegModel/Analyses/SVM/outputs/CarsoHuasteco/WBC'
CH.wbc.path <- file.path(CH.wbc.dir, dir(CH.wbc.dir, pattern = "WBC_mtp.tif$"))
CH.wbc <- terra::rast(CH.wbc.path)

#SVMsp masked ENMs
CH.SVMsp.dir <- '/home/erica/PhD/Mex_VegModel/Analyses/SVM/outputs/CarsoHuasteco/SVMsp'
CH.SVMsp.path <- file.path(CH.SVMsp.dir, dir(CH.SVMsp.dir, pattern = "SVMsp_mtp.tif$"))
CH.SVMsp <- terra::rast(CH.SVMsp.path)

#SVMhyb masked ENMs
CH.SVMhyb.dir <- '/home/erica/PhD/Mex_VegModel/Analyses/SVM/outputs/CarsoHuasteco/SVMhyb'
CH.SVMhyb.path <- file.path(CH.SVMhyb.dir, dir(CH.SVMhyb.dir, pattern = "SVMhyb_mtp.tif$"))
CH.SVMhyb <- terra::rast(CH.SVMhyb.path)

#Importing the 3-year MODIS VCF average raster for the Carso Huasteco
VCF <- terra::rast("/home/erica/PhD/Mex_VegModel/Analyses/ForestMask/data/MODIS_VCF/MODISVCF_2016-2018Avg.tif")

#Loading the Carso Huasteco shapefile matching it to the the VCF raster's spatial reference system
CH <- terra::vect("/home/erica/PhD/Mex_VegModel/Analyses/ENM/data/CarsoHuasteco.shp")

CHt <- terra::project(CH, crs(VCF))

#Cropping the VCF layer to the extent of the Carso Huasteco
VCF <- terra::crop(VCF, CHt)
VCF<- terra::mask(VCF, CHt)

############################################
# 3. PERCENT TREE COVER ENM THRESHOLDING 
############################################

#' @param enm a raster object of enm predictions
#' @param mask raster to mask the enm predictions. Must be same extent and projection.

#In WGS 84 (no projection) for visualization
FC.thresh <- function(enm, mask, out.dir){
  CF <- enm[[grep("^CF_", names(enm))]]
  PO <- enm[[grep("^PO_", names(enm))]]
  
  #Cloud forest threshold(>=60% tree cover)
  CFm <- mask >= 60
  CFm[CFm == 0]<- NA
  
  # Reprojecting to match extent & resolution of ENMs
  CFm <- terra::project(CFm, CF, method = "bilinear")
  CFm <- terra::crop(CFm, CF)
  
  #Pine-oak forest threshold (=>40% Forest cover <=60%)
  POm <- mask >= 40
  POm[POm == 0]<- NA
  POm <- terra::project(POm, PO, method = "bilinear")
  POm <- terra::crop(POm, PO)
  
  # Masking the vegetation ENMs
  CFmsk <- terra::mask(CF, CFm)
  POmsk <- terra::mask(PO, POm)
  
  #Export maksed ENM raster files
  terra::writeRaster(CFmsk, file.path(out.dir, paste0(names(CFmsk), "_FC.tif")))

  terra::writeRaster(POmsk, file.path(out.dir, paste0(names(POmsk), "_FC.tif")))
  
  #Creating a single raster stack of masked ENMs
  enm.msk <- c(CFmsk, POmsk)
  return(enm.msk)
}

#In sinusoidal for areal calculations
FC.thresh2 <- function(enm, mask, out.dir){
  CF <- enm[[grep("^CF_", names(enm))]]
  PO <- enm[[grep("^PO_", names(enm))]]
  
  #Subsetting the mask by percent tree cover threshold
  CFm <- mask >= 60
  CFm[CFm == 0]<- NA
  
  # Reprojecting to match extent & resolution of ENM
  CF <- terra::project(CF, CFm)
  
  POm <- mask >= 40 
  POm[POm == 0]<- NA
  PO <- terra::project(PO, POm)
  
  #Masking the vegetation ENMs
  CFmsk <- terra::mask(CF, CFm)
  POmsk <- terra::mask(PO, POm)
  
  #Export maksed ENM raster files
  terra::writeRaster(CFmsk, file.path(out.dir, paste0(names(CFmsk), "_FC.tif")))

  terra::writeRaster(POmsk, file.path(out.dir, paste0(names(POmsk), "_FC.tif")))
  
  #Creating a single raster stack of masked ENMs
  enm.msk <- c(CFmsk, POmsk)
  
  return(enm.msk)
}

############################################
# 3. MASKING ENMs WITH TREE COVER TRHESHOLDS
############################################

#creating a list with all predictions
ENM.list <- list(CH.enm, CH.wbc, CH.SVMsp, CH.SVMhyb)

#applying the forest cover threshold to all pre and post-processed predictions
ENM.th <- lapply(ENM.list, FC.thresh, mask = VCF, out.dir = "/home/erica/PhD/Mex_VegModel/Analyses/ForestMask/outputs/WGS84")


#applying the forest cover threshold to all predictions (in sinusoidal proj)
ENM.th2 <- lapply(ENM.list, FC.thresh2, mask = VCF, out.dir = "/home/erica/PhD/Mex_VegModel/Analyses/ForestMask/outputs/Sinusoidal")
