################################################################
####                     SVM: HYBRID                       #####
################################################################

##'Support vector machine cell classification to identify which vegetation type 
##'is more likely to occupy areas where individual vegetation ENM predictions 
##'overlap based on ENM suitability values & geographic coordinates



library(rgdal)
library(sp)
library(raster)
library(maskRangeR)

#########################################
# 1. LOADING OCC DATA & ENM PREDICTIONS 
#########################################

#Loading VCF 25% occurrence data for each vegetation types
CH.occ.dir <- '/home/erica/PhD/Mex_VegModel/Analyses/ENM/data/OccData/VCF_RandFilter/CarsoHuasteco'
CH.occ.paths <-file.path(CH.occ.dir, dir(CH.occ.dir, pattern = "VCF25.csv"))
CH.occs<- lapply(CH.occ.paths, read.csv, header = T)


#Loading the continuous ENM predictions
CH.enm.dir <- '/home/erica/PhD/Mex_VegModel/Analyses/ENM/Outputs/RasPreds/CarsoHuasteco'
CH.paths <- file.path(CH.enm.dir, dir(CH.enm.dir, pattern = "VCF25.tif$"))
CH.enm <- raster::stack(CH.paths)

#########################################
# 2. SVM SPATIAL-ENVIRONMENTAL CLASSIFIER
#########################################
##' Create a spatial-environmental SVM classifier layer from occurrence points and 
##' ENM suitability values
##' @param occs.list list occurrence of dataframes with at least 2 columns: 
##' 'latitude' & 'longitude'
##' @param enm RasterStack with ENM predictions
##' @param file.name file path & name of the desired output raster file. 

SVMhyb <- function(occs.list, enm, file.name){
  require(maskRangeR)
  
  CF.occs <- occs.list[[1]][intersect(names(occs.list[[1]]), c("longitude", "latitude"))]
  PO.occs <- occs.list[[2]][intersect(names(occs.list[[2]]), c("longitude", "latitude"))]
  SS.occs <- occs.list[[3]][intersect(names(occs.list[[3]]), c("longitude", "latitude"))]
  
  svm <- maskRangeR::rangeSVM(CF.occs, PO.occs, SS.occs, sdm = enm, nrep = 100)
  svm.pred <- maskRangeR::rangeSVM_predict(svm = svm, r = enm[[1]], sdm = enm)
  
  writeRaster(svm.pred, file.name, format = "GTiff", overwrite = T)
  return(svm.pred)
}

#Creating the SVMsp classifier
CH.SVMhyb.class <- SVMhyb(CH.occs, CH.enm, file.name = "outputs/CarsoHuasteco/SVMhyb/SVMhyb_trin_VCF25.tif")

#############################################################
# 3. ENM POST-PROCESSING WITH SPATIAL-ENVIRONMENTAL SVM MASK
#############################################################
##' Masking continuous ENM predictions with a spatial-environmental SVM classifier
##' @param enm RasterStack of ENM predictions used to create the SVMhyb classifier
##' @param SVMhyb.class spatial-environmental SVM classifier layer
##' @param out.dir file directory to save outputs

SVMhyb.mask <- function(enm, SVMhyb.class, out.dir){
  #Cloud forest
  CF <- subset(enm, which(grepl("CF_", names(enm)) == TRUE))
  PO <- subset(enm, which(grepl("PO_", names(enm)) == TRUE))
  SS <- subset(enm, which(grepl("SS_", names(enm)) == TRUE))
  
  #Cloud forest
  CF.msk <- SVMhyb.class
  CF.msk[CF.msk != 1] <- NA
  CF.msk <- mask(CF, CF.msk) 
  
  #Pine-Oak forest
  PO.msk <- SVMhyb.class
  PO.msk[PO.msk != 2] <- NA
  PO.msk <- mask(PO, PO.msk)
  
  #Submontane scrubland
  SS.msk <- SVMhyb.class
  SS.msk[SS.msk != 3] <- NA
  SS.msk <- mask(SS, SS.msk)
  
  enm.SVMhyb <- stack(CF.msk, PO.msk, SS.msk)
  names(enm.SVMhyb) <- paste0(names(enm), "_SVMhyb")
  
  #Exporting SVMsp-env masked ENMs to output directory 
  writeRaster(enm.SVMhyb, paste(out.dir, names(enm.SVMhyb), sep = "/"), 
              bylayer = T, format = "GTiff", overwrite = T)
  
  return(enm.SVMhyb)
}

#Post-processing the vegetation ENMs with the SVMsp-env mask 
CH.SVMhyb <- SVMhyb.mask(CH.enm, CH.SVMhyb.class, out.dir = "outputs/CarsoHuasteco/SVMhyb")
