################################################################
####                     SVM: SPATIAL                      #####
################################################################

##'Support vector machine cell classification to identify which vegetation type 
##'is more likely to occupy areas where individual vegetation ENM predictions 
##'overlap based on geographic coordinates of occurrence/background points alone

library(rgdal)
library(sp)
library(raster)
library(maskRangeR)

#########################################
# 1. LOADING OCC DATA & ENM PREDICTIONS 
#########################################

#Loading VCF 25% occurrence data for each vegetation type
CH.occ.dir <- '/home/erica/PhD/Mex_VegModel/Analyses/ENM/data/OccData/VCF_RandFilter/CarsoHuasteco'
CH.occ.paths <-file.path(CH.occ.dir, dir(CH.occ.dir, pattern = "VCF25.csv"))

CH.occs<- lapply(CH.occ.paths, read.csv, header = T)

# Loading the optimal continuous ENMs for each veg type
CH.enm.dir <- '/home/erica/PhD/Mex_VegModel/Analyses/ENM/Outputs/RasPreds/CarsoHuasteco'
CH.paths <- file.path(CH.enm.dir, dir(CH.enm.dir, pattern = "VCF25.tif$"))
CH.enm <- raster::stack(CH.paths)

######################################
# 2. SVM SPATIAL CLASSIFIER
######################################

##' Create a spatial SVM classifier layer from occurrence points
##' @param occs.list list occurrence of dataframes with at least 2 columns: 
##' 'latitude' & 'longitude'
##' @param enm RasterStack with ENM predictions
##' @param file.name file path & name of the desired output raster file. 

SVMsp <- function(occs.list, enm, file.name){
  require(maskRangeR)
  
  CF.occs <- occs.list[[1]][intersect(names(occs.list[[1]]), c("longitude", "latitude"))]
  PO.occs <- occs.list[[2]][intersect(names(occs.list[[2]]), c("longitude", "latitude"))]
  SS.occs <- occs.list[[3]][intersect(names(occs.list[[3]]), c("longitude", "latitude"))]
  
  svm <- maskRangeR::rangeSVM(CF.occs, PO.occs, SS.occs, nrep = 100)
  svm.pred <- maskRangeR::rangeSVM_predict(svm = svm, r = enm[[1]])
  
  writeRaster(svm.pred, file.name, format = "GTiff", overwrite = T)
  return(svm.pred)
}

#Creating the SVMsp classifier 
CH.SVMsp.class <- SVMsp(CH.occs, CH.enm, file.name = "outputs/CarsoHuasteco/SVMsp/SVMsp_trin_VCF25.tif")

################################################
# 3. ENM POST-PROCESSING WITH SPATIAL SVM MASK
################################################

##' Post-processing continuous ENM predictions with a spatial SVM mask
##' @param enm RasterStack of ENM predictions used to create the SVMsp classifier
##' @param SVMsp.class spatial SVM classifier layer
##' @param out.dir file directory to save outputs

SVMsp.mask <- function(enm, SVMsp.class, out.dir){
  #Cloud forest
  CF <- subset(enm, which(grepl("CF_", names(enm)) == TRUE))
  PO <- subset(enm, which(grepl("PO_", names(enm)) == TRUE))
  SS <- subset(enm, which(grepl("SS_", names(enm)) == TRUE))
  
  #Cloud forest
  CF.msk <- SVMsp.class
  CF.msk[CF.msk != 1] <- NA
  CF.msk <- mask(CF, CF.msk) 
  
  #Pine-Oak forest
  PO.msk <- SVMsp.class
  PO.msk[PO.msk != 2] <- NA
  PO.msk <- mask(PO, PO.msk)
  
  #Submontane scrubland
  SS.msk <- SVMsp.class
  SS.msk[SS.msk != 3] <- NA
  SS.msk <- mask(SS, SS.msk)
  
  enm.SVMsp <- stack(CF.msk, PO.msk, SS.msk)
  names(enm.SVMsp) <- paste0(names(enm), "_SVMsp")
  
  #Exporting SVMsp-masked ENMs to output directory 
  writeRaster(enm.SVMsp, paste(out.dir, names(enm.SVMsp), sep = "/"), 
              bylayer = T, format = "GTiff", overwrite = T)
  
  return(enm.SVMsp)
}

#Post-processing the vegetation ENMs with the SVMsp mask  
CH.SVMsp <- SVMsp.mask(CH.enm, CH.SVMsp.class, out.dir = "outputs/CarsoHuasteco/SVMsp")
