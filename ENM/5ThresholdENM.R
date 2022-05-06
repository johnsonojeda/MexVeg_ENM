##############################################################################
###                          ENM THRESHOLDING                             ####
##############################################################################

#' Applying MTP and 10 pct thresholds the continuous ENMs for each vegetation type 

#loading the necessary R packages
library(raster)
library(rgdal)
library(rgeos)
library(dplyr)

#######################################
## 1. UPLOADING ENM raster predictions
#######################################

#specifying file paths for Raster predictions 

#Continuous ENMs
CH.enm.dir <- 'Outputs/RasPreds/CarsoHuasteco'
CH.enm.path <- file.path(CH.enm.dir, dir(CH.enm.dir, pattern = "VCF25.tif$"))
CH.enm <- raster::stack(CH.enm.path)

#WBC masked ENMs
CH.wbc.dir <- '/home/erica/PhD/Mex_VegModel/Analyses/SVM/outputs/CarsoHuasteco/WBC'
CH.wbc.path <- file.path(CH.wbc.dir, dir(CH.wbc.dir, pattern = "WBC.tif$"))
CH.wbc <- raster::stack(CH.wbc.path)

#SVMsp masked ENMs
CH.SVMsp.dir <- '/home/erica/PhD/Mex_VegModel/Analyses/SVM/outputs/CarsoHuasteco/SVMsp'
CH.SVMsp.path <- file.path(CH.SVMsp.dir, dir(CH.SVMsp.dir, pattern = "SVMsp.tif$"))
CH.SVMsp <- raster::stack(CH.SVMsp.path)

#SVMhyb masked ENMs
CH.SVMhyb.dir <- '/home/erica/PhD/Mex_VegModel/Analyses/SVM/outputs/CarsoHuasteco/SVMhyb'
CH.SVMhyb.path <- file.path(CH.SVMhyb.dir, dir(CH.SVMhyb.dir, pattern = "SVMhyb.tif$"))
CH.SVMhyb <- raster::stack(CH.SVMhyb.path)

######################################
# 2. CREATING A TRHESHOLDING FUNCTION
######################################
#loading file with all calculated MTP and 10 pct thresholds per veg type/region

thresh <- as_tibble(read.csv("Outputs/ENMThresholds.csv"))

#Reclassifying so values under the calculated threshold for each veg type is NA
##' @param enm a raster object of enm predictions 
##' @param th.tbl a dataframe with mtp & 10pct threshold values for a region/veg type


ENM.thresh <- function(enm, th.tbl, out.dir) {
  CF <- subset(enm, which(grepl("CF_", names(enm)) == TRUE))
  PO <- subset(enm, which(grepl("PO_", names(enm)) == TRUE))
  SS <- subset(enm, which(grepl("SS_", names(enm)) == TRUE))
  
  mtp <- as_tibble(th.tbl) %>% select(Veg.type, MTP)
  
  pct10 <- as_tibble(th.tbl) %>% select(Veg.type, pct10)
  
  #MTP thresholds
  mtp.CF <- as.numeric(subset(mtp, Veg.type == "Cloud forest")[,2])
  mtp.PO <- as.numeric(subset(mtp, Veg.type == "PineOak forest")[,2])
  mtp.SS <- as.numeric(subset(mtp, Veg.type == "Submontane scrub")[,2])
    
  #10 Pct thresholds
  pct10.CF <- as.numeric(subset(pct10, Veg.type == "Cloud forest")[,2])
  pct10.PO <- as.numeric(subset(pct10, Veg.type == "PineOak forest")[,2])
  pct10.SS <- as.numeric(subset(pct10, Veg.type == "Submontane scrub")[,2])  
  
 #Thresholding the predictions
 # Mininum training presence 
  CF.mtp <- CF
  CF.mtp[CF.mtp < mtp.CF] <- NA
  
  PO.mtp <- PO
  PO.mtp[PO.mtp < mtp.PO] <- NA
  
  SS.mtp <- SS
  SS.mtp[SS.mtp < mtp.SS] <- NA

  # 10 percentile 
  CF.10pct <- CF
  CF.10pct[CF.10pct < pct10.CF] <- NA
  
  PO.10pct <- PO
  PO.10pct[PO.10pct < pct10.PO] <- NA
  
  SS.10pct <- SS
  SS.10pct[SS.10pct < pct10.SS] <- NA
  
  enm.mtp <- stack(CF.mtp, PO.mtp, SS.mtp)
  names(enm.mtp) <- paste0(names(enm.mtp), "_mtp")
  
  enm.10pct <- stack(CF.10pct, PO.10pct, SS.10pct)
  names(enm.10pct) <- paste0(names(enm.10pct), "_10pct")
  
  writeRaster(enm.mtp, paste0(out.dir, names(enm.mtp)), bylayer = T, format = "GTiff", 
                            overwrite = T) 
  writeRaster(enm.10pct, paste0(out.dir, names(enm.10pct)),bylayer = T, format = "GTiff", 
              overwrite = T) 
  
  return(list(enm.mtp, enm.10pct))
  
}

######################################
# 3. THRESHOLDING  ENMs
######################################

#Obtaining the thresholded rasters for the continuous ENMs 
CH.th<- ENM.thresh(CH.enm, th.tbl = thresh, out.dir = "Outputs/RasPreds/CarsoHuasteco/")

#Obtaining the thresholded rasters for the WBC-masked ENMs
CH.wbc.th <- ENM.thresh(CH.wbc, th.tbl = thresh, , 
                        out.dir = "/home/erica/PhD/Mex_VegModel/Analyses/SVM/outputs/CarsoHuasteco/WBC/")

#Obtaining the thresholded rasters for the SVMPsp-masked ENMs
CH.SVMsp.th <- ENM.thresh(CH.SVMsp, th.tbl = thresh, 
                        out.dir = "/home/erica/PhD/Mex_VegModel/Analyses/SVM/outputs/CarsoHuasteco/SVMsp/")
#Obtaining the thresholded rasters for the SVMPhyb-masked ENMs
CH.SVMhyb.th <- ENM.thresh(CH.SVMhyb, th.tbl = thresh, 
                          out.dir = "/home/erica/PhD/Mex_VegModel/Analyses/SVM/outputs/CarsoHuasteco/SVMhyb/")
