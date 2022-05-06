################################################################
####                     WINNER BY CELL                    #####
################################################################

##'Winner by cell classification to identify which vegetation type is more likely 
##'to occupy areas where individual vegetation ENM predictions overlap based on
##'ENM suitability values alone

library(raster)
library(rgdal)

######################################
# 1. LOADING ENM PREDICTIONS 
######################################

# Loading the optimal continuous ENMs for each veg type
CH.enm.dir <- '/home/erica/PhD/Mex_VegModel/Analyses/ENM/Outputs/RasPreds/CarsoHuasteco'
CH.paths <- file.path(CH.enm.dir, dir(CH.enm.dir, pattern = "VCF25.tif$"))

#creating a single raster stack with all three vegetation ENMs 
CH.enm <- raster::stack(CH.paths)

######################################
# 2. WBC CLASSIFIER 
######################################

##' Create a Winner-by-cell classifier layer from ENM raster stack 
##' @param enm raster stack with ENM predictions 
##' @param file.name path to which resulting classifier should be exported to 

WBC.class <- function(enm, file.name){
  CF <- subset(enm, which(grepl("CF_", names(enm)) == TRUE))
  PO <- subset(enm, which(grepl("PO_", names(enm)) == TRUE))
  SS <- subset(enm, which(grepl("SS_", names(enm)) == TRUE))
  
  #Cloud forest wins
  CF.msk  <- CF > PO & CF > SS
  CF.msk[CF.msk != 0] <- 1
  
  #Pine-oak forest wins
  PO.msk  <- PO > CF & PO > SS
  PO.msk[PO.msk != 0] <- 2
  
  #Scrubland wins
  SS.msk  <- SS > CF & SS > PO
  SS.msk[SS.msk != 0] <- 3
  
  #creating a raster stack
  WBC.msk <- stack(CF.msk, PO.msk, SS.msk)
  names(WBC.msk) <- paste0(names(enm), "_WBC")
  
  #combining into a single raster
  WBC.trin <- WBC.msk[[1]] + WBC.msk[[2]] + WBC.msk[[3]]
  WBC.trin[WBC.trin == 0]<- NA
  
  #exporting the WBC classifier raster
  writeRaster(WBC.trin, file.name, format = "GTiff", overwrite = T)
  return(WBC.trin)
}

#Creating the WBC classifier 
CH.WBC <- WBC.class(CH.enm, file.name = "outputs/CarsoHuasteco/WBC/WBCtrin_VCF25.tif")

############################################
# 3. ENM POST-PROCESSING WITH WBC MASK
############################################

##' Post-processing continuous ENM predictions with a Winner-by-cell mask
##' @param enm RasterStack of ENM predictions used to create the WBC mask
##' @param wbc.class WBC classifier layer
##' @param out.dir file directory to save outputs 

WBC.mask <- function(enm, wbc.class, out.dir){
  #Cloud forest
  CF <- subset(enm, which(grepl("CF_", names(enm)) == TRUE))
  PO <- subset(enm, which(grepl("PO_", names(enm)) == TRUE))
  SS <- subset(enm, which(grepl("SS_", names(enm)) == TRUE))
  
  #Cloud forest
  CF.msk <- wbc.class
  CF.msk[CF.msk != 1] <- NA
  CF.msk <- mask(CF, CF.msk) 
  
  #Pine-Oak forest
  PO.msk <- wbc.class
  PO.msk[PO.msk != 2] <- NA
  PO.msk <- mask(PO, PO.msk)
  
  #Submontane scrubland
  SS.msk <- wbc.class
  SS.msk[SS.msk != 3] <- NA
  SS.msk <- mask(SS, SS.msk)
  
  enm.wbc <- stack(CF.msk, PO.msk, SS.msk)
  names(enm.wbc) <- paste0(names(enm), "_WBC")
 
  #Exporting WBC-masked ENMs to output directory 
  writeRaster(enm.wbc, paste(out.dir, names(enm.wbc), sep = "/"), 
              bylayer = T, format = "GTiff", overwrite = T)
  
  return(enm.wbc)
}

#Post-processing the vegetation ENMs with the WBC mask 
CH.WBCmsk <- WBC.mask(CH.enm, CH.WBC, "outputs/CarsoHuasteco/WBC")
