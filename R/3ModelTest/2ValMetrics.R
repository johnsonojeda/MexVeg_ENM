##############################################################################
###                           ENM MODEL TESTING                            ###
##############################################################################
##' Calculating model accuracy metrics based on 75% withheld occurrence points
library(dplyr)
library(rgdal)
library(raster)
library(rgeos)

##################################
## 1. LOADING OCC DATA           
##################################
#Loading the 75% withheld occurrences for each veg type into a single dataframe
occ.dir <- "data"
occ.path <- file.path(occ.dir, dir(occ.dir, pattern = ".csv"))

occs <- lapply(occ.path, read.csv)
occs <- occs %>% bind_rows()
##################################
## 2. LOADING ENMs       
##################################

#Continuous ENMs
CH.enm.dir <- '/home/erica/PhD/Mex_VegModel/Analyses/ENM/Outputs/RasPreds/CarsoHuasteco'
CH.enm.path <- file.path(CH.enm.dir, dir(CH.enm.dir, pattern = "VCF25_mtp.tif$"))
CH.enm <- raster::stack(CH.enm.path)

#WBC ENMs
CH.wbc.dir <- '/home/erica/PhD/Mex_VegModel/Analyses/SVM/outputs/CarsoHuasteco/WBC'
CH.wbc.path <- file.path(CH.wbc.dir, dir(CH.wbc.dir, pattern = "WBC_mtp.tif$"))
CH.wbc <- raster::stack(CH.wbc.path)

#SVMsp ENMs
CH.SVMsp.dir <- '/home/erica/PhD/Mex_VegModel/Analyses/SVM/outputs/CarsoHuasteco/SVMsp'
CH.SVMsp.path <- file.path(CH.SVMsp.dir, dir(CH.SVMsp.dir, pattern = "SVMsp_mtp.tif$"))
CH.SVMsp <- raster::stack(CH.SVMsp.path)

#SVMhyb ENMs
CH.SVMhyb.dir <- '/home/erica/PhD/Mex_VegModel/Analyses/SVM/outputs/CarsoHuasteco/SVMhyb'
CH.SVMhyb.path <- file.path(CH.SVMhyb.dir, dir(CH.SVMhyb.dir, pattern = "SVMhyb_mtp.tif$"))
CH.SVMhyb <- raster::stack(CH.SVMhyb.path)

ENM.list <- list(CH.enm, CH.wbc, CH.SVMsp, CH.SVMhyb)

##########################
## 3. MODEL ERROR RATES
###########################
#' Calculating omission, commission, and overall misclassification error rates 
#' for each vegetation ENM based on withheld occurrences

#' @param enm. Raster object of ENM output
#' @param occs Data frame of occurrence records

OmCom.enm <- function(enm, occs){
  # Defining presences and absences for each veg type
  # Note: true presences of one veg type = true absences of the other two
  occs <- as_tibble(occs)
  CF.pres <- occs %>% filter(scientific_name == "Cloud_forest")
  PO.pres <- occs %>% filter(scientific_name == "PineOak_forest")
  SS.pres <- occs %>% filter(scientific_name == "Submontane_scrub")
  
  CF.abs <- occs %>% filter(scientific_name != "Cloud_forest")
  PO.abs <- occs %>% filter(scientific_name != "PineOak_forest")
  SS.abs <- occs %>% filter(scientific_name != "Submontane_scrub")

  #Determining which ENM raster layer belongs to which veg type
  CF <- subset(enm, which(grepl("CF_", names(enm)) == TRUE))
  PO <- subset(enm, which(grepl("PO_", names(enm)) == TRUE))
  SS <- subset(enm, which(grepl("SS_", names(enm)) == TRUE))
  
  #Calculating omission rates
  CF.pres <- raster::extract(CF, CF.pres[c('longitude', 'latitude')])
  CF.OR <- (sum(is.na(CF.pres)))/length(CF.pres)
  
  PO.pres <- raster::extract(PO, PO.pres[c('longitude', 'latitude')])
  PO.OR <- (sum(is.na(PO.pres)))/length(PO.pres)
  
  SS.pres <- raster::extract(SS, SS.pres[c('longitude', 'latitude')])
  SS.OR <- (sum(is.na(SS.pres)))/length(SS.pres)
  
  #Calculating commission rates
  CF.abs <- raster::extract(CF, CF.abs[c('longitude', 'latitude')])
  CF.CR <- (sum(!is.na(CF.abs)))/length(CF.abs)
  
  PO.abs <- raster::extract(PO, PO.abs[c('longitude', 'latitude')])
  PO.CR <- (sum(!is.na(PO.abs)))/length(PO.abs)
  
  SS.abs <- raster::extract(SS, SS.abs[c('longitude', 'latitude')])
  SS.CR <- (sum(!is.na(SS.abs)))/length(SS.abs)
  
  #Calculating overall classification error rate
  CF.ER <- (sum(is.na(CF.pres)) + sum(!is.na(CF.abs)))/(length(CF.pres) + length(CF.abs))            
  
  PO.ER <- (sum(is.na(PO.pres)) + sum(!is.na(PO.abs)))/(length(PO.pres) + length(PO.abs))  
  
  SS.ER <- (sum(is.na(SS.pres)) + sum(!is.na(SS.abs)))/(length(SS.pres) + length(SS.abs))            
  
  Om.rate <- c(CF.OR, PO.OR, SS.OR)
  Com.rate <- c(CF.CR, PO.CR, SS.CR)
  Overall.error <- c(CF.ER, PO.ER, SS.ER)
  Model.name <- names(enm)
  Veg.type <- unique(occs$scientific_name)
  
  results <- data.frame(Model.name, Veg.type, Om.rate, Com.rate, Overall.error)
  
  return(results)
}

#Calculating model error rates for each vegetation type
Veg.OmCom <- lapply(ENM.list, OmCom.enm, occs = occs)

#Concatenating into single dataframe and exporting outpu file
Veg.OmCom <- Veg.OmCom %>% bind_rows()%>% arrange(Veg.type)

write.csv(Veg.OmCom, "outputs/ModelValidation_INEGI.csv")

