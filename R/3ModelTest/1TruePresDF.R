##############################################################################
###                VEGETATION PRESENCE/ABSENCE POINTS                      ###
##############################################################################
##' Generating separate "true" presence and absence dataframes from the 75% 
##' withheld occurrence data for each vegetation type for model testing

library(dplyr)

##################################
## 1. LOADING OCC DATA           
##################################

#Loading the VCF-100% and VCF-25% occurrence points
VCF100.dir <- "/home/erica/PhD/Mex_VegModel/Analyses/ENM/data/OccData/VCF_thresholded/CarsoHuasteco"

VCF25.dir <- "/home/erica/PhD/Mex_VegModel/Analyses/ENM/data/OccData/VCF_RandFilter/CarsoHuasteco"

VCF100.occ <- lapply(paste(file.path(VCF100.dir), 
                           dir(VCF100.dir, pattern = "VCF.csv"), sep = "/"),
                           read.csv)

VCF25.occ <- lapply(paste(file.path(VCF25.dir), 
                           dir(VCF25.dir, pattern = "VCF25.csv"), sep = "/"),
                     read.csv)

##############################
## 2. OBTAIN TRUE PRESENCES
##############################
#obtaining the 75% of the VCF occurrences that were not used to build the ENMs
#these will constitute true presences for each veg type in model testing
val.pres <- mapply(anti_join, VCF100.occ, VCF25.occ, SIMPLIFY = FALSE)

#Exporting the resulting dataframes 
out.dir <- "data/"

lapply(val.pres, function(x){
  name = unique(x$scientific_name)
  df <- write.csv(x, file = paste0(out.dir, name, "_pres.csv"))
})


