##########################################################
###                      AREA DIFFERENCES             ####
##########################################################
##' Calculating the difference in suitable area for both forest types before and 
##' after applying a percent tree cover mask to pre- and post-processed ENMs. 

#library(raster)
#library(rgdal)
#library(rgeos)
library(dplyr)
library(tidyr)
library(terra)

#####################################
## 1. UPLOADING UNMAKSED/MASKED ENMs
#####################################
# Loading pre and post-processed ENMs

#Continuous ENMs
CH.enm.dir <- '/home/erica/PhD/Mex_VegModel/Analyses/ENM/Outputs/RasPreds/CarsoHuasteco'
CH.enm.path <- file.path(CH.enm.dir, dir(CH.enm.dir, pattern = "25_mtp.tif$"))
CH.enm <- terra::rast(CH.enm.path)[[1:2]]

#WBC masked ENMs
CH.wbc.dir <- '/home/erica/PhD/Mex_VegModel/Analyses/SVM/outputs/CarsoHuasteco/WBC'
CH.wbc.path <- file.path(CH.wbc.dir, dir(CH.wbc.dir, pattern = "WBC_mtp.tif$"))
CH.wbc <- terra::rast(CH.wbc.path)[[1:2]]

#SVMsp masked ENMs
CH.SVMsp.dir <- '/home/erica/PhD/Mex_VegModel/Analyses/SVM/outputs/CarsoHuasteco/SVMsp'
CH.SVMsp.path <- file.path(CH.SVMsp.dir, dir(CH.SVMsp.dir, pattern = "SVMsp_mtp.tif$"))
CH.SVMsp <- terra::rast(CH.SVMsp.path)[[1:2]]

#SVMhyb masked ENMs
CH.SVMhyb.dir <- '/home/erica/PhD/Mex_VegModel/Analyses/SVM/outputs/CarsoHuasteco/SVMhyb'
CH.SVMhyb.path <- file.path(CH.SVMhyb.dir, dir(CH.SVMhyb.dir, pattern = "SVMhyb_mtp.tif$"))
CH.SVMhyb <- terra::rast(CH.SVMhyb.path)[[1:2]]

#Loading the percent tree cover masked ENMs
FC.dir <- "/home/erica/PhD/Mex_VegModel/Analyses/ForestMask/outputs/Sinusoidal/"

FC.enm.path <- file.path(FC.dir, dir(FC.dir, pattern = "25_mtp_FC.tif"))
FC.enm <- terra::rast(FC.enm.path)

FC.wbc.path <- file.path(FC.dir, dir(FC.dir, pattern = "WBC_mtp_FC.tif"))
FC.wbc <- terra::rast(FC.wbc.path)

FC.SVMsp.path <- file.path(FC.dir, dir(FC.dir, pattern = "SVMsp_mtp_FC.tif"))
FC.SVMsp <- terra::rast(FC.SVMsp.path)

FC.SVMhyb.path <- file.path(FC.dir, dir(FC.dir, pattern = "SVMhyb_mtp_FC.tif"))
FC.SVMhyb <- terra::rast(FC.SVMhyb.path)

##################################################
## 2. CHANGING ENM PROJECTIONS TO CALCULATE AREAS
##################################################
#Setting the projection of the "pre-masked" ENMs to match the MODIS VCF reference system 
ENM.list <- list(CH.enm, CH.wbc, CH.SVMsp, CH.SVMhyb) 

ENM.list <- lapply(ENM.list, project, y = FC.enm)

ENM.list <- append(ENM.list, list(FC.enm, FC.wbc, FC.SVMsp, FC.SVMhyb))
  
##################################
## 3. AREA CALCULATIONS
##################################
ENM.area <- lapply(ENM.list, function (x){

  #Names of raster layers
  n <- names(x)

  #Calculate area of each cell in Km2
  cell_area <- terra::cellSize(x, unit = "km")

  #Calculate total area
  a <- sapply(seq_len(nlyr(x)), function(i){
    #Calculate area
    terra::global(!is.na(x[[i]]), 
    "sum", na.rm = T)[,1] * prod(res(test.enm))/1e6
  })
  #Create a dataframe
  b <- data.frame(
    Model.name = n, 
    Area = a
  )

   # Classify model and vegetation type
  b <- b %>%
    mutate(
      Model = case_when(
        grepl("25_mtp", Model.name) ~ "ENM",
        grepl("WBC", Model.name) ~ "WBC",
        grepl("SVMsp", Model.name) ~ "SVMsp",
        grepl("SVMhyb", Model.name) ~ "SVMsp-env"
      ),
      Veg.type = case_when(
        grepl("CF_", Model.name) ~ "Cloud_forest",
        grepl("PO_", Model.name) ~ "Pine Oak_forest",
        grepl("SS_", Model.name) ~ "Submontane_scrubland"
      )
    )

  return(b)
})

# Combine results into a single table
ENM.area.tbl <- ENM.area %>%
  dplyr::bind_rows() %>%
  dplyr::arrange(Veg.type) %>%
  dplyr::mutate(
    Mask = case_when(
      grepl("_FC", Model.name) ~ "After mask",
      TRUE ~ "Before mask"
    )
  ) %>%
  dplyr::select(-Model.name) %>%
  tidyr::pivot_wider(
    names_from = Mask,
    values_from = Area
  ) %>%
  mutate(
    Area.dif = `Before mask` - `After mask`,
    Pct.dif = (Area.dif / `Before mask`) * 100
  )

# Export results
write.csv(
  ENM.area.tbl,
  "/home/erica/PhD/Mex_VegModel/Analyses/ForestMask/outputs/Sinusoidal/ENM_areas.csv",
  row.names = FALSE
)
