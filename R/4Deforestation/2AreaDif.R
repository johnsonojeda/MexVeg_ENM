##########################################################
###                      AREA DIFFERENCES             ####
##########################################################
##' Calculating the difference in suitable area for both forest types before and 
##' after applying a percent tree cover mask to pre- and post-processed ENMs. 

library(raster)
library(rgdal)
library(rgeos)
library(dplyr)
library(tidyr)

#####################################
## 1. UPLOADING UNMAKSED/MASKED ENMs
#####################################
# Loading pre and post-processed ENMs

#Continuous ENMs
CH.enm.dir <- '/home/erica/PhD/Mex_VegModel/Analyses/ENM/Outputs/RasPreds/CarsoHuasteco'
CH.enm.path <- file.path(CH.enm.dir, dir(CH.enm.dir, pattern = "25_mtp.tif$"))
CH.enm <- raster::stack(CH.enm.path)[[1:2]]

#WBC masked ENMs
CH.wbc.dir <- '/home/erica/PhD/Mex_VegModel/Analyses/SVM/outputs/CarsoHuasteco/WBC'
CH.wbc.path <- file.path(CH.wbc.dir, dir(CH.wbc.dir, pattern = "WBC_mtp.tif$"))
CH.wbc <- raster::stack(CH.wbc.path)[[1:2]]

#SVMsp masked ENMs
CH.SVMsp.dir <- '/home/erica/PhD/Mex_VegModel/Analyses/SVM/outputs/CarsoHuasteco/SVMsp'
CH.SVMsp.path <- file.path(CH.SVMsp.dir, dir(CH.SVMsp.dir, pattern = "SVMsp_mtp.tif$"))
CH.SVMsp <- raster::stack(CH.SVMsp.path)[[1:2]]

#SVMhyb masked ENMs
CH.SVMhyb.dir <- '/home/erica/PhD/Mex_VegModel/Analyses/SVM/outputs/CarsoHuasteco/SVMhyb'
CH.SVMhyb.path <- file.path(CH.SVMhyb.dir, dir(CH.SVMhyb.dir, pattern = "SVMhyb_mtp.tif$"))
CH.SVMhyb <- raster::stack(CH.SVMhyb.path)[[1:2]]

#Loading the percent tree cover masked ENMs
FC.dir <- "outputs/Sinusoidal"

FC.enm.path <- file.path(FC.dir, dir(FC.dir, pattern = "25_mtp_FC.tif"))
FC.enm <- raster::stack(FC.enm.path)

FC.wbc.path <- file.path(FC.dir, dir(FC.dir, pattern = "WBC_mtp_FC.tif"))
FC.wbc <- raster::stack(FC.wbc.path)

FC.SVMsp.path <- file.path(FC.dir, dir(FC.dir, pattern = "SVMsp_mtp_FC.tif"))
FC.SVMsp <- raster::stack(FC.SVMsp.path)

FC.SVMhyb.path <- file.path(FC.dir, dir(FC.dir, pattern = "SVMhyb_mtp_FC.tif"))
FC.SVMhyb <- raster::stack(FC.SVMhyb.path)

##################################################
## 2. CHANGING ENM PROJECTIONS TO CALCULATE AREAS
##################################################
#Setting the projection of the "pre-masked" ENMs to match the MODIS VCF reference system 
ENM.list <- list(CH.enm, CH.wbc, CH.SVMsp, CH.SVMhyb) 

ENM.list <- lapply(ENM.list, projectRaster, to = FC.enm)

ENM.list <- append(ENM.list, list(FC.enm, FC.wbc, FC.SVMsp, FC.SVMhyb))
  
##################################
## 3. AREA CALCULATIONS
##################################
#Calculating the avaiable suitable area for each model per forest type
ENM.area <- lapply(ENM.list, function (x){
  require(dplyr)
  n <- names(x)
  a <- lapply(unstack(x) , function(y){
    (ncell(y[!is.na(y)])* res(y)[1] * res(y)[2])/(1000^2) #area in km2    
  }) 
  b <- data.frame(n, unlist(a))
  names(b) <- c("Model.name", "Area")
  
  b <- dplyr::as_tibble(b) %>% 
    mutate(Model = case_when(grepl("25_mtp", Model.name)~ "ENM", 
                             grepl("WBC", Model.name) ~ "WBC",
                             grepl("SVMsp", Model.name)~ "SVMsp",
                             grepl("SVMhyb", Model.name)~ "SVMsp-env")) %>%
    mutate(Veg.type = case_when(grepl("CF_", Model.name) ~ "Cloud_forest",
           grepl("PO_", Model.name) ~ "Pine Oak_forest", 
  grepl("SS_", Model.name) ~ "Submontane_scrubland"))
  
  return(b)
  print(b)
  })

#Exporting results as a single table

ENM.area.tbl <- ENM.area %>% bind_rows() %>% arrange(Veg.type) %>% 
  mutate(Mask = case_when(grepl("_FC", Model.name) ~ "After mask", 
                          TRUE ~ "Before mask"))%>% 
  select(-Model.name) %>% spread(key = "Mask", value = "Area") %>% 
  mutate(Area.dif = `Before mask`- `After mask`) %>%
  mutate(Pct.dif = (`Area.dif`/`Before mask`)*100)

write.csv(ENM.area.tbl, "outputs/Sinusoidal/ENM_areas.csv")
