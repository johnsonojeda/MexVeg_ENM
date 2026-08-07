##########################################################
###          PERCENT OVERLAP OTHER NATURAL LAND COVER ####
##########################################################
#' Estimating the % of overlap between predictions and other natural land cover

library(terra)
library(dplyr)

#########################################
## 1. UPLOADING ENMs & OTHER LAND COVER
#########################################
#Loading VCF data
VCF <- terra::rast("/home/erica/PhD/Mex_VegModel/Analyses/ForestMask/data/MODIS_VCF/MODISVCF_2016-2018Avg.tif")

# Loading unprocessed and post-processed ENMs

#Continuous ENMs
CH.enm.dir <- '/home/erica/PhD/Mex_VegModel/Analyses/ENM/Outputs/RasPreds/CarsoHuasteco'
CH.enm.path <- file.path(CH.enm.dir, dir(CH.enm.dir, pattern = "25_mtp.tif$"))
CH.enm <- terra::rast(CH.enm.path)
#CH.enm <- terra::project(CH.enm, crs(VCF))

#WBC masked ENMs
CH.wbc.dir <- '/home/erica/PhD/Mex_VegModel/Analyses/SVM/outputs/CarsoHuasteco/WBC'
CH.wbc.path <- file.path(CH.wbc.dir, dir(CH.wbc.dir, pattern = "WBC_mtp.tif$"))
CH.wbc <- terra::rast(CH.wbc.path)
#CH.wbc <- terra::project(CH.wbc, crs(VCF))

#SVMsp masked ENMs
CH.SVMsp.dir <- '/home/erica/PhD/Mex_VegModel/Analyses/SVM/outputs/CarsoHuasteco/SVMsp'
CH.SVMsp.path <- file.path(CH.SVMsp.dir, dir(CH.SVMsp.dir, pattern = "SVMsp_mtp.tif$"))
CH.SVMsp <- terra::rast(CH.SVMsp.path)
#CH.SVMsp <- terra::project(CH.SVMsp, crs(VCF))

#SVMhyb masked ENMs
CH.SVMhyb.dir <- '/home/erica/PhD/Mex_VegModel/Analyses/SVM/outputs/CarsoHuasteco/SVMhyb'
CH.SVMhyb.path <- file.path(CH.SVMhyb.dir, dir(CH.SVMhyb.dir, pattern = "SVMhyb_mtp.tif$"))
CH.SVMhyb <- terra::rast(CH.SVMhyb.path)
#CH.SVMhyb <- terra::project(CH.SVMhyb, crs(VCF))

#Creating a list of ENMs & projecting to sinusoidal
ENM.list <- list(CH.enm, CH.wbc, CH.SVMsp, CH.SVMhyb)
ENM.list <- lapply(ENM.list, project, y = VCF)

#Importing other land cover data to exclude
CH.Natural <- terra::vect("/home/erica/PhD/Mex_VegModel/Analyses/Data/CountryData/LandCover/OtherNatural_LandCover.shp")
CH.Natural <- terra::project(CH.Natural, crs(VCF))

############################################
## 2. % AREA OVERLAP WITH OTHER NATURAL VEG
############################################
#Calculating the percent area of pre and post-processed ENMs fallin within
#polygons identified as other natural vegetation not considered in this study

ENM.area.overlap <- lapply(ENM.list, function(x){
#Area of one raster cell (Km2)
cell.area <- prod(terra::res(x))/1e6

#Calculating pct overlap with other natural vegetation per raster
pct_overlap <- lapply(seq_len(terra::nlyr(x)), function(i){
  r <- x[[i]]
  
  # Total ENM area
    enm.area <- terra::global(
      !is.na(r),
      "sum",
      na.rm = TRUE
    )[1,1] * cell.area

    # Intersect with natural vegetation
    r.int <- terra::mask(
      terra::crop(r, CH.Natural),
      CH.Natural
    )

    # Area inside polygon
    inter.area <- terra::global(
      !is.na(r.int),
      "sum",
      na.rm = TRUE
    )[1,1] * cell.area

    data.frame(
      Model.name = names(x)[i],
      ENM.area = enm.area,
      Intersect.area = inter.area,
      Percent.overlap = 100 * inter.area / enm.area
    )
  })

  do.call(rbind, pct_overlap)
})

ENM.area.overlap <- dplyr::bind_rows(ENM.area.overlap) %>%
  dplyr::mutate(Veg.type = case_when(
    grepl("^CF", Model.name) ~ "Cloud_forest",
  grepl("^PO", Model.name) ~ "PineOak_forest",
grepl("^SS", Model.name)~ "Submontane_scrubland",
TRUE ~ NA_character_))%>%
  dplyr::arrange(Veg.type) %>%
  dplyr::relocate(Veg.type, .before = Model.name)

# Export results
write.csv(
  ENM.area.overlap,
  "/home/erica/PhD/Mex_VegModel/Analyses/ENM/Outputs/Overlap_OtherVeg.csv",
  row.names = FALSE
)

#'THE PERCENT OVERLAP IS FAIRLY THE SAME ACROSS PRE AND POST-PROCESSED MODELS, CHECK THE INTERSECT TO MAKE SURE 
#' MASKING IS BEING DONE RIGHT AND RECALCULATE IF NEEDED
