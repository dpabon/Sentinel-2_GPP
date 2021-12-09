#' Predicting Ecosystem Functional Properties using Sentinel-2
#'
#' This function use pre-trained regression trees to predict GPP.
#'
#' @param sentinel2_product The path to a Sentinel-2 L2A product downloaded from https://scihub.copernicus.eu/ or generated with sen2cor. (String)
#' @param model RF model to predict GPP
#' @param product Can be "gpp"
#'
#' @return
#'  An object of class raster if product = "gpp".
#' @export
#'
#' @examples
sentinel2.efps <- function(sentinel2_product, model, product = "gpp") {
  library(rgdal)
  library(raster)
  library(randomForest)
  source("~/dpabon/bin_EFPs_ML/EFPs_ML/VI_functions.R")
  
  images <- Sys.glob(file.path(sentinel2_product, "GRANULE", "*", "IMG_DATA", "*", "*"))
  length(images)
  bands_20m_for_forest <- c("B2_mean", "B3_mean", "B4_mean", "B5_mean", "B6_mean", "B7_mean", "B11_mean", "B12_mean", "B8A_mean")
  
  images.position <- c(9:17)
  
  bands20.raster <- stack(readGDAL(images[images.position[1]]))
  
  for (i in 2:length(images.position)){
    bands20.raster <- addLayer(bands20.raster, stack(readGDAL(images[images.position[i]])))
  }
  
  names(bands20.raster) <- bands_20m_for_forest
  
  # bands to resample
  
  B08 <- raster(readGDAL(images[5]))
  B08 <- resample(B08, bands20.raster$B2_mean, method = "ngb")
  bands20.raster <- addLayer(bands20.raster, B08)
  names(bands20.raster)[length(names(bands20.raster))] <- "B8_mean"
  rm(B08)
  gc()
  
  B01 <- raster(readGDAL(images[22]))
  B09 <- raster(readGDAL(images[29]))
  
  B01 <- resample(B01, bands20.raster$B2_mean, method = "ngb")
  B09 <- resample(B09, bands20.raster$B2_mean, method = "ngb")
  
  bands20.raster <- addLayer(bands20.raster, B01, B09)
  names(bands20.raster)[length(names(bands20.raster))-1] <- "B1_mean"
  names(bands20.raster)[length(names(bands20.raster))] <- "B9_mean"
  rm(B01)
  rm(B09)
  gc()
  
  names(bands20.raster)
  
  stack.vi <- all_vi(red = bands20.raster$B4_mean, blue = bands20.raster$B2_mean, green = bands20.raster$B3_mean, red_b4 = bands20.raster$B4_mean, red_b5 = bands20.raster$B5_mean, red_b6 = bands20.raster$B6_mean, red_b7 = bands20.raster$B7_mean, nir = bands20.raster$B8_mean, nir_b5 = bands20.raster$B5_mean, nir_b6 = bands20.raster$B6_mean, nir_b7 = bands20.raster$B7_mean, nir_b8 = bands20.raster$B8_mean)
  
  stack.vi <- addLayer(stack.vi, bands20.raster)
  
  rm(bands20.raster)
  gc()
  
  # applying quality scene classification
  
  csc <- raster(readGDAL(images[18]))
  
  stack.vi <- raster::mask(stack.vi, csc, maskvalue = 4, inverse = T)
  
  if(product == "gpp"){
    output.gpp <- raster::predict(stack.vi, model, na.rm = T, inf.rm = T)
    names(output.gpp) <- "GPP"
  }
  return(output.gpp)
}
