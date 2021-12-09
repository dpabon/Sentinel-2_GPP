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


# VI functions 

bi2 <- function(red, green, nir, red_factor = 1, green_factor = 1, nir_factor =1) {
  BI2 = sqrt( ( (red_factor * red * red_factor * red) + (green_factor * green * green_factor * green) + (nir_factor * nir * nir_factor * nir) ) / 3 )
  return(BI2)
}

bi <- function(red, green, red_factor = 1, green_factor = 1) {
  BI = sqrt( ( (red_factor * red * red_factor * red) + (green_factor * green * green_factor * green) ) / 2 )
  return(BI)
}

ci <- function(red, green, red_factor = 1, green_factor = 1) {
  CI = (red_factor * red - green_factor * green) / (red_factor * red + green_factor * green) 
  return(CI)
}

ri <- function(red, green, red_factor = 1, green_factor = 1) {
  RI = (red_factor * red * red_factor * red) / (green_factor * green * green_factor * green * green_factor * green)
  return(RI)
}

arvi <- function(red, blue, nir, red_factor = 1, blue_factor = 1, nir_factor = 1, gamma = 1) {
  rb = (red_factor * red) - gamma * (blue_factor * blue - red_factor * red)
  ARVI = (nir_factor * nir - rb) / (nir_factor * nir + rb)
  return(ARVI)
}

dvi <- function(red, nir, red_factor = 1, nir_factor = 1) {
  DVI = (nir_factor * nir - red_factor * red)
  return(DVI)
}

gemi <- function(red, nir, red_factor = 1, nir_factor = 1) {
  eta = (2 * (nir_factor * nir * nir_factor * nir - red_factor * red * red_factor * red) + 1.5 * nir_factor * nir + 0.5 * red_factor * red) / (nir_factor * nir + red_factor * red + 0.5)
  
  GEMI = eta * (1 - 0.25 * eta) - (red_factor * red - 0.125) / (1 - red_factor * red) 
  return(GEMI)
}

gndvi <- function(green, nir, green_factor = 1, nir_factor = 1) {
  GNDVI = (nir_factor * nir - green_factor * green) / (nir_factor * nir + green_factor * green)
  return(GNDVI)
}

ipvi <- function(red, nir, red_factor = 1, nir_factor = 1) {
  IPVI = (nir_factor * nir) / (nir_factor * nir + red_factor * red)
  return(IPVI)
}

ireci <- function(red_b4, red_b5, red_b6, nir_b7, red_b4_factor = 1, red_b5_factor = 1, red_b6_factor = 1, nir_b7_factor = 1){
  IRECI = (nir_b7_factor * nir_b7_factor - red_b4_factor * red_b4) / (red_b5_factor * red_b5 / red_b6_factor * red_b6)
  return(IRECI)
}

# ireci with B8 by BAA??

mcari <- function(red_b4, red_b5, green, red_b4_factor = 1, red_b5_factor = 1, green_factor = 1){
  MCARI = ((red_b5_factor * red_b5 - red_b4_factor * red_b4) - 0.2 * (red_b5_factor * red_b5 - green_factor * green)) * (red_b5_factor * red_b5 / red_b4_factor * red_b4)  
  return(MCARI)
}

msavi2 <- function(red, nir, red_factor = 1, nir_factor = 1) {
  MSAVI2 = (1/2) * ( 2 * nir_factor * nir + 1 - sqrt( ( 2 * nir_factor * nir + 1) * ( 2 * nir_factor * nir + 1) - 8 * (nir_factor * nir - red_factor * red) ) )
}

msavi <- function(red, nir, red_factor = 1, nir_factor = 1, slope = 0.5) {
  L = 1 - 2 * slope * ((nir_factor * nir - red_factor * red)/(nir_factor * nir + red_factor * red)) * (nir_factor * nir - slope * red_factor * red)
  MSAVI = (1 + L) * (nir_factor * nir - red_factor * red) / (nir_factor * nir + red_factor * red + L)
  return(MSAVI)
}

mtci <- function(red_b4, red_b5, nir_b6, red_b4_factor = 1, red_b5_factor = 1, nir_b6_factor = 1) {
  MTCI = (nir_b6_factor * nir_b6 - red_b5_factor * red_b5) / (red_b5_factor * red_b5 - red_b4_factor * red_b4)
  return(MTCI)
}

ndi45 <- function(red_b4, nir_b5, red_b4_factor = 1, nir_b5_factor = 1) {
  NDI45 = (nir_b5_factor * nir_b5 - red_b4_factor * red_b4) / (nir_b5_factor * nir_b5 + red_b4_factor * red_b4)
  return(NDI45)
}

pssra <- function(red_b4, nir_b7, red_b4_factor = 1, nir_b7_factor = 1) {
  PSSRa = (nir_b7_factor * nir_b7) / (red_b4_factor * red_b4)
  return(PSSRa)
}

pvi <- function(red_b4, nir_b8, red_b4_factor = 1, nir_b8_factor = 1, angle_soil_line = 45){
  PVI = sin(angle_soil_line) * nir_b8_factor * nir_b8 - cos(angle_soil_line) * red_b4_factor * red_b4
  return(PVI)
}

reip <- function(red_b4, red_b5, red_b6, nir_b7, red_b4_factor = 1, red_b5_factor =1, red_b6_factor = 1, nir_b7_factor = 1){
  REIP = 700 + 40 * ( (red_b4_factor * red_b4 + nir_b7_factor * nir_b7)/2 - red_b5_factor * red_b5) / (red_b6_factor * red_b6 - red_b5_factor * red_b5) 
  return(REIP)
}

rvi <- function(red_b4, nir_b8, red_b4_factor = 1, nir_b8_factor = 1) {
  RVI = (nir_b8_factor * nir_b8) / (red_b4_factor * red_b4) 
  return(RVI)
}


s2rep <- function(red_b4, red_b5, red_b6, nir_b7, red_b4_factor = 1, red_b5_factor = 1, red_b6_factor = 1, nir_b7_factor = 1){
  S2REP = 705 + 35 * ( (red_b4_factor * red_b4 + nir_b7_factor * nir_b7)/2 - red_b5_factor * red_b5) / (red_b6_factor * red_b6 - red_b5_factor * red_b5)
  return(S2REP)
}

savi <- function(red_b4, nir_b8, red_b4_factor =1, nir_b8_factor = 1, soil_correction = 0.5){
  SAVI = (1 + soil_correction) * (nir_b8_factor * nir_b8 - red_b4_factor * red_b4) / (nir_b8_factor * nir_b8 + red_b4_factor * red_b4 + soil_correction)
  return(SAVI)
}

tndvi <- function(red_b4, nir_b8, red_b4_factor = 1, nir_b8_factor = 1){
  TNDVI = sqrt( (nir_b8_factor * nir_b8 - red_b4_factor * red_b4) / (nir_b8_factor * nir_b8 + red_b4_factor * red_b4) + 0.5)
  return(TNDVI)
}

tsavi <- function(red_b4, nir_b8, red_b4_factor = 1, nir_b8_factor = 1, slope = 0.5, intercept = 0.5, adjustment = 0.08) {
  TSAVI = slope * (nir_b8_factor * nir_b8 - slope * red_b4_factor * red_b4 - intercept) / (slope * nir_b8_factor * nir_b8 + red_b4_factor * red_b4 - intercept * slope + adjustment * ( 1 + slope * slope ))
  return(TSAVI)
}


wdvi <- function(red_b4, nir_b8, red_b4_factor = 1, nir_b8_factor = 1, slope = 0.5) {
  WDVI = (nir_b8_factor * nir_b8 - slope * red_b4_factor * red_b4) 
  return(WDVI)
}

cir <- function(red_b5, red_b7, red_b5_factor = 1, red_b7_factor = 1) {
  CIR = ((red_b7_factor * red_b7) / (red_b5_factor * red_b5)) -1
  return(CIR)
}

cig <- function(green, red_b7, green_factor = 1, red_b7_factor = 1) {
  CIG = ((red_b7_factor * red_b7) / (green_factor * green)) -1
}

ndvi <- function(red_b4, nir_b8, red_b4_factor = 1, nir_b8_factor = 1) {
  NDVI = (nir_b8_factor * nir_b8 - red_b4_factor * red_b4) / (nir_b8_factor + nir_b8 - red_b4_factor * red_b4)
  return(NDVI)
}

nirv <- function(red_b4, nir_b8, red_b4_factor = 1, nir_b8_factor = 1) {
  NIRV = ((nir_b8_factor * nir_b8 - red_b4_factor * red_b4) / (nir_b8_factor + nir_b8 - red_b4_factor * red_b4)) * nir_b8_factor * nir_b8
}

all_vi <- function(red, blue, green, red_b4, red_b5, red_b6, red_b7, nir, nir_b5, nir_b6, nir_b7, nir_b8){
  r1 <- bi2(red = red, green = green, nir = nir)
  r2 <- bi(red = red, green = green)
  r3 <- ci(red = red, green = green)
  r4 <- ri(red = red, green = green)
  r5 <- arvi(red = red, blue = blue, nir = nir)
  r6 <- dvi(red = red, nir = nir)
  r7 <- gemi(red = red, nir = nir)
  r8 <- gndvi(green = green, nir = nir)
  r9 <- ipvi(red = red, nir = nir)
  r10 <- ireci(red_b4 = red_b4, red_b5 = red_b5, red_b6 = red_b6, nir_b7 = nir_b7)
  r11 <- mcari(red_b4 = red_b4, red_b5 = red_b5, green = green)
  r12 <- msavi2(red = red, nir = nir)
  r13 <- msavi(red = red, nir = nir)
  r14 <- mtci(red_b4 = red_b4, red_b5 = red_b5, nir_b6 = nir_b6)
  r15 <- ndi45(red_b4 = red_b4, nir_b5 = nir_b5)
  r16 <- pssra(red_b4 = red_b4, nir_b7 = nir_b7)
  r17 <- pvi(red_b4 = red_b4, nir_b8 = nir_b8)
  r18 <- reip(red_b4 = red_b4, red_b5 = red_b5, red_b6 = red_b6, nir_b7 = nir_b7)
  r19 <- rvi(red_b4 = red_b4, nir_b8 = nir_b8)
  r20 <- s2rep(red_b4 = red_b4, red_b5 = red_b5, red_b6 = red_b6, nir_b7 = nir_b7)
  r21 <- savi(red_b4 = red_b4, nir_b8 = nir_b8)
  r22 <- tndvi(red_b4 = red_b4, nir_b8 = nir_b8)
  r23 <- tsavi(red_b4 = red_b4, nir_b8 = nir_b8)
  r24 <- wdvi(red_b4 = red_b4, nir_b8 = nir_b8)
  r25 <- cir(red_b5 = red_b5, red_b7 = red_b7)
  r26 <- cig(green = green, red_b7 = red_b7)
  r27 <- ndvi(red_b4 = red_b4, nir_b8 = nir_b8)
  r28 <- nirv(red_b4 = red_b4, nir_b8 = nir_b8)
  output <- stack(r1,
                  r2,
                  r3,
                  r4,
                  r5,
                  r6,
                  r7,
                  r8,
                  r9,
                  r10,
                  r11,
                  r12,
                  r13,
                  r14,
                  r15,
                  r16,
                  r17,
                  r18,
                  r19,
                  r20,
                  r21,
                  r22,
                  r23,
                  r24,
                  r25,
                  r26,
                  r27,
                  r28)
  names(output) <- c("bi2_mean",
                     "bi_mean",
                     "ci_mean",
                     "ri_mean",
                     "arvi_mean",
                     "dvi_mean",
                     "gemi_mean",
                     "gndvi_mean",
                     "ipvi_mean",
                     "ireci_mean",
                     "mcari_mean",
                     "msavi2_mean",
                     "msavi_mean",
                     "mtci_mean",
                     "ndi45_mean",
                     "pssra_mean",
                     "pvi_mean",
                     "reip_mean",
                     "rvi_mean",
                     "s2rep_mean",
                     "savi_mean",
                     "tndvi_mean",
                     "tsavi_mean",
                     "wdvi_mean",
                     "cir_mean",
                     "cig_mean",
                     "ndvi_mean",
                     "nirv_mean")
  return(output)
}