# Upscaling GPP from Sentinel-2 L2A products
# created by : Daniel E. Pabon Moreno
# email: dpabon@bgc-jena.mpg.de


# Require 

library(rgdal)
library(raster)
library(randomForest)
library(RStoolbox)
library(terrainr)
library(ggplot2)

# loading function to predict GPP ----

source("final_models/sentinel-2_gpp_function.R")

# loading pre-trained model

smoter_GPP_model <- readRDS("final_models/final_model_rf_smoter.rds") ---

# locating L2A product

sentinel2.l2a <- "/Supplement_material_6/S2A_MSIL2A_20200623T103031_N0214_R108_T32ULU_20200623T142851.SAFE/"

# upscaling GPP

sentinel2.l2a.GPP <- sentinel2.efps(sentinel2_product = sentinel2.l2a, model = smoter_GPP_model,  product = "gpp")

# plotting the result

ggplot()+
  geom_raster(data = sentinel2.l2a.GPP, aes(x = x, y = y, fill = GPP)) +
  scale_fill_viridis(direction = -1) +
  coord_quickmap() +
  xlab(expression(paste("Longitude (",degree,")"))) + ylab(expression(paste("Latitude (",degree,")"))) +
  guides(fill = guide_colourbar(title = latex2exp::TeX('GPP $\\lbrack \\mu mol  CO_{2}  m^{-2}  s^{-1} \\rbrack$') , title.position = "right"))+
  theme_bw() +
  theme(legend.title = element_text(angle = 90, hjust = 0.5, vjust = 0.5))
