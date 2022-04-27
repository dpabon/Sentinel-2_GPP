library(ncdf4)
library(raster)
library(easyNCDF)
library(tictoc)
library(rgdal)
library(abind)
library(tidyverse)
library(outliers)

setwd("~/l1_work3/cube_f/")

outlier.var <- function(x){
  y <- x[!is.na(x)]
  y.good <- which(scores(y, type = "z", prob = 0.99) == F)
  y <- as.array(y[y.good])
  return(y)
}

useful.inf  <- function(x, y, mask){
  # x = a vector with useful information
  # y = an quality scene classification array where dim = lat, lon, time
  # mask = m.masked
  # new dimension with m.mask
  l.mask <- length(which(is.na(mask) == F))
  return((length(which(is.element(y, x))))  / l.mask)
}

sites <- read_csv("~/dpabon/data_EFPS_ML/summary_ec_sites.csv")

sites <- sites %>%
  filter(!igbp %in% c("URB", "WSA", "WAT"))

site_name <- c("Be-Bra",
               "BE-Lon",
               "BE-Vie",
               "CH-Aws",
               "CH-Cha",
               "CH-Dav",
               "CH-Fru",
               "CH-Lae",
               "CH-Oe2",
               "CZ-BK1",
               "CZ-Lnz",
               "CZ-RAJ",
               "CZ-Stn",
               "CZ-wet",
               "DE-Akm",
               "DE-Geb",
               "DE-Gri",
               "DE-Hai",
               "DE-HoH",
               "DE-Hte",
               "DE-Hzd",
               "DE-Kli",
               "DE-Obe",
               "DE-RuR",
               "DE-RuS",
               "DE-RuW",
               "DE-Tha",
               "DK-Sor",
               "Es-Abr",
               "ES-LM1",
               "ES-LM2",
               "FI-Hyy",
               "FI-Let",
               "FI-Sii",
               "FI-Var",
               "FR-EM2",
               "FR-Hes",
               "IT-BCi",
               "IT-Lsn",
               "IT-Tor",
               "NL-Loo",
               "RU-Fy2",
               "RU-Fyo",
               "SE-Deg",
               "SE-Htm",
               "SE-Lnn",
               "SE-Nor",
               "SE-Ros",
               "SE-Svb",
               "US-ARM",
               "US-Bar",
               "US-Ho1",
               "US-MMS",
               "US-Seg",
               "US-Ses",
               "US-UMB",
               "US-Vcm",
               "US-Wjs")
# First cube 

if(sites$Multiple_Tiles[1] == FALSE){
  cube <- nc_open(paste(site_name[1], "/", site_name[1], "_2015_2018.nc", sep = ""))
  #cube
  
  dates <- ncvar_get(nc = cube, varid = "time")
  
  dates <- as.Date(dates, origin = "1970-01-01")
  
  lat <- ncvar_get(cube, "lat")
  lon <- ncvar_get(cube, "lon")
  
  qc <- ncvar_get(cube, "quality_scene_classification")
  
  # masking 200 m radius around the EC tower
  
  template_raster <- raster(paste(site_name[1], "/", site_name[1], "_2015_2018.nc", sep = ""), level = 1, varname = "bi2", stopIfNotEqualSpaced = F)
  template_raster
  # shapefile with the area to cut
  mask_shape  <- readOGR('~/dpabon/data_EFPS_ML/buffer_100m/buffer_100m.shp')
  
  
  mask_shape <- mask_shape[mask_shape$site_name == site_name[1],]
  
  #mask_shape <- extent(mask_shape)
  
  r <- rasterize(mask_shape, template_raster)
  
  masked <- mask(x = template_raster, mask = r)
  
  m.masked <- as.matrix(masked)
  m.masked <- t(m.masked)
  # image(m.masked)
  for(i in 1:ncol(m.masked)){
    for(d in 1:nrow(m.masked)){
      if(is.na(m.masked[d,i]) == T) {
        qc[d,i,] <- NA
      }
    }
  }
  qc_sequence <- apply(X = qc, MARGIN = c(3), FUN = useful.inf, x = useful, mask = m.masked)
  
  n.images <- length(which(qc_sequence >= qc_filter))
  
  images.to.select <- which(qc_sequence >= qc_filter)
  dates.to.select <- dates[which(qc_sequence >= qc_filter)]
  
  qc.mask <- qc[,,images.to.select]
  
  #dim(qc.mask)
  
  # Open all the variables 
  
  variables_names <- names(cube$var)
  
  # removing quality scene classification
  variables_names <- variables_names[-26]
  
  variables <- NcOpen(file_path = paste(site_name[1], "/", site_name[1], "_2015_2018.nc", sep = ""))
  
  
  variables <- NcToArray(variables, vars_to_read = variables_names, dim_indices = list(time = images.to.select))
  
  
  for (t in 1:dim(variables)[4]) {
    for (v in 1:dim(variables)[1]) {
      variables[v,,,t][-which(is.element(qc.mask[,,t], useful))] <- NA
    }
  }
  
  
  # estimating other VIs
  
  temp <- array(data = NA, dim = c(var =4, lon = dim(variables)[2], lat = dim(variables)[3], time = dim(variables)[4]))
  
  for(i in 1:dim(variables)[4]) {
    for(x in 1:dim(variables)[3]){
      for (d in 1:dim(variables)[2]){
        # CI red (CIr) (B7/B5) -1
        temp[1,d,x,i] <- (variables[20,d,x,i] / variables[18,d,x,i]) - 1
        # CI green (CIg) (B7/B3) -1
        temp[2,d,x,i] <- (variables[20,d,x,i] / variables[16,d,x,i]) - 1
        # ndvi
        temp[3,d,x,i] <- ((variables[21,d,x,i] - variables[17,d,x,i]) / (variables[21,d,x,i] +  variables[17,d,x,i]))
        # nirv
        temp[4,d,x,i] <- ((variables[21,d,x,i] - variables[17,d,x,i]) / (variables[21,d,x,i] +  variables[17,d,x,i])) * variables[21,d,x,i]
      }
    }
  }
  
  variables <- abind(variables, temp, along = 1)
  
  variables_names <- c(variables_names, "cir", "cig", "ndvi", "nirv")
  
  # masking forest pixels
  
  # forest <- raster("~/dpabon/to_cluster/data/euroflux/land_cover_forest/TCD_2015_020m_eu_03035_d05_Full/TCD_2015_020m_eu_03035_d05_full.tif")
  
  # crs_o <- "+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext +no_defs"
  
  # forest_proj <- projectRaster(forest, crs_o)
  
  # plot(forest)
  
  # cropping the buffer area
  # extend (min_lon, max_lon, min_lat, max_lat)
  
  # min(lat)
  # max(lat)
  # b_forest <- extent(min(lon), max(lon), min(lat), max(lat))
  # plot(b_forest)
  # crop(forest, b_forest)
  
  # outlier detection 
  
  
  depured.variables <- apply(X = variables, MARGIN = c(1,4), outlier.var)
  
  
  # estimating sd and mean for each variable - time
  
  mean.variable <- matrix(data = NA, nrow = dim(depured.variables)[1], ncol = dim(depured.variables)[2])
  sd.variable <- matrix(data = NA, nrow = dim(depured.variables)[1], ncol = dim(depured.variables)[2])
  
  for (i in 1:dim(depured.variables)[1]){
    for (j in 1:dim(depured.variables)[2]){
      mean.variable[i,j] <- mean(unlist(depured.variables[i,j]))
      sd.variable[i,j] <- sd(unlist(depured.variables[i,j]))
    }
  }
  
  # Genereting a dataframe with the structure:
  # Site | Time | Var1.mean | Var1.sd | ... |
  
  output.info.ML <- data.frame(site = rep(site_name[1], length(dates.to.select)))
  
  output.info.ML[,"date"] <- dates.to.select
  
  
  for (i in 1:dim(mean.variable)[1]) {
    # mean
    output.info.ML[,paste(variables_names[i], "_mean", sep = "")] <- mean.variable[i,]
    # sd
    output.info.ML[,paste(variables_names[i], "_sd", sep = "")] <- sd.variable[i,]
  }
  
  
}else{
  for (tile in 1:2){
    cube <- nc_open(paste(site_name[1], "/", site_name[1], "_", tile, "_2015_2018.nc", sep = ""))
    #cube
    
    
    dates <- ncvar_get(nc = cube, varid = "time")
    
    dates <- as.Date(dates, origin = "1970-01-01")
    
    lat <- ncvar_get(cube, "lat")
    lon <- ncvar_get(cube, "lon")
    
    qc <- ncvar_get(cube, "quality_scene_classification")
    
    # masking 200 m radius around the EC tower
    
    template_raster <- raster(paste(site_name[1], "/", site_name[1], "_", tile, "_2015_2018.nc", sep = ""), level = 1, varname = "bi2", stopIfNotEqualSpaced = F)
    template_raster
    # shapefile with the area to cut
    mask_shape  <- readOGR('~/dpabon/data_EFPS_ML/buffer_100m/buffer_100m.shp')
    
    
    mask_shape <- mask_shape[mask_shape$site_name == site_name[1],]
    
    #mask_shape <- extent(mask_shape)
    
    r <- rasterize(mask_shape, template_raster)
    
    masked <- mask(x = template_raster, mask = r)
    
    m.masked <- as.matrix(masked)
    m.masked <- t(m.masked)
    # image(m.masked)
    for(i in 1:ncol(m.masked)){
      for(d in 1:nrow(m.masked)){
        if(is.na(m.masked[d,i]) == T) {
          qc[d,i,] <- NA
        }
      }
    }
    
    
    
    qc_sequence <- apply(X = qc, MARGIN = c(3), FUN = useful.inf, x = useful, mask = m.masked)
    
    n.images <- length(which(qc_sequence >= qc_filter))
    images.to.select <- which(qc_sequence >= qc_filter)
    dates.to.select <- dates[which(qc_sequence >= qc_filter)]
    
    qc.mask <- qc[,,images.to.select]
    
    #dim(qc.mask)
    
    # Open all the variables 
    
    variables_names <- names(cube$var)
    
    # removing quality scene classification
    variables_names <- variables_names[-26]
    
    variables <- NcOpen(file_path = paste(site_name[1], "/", site_name[1], "_", tile, "_2015_2018.nc", sep = ""))
    
    
    variables <- NcToArray(variables, vars_to_read = variables_names, dim_indices = list(time = images.to.select))
    
    
    # masking non useful values for all variables all times with NA
    
    for (t in 1:dim(variables)[4]) {
      for (v in 1:dim(variables)[1]) {
        variables[v,,,t][-which(is.element(qc.mask[,,t], useful))] <- NA
      }
    }
    
    
    # estimating other VIs
    
    temp <- array(data = NA, dim = c(var =4, lon = dim(variables)[2], lat = dim(variables)[3], time = dim(variables)[4]))
    
    for(i in 1:dim(variables)[4]) {
      for(x in 1:dim(variables)[3]){
        for (d in 1:dim(variables)[2]){
          # CI red (CIr) (B7/B5) -1
          temp[1,d,x,i] <- (variables[20,d,x,i] / variables[18,d,x,i]) - 1
          # CI green (CIg) (B7/B3) -1
          temp[2,d,x,i] <- (variables[20,d,x,i] / variables[16,d,x,i]) - 1
          # ndvi
          temp[3,d,x,i] <- ((variables[21,d,x,i] - variables[17,d,x,i]) / (variables[21,d,x,i] +  variables[17,d,x,i]))
          # nirv
          temp[4,d,x,i] <- ((variables[21,d,x,i] - variables[17,d,x,i]) / (variables[21,d,x,i] +  variables[17,d,x,i])) * variables[21,d,x,i]
        }
      }
    }
    
    variables <- abind(variables, temp, along = 1)
    
    variables_names <- c(variables_names, "cir", "cig", "ndvi", "nirv")
    
    # masking forest pixels
    
    # forest <- raster("~/dpabon/to_cluster/data/euroflux/land_cover_forest/TCD_2015_020m_eu_03035_d05_Full/TCD_2015_020m_eu_03035_d05_full.tif")
    
    # crs_o <- "+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext +no_defs"
    
    # forest_proj <- projectRaster(forest, crs_o)
    
    # plot(forest)
    
    # cropping the buffer area
    # extend (min_lon, max_lon, min_lat, max_lat)
    
    # min(lat)
    # max(lat)
    # b_forest <- extent(min(lon), max(lon), min(lat), max(lat))
    # plot(b_forest)
    # crop(forest, b_forest)
    
    # outlier detection 

    
    depured.variables <- apply(X = variables, MARGIN = c(1,4), outlier.var)
    
    
    # estimating sd and mean for each variable - time
    
    mean.variable <- matrix(data = NA, nrow = dim(depured.variables)[1], ncol = dim(depured.variables)[2])
    sd.variable <- matrix(data = NA, nrow = dim(depured.variables)[1], ncol = dim(depured.variables)[2])
    
    for (i in 1:dim(depured.variables)[1]){
      for (j in 1:dim(depured.variables)[2]){
        mean.variable[i,j] <- mean(unlist(depured.variables[i,j]))
        sd.variable[i,j] <- sd(unlist(depured.variables[i,j]))
      }
    }
    
    # Genereting a dataframe with the structure:
    if(tile ==1){
    # Site | Time | Var1.mean | Var1.sd | ... |
    
      output.info.ML <- data.frame(site = rep(site_name[1], length(dates.to.select)))
      
      output.info.ML[,"date"] <- dates.to.select
      
      
      for (i in 1:dim(mean.variable)[1]) {
        # mean
        output.info.ML[,paste(variables_names[i], "_mean", sep = "")] <- mean.variable[i,]
        # sd
        output.info.ML[,paste(variables_names[i], "_sd", sep = "")] <- sd.variable[i,]
      }
    }else{
      output.info.ML.temp <- data.frame(site = rep(sites$site_name[1], length(dates.to.select)))
      
      output.info.ML.temp[,"date"] <- dates.to.select
      
      
      for (i in 1:dim(mean.variable)[1]) {
        # mean
        output.info.ML.temp[,paste(variables_names[i], "_mean", sep = "")] <- mean.variable[i,]
        # sd
        output.info.ML.temp[,paste(variables_names[i], "_sd", sep = "")] <- sd.variable[i,]
      }
    }
  }
  output.info.ML <- rbind(output.info.ML, output.info.ML.temp)
}
nrow(output.info.ML)

  write.csv(output.info.ML, paste("~/dpabon/data_EFPS_ML/sentinel2_100m_info_v4/", sites$site_name[1], "_sentinel2_100m_info.csv", sep = "") , row.names = F)

# The rest of the sites indexing to the same data.frame

for(st in 23:nrow(sites)) {
  print(paste("extracting information for:", site_name[st]))
  
  if(sites$Multiple_Tiles[st] == FALSE){
    cube <- nc_open(paste(site_name[st], "/", site_name[st], "_2015_2018.nc", sep = ""))
    #cube
    
    dates <- ncvar_get(nc = cube, varid = "time")
    
    dates <- as.Date(dates, origin = "1970-01-01")
    
    lat <- ncvar_get(cube, "lat")
    lon <- ncvar_get(cube, "lon")
    
    qc <- ncvar_get(cube, "quality_scene_classification")
    
    # masking 200 m radius around the EC tower
    
    template_raster <- raster(paste(site_name[st], "/", site_name[st], "_2015_2018.nc", sep = ""), level = 1, varname = "bi2", stopIfNotEqualSpaced = F)
    template_raster
    # shapefile with the area to cut
    mask_shape  <- readOGR('~/dpabon/data_EFPS_ML/buffer_100m/buffer_100m.shp')
    
    
    mask_shape <- mask_shape[mask_shape$site_name == site_name[st],]
    
    #mask_shape <- extent(mask_shape)
    
    r <- rasterize(mask_shape, template_raster)
    
    masked <- mask(x = template_raster, mask = r)
    
    m.masked <- as.matrix(masked)
    m.masked <- t(m.masked)
    # image(m.masked)
    for(i in 1:ncol(m.masked)){
      for(d in 1:nrow(m.masked)){
        if(is.na(m.masked[d,i]) == T) {
          qc[d,i,] <- NA
        }
      }
    }
    qc_sequence <- apply(X = qc, MARGIN = c(3), FUN = useful.inf, x = useful, mask = m.masked)
    
    n.images <- length(which(qc_sequence >= qc_filter))
    if(n.images == 0){
      next
    }
    images.to.select <- which(qc_sequence >= qc_filter)
    dates.to.select <- dates[which(qc_sequence >= qc_filter)]
    
    qc.mask <- qc[,,images.to.select]
    
    #dim(qc.mask)
    
    # Open all the variables 
    
    variables_names <- names(cube$var)
    
    # removing quality scene classification
    variables_names <- variables_names[-26]
    
    variables <- NcOpen(file_path = paste(site_name[st], "/", site_name[st], "_2015_2018.nc", sep = ""))
    
    
    variables <- NcToArray(variables, vars_to_read = variables_names, dim_indices = list(time = images.to.select))
    
    
    # masking non useful values for all variables all times with NA
    
    
    for (t in 1:dim(variables)[4]) {
      for (v in 1:dim(variables)[1]) {
        variables[v,,,t][-which(is.element(qc.mask[,,t], useful))] <- NA
      }
    }
    
    
    # estimating other VIs
    
    temp <- array(data = NA, dim = c(var =4, lon = dim(variables)[2], lat = dim(variables)[3], time = dim(variables)[4]))
    
    for(i in 1:dim(variables)[4]) {
      for(x in 1:dim(variables)[3]){
        for (d in 1:dim(variables)[2]){
          # CI red (CIr) (B7/B5) -1
          temp[1,d,x,i] <- (variables[20,d,x,i] / variables[18,d,x,i]) - 1
          # CI green (CIg) (B7/B3) -1
          temp[2,d,x,i] <- (variables[20,d,x,i] / variables[16,d,x,i]) - 1
          # ndvi
          temp[3,d,x,i] <- ((variables[21,d,x,i] - variables[17,d,x,i]) / (variables[21,d,x,i] +  variables[17,d,x,i]))
          # nirv
          temp[4,d,x,i] <- ((variables[21,d,x,i] - variables[17,d,x,i]) / (variables[21,d,x,i] +  variables[17,d,x,i])) * variables[21,d,x,i]
        }
      }
    }
    
    variables <- abind(variables, temp, along = 1)
    
    variables_names <- c(variables_names, "cir", "cig", "ndvi", "nirv")
    
    # masking forest pixels
    
    # forest <- raster("~/dpabon/to_cluster/data/euroflux/land_cover_forest/TCD_2015_020m_eu_03035_d05_Full/TCD_2015_020m_eu_03035_d05_full.tif")
    
    # crs_o <- "+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext +no_defs"
    
    # forest_proj <- projectRaster(forest, crs_o)
    
    # plot(forest)
    
    # cropping the buffer area
    # extend (min_lon, max_lon, min_lat, max_lat)
    
    # min(lat)
    # max(lat)
    # b_forest <- extent(min(lon), max(lon), min(lat), max(lat))
    # plot(b_forest)
    # crop(forest, b_forest)
    
    # outlier detection 
    
    
    depured.variables <- apply(X = variables, MARGIN = c(1,4), outlier.var)
    
    
    # estimating sd and mean for each variable - time
    
    mean.variable <- matrix(data = NA, nrow = dim(depured.variables)[1], ncol = dim(depured.variables)[2])
    sd.variable <- matrix(data = NA, nrow = dim(depured.variables)[1], ncol = dim(depured.variables)[2])
    
    for (i in 1:dim(depured.variables)[1]){
      for (j in 1:dim(depured.variables)[2]){
        mean.variable[i,j] <- mean(unlist(depured.variables[i,j]))
        sd.variable[i,j] <- sd(unlist(depured.variables[i,j]))
      }
    }
    
    
    # Genereting a dataframe with the structure:
    # Site | Time | Var1.mean | Var1.sd | ... |
    
    output.info.ML <- data.frame(site = rep(sites$site_name[st], length(dates.to.select)))
    
    output.info.ML[,"date"] <- dates.to.select
    
    
    for (i in 1:dim(mean.variable)[1]) {
      # mean
      output.info.ML[,paste(variables_names[i], "_mean", sep = "")] <- mean.variable[i,]
      # sd
      output.info.ML[,paste(variables_names[i], "_sd", sep = "")] <- sd.variable[i,]
    }
    
    write.csv(output.info.ML, paste("~/dpabon/data_EFPS_ML/sentinel2_100m_info_v4/", sites$site_name[st], "_sentinel2_100m_info.csv", sep = "") , row.names = F)
    
  }else{
    for (tile in 1:2){
      cube <- nc_open(paste(site_name[st], "/", site_name[st], "_", tile, "_2015_2018.nc", sep = ""))
      #cube
      
      
      dates <- ncvar_get(nc = cube, varid = "time")
      
      dates <- as.Date(dates, origin = "1970-01-01")
      
      lat <- ncvar_get(cube, "lat")
      lon <- ncvar_get(cube, "lon")
      
      qc <- ncvar_get(cube, "quality_scene_classification")
      
      # masking 200 m radius around the EC tower
      
      template_raster <- raster(paste(site_name[st], "/", site_name[st], "_", tile, "_2015_2018.nc", sep = ""), level = 1, varname = "bi2", stopIfNotEqualSpaced = F)
      template_raster
      # shapefile with the area to cut
      mask_shape  <- readOGR('~/dpabon/data_EFPS_ML/buffer_100m/buffer_100m.shp')
      
      
      mask_shape <- mask_shape[mask_shape$site_name == site_name[st],]
      
      #mask_shape <- extent(mask_shape)
      
      r <- rasterize(mask_shape, template_raster)
      
      masked <- mask(x = template_raster, mask = r)
      
      m.masked <- as.matrix(masked)
      m.masked <- t(m.masked)
      # image(m.masked)
      for(i in 1:ncol(m.masked)){
        for(d in 1:nrow(m.masked)){
          if(is.na(m.masked[d,i]) == T) {
            qc[d,i,] <- NA
          }
        }
      }
      
      
      
      qc_sequence <- apply(X = qc, MARGIN = c(3), FUN = useful.inf, x = useful, mask = m.masked)
      
      n.images <- length(which(qc_sequence >= qc_filter))
      images.to.select <- which(qc_sequence >= qc_filter)
      dates.to.select <- dates[which(qc_sequence >= qc_filter)]
      
      
      qc.mask <- qc[,,images.to.select]
      
      #dim(qc.mask)
      
      # Open all the variables 
      
      variables_names <- names(cube$var)
      
      # removing quality scene classification
      variables_names <- variables_names[-26]
      
      variables <- NcOpen(file_path = paste(site_name[st], "/", site_name[st], "_", tile, "_2015_2018.nc", sep = ""))
      
      
      variables <- NcToArray(variables, vars_to_read = variables_names, dim_indices = list(time = images.to.select))
      
      
      # masking non useful values for all variables all times with NA
      
      
      for (t in 1:dim(variables)[4]) {
        for (v in 1:dim(variables)[1]) {
          variables[v,,,t][-which(is.element(qc.mask[,,t], useful))] <- NA
        }
      }
      
      # estimating other VIs
      
      temp <- array(data = NA, dim = c(var =4, lon = dim(variables)[2], lat = dim(variables)[3], time = dim(variables)[4]))
      
      for(i in 1:dim(variables)[4]) {
        for(x in 1:dim(variables)[3]){
          for (d in 1:dim(variables)[2]){
            # CI red (CIr) (B7/B5) -1
            temp[1,d,x,i] <- (variables[20,d,x,i] / variables[18,d,x,i]) - 1
            # CI green (CIg) (B7/B3) -1
            temp[2,d,x,i] <- (variables[20,d,x,i] / variables[16,d,x,i]) - 1
            # ndvi
            temp[3,d,x,i] <- ((variables[21,d,x,i] - variables[17,d,x,i]) / (variables[21,d,x,i] +  variables[17,d,x,i]))
            # nirv
            temp[4,d,x,i] <- ((variables[21,d,x,i] - variables[17,d,x,i]) / (variables[21,d,x,i] +  variables[17,d,x,i])) * variables[21,d,x,i]
          }
        }
      }
      
      variables <- abind(variables, temp, along = 1)
      
      variables_names <- c(variables_names, "cir", "cig", "ndvi", "nirv")
      
      # masking forest pixels
      
      # forest <- raster("~/dpabon/to_cluster/data/euroflux/land_cover_forest/TCD_2015_020m_eu_03035_d05_Full/TCD_2015_020m_eu_03035_d05_full.tif")
      
      # crs_o <- "+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext +no_defs"
      
      # forest_proj <- projectRaster(forest, crs_o)
      
      # plot(forest)
      
      # cropping the buffer area
      # extend (min_lon, max_lon, min_lat, max_lat)
      
      # min(lat)
      # max(lat)
      # b_forest <- extent(min(lon), max(lon), min(lat), max(lat))
      # plot(b_forest)
      # crop(forest, b_forest)
      
      # outlier detection 
      
      
      depured.variables <- apply(X = variables, MARGIN = c(1,4), outlier.var)
      
      
      # estimating sd and mean for each variable - time
      
      mean.variable <- matrix(data = NA, nrow = dim(depured.variables)[1], ncol = dim(depured.variables)[2])
      sd.variable <- matrix(data = NA, nrow = dim(depured.variables)[1], ncol = dim(depured.variables)[2])
      
      for (i in 1:dim(depured.variables)[1]){
        for (j in 1:dim(depured.variables)[2]){
          mean.variable[i,j] <- mean(unlist(depured.variables[i,j]))
          sd.variable[i,j] <- sd(unlist(depured.variables[i,j]))
        }
      }
      
      
      # Genereting a dataframe with the structure:
      # Site | Time | Var1.mean | Var1.sd | ... |
      
      if(tile ==1){
        # Site | Time | Var1.mean | Var1.sd | ... |
        
        output.info.ML <- data.frame(site = rep(sites$site_name[st], length(dates.to.select)))
        
        output.info.ML[,"date"] <- dates.to.select
        
        
        for (i in 1:dim(mean.variable)[1]) {
          # mean
          output.info.ML[,paste(variables_names[i], "_mean", sep = "")] <- mean.variable[i,]
          # sd
          output.info.ML[,paste(variables_names[i], "_sd", sep = "")] <- sd.variable[i,]
        }
      }else{
        output.info.ML.temp <- data.frame(site = rep(sites$site_name[st], length(dates.to.select)))
        
        output.info.ML.temp[,"date"] <- dates.to.select
        
        
        for (i in 1:dim(mean.variable)[1]) {
          # mean
          output.info.ML.temp[,paste(variables_names[i], "_mean", sep = "")] <- mean.variable[i,]
          # sd
          output.info.ML.temp[,paste(variables_names[i], "_sd", sep = "")] <- sd.variable[i,]
        }
      }
    }
    output.info.ML <- rbind(output.info.ML, output.info.ML.temp)
    write.csv(output.info.ML, paste("~/dpabon/data_EFPS_ML/sentinel2_100m_info_v4/", sites$site_name[st], "_sentinel2_100m_info.csv", sep = "") , row.names = F)
    rm(output.info.ML)
    rm(output.info.ML.temp)
  }
  print(paste("information extracted for:", site_name[st]))
}
