setwd("~/dpabon/to_cluster/")
library(ncdf4)
#library(tidyverse)
library(raster)
library("easyNCDF")
library(lubridate)
source("dates_extraction.R")
# sites names 

sites <- c("Be-Bra",
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
           "US-Wjs"
            )
sites


  
  
for(st in 1:length(sites)) {
  system2(command = "curl",args =  paste("-s -X POST $URL -d chat_id=$CHAT_ID -d text='Starting cube generation for: ", sites[st], "'",  sep = ""))

  dates <- dates_extraction(path = paste("~/l1_work3/cube/", sites[st],"/", sep = ""), pattern = "*.nc")
  files <- list.files(path = paste("~/l1_work3/cube/", sites[st], "/", sep = ""), pattern = "*.nc",full.names = T)

  nc <- nc_open(files[1])

  lat <- ncvar_get(nc, varid = "lat")
  lon <- ncvar_get(nc, varid = "lon")
  ncvar_get(nc, varid = "crs")
  var <- NcReadVarNames(files[1])

  var <- var[c(seq(2,27, by = 2),28:39,44,47,48,seq(75,106, by = 2))]

  test <- NcToArray(file_to_read = files[1], vars_to_read = var, unlist = T)

  print(dim(test))

  # defined dimensions
  londim <- ncdim_def("lon","degrees_east",as.double(lon))
  latdim <- ncdim_def("lat","degrees_north",as.double(lat))
  timedim <- ncdim_def("time","days of the year since 1970-01-01",as.numeric(dates))


  # defined variables

  output <- array(data = NA, dim = c(44, length(lon), length(lat), length(dates)))
  #, dimnames = list(as.character(dates), var,lon,lat)
  names(dim(output)) <- c("var", "lon", "lat", "time")

  metadata <- list(
    bi2 = list(name = "bi2"),
    bi = list(name = "bi"),
    ci = list(name = "ci"),
    ri = list(name = "ri"),
    arvi = list(name = "arvi"),
    dvi = list(name = "dvi"),
    gemi = list(name = "gemi"),
    gndvi = list(name = "gndvi"),
    ipvi = list(name = "ipvi"),
    ireci = list(name = "ireci"),
    mcari = list(name = "mcari"),
    msavi2 = list(name = "msavi2"),
    msavi = list(name = "msavi"),
    B1 = list(name = "B1"),
    B2 = list(name = "B2"),
    B3 = list(name = "B3"),
    B4 = list(name = "B4"),
    B5 = list(name = "B5"),
    B6 = list(name = "B6"),
    B7 = list(name = "B7"),
    B8 = list(name = "B8"),
    B8A = list(name = "B8A"),
    B9 = list(name = "B9"),
    B11 = list(name = "B11"),
    B12 = list(name = "B12"),
    quality_scene_classification = list(name = "quality_scene_classification"),
    sun_zenith = list(name = "sun_zenith"),
    sun_azimuth = list(name = "sun_azimuth"),
    mtci = list(name = "mtci"),
    ndi45 = list(name = "ndi45"),
    pssra = list(name = "pssra"),
    pvi = list(name = "pvi"),
    reip = list(name = "reip"),
    rvi = list(name = "rvi"),
    s2rep = list(name = "s2rep"),
    savi = list(name = "savi"),
    tndvi = list(name = "tndvi"),
    tsavi = list(name = "tsavi"),
    wdvi = list(name = "wdvi"),
    lai = list(name = "lai"),
    lai_cab = list(name = "lai_cab"),
    lai_cw = list(name = "lai_cw"),
    fapar = list(name = "fapar"),
    fcover = list(name = "fcover")
  )
  attr(output, 'variables') <- metadata


  dim(lon) <- length(lon)
  metadata <- list(lon = list(units = 'degrees_east'))
  attr(lon, 'variables') <- metadata
  names(dim(lon)) <- 'lon'


  dim(lat) <- length(lat)
  metadata <- list(lat = list(units = 'degrees_north'))
  attr(lat, 'variables') <- metadata
  names(dim(lat)) <- 'lat'

  time <- as.numeric(dates)
  dim(time) <- length(time)
  metadata <- list(time = list(units = 'day of the year since 1970-01-01'))
  attr(time, 'variables') <- metadata
  names(dim(time)) <- 'time'

  output[,,,1] <- test


  #image(output[1, 16,,])

  for (i in 2:length(dates)) {
    test <- NcToArray(file_to_read = files[i], vars_to_read = var, unlist = T)
    print(dim(test))
    if(all(dim(test) == dim(output[,,,i])) == FALSE) {
      print(paste("Problem with", files[i]))
    }else{
      output[,,,i] <- test
      print(paste(files[i], "added to the datacube"))
    }
  }

  dim(output)

  ArrayToNc(list(output,lon,lat, time), paste("~/l1_work3/cube_f/", sites[st], "/", sites[st], "_2015_2018.nc",
                                              sep = ""))

  system2(command = "curl",args =  paste("-s -X POST $URL -d chat_id=$CHAT_ID -d text='Cube generated for ", sites[st], "'",  sep = ""))
}

### Multiple tiles

sites_multiple <- read.csv("~/dpabon/data_EFPS_ML/sites_2_tiles.csv")

sites_multiple <- sites_multiple[-c(1,2,3,6,7,8,9,10,11,13),]

for(st in 1:nrow(sites_multiple)) {
  for(tile in 1:2){
    system2(command = "curl",args =  paste("-s -X POST $URL -d chat_id=$CHAT_ID -d text='Starting cube generation for: ", sites_multiple[st,1], "'",  sep = ""))
    dates <- dates_extraction(path = paste("~/l1_work3/cube/", sites_multiple[st,1],"/", sep = ""), pattern = paste("*",sites_multiple[st,tile+1], "*.nc", sep = ""))
    files <- list.files(path = paste("~/l1_work3/cube/", sites_multiple[st,1], "/", sep = ""), pattern = glob2rx(paste("*",sites_multiple[st,tile+1], "*.nc", sep = "")),full.names = T)
    
    nc <- nc_open(files[1])
    
    lat <- ncvar_get(nc, varid = "lat")
    lon <- ncvar_get(nc, varid = "lon")
    ncvar_get(nc, varid = "crs")
    var <- NcReadVarNames(files[1])
    
    var <- var[c(seq(2,27, by = 2),28:39,44,47,48,seq(75,106, by = 2))]
    
    test <- NcToArray(file_to_read = files[1], vars_to_read = var, unlist = T)
    
    print(dim(test))
    
    # defined dimensions
    londim <- ncdim_def("lon","degrees_east",as.double(lon)) 
    latdim <- ncdim_def("lat","degrees_north",as.double(lat)) 
    timedim <- ncdim_def("time","days of the year since 1970-01-01",as.numeric(dates))
    
    
    # defined variables
    
    output <- array(data = NA, dim = c(44, length(lon), length(lat), length(dates)))
    #, dimnames = list(as.character(dates), var,lon,lat)
    names(dim(output)) <- c("var", "lon", "lat", "time")
    
    metadata <- list(
      bi2 = list(name = "bi2"),
      bi = list(name = "bi"),
      ci = list(name = "ci"),
      ri = list(name = "ri"),
      arvi = list(name = "arvi"),
      dvi = list(name = "dvi"),
      gemi = list(name = "gemi"),
      gndvi = list(name = "gndvi"),
      ipvi = list(name = "ipvi"),
      ireci = list(name = "ireci"),
      mcari = list(name = "mcari"),
      msavi2 = list(name = "msavi2"),
      msavi = list(name = "msavi"),
      B1 = list(name = "B1"),
      B2 = list(name = "B2"),
      B3 = list(name = "B3"),
      B4 = list(name = "B4"),
      B5 = list(name = "B5"),
      B6 = list(name = "B6"),
      B7 = list(name = "B7"),
      B8 = list(name = "B8"),
      B8A = list(name = "B8A"),
      B9 = list(name = "B9"),
      B11 = list(name = "B11"),
      B12 = list(name = "B12"),
      quality_scene_classification = list(name = "quality_scene_classification"),
      sun_zenith = list(name = "sun_zenith"),
      sun_azimuth = list(name = "sun_azimuth"),
      mtci = list(name = "mtci"),
      ndi45 = list(name = "ndi45"),
      pssra = list(name = "pssra"),
      pvi = list(name = "pvi"),
      reip = list(name = "reip"),
      rvi = list(name = "rvi"),
      s2rep = list(name = "s2rep"),
      savi = list(name = "savi"),
      tndvi = list(name = "tndvi"),
      tsavi = list(name = "tsavi"),
      wdvi = list(name = "wdvi"),
      lai = list(name = "lai"),
      lai_cab = list(name = "lai_cab"),
      lai_cw = list(name = "lai_cw"),
      fapar = list(name = "fapar"),
      fcover = list(name = "fcover")
    )
    attr(output, 'variables') <- metadata
    
    
    dim(lon) <- length(lon)
    metadata <- list(lon = list(units = 'degrees_east'))
    attr(lon, 'variables') <- metadata
    names(dim(lon)) <- 'lon'
    
    
    dim(lat) <- length(lat)
    metadata <- list(lat = list(units = 'degrees_north'))
    attr(lat, 'variables') <- metadata
    names(dim(lat)) <- 'lat'
    
    time <- as.numeric(dates)
    dim(time) <- length(time)
    metadata <- list(time = list(units = 'day of the year since 1970-01-01'))
    attr(time, 'variables') <- metadata
    names(dim(time)) <- 'time'
    
    output[,,,1] <- test
    
    
    #image(output[1, 16,,])
    
    for (i in 2:length(dates)) {
      test <- NcToArray(file_to_read = files[i], vars_to_read = var, unlist = T)
      print(dim(test))
      if(all(dim(test) == dim(output[,,,i])) == FALSE) {
        print(paste("Problem with", files[i]))
      }else{
        output[,,,i] <- test
        print(paste(files[i], "added to the datacube"))
      }
    }
    
    dim(output)
    
    ArrayToNc(list(output,lon,lat, time), paste("~/l1_work3/cube_f/", sites_multiple[st,1], "/", sites_multiple[st,1], "_", tile, "_2015_2018.nc",
                                                sep = ""))
    
    system2(command = "curl",args =  paste("-s -X POST $URL -d chat_id=$CHAT_ID -d text='Cube generated for ", sites_multiple[st,1], "'",  sep = ""))
  }
  
}

