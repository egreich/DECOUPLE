### Script to download and save SoilGrid2.0 data for each site coordinate

# Load site data
df_sites <- read.csv("./data_formatted/sitedata.csv")

#################################### for a single point queries using direct webdav access
# Following the tutorial found here: https://git.wur.nl/isric/soilgrids/soilgrids.notebooks/-/blob/master/markdown/xy_info_from_R.md

# Load libraries
library(tidyverse)
library(rgdal) # Please note that rgdal will be retired by the end of 2023, plan transition to sf/stars/terra functions using GDAL and PROJ
library(gdalUtils)
library(sf)

# For sand fraction
for(i in 1:nrow(df_sites)){
  print(i) # so we know where the errors are
  
  # We define the variables for the soil property and layer of interest. See here (https://www.isric.org/explore/soilgrids/faq-soilgrids#What_do_the_filename_codes_mean) for the naming conventions
  voi = "sand" # variable of interest
  depth = "0-5cm"
  layer = "mean"
  
  voi_layer = paste(voi,depth,layer, sep="_") # layer of interest 
  
  # We set other variables necessary for the WCS call for all kinds of requests
  webdav_path = '/vsicurl?max_retry=3&retry_delay=1&list_dir=no&url=https://files.isric.org/soilgrids/latest/data/'
  
  # Read in point data and re-project to Homolosine
  
  data=df_sites
  
  data <- data %>%
    mutate(longitude = lon, latitude = lat, my_id = Site) %>%
    select(longitude, latitude, my_id)
  
  data <- data[i,] # select one row
  
  spdata=st_as_sf(data,coords = c("longitude", "latitude"), crs = 4326)
  
  igh='+proj=igh +lat_0=0 +lon_0=0 +datum=WGS84 +units=m +no_defs'
  spdata_igh=st_transform(spdata, igh)
  
  data_igh=data.frame(st_coordinates(spdata_igh),id=spdata_igh$my_id)
  
  # Use gdallocationinfo to get the values from the pixels
  
  fun_pixel_values=function(rowPX,data,VOI,VOI_LYR){
    as.numeric(
      gdallocationinfo(
        srcfile=paste0(webdav_path,"/",VOI,"/", VOI_LYR,".vrt"),
        x=data[rowPX,"X"],
        y=data[rowPX,"Y"],
        geoloc=TRUE,
        valonly=TRUE))
  }
  
  value_pixels=unlist(lapply(1:3,function(x){fun_pixel_values(x,data_igh,voi,voi_layer)}))
  
  df_sites$sand[i] <- value_pixels[1]/10 # Note: we need to divide by ten to get the units in %
}

# For clay fraction
for(i in 1:nrow(df_sites)){
  print(i) # so we know where the errors are
  
  # We define the variables for the soil property and layer of interest. See here (https://www.isric.org/explore/soilgrids/faq-soilgrids#What_do_the_filename_codes_mean) for the naming conventions
  voi = "clay" # variable of interest
  depth = "0-5cm"
  layer = "mean"
  
  voi_layer = paste(voi,depth,layer, sep="_") # layer of interest 
  
  # We set other variables necessary for the WCS call for all kinds of requests
  webdav_path = '/vsicurl?max_retry=3&retry_delay=1&list_dir=no&url=https://files.isric.org/soilgrids/latest/data/'
  
  # Read in point data and re-project to Homolosine
  
  data=df_sites
  
  data <- data %>%
    mutate(longitude = lon, latitude = lat, my_id = Site) %>%
    select(longitude, latitude, my_id)
  
  data <- data[i,] # select one row
  
  spdata=st_as_sf(data,coords = c("longitude", "latitude"), crs = 4326)
  
  igh='+proj=igh +lat_0=0 +lon_0=0 +datum=WGS84 +units=m +no_defs'
  spdata_igh=st_transform(spdata, igh)
  
  data_igh=data.frame(st_coordinates(spdata_igh),id=spdata_igh$my_id)
  
  # Use gdallocationinfo to get the values from the pixels
  
  fun_pixel_values=function(rowPX,data,VOI,VOI_LYR){
    as.numeric(
      gdallocationinfo(
        srcfile=paste0(webdav_path,"/",VOI,"/", VOI_LYR,".vrt"),
        x=data[rowPX,"X"],
        y=data[rowPX,"Y"],
        geoloc=TRUE,
        valonly=TRUE))
  }
  
  value_pixels=unlist(lapply(1:3,function(x){fun_pixel_values(x,data_igh,voi,voi_layer)}))
  
  df_sites$clay[i] <- value_pixels[1]/10  # Note: we need to divide by ten to get the units in %
}


write.csv(df_sites, file = "./data_misc/soildata.csv") # save extracted data


#################################### Use the following code after 2023, not using now because terra won't update for some reason

library(geodata)
library(terra)

#sample random coordinates using terra
f <- system.file("ex/elev.tif", package="terra")
r <- rast(f)
s <- spatSample(r, 10, "random", xy=TRUE, values=FALSE)

#set coordinate system to homosline crs as SoilGrids uses this
igh <- '+proj=igh +lat_0=0 +lon_0=0 +datum=WGS84 +units=m +no_defs'
sites <- vect(s, crs = igh) 
#You may include geom = (c"x", "y") in vect() if using a data.frame of coordinates

#download data using soil_world; giving an example of 1 variable
rpH <- soil_world(var="phh2o", depth=15, stat="mean", path="Rdata/") #I created an Rdata folder in my Rproj so I can organize the data better

#extract the SoilGrids data based on the coordinate
values_pH <- terra::extract(rpH, sites)
values_pH #check if any NA values? If so, may need to add na.rm = TRUE in extract()
sites <- cbind.data.frame(crds(sites), values_pH) #binding together data and the original coordinates.








