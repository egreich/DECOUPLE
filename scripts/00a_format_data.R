### Script to format data in a consistent way with README files
### input data: Qian's compiled EC-SIF data, ICOS fluxnet archive data, Jake Nelson's LAI and TEA estimates
### output data: csv dataframe for each site with consistent columns and units

### Call functions
library(tidyverse)
library(R.matlab) # to read in matlab files
library(lubridate) # for date functions like ymd_hms()
library(chron) # for fractional DOY to time functions
library(bigleaf) # for converting LE to ET
library(ncdf4) # for reading in netcdf files
library(RNetCDF) # for "time since" conversions
source("./scripts/functions.R") # for assign_block function



### Create new dataframes site-by-site
# For variable order and units, see README in data_formatted folder


#################################################### Laegeren

# Load Laegeren data

# Data from Jake

# TEA
dname <- "T_TEA_DT" 
ncin <- nc_open("./data_raw/CH-Lae/CH-Lae_TEA.nc")

# get time
time <- ncvar_get(ncin,"time")
tunits <- ncatt_get(ncin,"time","units")
nt <- dim(time)
nt

# store the data in an array/matrix
tmp_array <- ncvar_get(ncin,dname)
dlname <- ncatt_get(ncin,dname,"long_name")
dunits <- ncatt_get(ncin,dname,"units")
fillvalue <- ncatt_get(ncin,dname,"_FillValue")
dim(tmp_array)

# replace netCDF fill values with NA's
tmp_array[tmp_array==fillvalue$value] <- NA

# reshape the array into vector
tmp_vec_long <- as.vector(tmp_array)
length(tmp_vec_long)

# reshape the vector into a matrix
tmp_mat <- matrix(tmp_vec_long, nrow=nt, ncol=1)
dim(tmp_mat)
tmp_tea <- data.frame(cbind(time,tmp_mat))
tmp_tea <- tmp_tea %>%
  rename(T_TEA = V2)
normtime<- utcal.nc(tunits$value, tmp_tea$time)
tmp_tea <- cbind(tmp_tea, normtime)
tmp_tea$minute = tmp_tea$minute - 15 # make the minutes the TIME_START minutes
tmp_tea$TIMESTAMP <- ymd_hms(paste(tmp_tea$year,tmp_tea$month,tmp_tea$day,tmp_tea$hour,tmp_tea$minute,tmp_tea$second))


# LAI
dname <- "LAI" 
ncin <- nc_open("./data_raw/CH-Lae/CH-Lae_lai.nc")

# get time
time <- ncvar_get(ncin,"time")
tunits <- ncatt_get(ncin,"time","units")
nt <- dim(time)
nt

# store the data in an array/matrix
tmp_array <- ncvar_get(ncin,dname)
dlname <- ncatt_get(ncin,dname,"long_name")
dunits <- ncatt_get(ncin,dname,"units")
fillvalue <- ncatt_get(ncin,dname,"_FillValue")
dim(tmp_array)

# replace netCDF fill values with NA's
tmp_array[tmp_array==fillvalue$value] <- NA

# reshape the array into vector
tmp_vec_long <- as.vector(tmp_array)
length(tmp_vec_long)

# reshape the vector into a matrix
tmp_mat <- matrix(tmp_vec_long, nrow=nt, ncol=1)
dim(tmp_mat)
tmp_lai <- data.frame(cbind(time,tmp_mat))
tmp_lai <- tmp_lai %>%
  rename(LAI = V2)
normtime<- utcal.nc(tunits$value, tmp_lai$time)
tmp_lai <- cbind(tmp_lai, normtime)
tmp_lai$date <- ymd(paste(tmp_lai$year,tmp_lai$month, tmp_lai$day))
tmp_lai <- tmp_lai %>%
  select(date, LAI)


# FLUXNET data
options(scipen=999) # turn off scientific notation for dates
fluxlaedatin <- read.csv("./data_raw/CH-Lae/FLX_CH-Lae_FLUXNET2015_FULLSET_HH_2004-2020_beta-3.csv")
fluxlaedatin <- fluxlaedatin %>%
  filter(TIMESTAMP_START %in% c(201801010000:201901010000)) %>% # filter for 2018
  select(TIMESTAMP_START, starts_with("TS")) # just take what we need, because it's a large file
fluxlaedatin$TIMESTAMP <-as.POSIXct(as.character(fluxlaedatin$TIMESTAMP_START), format = "%Y%m%d%H%M", 
                                    origin = "2018-01-01")
# Data from Qian
laedatin <- read.csv("./data_raw/CH-Lae/CH-Lae_dat.csv")
laedatin <- laedatin[-1,] # get rid of unit row, as these are recorded in the README
laedatin$TIMESTAMP <-as.POSIXct(laedatin$TIMESTAMP, format = "%Y%m%d %H:%M:%S")

# Combine dataframes
laedatin <- left_join(laedatin, fluxlaedatin)
laedatin <- left_join(laedatin, tmp_tea)

# Reformat
laedat <- data.frame(
  Year = laedatin$Year,
  DOY = laedatin$DOY,
  TIMESTAMP = laedatin$TIMESTAMP,
  TA = laedatin$TA_1_1_1,
  SW_IN = laedatin$SWIN_1_1_1,
  LW_IN = laedatin$LWIN_1_1_1,
  VPD = laedatin$VPD,
  RH = laedatin$RH,
  PA = laedatin$air_pressure,
  P = as.numeric(laedatin$P_RAIN_1_1_1)*1000, # convert m to mm
  WS = laedatin$wind_speed,
  PAR = laedatin$PAR,
  APAR = laedatin$APAR,
  SW_OUT = laedatin$SWOUT_1_1_1,
  LW_OUT = laedatin$LWOUT_1_1_1,
  LE = laedatin$LE,
  ET = as.numeric(laedatin$LE)*0.0352512,
  H = laedatin$H,
  NEE = laedatin$NEE,
  GPP = laedatin$GPP,
  SWC_1 = laedatin$SWC_1_1_1,
  SWC_2 = laedatin$SWC_1_2_1,
  SWC_3 = laedatin$SWC_1_3_1,
  SWC_4 = laedatin$SWC_1_4_1,
  SWC_5 = laedatin$SWC_1_5_1,
  TS_1 = laedatin$TS_F_MDS_1,
  SIF_O2A = laedatin$SIF_O2A,
  SIF_O2B = laedatin$SIF_O2B,
  NDVI = laedatin$NDVI,
  PRI = laedatin$PRI,
  Ref_650_660 = laedatin$R650_660,
  Ref_660_670 = laedatin$R660_670,
  Ref_670_680 = laedatin$R670_680,
  Ref_700_710 = laedatin$R700_710,
  Ref_770_780 = laedatin$R770_780,
  Ref_780_790 = laedatin$R780_790,
  Gs = laedatin$Gs_ms,
  T_TEA = laedatin$T_TEA
)

# Join LAI by date
laedat$date <- as.Date(laedat$TIMESTAMP)
laedat <- left_join(laedat, tmp_lai, by="date")
laedat <- laedat %>%
  select(-c(date))

laedat[,4:ncol(laedat)] <- sapply(laedat[,4:ncol(laedat)], as.numeric)

laedat[laedat==-9999] <- NA # replace -9999 with NAs

write.csv(laedat, file = "data_formatted/CH-Lae/CH-Lae_dat.csv") # save data as csv file

#################################################### Crk

# Read in Crk data
crkdatin <- read.csv("./data_raw/KR-Crk/Crk_dat.csv")
crkdatin <- crkdatin[-1,] # get rid of unit row, as these are recorded in the README
crkdatin[,c(3:58)]<- lapply(crkdatin[,c(3:58)],as.numeric) # make sure certain columns are numeric
#parse_date_time(crkdatin$Timestamp, "abdHMSzY")
#crkdatin$Timestamp <- as.POSIXct(strptime(crkdatin$Timestamp, "%a/%b/%d %H:%M:%S %z %Y"))


# Read in the MATLAB file, as some columns were corrupted
crktemp <- readMat("./data_raw/KR-Crk/CRK2016.mat")
crktemp2 <- crktemp$datr
namelist <- c(rownames(crktemp2))
crktemp3 <- lapply(crktemp2, unlist, use.names=TRUE)
names(crktemp3) <- namelist
crktemp4 <- as.data.frame(crktemp3)

crk_ref <- as.data.frame(crktemp[["datr"]][[27]])

crkdat <- data.frame(
  Year = 2016,
  DOY = crkdatin$DOY,
  TIMESTAMP = crkdatin$Timestamp,
  TA = crkdatin$T_AIR.1.,
  SW_IN = crkdatin$RSDN.1.,
  LW_IN = crkdatin$RLDN.1.,
  VPD = crktemp4$VPD,
  RH = crktemp4$RH,
  PA = as.numeric(crkdatin$Press.1.)/10, # convert hPa to kPa
  P = crkdatin$PPT.1.,
  WS = crktemp4$wind,
  NETRAD = crkdatin$RNET.1.,
  PAR = crktemp4$PAR.tot,
  APAR = crktemp4$APAR,
  SW_OUT = crkdatin$RSUP.1.,
  LW_OUT = crkdatin$RLUP.1.,
  G = as.numeric(crkdatin$G.1.),
  LE = as.numeric(crkdatin$LE),
  ET = as.numeric(crkdatin$LE)*0.0352512,
  H = as.numeric(crkdatin$H),
  NEE = as.numeric(crkdatin$NEE_MDS),
  RECO = as.numeric(crkdatin$RE_MDS),
  GPP = as.numeric(crktemp4$GPP),
  SWC_1 = (crkdatin$SWC.1. + crkdatin$SWC.3.)/2, #0-10
  SWC_2 = (crkdatin$SWC.2. + crkdatin$SWC.4.)/2, #10-20
  TS_1 = (crkdatin$T_SOIL.2. + crkdatin$T_SOIL.5.)/2, #0-10
  TS_2 = (crkdatin$T_SOIL.3. + crkdatin$T_SOIL.6.)/2, #10-20
  SIF_O2A = crktemp4$SIF760/pi, # need to divide by pi to scale to sfm SIF from other sites
  #SIF_O2B,
  #NDVI,
  #PRI,
  Ref_1 = rowMeans(crk_ref, na.rm = T),
  Gs = crkdatin$Gs_ms,
  T_TEA = NA,
  LAI = NA
)

crkdat[crkdat==-9999] <- NA # replace -9999 with NAs

# Fix TIMESTAMP
crkdat$TIMESTAMP <- mdy_hm(paste(crkdat$TIMESTAMP))

write.csv(crkdat, file = "data_formatted/Crk/Crk_dat.csv") # save data as csv file


#################################################### Gebesee

# Data from Jake

# TEA
dname <- "T_TEA_DT" 
ncin <- nc_open("./data_raw/DE-Geb/DE-Geb_TEA.nc")

# get time
time <- ncvar_get(ncin,"time")
tunits <- ncatt_get(ncin,"time","units")
nt <- dim(time)
nt

# store the data in an array/matrix
tmp_array <- ncvar_get(ncin,dname)
dlname <- ncatt_get(ncin,dname,"long_name")
dunits <- ncatt_get(ncin,dname,"units")
fillvalue <- ncatt_get(ncin,dname,"_FillValue")
dim(tmp_array)

# replace netCDF fill values with NA's
tmp_array[tmp_array==fillvalue$value] <- NA

# reshape the array into vector
tmp_vec_long <- as.vector(tmp_array)
length(tmp_vec_long)

# reshape the vector into a matrix
tmp_mat <- matrix(tmp_vec_long, nrow=nt, ncol=1)
dim(tmp_mat)
tmp_tea <- data.frame(cbind(time,tmp_mat))
tmp_tea <- tmp_tea %>%
  rename(T_TEA = V2)
normtime<- utcal.nc(tunits$value, tmp_tea$time)
tmp_tea <- cbind(tmp_tea, normtime)
tmp_tea$minute = tmp_tea$minute - 15 # make the minutes the TIME_START minutes
tmp_tea$TIMESTAMP <- ymd_hms(paste(tmp_tea$year,tmp_tea$month,tmp_tea$day,tmp_tea$hour,tmp_tea$minute,tmp_tea$second))


# LAI
dname <- "LAI" 
ncin <- nc_open("./data_raw/DE-Geb/DE-Geb_lai.nc")

# get time
time <- ncvar_get(ncin,"time")
tunits <- ncatt_get(ncin,"time","units")
nt <- dim(time)
nt

# store the data in an array/matrix
tmp_array <- ncvar_get(ncin,dname)
dlname <- ncatt_get(ncin,dname,"long_name")
dunits <- ncatt_get(ncin,dname,"units")
fillvalue <- ncatt_get(ncin,dname,"_FillValue")
dim(tmp_array)

# replace netCDF fill values with NA's
tmp_array[tmp_array==fillvalue$value] <- NA

# reshape the array into vector
tmp_vec_long <- as.vector(tmp_array)
length(tmp_vec_long)

# reshape the vector into a matrix
tmp_mat <- matrix(tmp_vec_long, nrow=nt, ncol=1)
dim(tmp_mat)
tmp_lai <- data.frame(cbind(time,tmp_mat))
tmp_lai <- tmp_lai %>%
  rename(LAI = V2)
normtime<- utcal.nc(tunits$value, tmp_lai$time)
tmp_lai <- cbind(tmp_lai, normtime)
tmp_lai$date <- ymd(paste(tmp_lai$year,tmp_lai$month, tmp_lai$day))
tmp_lai <- tmp_lai %>%
  select(date, LAI)


# FLUXNET data
fluxgebdatin <- read.csv("./data_raw/DE-Geb/FLX_DE-Geb_FLUXNET2015_FULLSET_HH_2001-2020_beta-3.csv")
fluxgebdatin <- fluxgebdatin %>%
  filter(TIMESTAMP_START %in% c(201901010000:202001010000)) %>% # filter for 2019
  select(TIMESTAMP_START, starts_with("TS")) # just take what we need, because it's a large file
fluxgebdatin$TIMESTAMP <-as.POSIXct(as.character(fluxgebdatin$TIMESTAMP_START), format = "%Y%m%d%H%M", 
                                    origin = "2019-01-01")

# Data from Qian
gebdatin <- read.csv("./data_raw/DE-Geb/DE-Geb_dat.csv")
gebdatin$TIMESTAMP <- paste(gebdatin$DATE, gebdatin$TIMESTAMP)
gebdatin$TIMESTAMP <-as.POSIXct(gebdatin$TIMESTAMP, format = "%Y%m%d %H:%M:%S")

# Combine dataframes
gebdatin <- left_join(gebdatin, fluxgebdatin)
gebdatin <- left_join(gebdatin, tmp_tea)

gebdat <- data.frame(
  Year = gebdatin$Year,
  DOY = gebdatin$DOY,
  TIMESTAMP = gebdatin$TIMESTAMP,
  TA = gebdatin$TA_1_1_1,
  SW_IN = gebdatin$SWIN_1_1_1,
  LW_IN = gebdatin$LWIN_1_1_1,
  VPD = gebdatin$VPD,
  RH = gebdatin$RH,
  PA = gebdatin$air_pressure,
  P = gebdatin$P_RAIN_1_1_1,
  WS = gebdatin$wind_speed,
  PAR = gebdatin$PAR_eddy,
  APAR = gebdatin$APAR_flox,
  SW_OUT = gebdatin$SWOUT_1_1_1,
  LW_OUT = gebdatin$LWOUT_1_1_1,
  LE = gebdatin$LE,
  ET = gebdatin$LE*0.0352512,
  H = gebdatin$H,
  NEE = gebdatin$NEE,
  GPP = gebdatin$GPP,
  SWC_1 = gebdatin$SWC_1_1_1,
  SWC_2 = gebdatin$SWC_1_2_1,
  SWC_3 = gebdatin$SWC_1_3_1,
  SWC_4 = gebdatin$SWC_1_4_1,
  SWC_5 = gebdatin$SWC_1_5_1,
  TS_1 = gebdatin$TS_F_MDS_1,
  TS_2 = gebdatin$TS_F_MDS_2,
  TS_3 = gebdatin$TS_F_MDS_3,
  TS_4 = gebdatin$TS_F_MDS_4,
  SIF_O2A = gebdatin$SIF_O2A_sfm,
  SIF_O2B = gebdatin$SIF_O2B_sfm,
  NDVI = gebdatin$NDVI,
  PRI = gebdatin$PRI,
  Ref_650_660 = gebdatin$Ref_650_660_Reflectance,
  Ref_660_670 = gebdatin$Ref_660_670_Reflectance,
  Ref_670_680 = gebdatin$Ref_670_680_Reflectance,
  Ref_700_710 = gebdatin$Ref_700_710_Reflectance,
  Ref_770_780 = gebdatin$Ref_770_780_Reflectance,
  Ref_780_790 = gebdatin$Ref_780_790_Reflectance,
  Gs = gebdatin$Gs,
  T_TEA = gebdatin$T_TEA
)

# Join LAI by date
gebdat$date <- as.Date(gebdat$TIMESTAMP)
gebdat <- left_join(gebdat, tmp_lai, by="date")
gebdat <- gebdat %>%
  select(-c(date))

gebdat[,4:ncol(gebdat)] <- sapply(gebdat[,4:ncol(gebdat)], as.numeric)

gebdat[gebdat==-9999] <- NA # replace -9999 with NAs

write.csv(gebdat, file = "data_formatted/DE-geb/DE-Geb_dat.csv") # save data as csv file

#################################################### Leinefelde

# Data from Jake - note: this only goes up to 2012 so we'll not use it for now

# TEA
dname <- "T_TEA_DT" 
ncin <- nc_open("./data_raw/DE-Lnf/DE-Lnf_TEA.nc")

# get time
time <- ncvar_get(ncin,"time")
tunits <- ncatt_get(ncin,"time","units")
nt <- dim(time)
nt

# store the data in an array/matrix
tmp_array <- ncvar_get(ncin,dname)
dlname <- ncatt_get(ncin,dname,"long_name")
dunits <- ncatt_get(ncin,dname,"units")
fillvalue <- ncatt_get(ncin,dname,"_FillValue")
dim(tmp_array)

# replace netCDF fill values with NA's
tmp_array[tmp_array==fillvalue$value] <- NA

# reshape the array into vector
tmp_vec_long <- as.vector(tmp_array)
length(tmp_vec_long)

# reshape the vector into a matrix
tmp_mat <- matrix(tmp_vec_long, nrow=nt, ncol=1)
dim(tmp_mat)
tmp_tea <- data.frame(cbind(time,tmp_mat))
tmp_tea <- tmp_tea %>%
  rename(T_TEA = V2)
normtime<- utcal.nc(tunits$value, tmp_tea$time)
tmp_tea <- cbind(tmp_tea, normtime)
tmp_tea$minute = tmp_tea$minute - 15 # make the minutes the TIME_START minutes
tmp_tea$TIMESTAMP <- ymd_hms(paste(tmp_tea$year,tmp_tea$month,tmp_tea$day,tmp_tea$hour,tmp_tea$minute,tmp_tea$second))


# LAI
dname <- "LAI" 
ncin <- nc_open("./data_raw/DE-Lnf/DE-Lnf_lai.nc")

# get time
time <- ncvar_get(ncin,"time")
tunits <- ncatt_get(ncin,"time","units")
nt <- dim(time)
nt

# store the data in an array/matrix
tmp_array <- ncvar_get(ncin,dname)
dlname <- ncatt_get(ncin,dname,"long_name")
dunits <- ncatt_get(ncin,dname,"units")
fillvalue <- ncatt_get(ncin,dname,"_FillValue")
dim(tmp_array)

# replace netCDF fill values with NA's
tmp_array[tmp_array==fillvalue$value] <- NA

# reshape the array into vector
tmp_vec_long <- as.vector(tmp_array)
length(tmp_vec_long)

# reshape the vector into a matrix
tmp_mat <- matrix(tmp_vec_long, nrow=nt, ncol=1)
dim(tmp_mat)
tmp_lai <- data.frame(cbind(time,tmp_mat))
tmp_lai <- tmp_lai %>%
  rename(LAI = V2)
normtime<- utcal.nc(tunits$value, tmp_lai$time)
tmp_lai <- cbind(tmp_lai, normtime)
tmp_lai$date <- ymd(paste(tmp_lai$year,tmp_lai$month, tmp_lai$day))
tmp_lai <- tmp_lai %>%
  select(date, LAI)


# Load Leinefelde data
lnfdatin_full <- read.csv("./data_raw/DE-Lnf/DE-Lnf_dat_full.csv")
lnfdatin <- read.csv("./data_raw/DE-Lnf/DE-Lnf_dat.csv")

lnfdat_cropped <- data.frame(
  Year = 2019,
  DOY = lnfdatin$DOY,
  date = lnfdatin$datetime..UTC.,
  time = lnfdatin$time,
  #TIMESTAMP = lnfdatin$TIMESTAMP,
  TA = lnfdatin$Tair_fall,
  SW_IN = lnfdatin$SW_IN_SPN1,
  #LW_IN = lnfdatin,
  VPD = lnfdatin$VPD,
  RH = lnfdatin$RH,
  #PA = lnfdatin,
  #P = lnfdatin$P44,
  #WS = lnfdatin,
  #NETRAD
  PAR = lnfdatin$PPFD_IN,
  #APAR = lnfdatin,
  #SW_OUT = lnfdatin$SW,
  LW_OUT = lnfdatin$LW_OUT_CNR1,
  #CO2
  #G
  #LE = lnfdatin,
  ET = lnfdatin$ET,
  #H = lnfdatin,
  #NEE = lnfdatin,
  #RECO,
  GPP = lnfdatin$GPP_MR_f,
  SWC_1 = lnfdatin$SWCa008,
  SWC_2 = lnfdatin$SWCa016,
  SWC_3 = lnfdatin$SWCa032,
  SWC_4 = lnfdatin$SWCa062,
  #TS,
  SIF_O2A = lnfdatin$SIF_A_sfm..mW.m.2nm.1sr.1.,
  SIF_O2B = lnfdatin$SIF_B_sfm..mW.m.2nm.1sr.1.,
  NDVI = lnfdatin$NDVI,
  PRI = lnfdatin$PRI,
  Ref_750 = lnfdatin$Reflectance.750....,
  Ref_760 = lnfdatin$Reflectance.760....,
  Gs = lnfdatin$Gs_ms
)

#lnfdat_cropped$DOY <- format(round(lnfdat_cropped$DOY, 3), nsmall = 3) # crop to match other dataframe
#lnfdat_cropped$DOY <- as.character(lnfdat_cropped$DOY)


lnfdatin <- lnfdatin_full


lnfdat_full <- data.frame(
  Year = 2019,
  DOY = lnfdatin$DOY,
  date = lnfdatin$date,
  time = lnfdatin$time,
  #TIMESTAMP = lnfdatin$TIMESTAMP,
  TA = lnfdatin$Tair_fall,
  SW_IN = lnfdatin$SW_IN_SPN1,
  VPD = lnfdatin$VPD,
  RH = lnfdatin$RH,
  PA = lnfdatin$air_pressure/1000, # convert to kPa
  P = lnfdatin$P44,
  WS = lnfdatin$WS,
  PAR = lnfdatin$PPFD_IN,
  LW_OUT = lnfdatin$LW_OUT_CNR1,
  ET = lnfdatin$ET,
  GPP = lnfdatin$GPP_MR_f,
  SWC_1 = lnfdatin$SWCa008,
  SWC_2 = lnfdatin$SWCa016,
  SWC_3 = lnfdatin$SWCa032,
  SWC_4 = lnfdatin$SWCa062,
  TS_1 = lnfdatin$TSa008,
  TS_2 = lnfdatin$TSa016,
  TS_3 = lnfdatin$TSa032,
  TS_4 = lnfdatin$TSa064,
  Gs = lnfdatin$Gs_ms
)


# Fix dates
lnfdat_full$TIMESTAMP = mdy_hm(paste(lnfdat_full$date, lnfdat_full$time))
lnfdat_cropped$TIMESTAMP = mdy_hm(paste(lnfdat_cropped$date, lnfdat_cropped$time))

# Join to include met info from the full day
lnfdat <- left_join(lnfdat_full, lnfdat_cropped, by = c("TIMESTAMP"))

# Pick which columns to include
lnfdat <- data.frame(
  Year = 2019,
  DOY = lnfdat$DOY.x,
  TIMESTAMP = lnfdat$TIMESTAMP,
  TA = lnfdat$TA.x,
  SW_IN = lnfdat$SW_IN.x,
  VPD = lnfdat$VPD.x,
  RH = lnfdat$RH.x,
  P = lnfdat$P,
  WS = lnfdat$WS,
  PAR = lnfdat$PAR.x,
  LW_OUT = lnfdat$LW_OUT.x,
  ET = lnfdat$ET.x,
  GPP = lnfdat$GPP.x,
  SWC_1 = lnfdat$SWC_1.x,
  SWC_2 = lnfdat$SWC_2.x,
  SWC_3 = lnfdat$SWC_3.x,
  SWC_4 = lnfdat$SWC_4.x,
  TS_1 = lnfdat$TS_1,
  TS_2 = lnfdat$TS_2,
  TS_3 = lnfdat$TS_3,
  TS_4 = lnfdat$TS_4,
  SIF_O2A = lnfdat$SIF_O2A,
  SIF_O2B = lnfdat$SIF_O2B,
  NDVI = lnfdat$NDVI,
  PRI = lnfdat$PRI,
  Ref_750 = lnfdat$Ref_750,
  Ref_760 = lnfdat$Ref_760,
  Gs = lnfdat$Gs.x,
  T_TEA = NA,
  LAI = NA
)

lnfdat[lnfdat==-9999] <- NA # replace -9999 with NAs

write.csv(lnfdat, file = "data_formatted/DE-Lnf/DE-Lnf_dat.csv") # save data as csv file

#################################################### Yat

# Data from Jake

# TEA
dname <- "T_TEA_DT" 
ncin <- nc_open("./data_raw/IL-Yat/IL-Yat_TEA.nc")

# get time
time <- ncvar_get(ncin,"time")
tunits <- ncatt_get(ncin,"time","units")
nt <- dim(time)
nt

# store the data in an array/matrix
tmp_array <- ncvar_get(ncin,dname)
dlname <- ncatt_get(ncin,dname,"long_name")
dunits <- ncatt_get(ncin,dname,"units")
fillvalue <- ncatt_get(ncin,dname,"_FillValue")
dim(tmp_array)

# replace netCDF fill values with NA's
tmp_array[tmp_array==fillvalue$value] <- NA

# reshape the array into vector
tmp_vec_long <- as.vector(tmp_array)
length(tmp_vec_long)

# reshape the vector into a matrix
tmp_mat <- matrix(tmp_vec_long, nrow=nt, ncol=1)
dim(tmp_mat)
tmp_tea <- data.frame(cbind(time,tmp_mat))
tmp_tea <- tmp_tea %>%
  rename(T_TEA = V2)
normtime<- utcal.nc(tunits$value, tmp_tea$time)
tmp_tea <- cbind(tmp_tea, normtime)
tmp_tea$minute = tmp_tea$minute - 15 # make the minutes the TIME_START minutes
tmp_tea$TIMESTAMP <- ymd_hms(paste(tmp_tea$year,tmp_tea$month,tmp_tea$day,tmp_tea$hour,tmp_tea$minute,tmp_tea$second))


# LAI
dname <- "LAI" 
ncin <- nc_open("./data_raw/IL-Yat/IL-Yat_lai.nc")

# get time
time <- ncvar_get(ncin,"time")
tunits <- ncatt_get(ncin,"time","units")
nt <- dim(time)
nt

# store the data in an array/matrix
tmp_array <- ncvar_get(ncin,dname)
dlname <- ncatt_get(ncin,dname,"long_name")
dunits <- ncatt_get(ncin,dname,"units")
fillvalue <- ncatt_get(ncin,dname,"_FillValue")
dim(tmp_array)

# replace netCDF fill values with NA's
tmp_array[tmp_array==fillvalue$value] <- NA

# reshape the array into vector
tmp_vec_long <- as.vector(tmp_array)
length(tmp_vec_long)

# reshape the vector into a matrix
tmp_mat <- matrix(tmp_vec_long, nrow=nt, ncol=1)
dim(tmp_mat)
tmp_lai <- data.frame(cbind(time,tmp_mat))
tmp_lai <- tmp_lai %>%
  rename(LAI = V2)
normtime<- utcal.nc(tunits$value, tmp_lai$time)
tmp_lai <- cbind(tmp_lai, normtime)
tmp_lai$date <- ymd(paste(tmp_lai$year,tmp_lai$month, tmp_lai$day))
tmp_lai <- tmp_lai %>%
  select(date, LAI)


# Read in Yatir data from FLUXNET
fluxyatdatin <- read.csv("./data_raw/IL-Yat/FLX_IL-Yat_FLUXNET2015_FULLSET_HH_2000-2020_beta-3.csv")
fluxyatdatin <- fluxyatdatin %>%
  filter(TIMESTAMP_START %in% c(201701010000:201801010000)) %>% # filter for 2017
  select(TIMESTAMP_START, starts_with("TS")) # just take what we need, because it's a large file
fluxyatdatin$TIMESTAMP <-as.POSIXct(as.character(fluxyatdatin$TIMESTAMP_START), format = "%Y%m%d%H%M", 
                                    origin = "2017-01-01")

# Read in data from Qian
yatdatin_full <- read.csv("./data_raw/IL-Yat/IL-Yat_full.csv")
yatdatin <- read.csv("./data_raw/IL-Yat/IL-Yat_dat.csv") # cropped data
yatdatin_full <- yatdatin_full[-1,] # get rid of unit row, as these are recorded in the README

# Fix dates
yatdatin_full$date <- yatdatin_full$TIMESTAMP
yatdatin_full$time <- yatdatin_full$TIMESTAMP.2
yatdatin_full$TIMESTAMP = mdy_hm(paste(yatdatin_full$date, yatdatin_full$time))
yatdatin_full$DOY = difftime(yatdatin_full$TIMESTAMP,as.POSIXct(as.Date("2017-01-01 00:00")))
yatdatin_full$DOY<-gsub(" days","",yatdatin_full$DOY)
yatdatin_full$TIMESTAMP <-as.POSIXct(yatdatin_full$TIMESTAMP, format = "%Y%m%d %H:%M:%S")

# Combine dataframes
yatdatin_full <- left_join(yatdatin_full, fluxyatdatin)
yatdatin_full <- left_join(yatdatin_full, tmp_tea, by = "TIMESTAMP")

yatdat_full <- data.frame(
  Year = yatdatin_full$Year,
  DOY = yatdatin_full$DOY,
  TIMESTAMP = yatdatin_full$TIMESTAMP,
  TA = yatdatin_full$TA_1_1_1,
  SW_IN = yatdatin_full$SWIN_1_1_1,
  LW_IN = yatdatin_full$LWIN_1_1_1,
  VPD = as.numeric(yatdatin_full$VPD)/100, # Pa to kPa
  RH = yatdatin_full$RH,
  PA = yatdatin_full$air_pressure,
  P = yatdatin_full$P_RAIN_1_1_1,
  WS = yatdatin_full$wind_speed,
  PAR = yatdatin_full$PAR,
  APAR = yatdatin_full$APAR,
  SW_OUT = yatdatin_full$SWOUT_1_1_1,
  LW_OUT = yatdatin_full$LWOUT_1_1_1,
  LE = yatdatin_full$LE,
  ET = as.numeric(yatdatin_full$LE)*0.0352512,
  H = yatdatin_full$H,
  NEE = yatdatin_full$NEE,
  GPP = as.numeric(yatdatin_full$GPP)*-1, # fix sign inconcsistency
  SWC_1 = yatdatin_full$SWC_1_1_1,
  SWC_2 = yatdatin_full$SWC_1_2_1,
  SWC_3 = yatdatin_full$SWC_1_3_1,
  SWC_4 = yatdatin_full$SWC_1_4_1,
  TS_1 = yatdatin_full$TS_F_MDS_1,
  TS_2 = yatdatin_full$TS_F_MDS_2,
  T_TEA = yatdatin_full$T_TEA
  #SIF_O2A = yatdatin_full$SIF_O2A,
  #SIF_O2B = yatdatin_full$SIF_O2B,
  #NDVI = yatdatin_full$NDVI,
  #PRI = yatdatin_full$PRI#,
  #Ref_650_660 = yatdatin_full$Ref_650_660_Reflectance,
  #Ref_660_670 = yatdatin_full$Ref_660_670_Reflectance,
  #Ref_670_680 = yatdatin_full$Ref_670_680_Reflectance,
  #Ref_700_710 = yatdatin_full$Ref_700_710_Reflectance,
  #Ref_770_780 = yatdatin_full$Ref_770_780_Reflectance,
  #Ref_780_790 = yatdatin_full$Ref_780_790_Reflectance,
  #Gs = yatdatin_full$Gs_ms
)



yatdat_cropped <- data.frame(
  Year = 2017,
  DOY = yatdatin$doy.dayfract,
  TIMESTAMP = yatdatin$TIMESTAMP,
  SIF_O2A = yatdatin$SIF_A_sfm_corrected,
  SIF_O2B = yatdatin$SIF_B_sfm_corrected,
  NDVI = yatdatin$NDVI,
  PRI = yatdatin$PRI,
  Ref_687 = yatdatin$Reflected.687,
  Ref_760 = yatdatin$Reflected.760
  #Gs
)
# Fix dates
yatdat_cropped$time <- paste(yatdat_cropped$DOY %/% 1, times(yatdat_cropped$DOY %% 1))
yatdat_cropped$time <- sub("^\\S+\\s+", '', yatdat_cropped$time)
yatdat_cropped$TIMESTAMP = ymd_hms(paste(yatdat_cropped$TIMESTAMP, yatdat_cropped$time))

yatdat_cropped <- yatdat_cropped %>%
  select(-c(Year,DOY,time))


# Join to include met info from the full day
yatdat <- left_join(yatdat_full, yatdat_cropped, by = "TIMESTAMP")

# Reorder columns
yatdat <- yatdat %>%
  relocate(Year, DOY, TIMESTAMP)

# Join LAI by date
yatdat$date <- as.Date(yatdat$TIMESTAMP)
yatdat <- left_join(yatdat, tmp_lai, by="date")
yatdat <- yatdat %>%
  select(-c(date))

yatdat[yatdat==-9999] <- NA # replace -9999 with NAs

write.csv(yatdat, file = "data_formatted/IL-yat/IL-yat_dat.csv") # save data as csv file


##################################################### Jrs

# Load JRS data
jrsdatin <- read.csv("./data_raw/Jrs/Jrs_dat.csv")

jrsdatin <- jrsdatin[-1,] # get rid of unit row, as these are recorded in the README

jrsdat <- data.frame(
  Year = 2019,
  DOY = jrsdatin$DOY,
  #TIMESTAMP = NA,
  TA = (as.numeric(jrsdatin$TA_1_1_1) + as.numeric(jrsdatin$TA_1_2_1))/2,
  SW_IN = jrsdatin$SWIN_1_1_1,
  LW_IN = jrsdatin$LWIN_1_1_1,
  VPD = jrsdatin$VPD,
  RH = jrsdatin$RH,
  PA = as.numeric(jrsdatin$air_pressure)/100, # convert Pa to kPa
  P = as.numeric(jrsdatin$P_RAIN_1_1_1)*1000,
  WS = jrsdatin$wind_speed,
  #PAR = jrsdatin,
  SW_OUT = jrsdatin$SWOUT_1_1_1,
  LW_OUT = jrsdatin$LWOUT_1_1_1,
  LE = jrsdatin$LE,
  ET = jrsdatin$ET,
  H = jrsdatin$H,
  NEE = jrsdatin$NEE,
  GPP = jrsdatin$GPP_gapfilled,
  SWC_1 = jrsdatin$SWC_1_1_1,
  SWC_2 = jrsdatin$SWC_1_2_1,
  SWC_3 = jrsdatin$SWC_1_3_1,
  SWC_4 = jrsdatin$SWC_1_4_1,
  SWC_5 = jrsdatin$SWC_1_5_1,
  TS_1 = jrsdatin$TS_1_1_1,
  TS_2 = jrsdatin$TS_1_2_1,
  TS_3 = jrsdatin$TS_1_3_1,
  TS_4 = jrsdatin$TS_1_4_1,
  TS_5 = jrsdatin$TS_1_5_1,
  SIF_O2A = jrsdatin$SFM_A_linear,
  #SIF_O2B = jrsdatin$SIF_O2B,
  NDVI = jrsdatin$ndvi,
  PRI = jrsdatin$pri,
  Ref_740_750 = jrsdatin$Ref_740_750_Reflectance,
  Ref_750_760 = jrsdatin$Ref_750_760_Reflectance,
  Ref_760_770 = jrsdatin$Ref_760_770_Reflectance,
  Ref_770_780 = jrsdatin$Ref_770_780_Reflectance,
  Ref_780_790 = jrsdatin$Ref_780_790_Reflectance,
  Gs = jrsdatin$gs.ms,
  T_TEA = NA,
  LAI = NA
)

# Fix dates
jrsdat$time <- paste(as.numeric(jrsdat$DOY) %/% 1, times(as.numeric(jrsdat$DOY) %% 1))
jrsdat$time <- sub("^\\S+\\s+", '', jrsdat$time)
jrsdat$date <-as.Date(as.numeric(jrsdat$DOY), origin = "2019-01-01")
jrsdat$TIMESTAMP = ymd_hms(paste(jrsdat$date, jrsdat$time))

# Pick which columns to include
jrsdat <- jrsdat %>%
  select(-c(time,date)) %>%
  select(Year,DOY,TIMESTAMP, everything())

jrsdat[jrsdat==-9999] <- NA # replace -9999 with NAs

write.csv(jrsdat, file = "data_formatted/Jrs/Jrs_dat.csv") # save data as csv file



##################################################### Sq

# Load sq data
sqdatin17 <- read.csv("./data_raw/Sq/Sq_2017_dat.csv")
sqdatin17 <- sqdatin17%>%
     filter(!is.na(Year))
sqdatin18 <- read.csv("./data_raw/Sq/Sq_2018_dat.csv")
sqdatin18 <- sqdatin18%>%
  filter(!is.na(Year))
sqdatin19 <- read.csv("./data_raw/Sq/Sq_2019_dat.csv")
sqdatin19 <- sqdatin19%>%
  filter(!is.na(Year))
sqdatin <- rbind(sqdatin17, sqdatin18, sqdatin19)


sqdatin <- sqdatin[-1,] # get rid of unit row, as these are recorded in the README

sqdat <- data.frame(
  Year = sqdatin$Year,
  DOY = sqdatin$Doy,
  TIMESTAMP = NA,
  TA = sqdatin$TA_1_1_1,
  SW_IN = sqdatin$SWIN_1_1_1,
  LW_IN = sqdatin$LWIN_1_1_1,
  VPD = sqdatin$VPD,
  RH = sqdatin$RH,
  PA = sqdatin$air_pressure, # convert Pa to kPa
  P = as.numeric(sqdatin$P_RAIN_1_1_1)*1000,
  WS = sqdatin$wind_speed,
  PAR = sqdatin$PAR,
  APAR = sqdatin$APAR,
  SW_OUT = sqdatin$SWOUT_1_1_1,
  LW_OUT = sqdatin$LWOUT_1_1_1,
  LE = sqdatin$LE,
  ET = sqdatin$ET,
  H = sqdatin$H,
  NEE = sqdatin$NEE,
  GPP = sqdatin$GPP,
  SWC_1_1_1 = sqdatin$SWC_1_1_1,
  SWC_1_2_1 = sqdatin$SWC_1_2_1,
  SWC_1_3_1 = sqdatin$SWC_1_3_1,
  SIF_O2A = sqdatin$SIF_O2A,
  NDVI = sqdatin$NDVI,
  PRI = sqdatin$PRI,
  Ref_650_660 = sqdatin$Ref_650_660_Reflectance,
  Ref_660_670 = sqdatin$Ref_660_670_Reflectance,
  Ref_670_680 = sqdatin$Ref_670_680_Reflectance,
  Ref_700_710 = sqdatin$Ref_700_710_Reflectance,
  Ref_770_780 = sqdatin$Ref_770_780_Reflectance,
  Ref_780_790 = sqdatin$Ref_780_790_Reflectance,
  Gs = sqdatin$Gs_ms
)

# Fix dates
sqdat$time <- paste(as.numeric(sqdat$DOY) %/% 1, times(as.numeric(sqdat$DOY) %% 1))
sqdat$time <- sub("^\\S+\\s+", '', sqdat$time)
sqdat17 <- sqdat %>%
  filter(Year == "2017")
sqdat17$date <-as.Date(as.numeric(sqdat17$DOY), origin = "2017-01-01")
sqdat18 <- sqdat %>%
  filter(Year == "2018")
sqdat18$date <-as.Date(as.numeric(sqdat18$DOY), origin = "2018-01-01")
sqdat19 <- sqdat %>%
  filter(Year == "2019")
sqdat19$date <-as.Date(as.numeric(sqdat19$DOY), origin = "2019-01-01")
sqdat <- rbind(sqdat17, sqdat18, sqdat19)

sqdat$TIMESTAMP = ymd_hms(paste(sqdat$date, sqdat$time))
sqdat$TIMESTAMP <- round_date(sqdat$TIMESTAMP, "30 mins") # fix corrupted DOY column

# Pick which columns to include
sqdat <- sqdat %>%
  select(-c(time,date))

sqdat[sqdat==-9999] <- NA # replace -9999 with NAs

write.csv(sqdat, file = "data_formatted/sq/sq_dat.csv") # save data as csv file










