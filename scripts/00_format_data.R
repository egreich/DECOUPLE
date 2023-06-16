### Script to format data in a consistent way with README files
### input data: Qian's collected data
### output data: csv dataframe for each site with consistent columns and units

### Call functions
library(tidyverse)
library(R.matlab) # to read in matlab files
library(lubridate) # for ymd_hms() functions
library(chron) # for fractional DOY to time functions
source("./scripts/functions.R") # for assign_block function



### Create new dataframes site-by-site
# For variable order and units, see README in data_formatted folder


#################################################### Laegeren

# Load Laegeren data
laedatin <- read.csv("./data_raw/CH-Lae/CH-Lae_dat.csv")

laedatin <- laedatin[-1,] # get rid of unit row, as these are recorded in the README

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
  P = laedatin$P_RAIN_1_1_1,
  WS = laedatin$wind_speed,
  #NETRAD
  PAR = laedatin$PAR,
  APAR = laedatin$APAR,
  SW_OUT = laedatin$SWOUT_1_1_1,
  LW_OUT = laedatin$LWOUT_1_1_1,
  #CO2
  #G
  LE = laedatin$LE,
  ET = as.numeric(laedatin$LE)*0.0352512,
  H = laedatin$H,
  NEE = laedatin$NEE,
  #RECO,
  GPP = laedatin$GPP,
  SWC_1_1_1 = laedatin$SWC_1_1_1,
  SWC_1_2_1 = laedatin$SWC_1_2_1,
  SWC_1_3_1 = laedatin$SWC_1_3_1,
  SWC_1_4_1 = laedatin$SWC_1_4_1,
  SWC_1_5_1 = laedatin$SWC_1_5_1,
  #TS,
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
  Gs = laedatin$Gs_ms
)

laedat[laedat==-9999] <- NA # replace -9999 with NAs

write.csv(laedat, file = "data_formatted/CH-Lae/CH-Lae_dat.csv") # save data as csv file

#################################################### Crk

# Read in Crk data
crkdatin <- read.csv("./data_raw/Crk/Crk_dat.csv")
crkdatin <- crkdatin[-1,] # get rid of unit row, as these are recorded in the README
crkdatin[,c(3:58)]<- lapply(crkdatin[,c(3:58)],as.numeric) # make sure certain columns are numeric
#parse_date_time(crkdatin$Timestamp, "abdHMSzY")
#crkdatin$Timestamp <- as.POSIXct(strptime(crkdatin$Timestamp, "%a/%b/%d %H:%M:%S %z %Y"))


# Read in the MATLAB file, as some columns were corrupted
crktemp <- readMat("./data_raw/Crk/CRK2016.mat")
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
  #CO2
  G = as.numeric(crkdatin$G.1.),
  LE = as.numeric(crkdatin$LE),
  ET = as.numeric(crkdatin$LE)*0.0352512,
  H = as.numeric(crkdatin$H),
  NEE = as.numeric(crkdatin$NEE_MDS),
  RECO = as.numeric(crkdatin$RE_MDS),
  GPP = as.numeric(crktemp4$GPP),
  SWC_0_10 = (crkdatin$SWC.1. + crkdatin$SWC.3.)/2,
  SWC_10_20 = (crkdatin$SWC.2. + crkdatin$SWC.4.)/2,
  TS_0_10 = (crkdatin$T_SOIL.2. + crkdatin$T_SOIL.5.)/2,
  TS_10_20 = (crkdatin$T_SOIL.3. + crkdatin$T_SOIL.6.)/2,
  SIF_O2A = crktemp4$SIF760/pi, # need to divide by pi to scale to sfm SIF from other sites
  #SIF_O2B,
  #NDVI,
  #PRI,
  Ref_1 = rowMeans(crk_ref, na.rm = T),
  Gs = crkdatin$Gs_ms
)

crkdat[crkdat==-9999] <- NA # replace -9999 with NAs

# Fix TIMESTAMP
crkdat$TIMESTAMP <- mdy_hm(paste(crkdat$TIMESTAMP))

write.csv(crkdat, file = "data_formatted/Crk/Crk_dat.csv") # save data as csv file


#################################################### Gebesee
gebdatin <- read.csv("./data_raw/DE-Geb/DE-Geb_dat.csv")

gebdat <- data.frame(
  Year = gebdatin$Year,
  DOY = gebdatin$DOY,
  TIMESTAMP = paste(gebdatin$DATE, gebdatin$TIMESTAMP),
  TA = gebdatin$TA_1_1_1,
  SW_IN = gebdatin$SWIN_1_1_1,
  LW_IN = gebdatin$LWIN_1_1_1,
  VPD = gebdatin$VPD,
  RH = gebdatin$RH,
  PA = gebdatin$air_pressure,
  P = gebdatin$P_RAIN_1_1_1,
  WS = gebdatin$wind_speed,
  #NETRAD
  PAR = gebdatin$PAR_eddy,
  APAR = gebdatin$APAR_flox,
  SW_OUT = gebdatin$SWOUT_1_1_1,
  LW_OUT = gebdatin$LWOUT_1_1_1,
  #CO2
  #G
  LE = gebdatin$LE,
  ET = gebdatin$LE*0.0352512,
  H = gebdatin$H,
  NEE = gebdatin$NEE,
  #RECO,
  GPP = gebdatin$GPP,
  SWC_1_1_1 = gebdatin$SWC_1_1_1,
  SWC_1_2_1 = gebdatin$SWC_1_2_1,
  SWC_1_3_1 = gebdatin$SWC_1_3_1,
  SWC_1_4_1 = gebdatin$SWC_1_4_1,
  SWC_1_5_1 = gebdatin$SWC_1_5_1,
  #TS,
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
  Gs = gebdatin$Gs
)

gebdat[gebdat==-9999] <- NA # replace -9999 with NAs

write.csv(gebdat, file = "data_formatted/DE-geb/DE-Geb_dat.csv") # save data as csv file

#################################################### Leinefelde

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
  #P = lnfdatin,
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
  SWC_08 = lnfdatin$SWCa008,
  SWC_16 = lnfdatin$SWCa016,
  SWC_32 = lnfdatin$SWCa032,
  SWC_62 = lnfdatin$SWCa062,
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
  PAR = lnfdatin$PPFD_IN,
  LW_OUT = lnfdatin$LW_OUT_CNR1,
  ET = lnfdatin$ET,
  GPP = lnfdatin$GPP_MR_f,
  SWC_08 = lnfdatin$SWCa008,
  SWC_16 = lnfdatin$SWCa016,
  SWC_32 = lnfdatin$SWCa032,
  SWC_62 = lnfdatin$SWCa062,
  Gs = lnfdatin$Gs_ms
)


# Fix dates
lnfdat_full$TIMESTAMP = ymd_hm(paste(lnfdat_full$date, lnfdat_full$time))
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
  PAR = lnfdat$PAR.x,
  LW_OUT = lnfdat$LW_OUT.x,
  ET = lnfdat$ET.x,
  GPP = lnfdat$GPP.x,
  SWC_08 = lnfdat$SWC_08.x,
  SWC_16 = lnfdat$SWC_16.x,
  SWC_32 = lnfdat$SWC_32.x,
  SWC_62 = lnfdat$SWC_62.x,
  SIF_O2A = lnfdat$SIF_O2A,
  SIF_O2B = lnfdat$SIF_O2B,
  NDVI = lnfdat$NDVI,
  PRI = lnfdat$PRI,
  Ref_750 = lnfdat$Ref_750,
  Ref_760 = lnfdat$Ref_760,
  Gs = lnfdat$Gs.x
)

lnfdat[lnfdat==-9999] <- NA # replace -9999 with NAs

write.csv(lnfdat, file = "data_formatted/DE-Lnf/DE-Lnf_dat.csv") # save data as csv file


#################################################### Yat

# Read in Yatir data
yatdatin_full <- read.csv("./data_raw/IL-Yat/IL-Yat_full.csv")
yatdatin <- read.csv("./data_raw/IL-Yat/IL-Yat_dat.csv")

yatdatin_full <- yatdatin_full[-1,] # get rid of unit row, as these are recorded in the README

yatdat_full <- data.frame(
  Year = yatdatin_full$Year,
  DOY = NA,
  date = yatdatin_full$TIMESTAMP,
  time = yatdatin_full$TIMESTAMP.2,
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
  SWC_1_1_1 = yatdatin_full$SWC_1_1_1,
  SWC_1_2_1 = yatdatin_full$SWC_1_2_1,
  SWC_1_3_1 = yatdatin_full$SWC_1_3_1,
  SWC_1_4_1 = yatdatin_full$SWC_1_4_1,
  SWC_1_5_1 = yatdatin_full$SWC_1_5_1#,
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

# Fix dates
yatdat_full$TIMESTAMP = mdy_hm(paste(yatdat_full$date, yatdat_full$time))
yatdat_full$DOY = difftime(yatdat_full$TIMESTAMP,as.POSIXct(as.Date("2017-01-01 00:00")))
yatdat_full$DOY<-gsub(" days","",yatdat_full$DOY)


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
yatdat_full <- yatdat_full %>%
  select(-c(date,time))


# Join to include met info from the full day
yatdat <- left_join(yatdat_full, yatdat_cropped, by = "TIMESTAMP")

# Reorder columns
yatdat <- yatdat %>%
  relocate(Year, DOY, TIMESTAMP)

yatdat[yatdat==-9999] <- NA # replace -9999 with NAs

write.csv(yatdat, file = "data_formatted/IL-yat/IL-yat_dat.csv") # save data as csv file


##################################################### Jrs

# Load JRS data
jrsdatin <- read.csv("./data_raw/Jrs/Jrs_dat.csv")

jrsdatin <- jrsdatin[-1,] # get rid of unit row, as these are recorded in the README

jrsdat <- data.frame(
  Year = jrsdatin$Year,
  DOY = jrsdatin$DOY,
  TIMESTAMP = NA,
  TA = (as.numeric(jrsdatin$TA_1_1_1) + as.numeric(jrsdatin$TA_1_2_1))/2,
  SW_IN = jrsdatin$SWIN_1_1_1,
  LW_IN = jrsdatin$LWIN_1_1_1,
  VPD = jrsdatin$VPD,
  RH = jrsdatin$RH,
  PA = as.numeric(jrsdatin$air_pressure)/100, # convert Pa to kPa
  P = jrsdatin$P_RAIN_1_1_1,
  WS = jrsdatin$wind_speed,
  #NETRAD
  #PAR = jrsdatin,
  #APAR = jrsdatin,
  SW_OUT = jrsdatin$SWOUT_1_1_1,
  LW_OUT = jrsdatin$LWOUT_1_1_1,
  #CO2
  #G
  LE = jrsdatin$LE,
  ET = jrsdatin$ET,
  H = jrsdatin$H,
  NEE = jrsdatin$NEE,
  #RECO,
  GPP = jrsdatin$GPP_gapfilled,
  SWC_1_1_1 = jrsdatin$SWC_1_1_1,
  SWC_1_2_1 = jrsdatin$SWC_1_2_1,
  SWC_1_3_1 = jrsdatin$SWC_1_3_1,
  SWC_1_4_1 = jrsdatin$SWC_1_4_1,
  SWC_1_5_1 = jrsdatin$SWC_1_5_1,
  TS_1_1_1 = jrsdatin$TS_1_1_1,
  TS_1_2_1 = jrsdatin$TS_1_2_1,
  TS_1_3_1 = jrsdatin$TS_1_3_1,
  TS_1_4_1 = jrsdatin$TS_1_4_1,
  TS_1_5_1 = jrsdatin$TS_1_5_1,
  SIF_O2A = jrsdatin$SFM_A_linear,
  #SIF_O2B = jrsdatin$SIF_O2B,
  NDVI = jrsdatin$ndvi,
  PRI = jrsdatin$pri,
  Ref_740_750 = jrsdatin$Ref_740_750_Reflectance,
  Ref_750_760 = jrsdatin$Ref_750_760_Reflectance,
  Ref_760_770 = jrsdatin$Ref_760_770_Reflectance,
  Ref_770_780 = jrsdatin$Ref_770_780_Reflectance,
  Ref_780_790 = jrsdatin$Ref_780_790_Reflectance,
  Gs = jrsdatin$gs.ms
)

# Fix dates
jrsdat$time <- paste(as.numeric(jrsdat$DOY) %/% 1, times(as.numeric(jrsdat$DOY) %% 1))
jrsdat$time <- sub("^\\S+\\s+", '', jrsdat$time)
jrsdat$date <-as.Date(as.numeric(jrsdat$DOY), origin = "2019-01-01")
jrsdat$TIMESTAMP = ymd_hms(paste(jrsdat$date, jrsdat$time))

# Pick which columns to include
jrsdat <- jrsdat %>%
  select(-c(time,date))

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
  P = sqdatin$P_RAIN_1_1_1,
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










