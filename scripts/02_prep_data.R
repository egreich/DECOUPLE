##### Prep data to run decouple SAM model

### Call functions

library(tidyverse)
source("./scripts/functions.R") # for assign_block function


### Load data

# Leinefelde
lnfdatin <- read.csv("./data_input/DE-Lnf_dat.csv")

# Gebesee
gebdatin <- read.csv("./data_input/DE-Gebdata.csv")

### Modify data columns UNITS
WUE_GPP = gebdatin$GPP/gebdatin$ET
gebdat <- gebdatin %>%
  mutate(
    ET = LE*0.0352512,
    VPD = as.numeric(VPD),
    Tair = TA_1_1_1,
    PAR = PAR_flox_mol,
    SWC_shall = SWC_1_1_1,
    SWC_deep = SWC_1_4_1) %>%
  mutate(WUE_GPP = GPP/ET,
         WUE_SIF = SIF_O2A_sfm/ET)

###################### evap scenario
### Modify data columns
# To run the get_evap function we need columns of:
# PA, Tair, Z, Tsoil, ws, RH, fsand, fclay, fc, SWC_shall
# Take the high and low of the last three days and calculate the average temperature. 

# gebdat <- gebdatin %>%
#   mutate(ET = LE*0.0352512,
#          SWC_shall = SWC_1_1_1, 
#          Tair = TA_1_1_1,
#          Tsoil = TA_1_1_1*1.5,
#          PA = air_pressure,
#          ws = wind_speed,
#          Z = 3.1,
#          h = 	161.5,
#          fc = 0.1671617535,
#          fclay = .06286599798,
#          fsand = .7446068832)

# Low-level sanity filtering
# gebdat$ET <- ifelse(gebdat$ET<0, 0, gebdat$ET)
# # Fill small amounts of data
gebdat$Gs <- fill_small(gebdat, "Gs")
gebdat$VPD <- fill_small(gebdat, "VPD")
gebdat$PAR <- fill_small(gebdat, "PAR")
# 
# 
# # Need, PA, Tair, Z, Tsoil, ws, RH, fsand, fclay, fc, SWC_shall
# gebdat <- get_evap(gebdat)
# 
# # Low-level sanity filtering
# gebdat$E4.5 <- ifelse(gebdat$ET<gebdat$E4.5, gebdat$ET, gebdat$E4.5)
# gebdat$E3.5 <- ifelse(gebdat$ET<gebdat$E3.5 , gebdat$ET, gebdat$E3.5)
# 
# gebdat <- gebdat %>%
#   mutate(T = ET-E3.5)
# 
# ggplot(data = gebdat) +
#   geom_line(aes(x=c(1:nrow(gebdat)), y = ET, color = "ET")) +
#   geom_line(aes(x=c(1:nrow(gebdat)), y = E4.5, color = "E4.5")) +
#   geom_line(aes(x=c(1:nrow(gebdat)), y = E3.5, color = "E3.5")) +
#   geom_line(aes(x=c(1:nrow(gebdat)), y = T, color = "T")) +
#   theme_bw()


# Save out
save(gebdat, file = "data_clean/gebdat.Rdata")


p <- ggplot(data = gebdat) +
  geom_line(aes(x=c(1:nrow(gebdat)), y = ET, color = "ET")) +
  geom_line(aes(x=c(1:nrow(gebdat)), y = GPP, color = "GPP")) +
  xlim(1000,2000) +
  #geom_line(aes(x=c(1:nrow(gebdat)), y = SIF_O2B_sfm, color = "SIF_O2B")) +
  #geom_line(aes(x=c(1:nrow(gebdat)), y = SIF_O2A_sfm, color = "SIF_O2A")) +
  theme_bw()
p
