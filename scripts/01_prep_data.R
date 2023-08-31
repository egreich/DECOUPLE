### Script to filter data
### input data: formatted data from data_input
### output data: csv dataframe for each site with SIF filtered for quality

library(tidyverse)
library(lubridate)
library(gridExtra) # for grid.arrange for combining graphs
library(zoo) # for na.approx to linearly interpolate
source("./scripts/functions.R")

#################################################### Laegeren
# Load Laegeren data
dat <- read.csv("./data_formatted/CH-Lae/CH-Lae_dat.csv")
dat$TIMESTAMP <- ymd_hms(dat$TIMESTAMP) # for filtering
dat$date <- as.Date(dat$TIMESTAMP, "%Y-%m-%d") # for filtering

# For each day, take R2 of GPP and SIF, and GPP and T, and SIF and T
# low R2s will be considered "decoupling" days
r_df_list <- c()
all_days <- unique(dat$date)
for(i in c(1:length(all_days))){
  today <- all_days[i]
  df_temp <- dat %>%
    filter(date == today)
  reg <- tryCatch(lm(GPP ~ SIF_O2A, data = df_temp), error=function(err) NA)
  
  if(!is.na(reg)){
    r_df_list[[i]] <- data.frame(date = today, rsquared = summary(reg)$r.squared)
  }
}
r_df <- bind_rows(r_df_list)
lower_bound <- quantile(r_df$rsquared)[2]
r_df <- r_df %>%
  mutate(coupled = ifelse(rsquared<lower_bound, F, T))


dat <- dat %>%
  mutate(
    # decide what is SWC_shall and what is SWC_deep
    SWC_shall = SWC_1,
    SWC_deep = SWC_4,
    TS = TS_1,
    # define WUE to relate carbon to water fluxes
    WUE_GPP = GPP/ET,
    WUE_SIF = SIF_O2A/ET,
    WUE_GPP_TEA = GPP/T_TEA,
    WUE_SIF_TEA = SIF_O2A/T_TEA
  )

dat <- left_join(dat, r_df)

# Fill small gaps
dat$TA <- fill_small(dat, "TA")
dat$VPD <- fill_small(dat, "VPD")
dat$PAR <- fill_small(dat, "PAR")
#dat$SWC_shall <- fill_small(dat, "SWC_shall")
#dat$SWC_deep <- fill_small(dat, "SWC_deep")

save(dat, file = "data_clean/laedat.RData") # save prepped data as R file, this will be the model input

#################################################### Crk
# Read in Crk data
dat <- read.csv("./data_formatted/KR-Crk/Crk_dat.csv")
dat$date <- as.Date(dat$TIMESTAMP, "%Y-%m-%d") # for filtering
dat$TIMESTAMP <- ymd_hms(dat$TIMESTAMP) # for filtering


# For each day, take R2 of GPP and SIF, and GPP and T, and SIF and T
# low R2s will be considered "decoupling" days
r_df_list <- c()
all_days <- unique(dat$date)
for(i in c(1:length(all_days))){
  today <- all_days[i]
  df_temp <- dat %>%
    filter(date == today)
  reg <- tryCatch(lm(GPP ~ SIF_O2A, data = df_temp), error=function(err) NA)
  
  if(!is.na(reg)){
    r_df_list[[i]] <- data.frame(date = today, rsquared = summary(reg)$r.squared)
  }
}
r_df <- bind_rows(r_df_list)
lower_bound <- quantile(r_df$rsquared)[2]
r_df <- r_df %>%
  mutate(coupled = ifelse(rsquared<lower_bound, F, T))


dat <- dat %>%
  mutate(
    # decide what is SWC_shall and what is SWC_deep
    SWC_shall = SWC_1,
    SWC_deep = SWC_2,
    TS = TS_1,
    # define WUE to relate carbon to water fluxes
    WUE_GPP = GPP/ET,
    WUE_SIF = SIF_O2A/ET,
    WUE_GPP_TEA = GPP/T_TEA,
    WUE_SIF_TEA = SIF_O2A/T_TEA
  )

dat <- left_join(dat, r_df)

# Fill small gaps
dat$TA <- fill_small(dat, "TA")
dat$VPD <- fill_small(dat, "VPD")
dat$PAR <- fill_small(dat, "PAR")
#dat$SWC_shall <- fill_small(dat, "SWC_shall")
#dat$SWC_deep <- fill_small(dat, "SWC_deep")


save(dat, file = "data_clean/crkdat.RData") # save prepped data as R file, this will be the model input


#################################################### Gebesee
# Load Gebesee data
dat <- read.csv("./data_formatted/DE-Geb/DE-Geb_dat.csv")
dat$TIMESTAMP <- ymd_hms(dat$TIMESTAMP) # for filtering
dat$date <- as.Date(dat$TIMESTAMP, "%Y%m%d") # for filtering


# For each day, take R2 of GPP and SIF, and GPP and T, and SIF and T
# low R2s will be considered "decoupling" days
r_df_list <- c()
all_days <- unique(dat$date)
for(i in c(1:length(all_days))){
  today <- all_days[i]
  df_temp <- dat %>%
    filter(date == today)
  reg <- tryCatch(lm(GPP ~ SIF_O2A, data = df_temp), error=function(err) NA)
  
  if(!is.na(reg)){
    r_df_list[[i]] <- data.frame(date = today, rsquared = summary(reg)$r.squared)
  }
}
r_df <- bind_rows(r_df_list)
lower_bound <- quantile(r_df$rsquared)[2]
r_df <- r_df %>%
  mutate(coupled = ifelse(rsquared<lower_bound, F, T))


dat <- dat %>%
  mutate(
    # decide what is SWC_shall and what is SWC_deep
    SWC_shall = SWC_1,
    SWC_deep = SWC_4,
    TS = TS_1,
    # define WUE to relate carbon to water fluxes
    WUE_GPP = GPP/ET,
    WUE_SIF = SIF_O2A/ET,
    WUE_GPP_TEA = GPP/T_TEA,
    WUE_SIF_TEA = SIF_O2A/T_TEA
  )

dat <- left_join(dat, r_df)

# Fill small gaps
dat$TA <- fill_small(dat, "TA")
dat$VPD <- fill_small(dat, "VPD")
dat$PAR <- fill_small(dat, "PAR")
#dat$SWC_shall <- fill_small(dat, "SWC_shall")
#dat$SWC_deep <- fill_small(dat, "SWC_deep")


save(dat, file = "data_clean/gebdat.RData") # save prepped data as R file, this will be the model input


#################################################### Leinefelde
# Load Leinefelde data
dat <- read.csv("./data_formatted/DE-Lnf/DE-Lnf_dat.csv")
dat$date <- as.Date(dat$TIMESTAMP, "%Y-%m-%d") # for filtering
dat$TIMESTAMP <- ymd_hms(dat$TIMESTAMP) # for filtering

# For each day, take R2 of GPP and SIF, and GPP and T, and SIF and T
# low R2s will be considered "decoupling" days
r_df_list <- c()
all_days <- unique(dat$date)
for(i in c(1:length(all_days))){
  today <- all_days[i]
  df_temp <- dat %>%
    filter(date == today)
  reg <- tryCatch(lm(GPP ~ SIF_O2A, data = df_temp), error=function(err) NA)
  
  if(!is.na(reg)){
    r_df_list[[i]] <- data.frame(date = today, rsquared = summary(reg)$r.squared)
  }
}
r_df <- bind_rows(r_df_list)
lower_bound <- quantile(r_df$rsquared)[2]
r_df <- r_df %>%
  mutate(coupled = ifelse(rsquared<lower_bound, F, T))


dat <- dat %>%
  mutate(
    # decide what is SWC_shall and what is SWC_deep
    SWC_shall = SWC_1,
    SWC_deep = SWC_2,
    TS = TS_1,
    # define WUE to relate carbon to water fluxes
    WUE_GPP = GPP/ET,
    WUE_SIF = SIF_O2A/ET,
    WUE_GPP_TEA = GPP/T_TEA,
    WUE_SIF_TEA = SIF_O2A/T_TEA
  )

dat <- left_join(dat, r_df)

# Fill small gaps
dat$TA <- fill_small(dat, "TA")
dat$VPD <- fill_small(dat, "VPD")
#dat$SWC_shall <- fill_small(dat, "SWC_shall")
#dat$SWC_deep <- fill_small(dat, "SWC_deep")


save(dat, file = "data_clean/lnfdat.RData") # save prepped data as R file, this will be the model input



#################################################### Yat
# Read in Yatir data
dat <- read.csv("./data_formatted/IL-Yat/IL-Yat_dat.csv")
dat$date <- as.Date(dat$TIMESTAMP, "%Y-%m-%d") # for filtering
dat$TIMESTAMP <- ymd_hms(dat$TIMESTAMP) # for filtering


# For each day, take R2 of GPP and SIF, and GPP and T, and SIF and T
# low R2s will be considered "decoupling" days
r_df_list <- c()
all_days <- unique(dat$date)
for(i in c(1:length(all_days))){
  today <- all_days[i]
  df_temp <- dat %>%
    filter(date == today)
  
  reg <- tryCatch(lm(GPP ~ SIF_O2A, data = df_temp), error=function(err) NA)
  
  if(!is.na(reg)){
    r_df_list[[i]] <- data.frame(date = today, rsquared = summary(reg)$r.squared)
  }
}
r_df <- bind_rows(r_df_list)
lower_bound <- quantile(r_df$rsquared)[2]
r_df <- r_df %>%
  mutate(coupled = ifelse(rsquared<lower_bound, F, T))


dat <- dat %>%
  mutate(
    # decide what is SWC_shall and what is SWC_deep
    SWC_shall = SWC_1,
    SWC_deep = SWC_4,
    TS = TS_1,
    # define WUE to relate carbon to water fluxes
    WUE_GPP = GPP/ET,
    WUE_SIF = SIF_O2A/ET,
    WUE_GPP_TEA = GPP/T_TEA,
    WUE_SIF_TEA = SIF_O2A/T_TEA
  )

dat <- left_join(dat, r_df)

# Fill small gaps
dat$TA <- fill_small(dat, "TA")
dat$VPD <- fill_small(dat, "VPD")
#dat$SWC_shall <- fill_small(dat, "SWC_shall")
#dat$SWC_deep <- fill_small(dat, "SWC_deep")

save(dat, file = "data_clean/yatdat.RData") # save prepped data as R file, this will be the model input


##################################################### Jrs
# Load JRS data
dat <- read.csv("./data_formatted/Jrs/Jrs_dat.csv")
dat$date <- as.Date(dat$TIMESTAMP, "%Y-%m-%d") # for filtering
dat$TIMESTAMP <- ymd_hms(dat$TIMESTAMP) # for filtering



# For each day, take R2 of GPP and SIF
# low R2s will be considered "decoupling" days
r_df_list <- c()
all_days <- unique(dat$date)
for(i in c(1:length(all_days))){
  today <- all_days[i]
  df_temp <- dat %>%
    filter(date == today)
  reg <- tryCatch(lm(GPP ~ SIF_O2A, data = df_temp), error=function(err) NA)
  
  if(!is.na(reg)){
    r_df_list[[i]] <- data.frame(date = today, rsquared = summary(reg)$r.squared)
  }
}
r_df <- bind_rows(r_df_list)
lower_bound <- quantile(r_df$rsquared)[2]
r_df <- r_df %>%
  mutate(coupled = ifelse(rsquared<lower_bound, F, T))


dat <- dat %>%
  mutate(
    # decide what is SWC_shall and what is SWC_deep
    SWC_shall = SWC_1,
    SWC_deep = SWC_4,
    TS = TS_1,
    # define WUE to relate carbon to water fluxes
    WUE_GPP = GPP/ET,
    WUE_SIF = SIF_O2A/ET,
    WUE_GPP_TEA = GPP/T_TEA,
    WUE_SIF_TEA = SIF_O2A/T_TEA
  )

dat <- left_join(dat, r_df)

# Fill small gaps
dat$TA <- fill_small(dat, "TA")
dat$VPD <- fill_small(dat, "VPD")
#dat$SWC_shall <- fill_small(dat, "SWC_shall")
#dat$SWC_deep <- fill_small(dat, "SWC_deep")

save(dat, file = "data_clean/jrsdat.RData") # save prepped data as R file, this will be the model input


##################################################### Sq
# Load sq data
dat <- read.csv("./data_formatted/Sq/Sq_dat.csv")
dat$date <- as.Date(dat$TIMESTAMP, "%Y-%m-%d") # for filtering
dat$TIMESTAMP <- ymd_hms(dat$TIMESTAMP) # for filtering
roi <- dat$Ref_770_780 # reflectance of interest
doy_win <- c(220,260)# window to look at

# Look at SIF and Reflectance
# p1 <- ggplot(data = dat) +
#   geom_line(aes(x=DOY, y=SIF_O2A))+
#   xlim(doy_win)+
#   theme_bw()
# 
# p2 <- ggplot(data = dat) +
#   geom_line(aes(x=DOY, y=Ref_770_780))+
#   xlim(doy_win)+
#   theme_bw()
# 
# grid.arrange(p1,p2, nrow=2)

# Code for smoothing SIF, if we want to do that
# # preserve location of the true NAs in SIF data
# for(i in c(1:nrow(dat))){
#   if(is.na(dat$SIF_O2A[i])){
#     dat$SIF_O2A_NA[i] <- T
#   }else{
#     dat$SIF_O2A_NA[i] <- F
#   }
# }
# upper_bound <- quantile(roi, na.rm = T)[4] # high reflectance = more clouds
# dat$SIF_O2A <- ifelse(roi > upper_bound, NA, dat$SIF_O2A)
# dat$SIF_O2A <- na.approx(dat$SIF_O2A, na.rm = FALSE)
# # maintain the true NAs (night)
# for(i in c(1:nrow(dat))){
#   if(dat$SIF_O2A_NA[i] == T){
#     dat$SIF_O2A[i] <- NA
#   }
# }

# Look at GPP and SIF
# p1 <- ggplot(data = dat) +
#   geom_line(aes(x=DOY, y=SIF_O2A))+
#   xlim(doy_win)+
#   theme_bw()
# 
# p2 <- ggplot(data = dat) +
#   geom_line(aes(x=DOY, y=GPP))+
#   xlim(doy_win)+
#   theme_bw()
# 
# grid.arrange(p1,p2, nrow=2)


# For each day, take R2 of GPP and SIF, and GPP and T, and SIF and T
# low R2s will be considered "decoupling" days
r_df_list <- c()
all_days <- unique(dat$date)
for(i in c(1:length(all_days))){
  today <- all_days[i]
  df_temp <- dat %>%
    filter(date == today)
  reg <- tryCatch(lm(GPP ~ SIF_O2A, data = df_temp), error=function(err) NA)
  
  if(!is.na(reg)){
    r_df_list[[i]] <- data.frame(date = today, rsquared = summary(reg)$r.squared)
  }

}
r_df <- bind_rows(r_df_list)
lower_bound <- quantile(r_df$rsquared, na.rm=T)[2]
r_df <- r_df %>%
  mutate(coupled = ifelse(rsquared<lower_bound, F, T))


dat <- dat %>%
  mutate(
    # decide what is SWC_shall and what is SWC_deep
    SWC_shall = SWC_1_1_1,
    SWC_deep = SWC_1_3_1,
    # define WUE to relate carbon to water fluxes
    WUE_GPP = GPP/ET,
    WUE_SIF = SIF_O2A/ET
  )

dat <- left_join(dat, r_df)

# Fill small gaps
dat$TA <- fill_small(dat, "TA")
dat$VPD <- fill_small(dat, "VPD")
dat$PAR <- fill_small(dat, "PAR")
#dat$SWC_shall <- fill_small(dat, "SWC_shall")
#dat$SWC_deep <- fill_small(dat, "SWC_deep")


# fill SIF gaps in the middle of days
# all_days <- unique(dat$date)
# for(i in c(1:length(all_days))){
#   # get one day
#   today <- all_days[i]
#   
#   # get indices from 7:30 to 16:00 for each day
#   index_range <- dat %>%
#     rowid_to_column("ind") %>%
#     filter(date == today) %>%
#     mutate(time = as.numeric(format(TIMESTAMP, "%H%M"))) %>%
#     filter(between(time, 730, 1600)) %>%
#     select(ind)
#   if(nrow(index_range)<1){
#     next
#   }
#   first <- index_range[1,1]
#   last <- index_range[nrow(index_range),1]
#   
#   # filter a temp df from 7:30 to 16:00
#   df_temp <- dat %>%
#     filter(date == today) %>%
#     mutate(time = as.numeric(format(TIMESTAMP, "%H%M"))) %>%
#     filter(between(time, 730, 1600)) %>%
#     select(-time)
#   
#   df_temp$SIF_O2A <- na.approx(df_temp$SIF_O2A, na.rm = FALSE) # fill daytime gaps
#   dat$SIF_O2A[first:last] <- df_temp$SIF_O2A # add back into main dataframe
# }

save(dat, file = "data_clean/sqdat.RData") # save prepped data as R file, this will be the model input





