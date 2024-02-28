### Script to filter data
### input data: formatted data from data_input
### output data: csv dataframe for each site with SIF filtered for quality

library(tidyverse)
library(lubridate)
library(gridExtra) # for grid.arrange for combining graphs
library(zoo) # for na.approx to linearly interpolate
library(cowplot)
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
  reg <- tryCatch(lm(GPP ~ T_TEA, data = df_temp), error=function(err) NA)
  
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

##################### temp pasted from yatir
# temp - note that this shows decoupling
p <- dat %>%
  ggplot() +
  geom_line(aes(x= TIMESTAMP, y = rsquared)) +
  geom_hline(yintercept=lower_bound, linetype = "dashed", color = "red") +
  theme_bw()
p
path_out = "./plots/11-15/" # set save path
ggsave2("p_lae_rsquared.png", plot = p, path = path_out, width = 8, height = 6)

hist(dat$rsquared)

p <- dat %>%
  ggplot() +
  geom_line(aes(x= DOY, y = GPP/100, color = "GPP")) +
  geom_line(aes(x= DOY, y = SIF_O2B, color = "SIF_O2B")) +
  geom_line(aes(x= DOY, y = T_TEA, color = "T_TEA")) +
  xlim(c(230,240)) +
  theme_bw() +
  ggtitle("CH-Lae MIxed Deciduous Mountain Forest") +
  theme(legend.title = element_blank())
p 
path_out = "./plots/11-15/" # set save path
ggsave2("p_lae_decoup.png", plot = p, path = path_out, width = 8, height = 6)
############## temp pasted from yatir

# temp - note that this shows decoupling

p1 <- dat %>%
  ggplot() +
  geom_point(aes(x= TA, y = GPP, color = DOY), size = .5, alpha = .3) +
  stat_smooth(aes(x= TA, y = GPP), method = 'loess', span = 0.5, size = 2, color = "black") +
  theme_bw()

p2 <- dat %>%
  ggplot() +
  geom_point(aes(x= TA, y = T_TEA, color = DOY), size = .5, alpha = .3) +
  stat_smooth(aes(x= TA, y = T_TEA), method = 'loess', span = 0.5, size = 2, color = "black") +
  theme_bw()

p <- grid.arrange(p1,p2, nrow=2)

path_out = "./plots/11-15/" # set save path
ggsave2("p_lae_fig3.png", plot = p, path = path_out, width = 6, height = 8)


dat %>%
  ggplot() +
  geom_point(aes(x= TIMESTAMP, y = rsquared, color = "rsquared"), alpha=0.5, size=.5) +
  stat_smooth(aes(x= TIMESTAMP, y = rsquared, color = "rsquared"), method = 'loess', span = 0.01, size = .4) +
  geom_line(aes(x= TIMESTAMP, y = LW_OUT/10000, color = "LW_OUT/10000")) +
  geom_line(aes(x= TIMESTAMP, y = TA/1000, color = "TA/1000")) +
  geom_hline(yintercept=lower_bound, linetype = "dashed", color = "red") +
  theme_bw()

hist(dat$rsquared)


p1 <- dat %>%
  ggplot() +
  geom_point(aes(x= DOY, y = GPP/100, color = "GPP/100"), size = .5) +
  geom_point(aes(x= DOY, y = T_TEA, color = "T_TEA"), size = .5) +
  stat_smooth(aes(x= DOY, y = GPP/100, color = "GPP/100"), method = 'loess', span = 0.1, size = .5) +
  stat_smooth(aes(x= DOY, y = T_TEA, color = "T_TEA"), method = 'loess', span = 0.1, size = .5) +
  xlim(c(200,230)) +
  theme_bw()

p1 <- dat %>%
  ggplot() +
  geom_point(aes(x= DOY, y = SIF_O2B, color = "SIF_O2B"), size = .5) +
  geom_point(aes(x= DOY, y = T_TEA, color = "T_TEA"), size = .5) +
  stat_smooth(aes(x= DOY, y = SIF_O2B, color = "SIF_O2B"), method = 'loess', span = 0.1, size = .5) +
  stat_smooth(aes(x= DOY, y = T_TEA, color = "T_TEA"), method = 'loess', span = 0.1, size = .5) +
  xlim(c(200,230)) +
  theme_bw()

p2 <- dat %>%
  ggplot() +
  geom_line(aes(x= DOY, y = GPP/100, color = "GPP")) +
  geom_line(aes(x= DOY, y = SIF_O2B, color = "SIF_O2B")) +
  geom_line(aes(x= DOY, y = T_TEA, color = "T_TEA")) +
  xlim(c(125,150)) +
  theme_bw()

grid.arrange(p1,p2, nrow=2)


library(cleaRskyQuantileRegression)
dat$Rsdpot_12 = calc_PotRadiation_CosineResponsePower(doy = dat$DOY, hour = dat$DOY,
                                                          latDeg = 51.3282,
                                                          longDeg = 10.3678,
                                                          timeZone = 0, isCorrectSolartime = TRUE,
                                                          cosineResponsePower = 1.2 )

library(latticeExtra)
lab_Rsdpot12_name = expression(bold("Potential Shortwave Radiation ")(W*m^-2))
lab_Rsd      = expression(bold("Observed Shortwave Radiation ")(W*m^-2))
xyplot(SW_IN ~  Rsdpot_12, data = dat,
       xlab = list(lab_Rsdpot12_name,cex = 1.3), ylab = list(label=lab_Rsd, cex=1.3),
       type = c("p","g"), pch = ".", cex = 3,
       # main = paste0(dtyrmon[SiteCode == sico, unique(SiteName)]," ,Quantile Regression, tau = 0.90" ),
       panel = function(x, y, ...) {
         panel.xyplot(x, y, ...)
         panel.abline(rq(y~x, tau = 0.85), col=1, lwd = 2)
         # panel.abline(0,1, col=1, lty = 2, lwd = 2)
         panel.ablineq(0,1, col=1, lty = 2, lwd = 2, at = 0.89, label = "1:1", rotate = TRUE, fontfamily = "sans", cex = 1.5,pos = 3)
         panel.ablineq(rq(y~x, tau = 0.85), col=1, at = 0.7, rotate =TRUE,
                       pos = 3, label = "Slope of 85% Quantile", cex = 1.5, font = "Helvetica", fontface = 2 )
         panel.text(700,-20, "half hourly data", fontfamily = "Helvetica")
       },
       grid = TRUE)

test <- calc_ClearSky_QuantileRegression(dat$SW_IN, dat$Rsdpot_12)


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
#summary(lm(SIF_O2B ~ T_TEA, data = dat))
#summary(lm(GPP ~ T_TEA, data = dat))
r_df_list <- c()
all_days <- unique(dat$date)
for(i in c(1:length(all_days))){
  today <- all_days[i]
  df_temp <- dat %>%
    filter(date == today)
  
  reg <- tryCatch(lm(GPP ~ T_TEA, data = df_temp), error=function(err) NA)
  #reg2 <- tryCatch(lm(SIF_O2B ~ T_TEA, data = df_temp), error=function(err) NA)
  
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

# temp - note that this shows decoupling
p <- dat %>%
  ggplot() +
  geom_line(aes(x= TIMESTAMP, y = rsquared)) +
  geom_hline(yintercept=lower_bound, linetype = "dashed", color = "red") +
  theme_bw()
p
path_out = "./plots/11-15/" # set save path
ggsave2("p_yat_rsquared.png", plot = p, path = path_out, width = 8, height = 6)

hist(dat$rsquared)

p <- dat %>%
  ggplot() +
  geom_line(aes(x= DOY, y = GPP/100, color = "GPP")) +
  geom_line(aes(x= DOY, y = SIF_O2B, color = "SIF_O2B")) +
  geom_line(aes(x= DOY, y = T_TEA, color = "T_TEA")) +
  xlim(c(110,115)) +
  theme_bw() +
  ggtitle("IL-Yatir Mediterranean Evergreen Coniferous Forest") +
  theme(legend.title = element_blank())
 p 
  path_out = "./plots/11-15/" # set save path
  ggsave2("p_yat_decoup.png", plot = p, path = path_out, width = 8, height = 6)
  
  p1 <- dat %>%
    ggplot() +
    geom_point(aes(x= TA, y = GPP, color = DOY), size = .5, alpha = .3) +
    stat_smooth(aes(x= TA, y = GPP), method = 'loess', span = 0.5, size = 2, color = "black") +
    theme_bw()
  
  p2 <- dat %>%
    ggplot() +
    geom_point(aes(x= TA, y = T_TEA, color = DOY), size = .5, alpha = .3) +
    stat_smooth(aes(x= TA, y = T_TEA), method = 'loess', span = 0.5, size = 2, color = "black") +
    theme_bw()
  
  p <- grid.arrange(p1,p2, nrow=2)
  
  path_out = "./plots/11-15/" # set save path
  ggsave2("p_yat_fig3.png", plot = p, path = path_out, width = 6, height = 8)


dat$Rsdpot_12 = calc_PotRadiation_CosineResponsePower(doy = dat$DOY, hour = dat$DOY,
                                                      latDeg = 31.345306,
                                                      longDeg = 35.051861,
                                                      timeZone = 0, isCorrectSolartime = TRUE,
                                                      cosineResponsePower = 1.2 )

library(latticeExtra)
lab_Rsdpot12_name = expression(bold("Potential Shortwave Radiation ")(W*m^-2))
lab_Rsd      = expression(bold("Observed Shortwave Radiation ")(W*m^-2))
xyplot(SW_IN ~  Rsdpot_12, data = dat,
       xlab = list(lab_Rsdpot12_name,cex = 1.3), ylab = list(label=lab_Rsd, cex=1.3),
       type = c("p","g"), pch = ".", cex = 3,
       # main = paste0(dtyrmon[SiteCode == sico, unique(SiteName)]," ,Quantile Regression, tau = 0.90" ),
       panel = function(x, y, ...) {
         panel.xyplot(x, y, ...)
         panel.abline(rq(y~x, tau = 0.85), col=1, lwd = 2)
         # panel.abline(0,1, col=1, lty = 2, lwd = 2)
         panel.ablineq(0,1, col=1, lty = 2, lwd = 2, at = 0.89, label = "1:1", rotate = TRUE, fontfamily = "sans", cex = 1.5,pos = 3)
         panel.ablineq(rq(y~x, tau = 0.85), col=1, at = 0.7, rotate =TRUE,
                       pos = 3, label = "Slope of 85% Quantile", cex = 1.5, font = "Helvetica", fontface = 2 )
         panel.text(700,-20, "half hourly data", fontfamily = "Helvetica")
       },
       grid = TRUE)

test <- calc_ClearSky_QuantileRegression(dat$SW_IN, dat$Rsdpot_12)


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





