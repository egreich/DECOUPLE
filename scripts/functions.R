### Helper functions
# Many of these functions are on github and are edited here in minor ways here
# Functions by Michael Fell (PostJAGS package): https://github.com/fellmk/PostJAGS
# Functions by Emma Reich and Biplabendu (Billu) Das (coda4dummies package): https://github.com/egreich/coda4dummies

library(purrr)
# Create a "not in" function using negate from the purrr package
`%nin%` <- negate(`%in%`)


# conversion factor for LE to ET data to get at the half hour, hour, or daily scale
convert.LE.to.ET <- function(scale = "halfhour", TA, LE.in){
  TA <- as.numeric(TA)
  LE.in <- as.numeric(LE.in)
  conv.fact.time <- ifelse(scale=="hour", 3600, 1800) # if we need to convert seconds to hours or half hours
  conv.fact.time <- ifelse(scale=="day", 86400, conv.fact.time) # if we need to convert seconds to days
  conv.fact <- conv.fact.time/((2.501 - 0.00237*TA)*(10^6))
  ET.out <- LE.in * conv.fact
  return(ET.out)
}

# Function to gap-fill variables using average of previous 7 timesteps in data
fill_small <- function(site, column_to_fill){
  ws_list <- site[,column_to_fill]
  prev_ws <- NA
  last_week <- NA
  
  for (i in 1:nrow(site)){
    
    if(i<10) next # skip 1st 10 iterations and go to next iteration
    
    prev_ws <- i-1
    last_week <- i-6
    
    if (is.na(ws_list[i])){
      ws_list[i] = mean(ws_list[last_week:prev_ws], na.rm = T)
    } else {
      ws_list[i] = ws_list[i]
    }
  }
  
  return(ws_list)
  
}


# Function to linear interpolate a variable if there is no P
# Linear interpolating, except when P is > 0 over the gap-period
fill_linear <- function(site, S){
  S_list <- S
  last_P <- NA
  prev_S <- NA
  step1 <- NA
  step2 <- NA
  step3 <- NA
  step4 <- NA
  step5 <- NA
  step6 <- NA
  step7 <- NA
  step8 <- NA
  step9 <- NA
  step10 <- NA
  step11 <- NA
  step12 <- NA
  step13 <- NA
  step14 <- NA
  step15 <- NA
  step16 <- NA
  step17 <- NA
  step18 <- NA
  step19 <- NA
  step20 <- NA
  step21 <- NA
  step22 <- NA
  step1P <- NA
  step2P <- NA
  step3P <- NA
  step4P <- NA
  step5P <- NA
  step6P <- NA
  step7P <- NA
  step8P <- NA
  step9P <- NA
  step10P <- NA
  step11P <- NA
  step12P <- NA
  step13P <- NA
  step14P <- NA
  step15P <- NA
  step16P <- NA
  step17P <- NA
  step18P <- NA
  step19P <- NA
  step20P <- NA
  step21P <- NA
  step22P <- NA
  
  for (i in 1:nrow(site)){
    
    if(is.na(site$P[i])){ # if precip is NA, set to 1 to just ignore those blocks
      site$P[i]=1
    }
    
    if(i<10) next # skip 1st 10 iterations and go to next iteration
    
    last_P <- i-7
    prev_S <- i-1
    step1 <- i+1
    step2 <- i+2
    step3 <- i+3
    step4 <- i+4
    step5 <- i+5
    step6 <- i+6
    step7 <- i+7
    step8 <- i+8
    step9 <- i+9
    step10 <- i+10
    step11 <- i+11
    step12 <- i+12
    step13 <- i+13
    step14 <- i+14
    step15 <- i+15
    step16 <- i+16
    step17 <- i+17
    step18 <- i+18
    step19 <- i+19
    step20 <- i+20
    step21 <- i+21
    step22 <- i+22
    
    step1P <- i+1-7
    step2P <- i+2-7
    step3P <- i+3-7
    step4P <- i+4-7
    step5P <- i+5-7
    step6P <- i+6-7
    step7P <- i+7-7
    step8P <- i+8-7
    step9P <- i+9-7
    step10P <- i+10-7
    step11P <- i+11-7
    step12P <- i+12-7
    step13P <- i+13-7
    step14P <- i+14-7
    step15P <- i+15-7
    step16P <- i+16-7
    step17P <- i+17-7
    step18P <- i+18-7
    step19P <- i+19-7
    step20P <- i+20-7
    step21P <- i+21-7
    step22P <- i+22-7
    
    
    if(!is.na(S_list[prev_S])){
      if (is.na(S_list[i])){
        if (site$P[i] == 0){
          if(!is.na(S_list[step1])){
            if(mean(site$P[i:step1], na.rm=T) == 0){
              S_list[i] = seq(S_list[prev_S], S_list[step1], length.out=3)[2]
            } else{
              S_list[i] = S_list[i]}
          } else if(!is.na(S_list[step2])){
            if(mean(site$P[i:step2], na.rm=T) == 0){
              S_list[i:step1] = seq(S_list[prev_S], S_list[step2], length.out=4)[2:3]
            } else{
              S_list[i:step1] = S_list[i:step1]}
          } else if(!is.na(S_list[step3])){
            if(mean(site$P[i:step3], na.rm=T) == 0){
              S_list[i:step2] = seq(S_list[prev_S], S_list[step3], length.out=5)[2:4]
            } else{
              S_list[i:step2] = S_list[i:step2]}
          } else if(!is.na(S_list[step4])){
            if(mean(site$P[i:step4], na.rm=T) == 0){
              S_list[i:step3] = seq(S_list[prev_S], S_list[step4], length.out=6)[2:5]
            } else{
              S_list[i:step3] = S_list[i:step3]}
          } else if(!is.na(S_list[step5])){
            if(mean(site$P[i:step5], na.rm=T) == 0){
              S_list[i:step4] = seq(S_list[prev_S], S_list[step5], length.out=7)[2:6]
            } else{
              S_list[i:step4] = S_list[i:step4]}
          } else if(!is.na(S_list[step6])){
            if(mean(site$P[i:step6], na.rm=T) == 0){
              S_list[i:step5] = seq(S_list[prev_S], S_list[step6], length.out=8)[2:7]
            } else{
              S_list[i:step5] = S_list[i:step5]}
          } else if(!is.na(S_list[step7])){
            if(mean(site$P[i:step7], na.rm=T) == 0){
              S_list[i:step6] = seq(S_list[prev_S], S_list[step7], length.out=9)[2:8]
            } else{
              S_list[i:step6] = S_list[i:step6]}
          } else if(!is.na(S_list[step8])){
            if(mean(site$P[i:step8], na.rm=T) == 0){
              S_list[i:step7] = seq(S_list[prev_S], S_list[step8], length.out=10)[2:9]
            } else{
              S_list[i:step7] = S_list[i:step7]}
          } else if(!is.na(S_list[step9])){
            if(mean(site$P[i:step9], na.rm=T) == 0){
              S_list[i:step8] = seq(S_list[prev_S], S_list[step9], length.out=11)[2:10]
            } else{
              S_list[i:step8] = S_list[i:step8]}
          } else if(!is.na(S_list[step10])){
            if(mean(site$P[i:step10], na.rm=T) == 0){
              S_list[i:step9] = seq(S_list[prev_S], S_list[step10], length.out=12)[2:11]
            } else{
              S_list[i:step9] = S_list[i:step9]}
          } else if(!is.na(S_list[step11])){
            if(mean(site$P[i:step11], na.rm=T) == 0){
              S_list[i:step10] = seq(S_list[prev_S], S_list[step11], length.out=13)[2:12]
            } else{
              S_list[i:step10] = S_list[i:step10]}
          } else if(!is.na(S_list[step12])){
            if(mean(site$P[i:step12], na.rm=T) == 0){
              S_list[i:step11] = seq(S_list[prev_S], S_list[step12], length.out=14)[2:13]
            } else{
              S_list[i:step11] = S_list[i:step11]}
          } else if(!is.na(S_list[step13])){
            if(mean(site$P[i:step13], na.rm=T) == 0){
              S_list[i:step12] = seq(S_list[prev_S], S_list[step13], length.out=15)[2:14]
            } else{
              S_list[i:step12] = S_list[i:step12]}
          } else if(!is.na(S_list[step14])){
            if(mean(site$P[i:step14], na.rm=T) == 0){
              S_list[i:step13] = seq(S_list[prev_S], S_list[step14], length.out=16)[2:15]
            } else{
              S_list[i:step13] = S_list[i:step13]}
          } else if(!is.na(S_list[step15])){
            if(mean(site$P[i:step15], na.rm=T) == 0){
              S_list[i:step14] = seq(S_list[prev_S], S_list[step15], length.out=17)[2:16]
            } else{
              S_list[i:step14] = S_list[i:step14]}
          } else if(!is.na(S_list[step16])){
            if(mean(site$P[i:step16], na.rm=T) == 0){
              S_list[i:step15] = seq(S_list[prev_S], S_list[step16], length.out=18)[2:17]
            } else{
              S_list[i:step15] = S_list[i:step15]}
          } else if(!is.na(S_list[step17])){
            if(mean(site$P[i:step17], na.rm=T) == 0){
              S_list[i:step16] = seq(S_list[prev_S], S_list[step17], length.out=19)[2:18]
            } else{
              S_list[i:step16] = S_list[i:step16]}
          } else if(!is.na(S_list[step18])){
            if(mean(site$P[i:step18], na.rm=T) == 0){
              S_list[i:step17] = seq(S_list[prev_S], S_list[step18], length.out=20)[2:19]
            } else{
              S_list[i:step17] = S_list[i:step17]}
          } else if(!is.na(S_list[step19])){
            if(mean(site$P[i:step19], na.rm=T) == 0){
              S_list[i:step18] = seq(S_list[prev_S], S_list[step19], length.out=21)[2:20]
            } else{
              S_list[i:step18] = S_list[i:step18]}
          } else if(!is.na(S_list[step20])){
            if(mean(site$P[i:step20], na.rm=T) == 0){
              S_list[i:step19] = seq(S_list[prev_S], S_list[step20], length.out=22)[2:21]
            } else{
              S_list[i:step19] = S_list[i:step19]}
          } else if(!is.na(S_list[step21])){
            if(mean(site$P[i:step21], na.rm=T) == 0){
              S_list[i:step20] = seq(S_list[prev_S], S_list[step21], length.out=23)[2:22]
            } else{
              S_list[i:step20] = S_list[i:step20]}
          } else if(!is.na(S_list[step22])){
            if(mean(site$P[i:step22], na.rm=T) == 0){
              S_list[i:step21] =  seq(S_list[prev_S], S_list[step22], length.out=24)[2:23]
            } else{
              S_list[i:step21] = S_list[i:step21]}
          }
        }else {
          S_list[i] = S_list[i]
        }
      }
    }
    
    
    if(!is.na(S_list[prev_S])){
      if (is.na(S_list[i])){
        if (mean(site$P[i], na.rm=T) == 0){
          S_list[i] <- S_list[prev_S]
        }else {
          S_list[i] = S_list[i]
        }
      }
    }
    
    
    
  }
  
  S_GF <- S_list
  return(S_GF)
}


# Function to calculate evaporation, based on descriptions of methods in Merlin et al. 2016
# To run the get_evap function we need columns of:
# PA, TA, TS, WS, RH, SWC_shall
# scalar inputs: Z, h, fc, fclay, fsand
# Take the high and low of the last three days and calculate the average temperature. 
get_evap <- function(site, Z, h, fc, fclay, fsand, conv.fact.time) {
  
  # Define constants for calculating vapor pressure. Constants from a Vaisala publication .
  # The maximum error using the constants below is 0.083%
  A  = 6.116441
  m  = 7.591386
  Tn = 240.7263
  
  # More Constants
  g  = 9.80665     # acceleration due to gravity in m/s^2
  Rs = 287.058 # specific gas constant for dry air in J/(kg*K)
  Ma = 28.9634 # molar mass of dry air in g/mol
  Mw = 18.01528 # molar mass of water in g/mol
  rmw <- Mw/Ma
  pi = 3.14159265359
  k <- 0.40 # the von Karman constant
  Rwv = 461.52 # specific gas constant for water vapor in J/(kg*K)
  
  # Calculate density of air
  #Pa <- (p0*(1-(L*site$h)/T0)**((g*M)/(R0*L))) # atmospheric pressure in Pa
  Pa <- site$PA*1000 # convert atmospheric pressure to Pa
  pair <- Pa/(Rs * (site$TA + 273.15)) # density of air, calculated using the ideal gas law in (kg/m^3)
  
  # Calculate psychrometric contstant
  Cp <- 1000 # specific heat of air J/kg * K(or degree C)
  lambda <- 22.6 * 10^5 # latent heat of water vaporization J/kg
  gamma <- (Cp*Pa)/(rmw*lambda) #psychrometric constant in Pa/K
  
  # Calculate aerodynamic resistance to heat transfer
  # 0.001 is Z0m (m), the momentum soil roughness, [Yang et al., 2008; Stefan et al., 2015]
  #Z <- site$Z # reference height (m)
  Ri <- (5*g*Z*(site$TS-site$TA))/(site$TA*(site$WS**2)) # Richardson number, represents the importance of free versus forced convection
  # a coefficient set to 0.75 in unstable conditions (T > Ta) and to 2 in stable conditions (T < Ta)
  eta_unstable <- 0.75
  #eta_stable <- 2
  rah0 <- (1/((k**2)*site$WS))*((log(Z/0.001))**2) # s/m
  rah_unstable <- rah0/((1 + Ri)^eta_unstable) # aerodynamic resistance to heat transfer s/m, estimated as in Choudhury et al. [1986]
  #rah_stable <- rah0/((1 + Ri)**eta_stable)
  rah_stable <- ((1 - 0.75*Ri)*((log(Z/0.001))**2))/((k**2)*site$WS)
  rah <- ifelse( Ri > 0, rah_unstable, rah_stable)
  
  # Vapor pressure calculations
  ea1 <- A * 10**((m*site$TA)/(site$TA+Tn)) # Temperature from data; A, m, and Tn are constants. Saturated.
  ea  <- ea1 * site$RH/100 # air vapor pressure, should be in Pa
  es <- A * 10**((m*site$TS)/(site$TS +Tn)) # Saturated vapor pressure in soil
  
  # Snow vapor pressure
  #Tsnow = (site$Tsoil + site$Tair)/2
  #es_snow =  A * 10**((m*Tsnow)/(Tsnow +Tn)) # Saturated vapor pressure in snow
  
  # calculate parameters for alternative alpha and Bowen
  psisat <- -10*exp(1.88-1.31*fsand) # parameterized air entry pressure, in mm of water
  bch <- 2.91 + 15.9*fclay # The Clapp and Hornberger parameter est. as in Cosby et al. [1984]
  sres <- 0.15*fclay # residual soil moisture
  ssat <- 0.489 - (0.126*fsand)
  rssmin <- 50 # minimum soil resistance s/m
  psi <- psisat * (site$SWC_shall/ssat)**(-bch)
  
  # calculate soil resistance
  rss1 <- ((fc - sres)/(site$SWC_shall - sres)) * rssmin # soil resistance for S > sres
  rss2 <- exp(8.206 - 4.255*(site$SWC_shall/fc)) # another soil resistance derivation
  rss <- ifelse(site$SWC_shall>sres, rss1, rss2) # the resistance to the diffusion of vapor in large soil pores
  
  # Calculate alpha and Bowen ratios
  #alpha1 <- ifelse(site$SWC_shall > site$fc, 1, 0.5 - 0.5*cos((site$SWC_shall * pi)/site$fc)) # wetness function
  alpha2 <- exp((psi*g)/(10000 * Rwv * site$TS)) # Kelvin equation
  bowen1 <- ifelse(site$SWC_shall > fc, 1, (0.5 - 0.5*cos((site$SWC_shall * pi)/fc))**2) # wetness function
  
  # This will calculate LE in W/m^2, so let's convert to mm (i.e. LE to E)
  # 1 Watt /m2 = 0.0864 MJ /m2/day
  # 1 MJ /m2/day  = 0.408 mm /day
  # conversion factor from bigleaf:
  conv.fact <- conv.fact.time/((2.501 - 0.00237*site$TA)*10^6)
  E4.5 <- (bowen1 * ((pair * Cp)/gamma) * ((alpha2*(es - ea))/rah)) * conv.fact# 0.0864 * 0.408 # from 4.5 CLM
  E3.5 <- (((pair * Cp)/gamma) * ((alpha2*(es - ea))/(rah+rss))) * conv.fact #0.0864 * 0.408 # from 3.5 CLM
  
  # Calculate intercepted Eint using LAI from MODIS
  # intercepted E (q_intr) in kg m-2 s-1
  # Eint = (site$P)*(1 - exp(-0.5*(site$LAI_mod)))
  
  site_E <- site %>%
    mutate(E4.5 = E4.5, E3.5 = E3.5, pair = pair, Bowen = bowen1, alpha = alpha2, gamma = gamma, Ri = Ri, rah_unstable = rah_unstable, rah = rah, rss = rss, es = es, ea = ea, Pa = Pa)
  return(site_E)
}

# function that organizes jagsUI output as a dataframe
# part of the coda4dummies package
# intended to work with Michael Fell's "coda.fast" function
# Inputs: 
# jagsui: jagsUI object
# dim: max dimension of parameters (between 1-2)
# type: rjags or jagsUI
# Output: dataframe
# dumsum <- function(jagsobj, dim, type){
#   
#   if(type %nin% c("rjags", "jagUI")){
#     paste("Please indicate whether this is a rjags or jagsUI samples object")
#   }
#   
#   
#   if(type == "jagsUI"){
#     jagsui <- jagsobj
#     jm_coda <- jagsobj$samples # convert to coda form to work with postjags functions
#     
#     # Organize the coda object as a dataframe
#     df_sum <- coda.fast(jm_coda)
#     df_sum <- rownames_to_column(df_sum, "var")
#     df_sum <- df_sum %>% # make index column
#       mutate(ID = sub('.*\\[(.*)\\]', '\\1', df_sum$var))
#     df_sum <- df_sum %>% # separate index column into 1st and 2nd dimension
#       mutate(ID1 = sub('(.*)\\,.*', '\\1', df_sum$ID),
#              ID2 = sub('.*\\,(.*)', '\\1', df_sum$ID))
#     df_sum$ID2 <- ifelse(!grepl(',', df_sum$ID), NA, df_sum$ID2) # get rid of ID2 if there's no 2nd dimension
#     df_sum$ID1 <- ifelse(grepl('[[:alpha:]]', df_sum$ID), 1, df_sum$ID1) # make ID1=1 if there is only 1 instance
#     df_sum <- df_sum %>% 
#       mutate(var = sub('(.*)\\[.*', '\\1', df_sum$var)) # get rid of numbers in var col
#     df_mod <- df_sum %>% 
#       select("var","ID1","ID2","mean","median","pc2.5","pc97.5") %>% #reorder columns, drop ID
#       mutate(overlap0 = do.call(c, jagsui$overlap0), gel = do.call(c, jagsui$Rhat))
#     
#     return(df_mod)
#   }
#   
#   if(type == "rjags"){
#     jm_coda <- jagsobj
#     
#     # Organize the coda object as a dataframe
#     df_sum <- coda.fast(jm_coda)
#     df_sum <- rownames_to_column(df_sum, "var")
#     df_sum <- df_sum %>% # make index column
#       mutate(ID = sub('.*\\[(.*)\\]', '\\1', df_sum$var))
#     df_sum <- df_sum %>% # separate index column into 1st and 2nd dimension
#       mutate(ID1 = sub('(.*)\\,.*', '\\1', df_sum$ID),
#              ID2 = sub('.*\\,(.*)', '\\1', df_sum$ID))
#     df_sum$ID2 <- ifelse(!grepl(',', df_sum$ID), NA, df_sum$ID2) # get rid of ID2 if there's no 2nd dimension
#     df_sum$ID1 <- ifelse(grepl('[[:alpha:]]', df_sum$ID), 1, df_sum$ID1) # make ID1=1 if there is only 1 instance
#     df_sum <- df_sum %>% 
#       mutate(var = sub('(.*)\\[.*', '\\1', df_sum$var)) # get rid of numbers in var col
#     df_mod <- df_sum %>% 
#       select("var","ID1","ID2","mean","median","pc2.5","pc97.5") #%>% #reorder columns, drop ID
#       #mutate(overlap0 = do.call(c, jagsui$overlap0), gel = do.call(c, jagsui$Rhat))
#     
#     return(df_mod)
#   }
#   
# }
dumsum <- function(jagsobj, type){
  
  # Create a "not in" function using negate from the purrr package
  `%nin%` <- purrr::negate(`%in%`)
  
  if(type %nin% c("rjags", "jagsUI")){
    paste("Please indicate whether this is a rjags or jagsUI samples object")
  }
  
  
  if(type == "jagsUI"){
    
    jagsui <- jagsobj
    jm_coda <- jagsui$samples # convert to coda form to work with postjags functions
    
    # Organize the coda object as a dataframe
    df_sum <- coda.fast(jm_coda)
    df_sum <- tibble::rownames_to_column(df_sum, "var")
    df_sum <- df_sum %>% # make index column
      mutate(ID = sub('.*\\[(.*)\\]', '\\1', df_sum$var))
    df_sum$ID <- ifelse(grepl('[[:alpha:]]', df_sum$ID), 1, df_sum$ID) # make ID=1 if there is only 1 instance
    
    # make a lists of list of indices
    IDlist <- strsplit(df_sum$ID, ",") #temp
    
    # get number of dims in ID
    counter <- 1
    for(i in 1:length(IDlist)){
      counter <- ifelse(length(IDlist[[i]])>counter, length(IDlist[[i]]), counter)
    }
    
    # create a character vector of column names based on max dim
    new_columns <- list()
    for(i in 1:counter){
      new_columns[[i]] <- paste("ID",i, sep="")
    }
    new_columns <- as.character(new_columns)
    
    # for each dimension, create a new column with the correct ID
    df_sum <- df_sum %>%
      tidyr::separate(ID,new_columns,sep=",")
    
    df_sum <- df_sum %>%
      mutate(var = sub('(.*)\\[.*', '\\1', df_sum$var)) # get rid of numbers in var col
    
    if(is.null(jagsui$overlap0)){
      paste("Warning: jagsUI did not calculate statistics")
      df_mod <- df_sum %>%
        select("var",starts_with("ID"),"mean","median","pc2.5","pc97.5") %>% #reorder columns, drop ID
        mutate(overlap0 = NA, gel = NA)
    }else{
      df_mod <- df_sum %>%
        select("var",starts_with("ID"),"mean","median","pc2.5","pc97.5") %>% #reorder columns, drop ID
        mutate(overlap0 = do.call(c, jagsui$overlap0), gel = do.call(c, jagsui$Rhat))
    }
    
    df_mod[,2:(counter+1)] <- lapply(2:(counter+1), function(x) as.numeric(df_mod[[x]])) # make appropraite columns numeric
    
    return(df_mod)
  }
  
  if(type == "rjags"){
    jm_coda <- jagsobj
    
    # Organize the coda object as a dataframe
    df_sum <- coda.fast(jm_coda)
    df_sum <- tibble::rownames_to_column(df_sum, "var")
    df_sum <- df_sum %>% # make index column
      mutate(ID = sub('.*\\[(.*)\\]', '\\1', df_sum$var))
    df_sum$ID <- ifelse(grepl('[[:alpha:]]', df_sum$ID), 1, df_sum$ID) # make ID=1 if there is only 1 instance
    
    # make a lists of list of indices
    IDlist <- strsplit(df_sum$ID, ",") #temp
    
    # get number of dims in ID
    counter <- 1
    for(i in 1:length(IDlist)){
      counter <- ifelse(length(IDlist[[i]])>counter, length(IDlist[[i]]), counter)
    }
    
    # create a character vector of column names based on max dim
    new_columns <- list()
    for(i in 1:counter){
      new_columns[[i]] <- paste("ID",i, sep="")
    }
    new_columns <- as.character(new_columns)
    
    # for each dimension, create a new column with the correct ID
    df_sum <- df_sum %>%
      tidyr::separate(ID,new_columns,sep=",")
    
    df_sum <- df_sum %>%
      mutate(var = sub('(.*)\\[.*', '\\1', df_sum$var)) # get rid of numbers in var col
    
    df_mod <- df_sum %>%
      select("var",starts_with("ID"),"mean","median","pc2.5","pc97.5") #%>% #reorder columns, drop ID
    #mutate(overlap0 = do.call(c, jagsui$overlap0), gel = do.call(c, jagsui$Rhat))
    
    df_mod[,2:(counter+1)] <- lapply(2:(counter+1), function(x) as.numeric(df_mod[[x]])) # make appropraite columns numeric
    
    return(df_mod)
  }
  
}

lowdevrestart <- function(saved_state, vary_by = 10){
  
  initlow <- saved_state[[3]] # initlow is just the lowest dev chain number
  
  # take chain with lowest deviance, and make remaining chains vary around it
  saved_state[[2]][[1]] = saved_state[[2]][[initlow]] # Best (low dev) initials for chain 1
  saved_state[[2]][[2]] = lapply(saved_state[[2]][[initlow]],"*",vary_by)
  saved_state[[2]][[3]] = lapply(saved_state[[2]][[initlow]],"/",vary_by)
  
  return(saved_state)
  
}

findlowdev <- function(codaobj){
  
  if(length(codaobj)>3){
    codaobj <- codaobj$samples
  }
  
  jm_coda <- codaobj
  # Save inits based on chains with lowest deviance
  dev_col <- which(colnames(jm_coda[[1]]) == "deviance")
  dev1<- mean(jm_coda[[1]][,dev_col])
  dev2<- mean(jm_coda[[2]][,dev_col])
  dev3<- mean(jm_coda[[3]][,dev_col])
  dev_min <- min(dev1, dev2, dev3)
  if(dev1 == dev_min){
    devin = 1
  } else if(dev2 == dev_min){
    devin = 2
  } else if(dev3 == dev_min){
    devin = 3
  }
  
  initlow <- devin
  print(paste("chain with lowest deviance: ", initlow, sep=""))
  
  
  return(initlow) #returns the number of the lowest deviance chain
  
}

keepvars <- function(codaobj, to_keep, paramlist, type){
  
  if(!is.null(codaobj$samples)){
    codaobj <- codaobj$samples
  }
  
  # Create a "not in" function using negate from the purrr package
  `%nin%` <- purrr::negate(`%in%`)
  
  remove_vars <- get_remove_index(to_keep, paramlist, type)
  
  newinits <- initfind(codaobj, OpenBUGS = FALSE)
  saved_state <- removevars(initsin = newinits,
                            variables = remove_vars)
  
  initlow <- findlowdev(codaobj) # find the lowest deviance chain
  saved_state[[3]] <- initlow
  names(saved_state[[3]]) <- "lowdevchain"
  
  return(saved_state)
  
}



# function that connects date or doy timeseries to coda output, for use after dumsum
# part of the coda4dummies package
# Inputs: 
# dfobj: df you want to connect the dates to
# datevect: vector of dates or doys
# datename: what you want the date column to be called
# identifier: the optional identifier column name from the dfobj dataframe, how the dates will connect if the variable name is the same with a different index
# varlist: variable names in the dfobj you want to connect dates to
# Output: dataframe

#temp
# dfobj = df_yat_WUE_SIF
# varlist = c("WUE.pred")
# datevect = YIN$DOY
# datename = "DOY"
# identifier = NULL
dateconnect <- function(dfobj, datevect, datename = "date", identifier = NULL, varlist){
  
  dflistj <- list() # create an empty list for output
  for(j in 1:length(varlist)){
    dfobj2 <- dfobj %>%
      filter(var == varlist[j])
    
    if(!is.null(identifier)){ # if the variables we want to connect to dates have the same names but are indexed
      
      IDlist <- unique(as.vector(dfobj2[,identifier]))
      endloop <- length(IDlist)
      
      dflisti <- list() # create an empty list
      for( i in c(1:endloop)){
        
        dfobj3 <- dfobj2 %>%
          filter(ID2==as.numeric(IDlist[i]))
        
        #dfobj3 <- dfobj3[-c(1,2),] # temp for testing because I changed the processing btwn runs to include PAR
        
        dfobj3[,datename] <- datevect
        dflisti[[i]] <- dfobj3
      }
      dflistj[[j]] <- bind_rows(dflisti)
      
    }else{ # if the variables we want to connect to dates are just named differently
      dfobj3 <- dfobj2
      dfobj3[,datename] <- datevect
      dflistj[[j]] <- dfobj3
    }
  }
  
  dfout <- bind_rows(dflistj)
  
  if(datename %in% colnames(dfobj)){ # if the column name already exists

    dfout[,"tempcol"] <- dfout[,datename] # define a temp column name
    dfout <- dfout[,!names(dfout)==datename] # remove the datename column, or else it messes with joining
    
    dfoutout <- left_join(dfobj, dfout)
    
    # when the datename column has an NA, fill it with the tempcol
    for(k in 1:nrow(dfoutout)){
      if(is.na(dfoutout[k,datename])){
        dfoutout[k,datename] <- dfoutout$tempcol[k]
      }
    }
    
    dfoutout <- dfoutout %>% # delete the temp column
      select(-tempcol)
    
  }else{
    dfoutout <- left_join(dfobj, dfout)
  }
  
  return(dfoutout)

}



# function that finds the index number for variables from rjags or jagsUI output
# intended to work with Michael Fell's "removevars" function
# Inputs: 
  # to_keep: string of variable names
  # list: list of all parameters tracked
  # type: rjags or jagsUI
# Output: list of index values for all variables NOT included in the input
get_remove_index <- function(to_keep, list, type){
  
  if(type %nin% c("rjags", "jagUI")){
    paste("Please indicate whether this is a rjags or jagsUI samples object")
  }
  
  if(type == "rjags"){
    list <- list[list != "deviance"] # remove deviance
    list <- sort(list, method = "radix")
    out_list <- c()
    for(j in c(1:length(list))){
      if(list[j] %in% to_keep){
        out_list[j] = NA
      } else{
        out_list[j] = j
      }
    }
    out_list <- out_list[!is.na(out_list)]
    out_list
  }
  
  if(type == "jagsUI"){
    list <- list[list != "deviance"] # remove deviance
    out_list <- c()
    for(j in c(1:length(list))){
      if(list[j] %in% to_keep){
        out_list[j] = NA
      } else{
        out_list[j] = j
      }
    }
    out_list <- out_list[!is.na(out_list)]
    out_list
  }
  
}


###############################################################################
# A function that does pairwise comparisons among defined outputs from JAGS
# written by Biplabendu (Billu) Das 
# https://biplabendu.github.io/homepage/
# @another_ant_guy
###############################################################################

# requires multcompView, pacman
library(multcompView)
library(pacman)
bd.check_similarity <- function(data, params, which.params, m, l, u, anchor=NULL, ids) {
  pacman::p_load(tidyverse)
  foo <- data;
  if(!is.null(anchor)) {
    all.letters <- list()
    for(i in 1:length(which.params)) {
      bar1 <-
        foo %>% 
        select({{params}},{{m}},{{l}},{{u}},{{anchor}},{{ids}})
      bar1 <- bar1[which(bar1[1]==which.params[i]),]
      colnames(bar1) <- c("p","m","l","u","anchor","id")
      bar2 <-
        bar1 %>% 
        inner_join(bar1, by=c("anchor")) %>% 
        mutate(c = ifelse(m.x==m.y, "ns",
                          ifelse(m.y>=l.x & m.y<=u.x | m.x>=l.y & m.x<=u.y, 
                                 "ns","sig"))) %>% 
        group_by(anchor, id.x, id.y) %>% 
        mutate(ids = paste(sort(c(id.x,id.y)),collapse = "-")) %>% 
        ungroup() %>% 
        select(anchor, ids , c) %>% 
        distinct() %>% 
        mutate(id2=ids) %>% 
        separate(id2, c("id1","id2"), "-")
      
      which.anchor <- unique(bar1$anchor)
      which.ids <- unique(bar1$id)
      l.letters <- list()
      for (j in 1:length(which.anchor)) {
        dif3 <- bar2 %>% filter(anchor==which.anchor[j]) %>% mutate(c=ifelse(c=="ns",F,T)) %>% pull(c)
        names(dif3) <- bar2 %>% filter(anchor==which.anchor[j]) %>% pull(ids)
        dif3L <- multcompView::multcompLetters(dif3)
        l.letters[[j]] <-
          data.frame(id = which.ids,
                     anchor = which.anchor[j],
                     params=which.params[i],
                     letters = dif3L$Letters)
        colnames(l.letters[[j]]) <- c(ids, anchor, params, "letters")
      }
      all.letters[[i]] <- do.call(rbind,l.letters)
    }
    return(do.call(rbind,all.letters))
    
  } else {
    writeLines("Please specify an anchor column.")
  }  
}


###############################################################################
# Written for work on the Multicomp project
# Updated by Michael Fell 9/10/2018
#   Added an option for an arbitrary function
#   Added more informative error messages
###############################################################################

# A function to summarize output from a JAGS or OpenBUGS model.
coda.fast <- function(coda=NULL, thin=1, FUN=NULL, colname = "optfun", ...){
  
  if(is.null(coda)){
    message("No coda object provided. Summarizing nothing is too philosophical")
    message("a task for this function.")
    stop()
  }
  
  # Get the number of chains
  chains <- length(coda)
  
  codal <- length(coda[[1]][,1])
  
  # Combine chains
  Ftab <- numeric()
  for(i in 1:chains){
    Ftab <- rbind(Ftab, coda[[i]][(0:(codal-1))%%thin==0,])
  }
  
  # mean, sd, 95% CrI table
  pred <- matrix(nrow=dim(Ftab)[2], ncol=5)
  colnames(pred)<-c("mean", "median", "sd","pc2.5","pc97.5")
  
  # Fill table with stats
  pred[,1] <- colMeans(Ftab) #calculate predicted mean RW's
  pred[,2] <- apply(X=Ftab, MARGIN=2, FUN=median, na.rm=TRUE)
  pred[,3] <- apply(X=Ftab, MARGIN=2, FUN=sd,na.rm=TRUE) #calculate sd, apply column wise
  pred[,4] <- apply(X=Ftab, MARGIN=2, FUN=quantile,probs=c(0.025),na.rm=TRUE) 
  pred[,5] <- apply(X=Ftab, MARGIN=2, FUN=quantile,probs=c(0.975),na.rm=TRUE)
  
  pred <- data.frame(pred)
  if(length(rownames(pred)) == length(colnames(coda[[1]]))){
    rownames(pred) <- colnames(coda[[1]])
  }else{
    message("Could not determine row (variable) names from coda.")
  }
  
  # Optional Function
  if(!is.null(FUN))
  {
    placeholder <- tryCatch(
      {
        out <- apply(X=Ftab, MARGIN=2, FUN=FUN, na.rm=TRUE, ...)
        out <- as.matrix(out)
        if(ncol(out) == nrow(pred)){
          out <- t(out)
        }
        
        pred <- cbind(pred, out)
        colnames(pred) <- c("mean", "median", "sd","pc2.5","pc97.5", colname)
      },
      error=function(cond){
        message(paste0("A problem led to an error executing the optional function."))
        message("The result without the added function will be returned.")
        message("Here is the original error:")
        message(cond)
      },
      warning=function(cond){
        message("A warning occurred executing the optional function.")
        message("The result without the added function will be returned.")
        message("Here is the original warning:")
        message(cond)
      },
      finally={
        return(pred)
      }
    )
  }
  
  # Return the summary values
  return(pred)
}

# A function to find initial values for a JAGS or OpenBUGS model.
# Output:
# The output from this function is a list containing two elements. The first
# contains the names of the variables and their indicies. These are useful 
# when using removevars to remove variables that don't need initial values
# in JAGS. The second element contains a list of initial values (this is a 
# list of lists).

initfind <- function(coda, iteration=0, OpenBUGS=FALSE){
  mcmcin <- coda # TODO change all mcmcin to coda in the future MKF 11/27/18
  # If mcmc.list convert to mcmc
  if(is.mcmc.list(mcmcin)==TRUE){
    mcmcin <- mcmc(data=mcmcin, thin=1)
  }
  
  # Get the number of chains
  n.chains <- length(mcmcin)
  
  # get variable names from a list
  var.names <- colnames(mcmcin[[1]])
  var.dims <- dim(mcmcin[[1]])
  if(iteration==0){
    iteration <- var.dims[1]
  }
  
  if(sum(var.names=="deviance")>0){
    var.names <- var.names[-which(var.names=="deviance")]
    var.dims[2] <- var.dims[2]-1 # correct for removing deviance
  }
  
  # Get names and dimension of each variable since the output is a table
  var.names2 <- apply(X=as.matrix(var.names), MARGIN=c(1), FUN=strsplit, split="\\x5B", perl=TRUE)
  var.names2 <- lapply(X=var.names2, FUN=unlist)
  var.names2 <- lapply(X=var.names2, FUN=gsub, pattern="\\x5D", replacement="", perl=TRUE)
  
  # Create a table of names and dimensions
  # Column 1 is the variable me column 2 has the dimensions
  var.info <- matrix(nrow=length(var.names), ncol=3)
  for(i in 1:length(var.names2)){
    if(length(var.names2[[i]]) > 1){
      var.info[i,] <- c(var.names2[[i]], var.names[i])
    }else if(length(var.names2[[i]]) == 1){
      var.info[i,] <- c(var.names2[[i]], 1, var.names[i])
      #print(i)
      #print(var.names2[[i]])
    }else{
      stop("A variable name has incorrect dimensions for parsing.") 
    }
  }
  
  # Get variable names
  unique.names <- unique(var.info[,1])
  initsoutall <- list()
  
  
  for(k in 1:n.chains){
    initsout <- list()
    for(i in 1:length(unique.names)){
      sel <- which(var.info[,1]==unique.names[i])
      #sel2 <- grep(pattern=paste0("^",unique.names[i],"\\x5B"), x=var.names)
      
      # Make sure the above selections worked
      #if(length(sel) != length(sel2)){
      #  stop("Error matching unique variable names with MCMC output")  
      #}
      name.sel <- var.info[sel,3]
      index <- apply(X=as.matrix(var.info[sel,2]), MARGIN=1, FUN=strsplit, split=",", perl=TRUE)
      index <- lapply(X=index, FUN=unlist)
      index <- matrix(data=as.numeric(unlist(index)), nrow=length(index), ncol=length(index[[1]]), byrow=TRUE)
      
      # There are possibly easier ways to do this but they make more assumptions
      dims <- as.numeric(apply(X=index, MARGIN=2, FUN=max))
      variable <- array(data=NA, dim=dims)
      
      # Fill the new variable with the correct values
      for(j in 1:dim(index)[1]){
        # The output into mcmc objects lists names in the order R stacks them
        # in arrays so the single index for the variable references the 
        # correct array location.
        variable[j] <- mcmcin[[k]][iteration, which(colnames(mcmcin[[k]])==name.sel[j])]
      }
      
      # Use dims to produce a new array to store the data
      initsout[[i]] <- variable
    } # End of variable loop
    names(initsout) <- unique.names
    initsoutall[[k]] <- initsout
  } # End of chain loop
  
  listout <- list(unique.names, initsoutall)
  names(listout) <- c("variables", "initials")
  
  # Account for OpenBUGS by outputing 1 dimensional arrays as vectors.
  if(OpenBUGS==TRUE){
    for(i in 1:n.chains){
      for(j in 1:length(listout[[2]][[i]])){
        if(length(dim(listout[[2]][[i]][[j]]))==1){
          listout[[2]][[i]][[j]] <- as.vector(listout[[2]][[i]][[j]])
        }
      }
    }
  }
  
  return(listout)
} # End of function

###############################################################################
#
# Removes specific variables from the initial values
#
###############################################################################

# A function to remove variables that don't need initial values in JAGS.

removevars <- function(initsin, variables){
  n.chains <- length(initsin[[2]])
  n.vars <- 1:length(initsin[[1]])
  n.vars <- n.vars[-variables]
  
  var.names <- initsin[[1]][n.vars]
  
  new.inits <- list()
  for(i in 1:n.chains){
    chain.inits <- list()
    for(k in 1:length(n.vars)){
      chain.inits[[k]] <- initsin[[2]][[i]][[n.vars[k]]] 
    } # End of k loop
    names(chain.inits) <- var.names
    new.inits[[i]] <- chain.inits
  } # End of i loop
  
  output <- list(var.names, new.inits)
  names(output) <- c("variables", "initials")
  
  return(output)
  
} # End of function
