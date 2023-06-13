
library(purrr)
# Create a "not in" function using negate from the purrr package
`%nin%` <- negate(`%in%`)



# function that finds the index of variables to remove
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
#
# A function that does pairwise comparisons among defined outputs from JAGS
#
###############################################################################

# A function that does pairwise comparisons among defined outputs from JAGS
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




# Function to calculate evaporation, based on descriptions of methods in Merlin et al. 2016
# To run the get_evap function we need columns of:
# PA, Tair, Z, Tsoil, ws, RH, fsand, fclay, fc, SWC_shall
# Take the high and low of the last three days and calculate the average temperature. 
get_evap <- function(site) {
  
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
  pair <- Pa/(Rs * (site$Tair + 273.15)) # density of air, calculated using the ideal gas law in (kg/m^3)
  
  # Calculate psychrometric contstant
  Cp <- 1000 # specific heat of air J/kg * K(or degree C)
  lambda <- 22.6 * 10^5 # latent heat of water vaporization J/kg
  gamma <- (Cp*Pa)/(rmw*lambda) #psychrometric constant in Pa/K
  
  # Calculate aerodynamic resistance to heat transfer
  # 0.001 is Z0m (m), the momentum soil roughness, [Yang et al., 2008; Stefan et al., 2015]
  Z <- site$Z # reference height (m)
  Ri <- (5*g*Z*(site$Tsoil-site$Tair))/(site$Tair*(site$ws**2)) # Richardson number, represents the importance of free versus forced convection
  # a coefficient set to 0.75 in unstable conditions (T > Ta) and to 2 in stable conditions (T < Ta)
  eta_unstable <- 0.75
  #eta_stable <- 2
  rah0 <- (1/((k**2)*site$ws))*((log(Z/0.001))**2) # s/m
  rah_unstable <- rah0/((1 + Ri)**eta_unstable) # aerodynamic resistance to heat transfer s/m, estimated as in Choudhury et al. [1986]
  #rah_stable <- rah0/((1 + Ri)**eta_stable)
  rah_stable <- ((1 - 0.75*Ri)*((log(Z/0.001))**2))/((k**2)*site$ws)
  rah <- ifelse( Ri > 0, rah_unstable, rah_stable)
  
  # Vapor pressure calculations
  ea1 <- A * 10**((m*site$Tair)/(site$Tair+Tn)) # Temperature from data; A, m, and Tn are constants. Saturated.
  ea  <- ea1 * site$RH/100 # air vapor pressure, should be in Pa
  es <- A * 10**((m*site$Tsoil)/(site$Tsoil +Tn)) # Saturated vapor pressure in soil
  
  # Snow vapor pressure
  #Tsnow = (site$Tsoil + site$Tair)/2
  #es_snow =  A * 10**((m*Tsnow)/(Tsnow +Tn)) # Saturated vapor pressure in snow
  
  # calculate parameters for alternative alpha and Bowen
  psisat <- -10*exp(1.88-1.31*site$fsand) # parameterized air entry pressure, in mm of water
  bch <- 2.91 + 15.9*site$fclay # The Clapp and Hornberger parameter est. as in Cosby et al. [1984]
  sres <- 0.15*site$fclay # residual soil moisture
  ssat <- 0.489 - (0.126*site$fsand)
  rssmin <- 50 # minimum soil resistance s/m
  psi <- psisat * (site$SWC_shall/ssat)**(-bch)
  
  # calculate soil resistance
  rss1 <- ((site$fc - sres)/(site$SWC_shall - sres)) * rssmin # soil resistance for S > sres
  rss2 <- exp(8.206 - 4.255*(site$SWC_shall/site$fc)) # another soil resistance derivation
  rss <- ifelse(site$SWC_shall>sres, rss1, rss2) # the resistance to the diffusion of vapor in large soil pores
  
  # Calculate alpha and Bowen ratios
  #alpha1 <- ifelse(site$SWC_shall > site$fc, 1, 0.5 - 0.5*cos((site$SWC_shall * pi)/site$fc)) # wetness function
  alpha2 <- exp((psi*g)/(10000 * Rwv * site$Tsoil)) # Kelvin equation
  bowen1 <- ifelse(site$SWC_shall > site$fc, 1, (0.5 - 0.5*cos((site$SWC_shall * pi)/site$fc))**2) # wetness function
  
  # This will calculate LE in W/m^2, so let's convert to mm (i.e. LE to E)
  # 1 Watt /m2 = 0.0864 MJ /m2/day
  # 1 MJ /m2/day  = 0.408 mm /day
  E4.5 <- (bowen1 * ((pair * Cp)/gamma) * ((alpha2*(es - ea))/rah)) * 0.0864 * 0.408 # from 4.5 CLM
  E3.5 <- (((pair * Cp)/gamma) * ((alpha2*(es - ea))/(rah+rss))) * 0.0864 * 0.408 # from 3.5 CLM

  # Calculate intercepted Eint using LAI from MODIS
  # intercepted E (q_intr) in kg m-2 s-1
  # Eint = (site$P)*(1 - exp(-0.5*(site$LAI_mod)))
  
  site_E <- site %>%
    mutate(E4.5 = E4.5, E3.5 = E3.5, pair = pair, Bowen = bowen1, alpha = alpha2, gamma = gamma, Ri = Ri, rah_unstable = rah_unstable, rah = rah, rss = rss, es = es, ea = ea, Pa = Pa)
  
  

  
  return(site_E)
}