#### Function to run JAGS model depending on site and response variable

run_mod2 <- function(dataIN, varname, sitename, Z, h, fc, fclay, fsand, newinits = F, overwrite = F, lowdev = F, post_only = F){
  
  # Temp for testing, comment out before using the function
  if(test==T){
    dataIN = dat
    overwrite = F
    lowdev = F
    newinits = T
  }
  
  # Since this is the DEPART model
  if(varname == "WUE_GPP"){ # varname starts at 4 to be consistent with model 1
    photname = "GPP"
  } else if(varname == "WUE_SIF"){
    photname = "SIF_O2A"
  }
  
  # Define filenames
  modelname <- "./models/model2.R"
  initfilename <- paste("./output_model2/inits/", sitename, "/inits_", varname, ".RData", sep = "")
  initlowfilename <- paste("./output_model1/inits/", sitename, "/initlow_", varname, ".RData", sep = "")
  jm_codafilename <- paste("./output_model2/coda/", sitename, "/jm_coda_", varname,".RData", sep = "")
  mcmcfoldername <- paste("./output_model2/convergence/", sitename, "/", varname, sep = "")
  mcmcfilename <- paste("MCMC_", varname, sep = "")
  dffilename <- paste("./output_model2/df/df_", sitename, "_", varname, ".csv", sep = "")
  
  # Create necessary directories if they do not already exist
  if(!dir.exists(paste("output_model2/inits", sep = ""))){ dir.create(paste("output_model2/inits", sep = ""))}
  if(!dir.exists(paste("output_model2/coda", sep = ""))){ dir.create(paste("output_model2/coda", sep = ""))}
  if(!dir.exists(paste("output_model2/convergence", sep = ""))){ dir.create(paste("output_model2/convergence", sep = ""))}
  if(!dir.exists(paste("output_model2/inits/", sitename, sep = ""))){ dir.create(paste("output_model2/inits/", sitename, sep = ""))}
  if(!dir.exists(paste("output_model2/coda/", sitename, sep = ""))){ dir.create(paste("output_model2/coda/", sitename, sep = ""))}
  if(!dir.exists(paste("output_model2/convergence/", sitename, sep = ""))){ dir.create(paste("output_model2/convergence/", sitename, sep = ""))}
  if(!dir.exists(mcmcfoldername)){ dir.create(mcmcfoldername)}
  if(!dir.exists(paste("output_model2/df", sep = ""))){ dir.create(paste("output_model2/df", sep = ""))}

  #####################################################################
  #Part 1: Model setup
  
  # If we are missing env variables we don't want to gapfill 5 timesteps into the past,
  # make ET NA to skip over those periods
  # we will filter out the NAs in ET later in our YIN df to skip over those rows in our weight calculations
  ntimestep = 8 # look at C1 and C2 to see the timesteps
  for(i in c(1:nrow(dataIN))){
    if(i < (ntimestep+2)){
      next
    }
    for(j in c(0:(ntimestep+1))){
      if(is.na(dataIN$SWC_shall[i-j])){
        dataIN[i,photname] = NA
      }
      if(is.na(dataIN$SWC_deep[i-j])){
        dataIN[i,photname] = NA
      }
      if(is.na(dataIN$PAR[i-j])){
        dataIN[i,photname] = NA
      }
    }
    
  }
  
  # If LAI is NA, have the model ignore it, since we've already gap-filled what is possible
  dataIN$LAI <- ifelse(is.na(dataIN$LAI), 0, dataIN$LAI)
  
  # Make the response variable file of interest. This file includes 
  # growing-season response (Y) variables of interest along with indices 
  # that link the Y variables back to covariates of interest.
  # NOTE: Kym did this to only test growing season response variables (4/1 to 10/31 for each year)
  YIN = dataIN %>%
    rowid_to_column("ind") %>% # ind: index to link the growing season Y variables back to the appropriate row in the covariate data set
    mutate(time = as.numeric(format(TIMESTAMP, "%H%M"))) %>%
    filter(between(time, 730, 1600)) %>% # filter from 7:30 am to 4:00 pm everyday
    filter(!is.na(SIF_O2A)) %>% # filter any gaps for SIF out that we don't want to fill
    filter(!is.na(ET)) %>%
    filter(!is.na(GPP)) %>%
    filter(!is.na(TA)) %>% # filter any weird start/end gaps for envir variables, this won't cause gaps in the middle of the day, since we already gap-filled
    filter(!is.na(VPD)) %>%
    filter(!is.na(SWC_shall)) %>%
    filter(!is.na(PAR)) %>%
    filter(!is.na(WS)) %>% # filter away any evap equation-dependent variables we are missing or do not want to gapfill
    filter(!is.na(Ri)) %>%
    filter(!is.na(ea)) %>%
    filter(!is.na(es))
    
  # Change Y to the column for the response variable of interest
  Phot  = YIN[,photname]
  Yday = YIN$ind
  
  # jIND file provides indices to calculate interactions between covariates
  # Basically a matrix version of X2, defined in the first initialization later in the script
  # key: 1 VPD # 2 Tair # 3 Sshall # 4 Sdeep
  # X1a = cbind(X1[,1]^2, X1[,2]^2)
  # X2 = X1[,1]*X1[,2], X1[,1]*X1[,3], X1[,1]*X1[,4], X1[,2]*X1[,3], X1[,2]*X1[,4], X1[,3]*X1[,4], X1[,5]*X1[,1], X1[,5]*X1[,2], X1[,5]*X1[,3], X1[,5]*X1[,4]
  jIND <- data.frame(j = c(1:10),
                     ID1 = c(1,1,1,2,2,3,5,5,5,5),
                     ID2 = c(2,3,4,3,4,4,1,2,3,4))
  
  # an example from a previous SAM version:
  # X2  = cbind(X1[,1]*X1[,2], X1[,1]*X1[,3], X1[,1]*X1[,5], X1[,1]*X1[,6], X1[,2]*X1[,3], X1[,2]*X1[,5], X1[,2]*X1[,6], 
  # X1[,3]*X1[,5], X1[,3]*X1[,6], X1[,5]*X1[,6])
  # jIND <- data.frame(j = c(1:10),
  #                    ID1 = c(1,1,1,1,2,2,2,3,3,5),
  #                    ID2 = c(2,3,5,6,3,5,6,5,6,6))
  
  # Choose the starting index. This is an index for a row in the Y data file. 
  # The value in the indexed row should be greater than 1 to accommodate 
  # calculation of antecedent values.
  Nstart = 1
  # Choose the ending index. 
  Nend   = nrow(YIN)
  
  # Prepare data for JAGS -- covariates are scaled
  data = list(
    # ET partitioning (DEPART) model data
    Tsoil = YIN$TS,
    S = YIN$SWC_shall, # unscaled SWC for evap
    S.min = min(YIN$SWC_shall), # minimum SWC for priors, alternative for some if else statements
    ET = YIN$ET,
    Phot = Phot,
    ws = YIN$WS,
    conv.fact = 0.0864 * 0.408,
    rho = YIN$pair,
    Ri = YIN$Ri,
    rah_unstable = YIN$rah, # in the data, rah is just rah_unstable when appropriate. We are letting rah_stable vary (stochastic) so that's why we read this in.
    Cp = 1000,
    pi = 3.14159265359,
    Rwv = 461.52, # specific gas constant for water vapor in J/(kg*K)
    sres = 0.15*fclay, # residual soil moisture
    ssat = 0.489 - (0.126*fsand),
    psisat = -10*exp(1.88-1.31*fsand), # parameterized air entry pressure, in mm of water
    rssmin = 50, # minimum soil resistance in s/m
    g  = 9.80665,     # acceleration due to gravity in m/s^2
    gamma = YIN$gamma,
    #rah = dataIN$rah,
    e.sat = YIN$es,
    e.a = YIN$ea,
    Z = Z, # reference height (m)
    #fclay = dataIN$fclay[1],
    fc = fc,
    bch = 2.91 + 15.9*fclay, # The Clapp and Hornberger parameter est. as in Cosby et al. [1984]
    LAI = YIN$LAI,
    P = YIN$P,
    
              
              # SAM model data
              Nstart = Nstart,
              Nend = Nend,
              Nlag = 9,
              Nparms = 5, # Nparms is the number of driving variables included to calculate main effects
              Yday = YIN$ind, # Choose column in YIN that provides indices linking response variables with covariates
              ID1 = jIND[,2], 
              ID2 = jIND[,3],
              jlength = nrow(jIND),
              VPD = as.vector(scale(as.numeric(dataIN$VPD),center=TRUE,scale=TRUE)), # scale function takes vector of values, centers and scales by SD
              Tair = as.vector(scale(as.numeric(dataIN$TA),center=TRUE,scale=TRUE)),
              Sshall = as.vector(scale(as.numeric(dataIN$SWC_shall),center=TRUE,scale=TRUE)),
              Sdeep = as.vector(scale(as.numeric(dataIN$SWC_deep),center=TRUE,scale=TRUE)),
              PAR = as.vector(scale(as.numeric(dataIN$PAR),center=TRUE,scale=TRUE)),
              # covariate timesteps into the past
              # this code is flexible in case we want to combine timesteps,
              # or average over multiple timesteps
              C1 = c(0, 1, 2, 3, 4, 5, 6, 7, 8), #stop times for covariates
              C2 = c(0, 1, 2, 3, 4, 5, 6, 7, 8) #start times for covariates
              )

  
  # Get initials for WUE precision based on ET
  tau.Y = 1/(sd(YIN$ET)**2)
  
  # Initial values are estimated using a linear model. As in the data list (above), 
  # covariates are centered and standardized. Replace name of covariate in quotes
  # with the appropriate column number in the dataIN file.
  X1  = cbind(data$VPD[Yday[Nstart:Nend]], #1
              data$Tair[Yday[Nstart:Nend]], #2
              data$Sshall[Yday[Nstart:Nend]], #3
              data$Sdeep[Yday[Nstart:Nend]], #4
              data$PAR[Yday[Nstart:Nend]] #5
  ) 
  
  # Notes: The code below is indexed numerically, which you will have to pay attention to as you change covariates of interest
  # Squared terms calculated for VPD and Tair
  X1a = cbind(X1[,1]^2, X1[,2]^2) 
  # Put all covariates together;
  # Interactions incorporated into linear model used to estimate initial values
  X2  = cbind(X1[,1]*X1[,2], X1[,1]*X1[,3], X1[,1]*X1[,4], X1[,2]*X1[,3], X1[,2]*X1[,4], 
              X1[,3]*X1[,4], X1[,5]*X1[,1], X1[,5]*X1[,2], X1[,5]*X1[,3], X1[,5]*X1[,4])
  # Fit simple linear model
  WUE.Y <- YIN$GPP/YIN$ET
  fit <- lm(WUE.Y ~ X1[,1] + X1[,2] + X1[,3] + X1[,4] + X1[,5] + # main effects
              X1a[,1] + X1a[,2] + # squared
              X2[,1] + X2[,2] + X2[,3] + X2[,4] + X2[,5] + X2[,6] + X2[,7] + X2[,8] + X2[,9] + X2[,10]) # interactions
  # Extract coefficient estimates:
  beta0  = fit$coefficients[1] # the intercept
  beta1  = fit$coefficients[2:6] # main effects
  beta1a = fit$coefficients[7:8] # squared effects
  beta2 = fit$coefficients[9:18] # interactive effects
  
  # Create initials based on the above estimates:
  inits = list(list(beta0 = beta0, beta1 = beta1, beta1a = beta1a, beta2 = beta2, tau.Y = tau.Y), #pink
               list(beta0 = beta0/10, beta1 = beta1/10, beta1a = beta1a/10, beta2 = beta2/10, tau.Y =  tau.Y/5), #green
               list(beta0 = beta0*10, beta1 = beta1*10, beta1a = beta1a*10, beta2 = beta2*10, tau.Y =  tau.Y*5)) #blue
  
  # Initial values: from saved state (if we already ran the model, we want to start where we left off)
  if(newinits==F){
    if(file.exists(initfilename)){
      load(initfilename)
      
      if(lowdev == T){
        if(file.exists(initlowfilename)){
          load(initlowfilename) # initlow object is just the lowest dev chain number, 1,2, or 3
          
          # take chain with lowest deviance, and make remaining chains vary around it
          saved_state[[2]][[1]] = saved_state[[2]][[initlow]] # Best (low dev) initials for chain 1
          saved_state[[2]][[2]] = lapply(saved_state[[2]][[initlow]],"*",2)
          saved_state[[2]][[3]] = lapply(saved_state[[2]][[initlow]],"/",2)
        }
      }
      
      initslist <- saved_state[[2]]
      
    }else if(!file.exists(initfilename)){
      initslist <- inits
    }
  }else if(newinits==T){
    initslist <- inits
  }
  
  #####################################################################
  # Part 2: Run JAGS Model (jagsui initializes and runs in one step)
  
  # parameters to track
  params = c("deviance", # deviance
             "WUE.pred", "T.ratio", "T.pred", "ET.pred", # variables associated with DEPART
             "beta0","beta1","beta1a","beta2", # intercept, main effects, squared effects, interactive effects
             "dYdX", # net sensitivities
             "tau.Y", # variance of Y
             "wV","wT","wSs","wSd","wPAR", # importance weights
             "R2", "Dsum") # model fit
  
  if(post_only==F){ # if we want to run the model
    
    # Run model with jagsui package
    start<-proc.time() # keep track of run time
    jagsui <- jags(data = data,
                   inits = initslist,
                   model.file = modelname,
                   parameters.to.save = params,
                   n.chains = 3,
                   n.adapt = 500,
                   n.thin = ifelse(test==F,10,1),
                   n.iter = ifelse(test==F,70000,200),
                   parallel = ifelse(test==F,TRUE, FALSE))
    end<-proc.time()
    elapsed<- (end-start)/60
    print("jagsui done running; minutes to completion:")
    print(elapsed[3])
    
    if(overwrite==T){
      save(jagsui, file = jm_codafilename) # save out whole jagsui object for local
    }
    
  } # end post_only==F
  
  if(post_only==T){ # if we just want to load the model to summarize it differently
    load(jm_codafilename) # named jagsui
  }
  
  jm_coda <- jagsui$samples # convert to coda form to work with postjags functions
  
  df_mod <- dumsum(jagsui, type = "jagsUI") # organize the coda object as a dataframe
  
  # Organize the coda object as a dataframe
  # df_sum <- coda.fast(jm_coda)
  # df_sum <- rownames_to_column(df_sum, "var")
  # df_sum <- df_sum %>% # make index column
  #   mutate(ID = sub('.*\\[(.*)\\]', '\\1', df_sum$var))
  # df_sum <- df_sum %>% # separate index column into 1st and 2nd dimension
  #   mutate(ID1 = sub('(.*)\\,.*', '\\1', df_sum$ID),
  #          ID2 = sub('.*\\,(.*)', '\\1', df_sum$ID))
  # df_sum$ID2 <- ifelse(!grepl(',', df_sum$ID), NA, df_sum$ID2) # get rid of ID2 if there's no 2nd dimension
  # df_sum$ID1 <- ifelse(grepl('[[:alpha:]]', df_sum$ID), 1, df_sum$ID1) # make ID1=1 if there is only 1 instance
  # df_sum <- df_sum %>% 
  #   mutate(var = sub('(.*)\\[.*', '\\1', df_sum$var)) # get rid of numbers in var col
  # df_mod <- df_sum %>% 
  #   select("var","ID1","ID2","mean","median","pc2.5","pc97.5") %>% #reorder columns, drop ID
  #   mutate(overlap0 = do.call(c, jagsui$overlap0), gel = do.call(c, jagsui$Rhat))
  
  df_mod <- dateconnect(dfobj = df_mod, datevect = YIN$DOY, datename = "DOY", identifier = "ID2", varlist = c("dYdX"))
  df_mod <- dateconnect(dfobj = df_mod, datevect = YIN$TIMESTAMP, datename = "TIMESTAMP", identifier = "ID2", varlist = c("dYdX"))
  df_mod <- dateconnect(dfobj = df_mod, datevect = YIN$DOY, datename = "DOY", identifier = NULL, varlist = c("WUE.pred", "T.ratio", "T.pred", "ET.pred"))
  df_mod <- dateconnect(dfobj = df_mod, datevect = YIN$TIMESTAMP, datename = "TIMESTAMP", identifier = NULL, varlist = c("WUE.pred", "T.ratio", "T.pred", "ET.pred"))
  
  write.csv(df_mod, file = dffilename)
  
  #####################################################################
  # Part 3: Check diagnostics
  
  mcmcplot(jm_coda, parms = c("deviance", # deviance
                              "WUE.pred", "T.ratio", "T.pred", "ET.pred", # variables associated with DEPART
                              "beta0","beta1","beta1a","beta2", # intercept, main effects, squared effects, interactive effects
                              "tau.Y", # variance of Y
                              "wV","wT","wSs","wSd","wPAR", # importance weights
                              "R2"), # model fit
           random = 15, # so we'll only save 15 random plots for the longer variables to save space
           dir = mcmcfoldername,
           filename = mcmcfilename)
  
  #####################################################################
  # Part 5: Save inits for future runs
  
  # Save state
  
  if(lowdev == T){
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
    save(initlow, file = initlowfilename)
  }
  
  # inits to save
  init_names = names(initslist[[1]])
  
  # use get_remove_index function to find which variables to remove
  remove_vars = get_remove_index(init_names, params, type="jagsUI")
  
  newinits <- initfind(jm_coda, OpenBUGS = FALSE)
  newinits[[1]]
  saved_state <- removevars(initsin = newinits, 
                            variables = remove_vars)
  
  # if(lowdev == T){
  #   # saved chain with lowest deviance, and make remaining chains vary around it
  #   saved_state[[2]][[1]] = saved_state[[2]][[devin]] # Best (low dev) initials for chain 1
  #   saved_state[[2]][[2]] = lapply(saved_state[[2]][[devin]],"*",2)
  #   saved_state[[2]][[3]] = lapply(saved_state[[2]][[devin]],"/",2)
  # }
  
  #saved_state[[1]]
  
  if(overwrite==T){
    save(saved_state, file = initfilename) #for local
  }
  
  # If converged, run and save replicated data. I'm not running this, because I'm calculating Bayesian R2 inside the model code
  # jm_rep <- update(jagsui, parameters.to.save = "Y.rep",
  #                  n.iter = 15000, n.thin = 5)
  # 
  # if(overwrite==T){
  #   save(jm_rep, file = jm_repfilename) #for local
  # }
  # 
  
  
  
  
} # end of function
