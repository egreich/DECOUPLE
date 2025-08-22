#### Function to run JAGS model depending on site and response variable

run_modDEPARTSAM <- function(dataIN, varname, sitename, Z, h, fc, fclay, fsand, newinits = F, overwrite = F, lowdev = F, post_only = F){
  
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
  modelname <- "./models/modelDEPARTSAM.R"
  initfilename <- paste("./output_modelDEPARTSAM/inits/", sitename, "/inits_", varname,"_hour", ".RData", sep = "")
  jm_codafilename <- paste("./output_modelDEPARTSAM/coda/", sitename, "/jm_coda_", varname,".RData", sep = "")
  mcmcfoldername <- paste("./output_modelDEPARTSAM/convergence/", sitename, "/", varname, sep = "")
  mcmcfilename <- paste("MCMC_", varname, sep = "")
  dffilename <- paste("./output_modelDEPARTSAM/df/df_", sitename, "_", varname, ".csv", sep = "")
  
  # Create necessary directories if they do not already exist
  if(!dir.exists(paste("output_modelDEPARTSAM", sep = ""))){ dir.create(paste("output_modelDEPARTSAM", sep = ""))}
  if(!dir.exists(paste("output_modelDEPARTSAM/inits", sep = ""))){ dir.create(paste("output_modelDEPARTSAM/inits", sep = ""))}
  if(!dir.exists(paste("output_modelDEPARTSAM/coda", sep = ""))){ dir.create(paste("output_modelDEPARTSAM/coda", sep = ""))}
  if(!dir.exists(paste("output_modelDEPARTSAM/convergence", sep = ""))){ dir.create(paste("output_modelDEPARTSAM/convergence", sep = ""))}
  if(!dir.exists(paste("output_modelDEPARTSAM/inits/", sitename, sep = ""))){ dir.create(paste("output_modelDEPARTSAM/inits/", sitename, sep = ""))}
  if(!dir.exists(paste("output_modelDEPARTSAM/coda/", sitename, sep = ""))){ dir.create(paste("output_modelDEPARTSAM/coda/", sitename, sep = ""))}
  if(!dir.exists(paste("output_modelDEPARTSAM/convergence/", sitename, sep = ""))){ dir.create(paste("output_modelDEPARTSAM/convergence/", sitename, sep = ""))}
  if(!dir.exists(mcmcfoldername)){ dir.create(mcmcfoldername)}
  if(!dir.exists(paste("output_modelDEPARTSAM/df", sep = ""))){ dir.create(paste("output_modelDEPARTSAM/df", sep = ""))}

  #####################################################################
  #Part 1: Model setup
  
  # convert umol C to g C for model
  dataIN$GPP <- umolCO2.to.gC(dataIN$GPP, constants = bigleaf.constants())
  
  # If we are missing env variables we don't want to gapfill 5 timesteps into the past,
  # make ET NA to skip over those periods
  # we will filter out the NAs in ET later in our YIN df to skip over those rows in our weight calculations
  ntimestep = 7 # look at C1 and C2 to see the timesteps
  
  C1 = c(0, 1, 2, 3, 4, 5, 6) #stop times for covariates
  C2 = c(0, 1, 2, 3, 4, 5, 6) #start times for covariates
  T1 = c(0, 24, 47, 71, 95, 119, 143) #stop times for covariates over longer timescales
  T2 = c(0, 1, 25, 48, 72, 96, 120) #start times for covariates over longer timescales
  
  for(i in c(1:nrow(dataIN))){
    
    # if GPP or ET is less than 0, assume it is 0 for the purpose of the log scale
    #dataIN$GPP[i] <- ifelse(dataIN$GPP[i]<0, 0.0000001, dataIN$GPP[i])
    #dataIN$ET[i] <- ifelse(dataIN$ET[i]<0, 0.0000001, dataIN$ET[i])
    dataIN$GPP[i] <- ifelse(dataIN$GPP[i]<0, 0, dataIN$GPP[i])
    #dataIN$ET[i] <- ifelse(dataIN$ET[i]<0, 0, dataIN$ET[i])
    
    # ignore time periods where we don't have enough data and we don't want to gapfill
    for(j in C2){
      if(i <= j){
        next
      }
      if(is.na(dataIN$PAR[i-j])){
        dataIN[i,photname] = NA
      }
      if(is.na(dataIN$VPD[i-j])){
        dataIN[i,photname] = NA
      }
    }
    for(j in 1:ntimestep){
      if(i <= T2[j]){
        next
      }
      if(i <= T1[j]){
        next
      }
      if(is.na(mean(dataIN$TA[(i-T1[j]):(i-T2[j])]))){
        dataIN[i,photname] = NA
      }
      if(is.na(mean(dataIN$SWC_shall[(i-T1[j]):(i-T2[j])]))){
        dataIN[i,photname] = NA
      }
      if(is.na(mean(dataIN$SWC_deep[(i-T1[j]):(i-T2[j])]))){
        dataIN[i,photname] = NA
      }
    }
  }
  
  # ntimestep = 6 # look at C1 and C2 to see the timesteps
  # for(i in c(1:nrow(dataIN))){
  #   if(i < (ntimestep+2)){
  #     next
  #   }
  #   for(j in c(0:(ntimestep+1))){
  #     if(is.na(dataIN$SWC_shall[i-j])){
  #       dataIN[i,photname] = NA
  #     }
  #     if(is.na(dataIN$SWC_deep[i-j])){
  #       dataIN[i,photname] = NA
  #     }
  #     if(is.na(dataIN$PAR[i-j])){
  #       dataIN[i,photname] = NA
  #     }
  #   }
  #   
  # }
  
  # If LAI is NA, have the model ignore it, since we've already gap-filled what is possible
  dataIN$LAI <- ifelse(is.na(dataIN$LAI), 0, dataIN$LAI)
  #dataIN$LAI <- 0 # have the model ignore LAI
  
  # Make the response variable file of interest. This file includes 
  # growing-season response (Y) variables of interest along with indices 
  # that link the Y variables back to covariates of interest.
  YIN = dataIN %>%
      rowid_to_column("ind") %>% # ind: index to link the growing season Y variables back to the appropriate row in the covariate data set
      mutate(time = as.numeric(format(TIMESTAMP, "%H%M"))) %>%
      filter(between(time, 700, 1600)) # filter from 7:30 am to 4:00 pm everyday
  
  YIN = YIN %>%
    filter(!is.na(SIF_O2A)) %>% # filter any gaps for SIF out that we don't want to fill
    filter(!is.na(ET)) %>%
    filter(!is.na(GPP)) %>%
    filter(!is.na(TA)) %>% # filter any weird start/end gaps for envir variables, this won't cause gaps in the middle of the day, since we already gap-filled
    filter(!is.na(VPD)) %>%
    filter(!is.na(SWC_shall)) %>%
    filter(!is.na(SWC_deep)) %>%
    filter(!is.na(PAR)) %>%
    filter(!is.na(WS)) %>% # filter away any evap equation-dependent variables we are missing or do not want to gapfill
    filter(!is.na(Ri)) %>%
    filter(!is.na(ea)) %>%
    filter(!is.na(es))
  
  # jIND file provides indices to calculate interactions between covariates
  # Basically a matrix version of X2, defined in the first initialization later in the script
  # key: 1 VPD # 2 Tair # 3 Sshall # 4 Sdeep # 5 PAR
  # X1a = cbind(X1[,1]^2, X1[,2]^2)
  # X2  = cbind(X1[,1]*X1[,2], X1[,1]*X1[,4], X1[,2]*X1[,4], X1[,5]*X1[,1], X1[,5]*X1[,2], X1[,5]*X1[,4])
  jIND <- data.frame(j = c(1:6),
                     ID1 = c(1,1,2,5,5,5),
                     ID2 = c(2,4,4,1,2,4))
  
  # Choose the starting index. This is an index for a row in the Y data file. 
  # The value in the indexed row should be greater than 1 to accommodate 
  # calculation of antecedent values.
  Nstart = 1
  maxtimeback = max(T1) + 1
  for( i in 1:maxtimeback){ # if we don't have enough space at the start of YIN, move down
    nextrow <- i+1
    Nstart = ifelse(YIN$ind[i]<maxtimeback, nextrow, Nstart)
  }
  if(Nstart>1){
    YIN <- YIN[Nstart:nrow(YIN),] # trim YIN to reflect new Nstart
    Nstart = 1 # make Nstart 1 again
  }
  
  # Change Y to the column for the response variable of interest
  Phot  = pull(YIN, photname) #YIN[,photname]
  Yday = YIN$ind
  
  # Choose the ending index. 
  Nend   = nrow(YIN)
  
  
  conv.fact.time <- 3600 # if we need to convert seconds to hours
  #conv.fact.time <- ifelse(scale=="day", 86400, conv.fact.time) # if we need to convert seconds to days
  
  
  # Prepare data for JAGS -- covariates are scaled
  data = list(
    # ET partitioning (DEPART) model data
    Tsoil = YIN$TS,
    S = YIN$SWC_shall, # unscaled SWC for evap
    S.min = min(YIN$SWC_shall), # minimum SWC for priors, alternative for some if else statements
    ET = YIN$ET,
    ET.max = max(YIN$ET),
    Phot = as.vector(Phot),
    ws = YIN$WS,
    conv.fact = 1/((2.501 - 0.00237*YIN$TA)*(10^6)), # 1/the Latent heat of vaporization (J kg-1) to get ET in terms of kg m-2 s-1 from W m-2. 
    conv.fact.time = conv.fact.time,
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
              Nlag = 7,
              Nparms = 5, # Nparms is the number of driving variables included to calculate main effects
              Yday = Yday, # Choose column in YIN that provides indices linking response variables with covariates
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
              C1 = C1, #stop times for covariates
              C2 = C2, #start times for covariates
              T1 = T1, #stop times for covariates over longer timescales
              T2 = T2 #start times for covariates over longer timescales
              )

  
  # Get initials for precision based on ET
  tau.Y = 1/(sd(YIN$ET)**2)
  
  # Initial values are estimated using a linear model. As in the data list (above), 
  # covariates are centered and standardized. Replace name of covariate in quotes
  # with the appropriate column number in the dataIN file.
  X1 = cbind(data$VPD[Yday[Nstart:Nend]], #1
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
  X2  = cbind(X1[,1]*X1[,2], X1[,1]*X1[,4], X1[,2]*X1[,4], 
              X1[,5]*X1[,1], X1[,5]*X1[,2], X1[,5]*X1[,4])
  # Fit simple linear model
  Y <- YIN$GPP/YIN$ET
  Y <- ifelse(is.infinite(Y), NA, Y) # if WUE is infinite, ignore it
  fit <- lm(Y ~ X1[,1] + X1[,2] + X1[,3] + X1[,4] + X1[,5] + # main effects
              X1a[,1] + X1a[,2] + # squared
              X2[,1] + X2[,2] + X2[,3] + X2[,4] + X2[,5] + X2[,6]) # interactions
  # Extract coefficient estimates:
  beta0  = fit$coefficients[1] # the intercept
  beta1  = fit$coefficients[2:6] # main effects
  beta1a = fit$coefficients[7:8] # squared effects
  beta2 = fit$coefficients[9:14] # interactive effects
  
  # Create initials based on the above estimates:
  inits = list(list(beta0 = beta0, beta1 = beta1, beta1a = beta1a, beta2 = beta2, tau.Y = 50), #pink
               list(beta0 = beta0/10, beta1 = beta1/10, beta1a = beta1a/10, beta2 = beta2/10, tau.Y = 50), #green
               list(beta0 = beta0*10, beta1 = beta1*10, beta1a = beta1a*10, beta2 = beta2*10, tau.Y = 50)) #blue
  
  # Initial values: from saved state (if we already ran the model, we want to start where we left off)
  if(newinits==F){
    if(file.exists(initfilename)){
      load(initfilename)
      
      if(lowdev == T){
        if(length(saved_state)==3){ # temp
          saved_state <- lowdevrestart(saved_state, vary_by = 2)
        }
      }
      
      initslist <- saved_state[[2]]
      
      # temp
      # if(length(saved_state[[2]][[1]][["beta2"]])>6){
      #   inits2 <- saved_state[[2]]
      #   inits2[[1]][["beta2"]] <- inits2[[1]][["beta2"]][1:6]
      #   inits2[[2]][["beta2"]] <- inits2[[2]][["beta2"]][1:6]
      #   inits2[[3]][["beta2"]] <- inits2[[3]][["beta2"]][1:6]
      #   
      #   initslist <- inits2
      # }
      
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
             "beta0_p_temp", "beta1_p_temp", "beta1a_p_temp", "beta2_p_temp", # for p-values
             "dYdX", # net sensitivities
             "tau.Y", # variance of Y
             "wV","wT","wSs","wSd","wPAR", # importance weights
             "VPDant", "TAant", "Sshall_ant", "Sdeep_ant", "PAR_ant", # antecedent terms
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
  
  df_mod <- dateconnect(dfobj = df_mod, datevect = YIN$DOY, datename = "DOY", identifier = "ID2", varlist = c("dYdX"))
  df_mod <- dateconnect(dfobj = df_mod, datevect = YIN$TIMESTAMP, datename = "TIMESTAMP", identifier = "ID2", varlist = c("dYdX"))
  df_mod <- dateconnect(dfobj = df_mod, datevect = YIN$DOY, datename = "DOY", identifier = NULL, varlist = c("WUE.pred", "T.ratio", "T.pred", "ET.pred", "VPDant", "TAant", "Sshall_ant", "Sdeep_ant", "PAR_ant"))
  df_mod <- dateconnect(dfobj = df_mod, datevect = YIN$TIMESTAMP, datename = "TIMESTAMP", identifier = NULL, varlist = c("WUE.pred", "T.ratio", "T.pred", "ET.pred", "VPDant", "TAant", "Sshall_ant", "Sdeep_ant", "PAR_ant"))
  
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
  
  # inits to save
  init_names = names(initslist[[1]])
  
  # create a saved_state object with initials for next run
  # saved_state[[3]] is the chain number with lowest deviance
  saved_state <- keepvars(codaobj = jm_coda, to_keep = init_names, paramlist = params, type="jagsUI")
  
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
