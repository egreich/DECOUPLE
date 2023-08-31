#### Function to run JAGS model depending on site and response variable

run_mod1 <- function(dataIN, varname, sitename, newinits = F, overwrite = F, lowdev = F, post_only = F){
  
  # Temp for testing, comment out before using the function
  if(test==T){
    dataIN = dat
    overwrite = F
    newinits = T
    lowdev = F
  }
  
  # Define filenames
  modelname <- "./models/model1.R"
  initfilename <- paste("./output_model1/inits/", sitename, "/inits_", varname, ".RData", sep = "")
  initlowfilename <- paste("./output_model1/inits/", sitename, "/initlow_", varname, ".RData", sep = "")
  jm_codafilename <- paste("./output_model1/coda/", sitename, "/jm_coda_", varname,".RData", sep = "")
  mcmcfoldername <- paste("./output_model1/convergence/", sitename, "/", varname, sep = "")
  mcmcfilename <- paste("MCMC_", varname, sep = "")
  dffilename <- paste("./output_model1/df/df_", sitename, "_", varname, ".csv", sep = "")
  
  # Create necessary directories if they do not already exist
  if(!dir.exists(paste("output_model1/inits", sep = ""))){ dir.create(paste("output_model1/inits", sep = ""))}
  if(!dir.exists(paste("output_model1/coda", sep = ""))){ dir.create(paste("output_model1/coda", sep = ""))}
  if(!dir.exists(paste("output_model1/convergence", sep = ""))){ dir.create(paste("output_model1/convergence", sep = ""))}
  if(!dir.exists(paste("output_model1/inits/", sitename, sep = ""))){ dir.create(paste("output_model1/inits/", sitename, sep = ""))}
  if(!dir.exists(paste("output_model1/coda/", sitename, sep = ""))){ dir.create(paste("output_model1/coda/", sitename, sep = ""))}
  if(!dir.exists(paste("output_model1/convergence/", sitename, sep = ""))){ dir.create(paste("output_model1/convergence/", sitename, sep = ""))}
  if(!dir.exists(mcmcfoldername)){ dir.create(mcmcfoldername)}
  if(!dir.exists(paste("output_model1/df", sep = ""))){ dir.create(paste("output_model1/df", sep = ""))}

  #####################################################################
  #Part 1: Model setup
  
  # If we are missing env variables we don't want to gapfill Nlag timesteps into the past,
  # make our response varibale NA to skip over those periods
  ntimestep = 8 # look at C1 and C2 to see the timesteps
  for(i in c(1:nrow(dataIN))){
    if(i < (ntimestep+2)){
      next
    }
    for(j in c(0:(ntimestep+1))){
      if(is.na(dataIN$SWC_shall[i-j])){
        dataIN[i,varname] = NA
      }
      if(is.na(dataIN$SWC_deep[i-j])){
        dataIN[i,varname] = NA
      }
      if(is.na(dataIN$PAR[i-j])){
        dataIN[i,varname] = NA
      }
    }
  }
  
  # check for NAs
  which(is.na(dataIN$PAR))
  
  # If WUE is Inf, have the model ignore it, since 0 T doesn't really make since in this context
  dataIN[,varname] <- ifelse(is.infinite(dataIN[,varname]), NA, dataIN[,varname])
  
  # Make the response variable file of interest. This file includes 
  # growing-season response (Y) variables of interest along with indices 
  # that link the Y variables back to covariates of interest.
  # NOTE: Kym did this to only test growing season response variables (4/1 to 10/31 for each year)
  YIN = dataIN %>%
    rowid_to_column("ind") %>% # ind: index to link the growing season Y variables back to the appropriate row in the covariate data set
    mutate(time = as.numeric(format(TIMESTAMP, "%H%M"))) %>%
    filter(between(time, 730, 1600)) %>% # filter from 7:30 am to 4:00 pm everyday
    filter(!is.na(TA)) %>% # filter any weird start/end gaps for envir variables, this won't cause gaps in the middle of the day, since we already gap-filled
    filter(!is.na(VPD)) %>%
    filter(!is.na(SWC_shall)) %>%
    filter(!is.na(SWC_deep)) %>%
    filter(!is.na(PAR)) %>%
    drop_na(all_of(varname)) %>% # Filter out NAs in response variables
    select(ind, all_of(varname), DOY, TIMESTAMP) # select variables we will be using as Y variables and also index and time columns
  
  # jIND file provides indices to calculate interactions between covariates
  # Basically a matrix version of X2, defined in the first initialization later in the script
  # key: 1 VPD # 2 Tair # 3 Sshall # 4 Sdeep # 5 PAR
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
  lagplusone = ntimestep+1
  for( i in 1:lagplusone){ # if we don't have enough space at the start of YIN, move down
    nextrow <- i+1
    Nstart = ifelse(YIN$ind[i]<lagplusone, nextrow, Nstart)
  }
  if(Nstart>1){
    YIN <- YIN[Nstart:nrow(YIN),] # trim YIN to reflect new Nstart
    Nstart = 1 # make Nstart 1 again
  }
  
  # Change Y to the column for the response variable of interest
  Y  = YIN[,varname]
  Yday = YIN$ind
  
  # Choose the ending index. 
  Nend   = nrow(YIN)
  
  # Prepare data for JAGS -- covariates are scaled
  data = list(
              Nstart = Nstart,
              Nend = Nend,
              Nlag = 9,
              Nparms = 5, # Nparms is the number of driving variables included to calculate main effects
              Y = Y,
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
              C1 = c(0, 1, 2, 3, 4, 5, 6, 7, 8), #stop times for covariates
              C2 = c(0, 1, 2, 3, 4, 5, 6, 7, 8) #start times for covariates
              )

  
  # Get initials for WUE precision based on the log of previous model versions
  sig.Y = (sd(Y)**2)
  sig.Y = if(varname=="WUE_GPP" & sitename=="lnf"){ # The Geb site has some weird values, so basing initials off of the standard dev is not good enough
    sig.Y = 8 # This is based on the initial value we got for GPP for the Geb site
  }
  
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
  fit <- lm(Y[Nstart:Nend] ~ X1[,1] + X1[,2] + X1[,3] + X1[,4] + X1[,5] + # main effects
              X1a[,1] + X1a[,2] + # squared
              X2[,1] + X2[,2] + X2[,3] + X2[,4] + X2[,5] + X2[,6] + X2[,7] + X2[,8] + X2[,9] + X2[,10]) # interactions
  # Extract coefficient estimates:
  beta0  = fit$coefficients[1] # the intercept
  beta1  = fit$coefficients[2:6] # main effects
  beta1a = fit$coefficients[7:8] # squared effects
  beta2 = fit$coefficients[9:18] # interactive effects
  
  # Create initials based on the above estimates:
  inits = list(list(beta0 = beta0, beta1 = beta1, beta1a = beta1a, beta2 = beta2, sig.Y = sig.Y), #pink
               list(beta0 = beta0/10, beta1 = beta1/10, beta1a = beta1a/10, beta2 = beta2/10, sig.Y =  sig.Y/5), #green
               list(beta0 = beta0*10, beta1 = beta1*10, beta1a = beta1a*10, beta2 = beta2*10, sig.Y =  sig.Y*5)) #blue
  
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
             "beta0","beta1","beta1a","beta2", # intercept, main effects, squared effects, interactive effects
             "dYdX", # net sensitivites
             "sig.Y", # variance of Y
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
                   n.thin = ifelse(test==F,3,1),
                   n.iter = ifelse(test==F,20000,200),
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
            
  # Connect dates to dataframe by variable
  df_mod <- dateconnect(dfobj = df_mod, datevect = YIN$DOY, datename = "DOY", identifier = "ID2", varlist = c("dYdX"))
  df_mod <- dateconnect(dfobj = df_mod, datevect = YIN$TIMESTAMP, datename = "TIMESTAMP", identifier = "ID2", varlist = c("dYdX"))
  
  write.csv(df_mod, file = dffilename)
  
  
  #####################################################################
  # Part 3: Check diagnostics
  
  mcmcplot(jm_coda, parms = c("deviance", # deviance
                              "beta0","beta1","beta1a","beta2", # intercept, main effects, squared effects, interactive effects
                              "sig.Y", # variance of Y
                              "wV","wT","wSs","wSd","wPAR", # importance weights
                              "R2"), # model fit
           dir = mcmcfoldername,
           filename = mcmcfilename)
  
  #####################################################################
  # Part 5: Save inits for future runs
  
  # Save state
  
  if(lowdev == T){
    # Save inits based on chains with lowest deviance
    # dev_col <- which(colnames(jm_coda[[1]]) == "deviance")
    # dev1<- mean(jm_coda[[1]][,dev_col])
    # dev2<- mean(jm_coda[[2]][,dev_col])
    # dev3<- mean(jm_coda[[3]][,dev_col])
    # dev_min <- min(dev1, dev2, dev3)
    # if(dev1 == dev_min){
    #   devin = 1
    # } else if(dev2 == dev_min){
    #   devin = 2
    # } else if(dev3 == dev_min){
    #   devin = 3
    # }
    # initlow <- devin
    
    initlow <- findlowdev(jm_coda)
    #print(paste("chain with lowest deviance: ", initlow, sep=""))
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
