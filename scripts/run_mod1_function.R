#### Function to run JAGS model depending on site and response variable

run_mod1 <- function(dataIN, varname, sitename, scale, newinits = F, overwrite = F, lowdev = F){
  
  # Temp for testing, comment out before using the function
  if(test==T){
    dataIN = dat
    overwrite = F
    newinits = T
    lowdev = F
  }
  
  # Define filenames
  modelname <- "./models/model1.R"
  initfilename <- paste("./output_model1/inits/", sitename, "/inits_", varname, "_",scale, ".RData", sep = "")
  jm_codafilename <- paste("./output_model1/coda/", sitename, "/jm_coda_", varname, "_",scale,".RData", sep = "")
  mcmcfoldername <- paste("./output_model1/convergence/", sitename, "/", varname, "_",scale, sep = "")
  mcmcfilename <- paste("MCMC_", varname, "_",scale, sep = "")
  dffilename <- paste("./output_model1/df/df_", sitename, "_", varname, "_",scale, ".csv", sep = "")
  
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
  ntimestep = 7 # look at C1 and C2 to see the timesteps
  
  C1 = c(0, 1, 2, 3, 4, 5, 6) #stop times for covariates
  C2 = c(0, 1, 2, 3, 4, 5, 6) #start times for covariates
  T1 = c(0, 48, 96, 144, 192, 240, 288) #stop times for covariates over longer timescales
  T2 = c(0, 1, 49, 97, 145, 193, 241) #start times for covariates over longer timescales
  
  if(scale=="hour"){
    T1 = c(0, 24, 47, 71, 95, 119, 143) #stop times for covariates over longer timescales
    T2 = c(0, 1, 25, 48, 72, 96, 120) #start times for covariates over longer timescales
  }
  
  for(i in c(1:nrow(dataIN))){
   
    # ignore time periods where we don't have enough data and we don't want to gapfill
    for(j in C2){
      if(i <= j){
        next
      }
      if(is.na(dataIN$PAR[i-j])){
        dataIN[i,varname] = NA
      }
      if(is.na(dataIN$VPD[i-j])){
        dataIN[i,varname] = NA
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
        dataIN[i,varname] = NA
      }
      if(is.na(mean(dataIN$SWC_shall[(i-T1[j]):(i-T2[j])]))){
        dataIN[i,varname] = NA
      }
      if(is.na(mean(dataIN$SWC_deep[(i-T1[j]):(i-T2[j])]))){
        dataIN[i,varname] = NA
      }
    }
  }
  
  # check for NAs
  #which(is.na(dataIN$PAR))
  
  # If WUE is Inf, have the model ignore it, since 0 T doesn't really make since in this context
  dataIN[,varname] <- ifelse(is.infinite(pull(dataIN, varname)), NA, pull(dataIN, varname))
  
  # Make the response variable file of interest. This file includes 
  # growing-season response (Y) variables of interest along with indices 
  # that link the Y variables back to covariates of interest.
  if(scale == "halfhour"){
    YIN = dataIN %>%
      rowid_to_column("ind") %>% # ind: index to link the growing season Y variables back to the appropriate row in the covariate data set
      mutate(time = as.numeric(format(TIMESTAMP, "%H%M"))) %>%
      filter(between(time, 730, 1600)) # filter from 7:30 am to 4:00 pm everyday
  } else if(scale == "hour"){
    YIN = dataIN %>%
      rowid_to_column("ind") %>% # ind: index to link the growing season Y variables back to the appropriate row in the covariate data set
      mutate(time = as.numeric(format(TIMESTAMP, "%H%M"))) %>%
      filter(between(time, 700, 1600)) # filter from 7:30 am to 4:00 pm everyday
    }else if(scale == "daily"){
    YIN = dataIN %>%
      rowid_to_column("ind") # ind: index to link the growing season Y variables back to the appropriate row in the covariate data set
  }
  
  YIN = YIN %>%
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
  # X2  = cbind(X1[,1]*X1[,2], X1[,1]*X1[,4], X1[,2]*X1[,4], X1[,5]*X1[,1], X1[,5]*X1[,2], X1[,5]*X1[,4])
  jIND <- data.frame(j = c(1:6),
                     ID1 = c(1,1,2,5,5,5),
                     ID2 = c(2,4,4,1,2,4))
  
  # an example from a previous SAM version:
  # X2  = cbind(X1[,1]*X1[,2], X1[,1]*X1[,3], X1[,1]*X1[,5], X1[,1]*X1[,6], X1[,2]*X1[,3], X1[,2]*X1[,5], X1[,2]*X1[,6], 
  # X1[,3]*X1[,5], X1[,3]*X1[,6], X1[,5]*X1[,6])
  # jIND <- data.frame(j = c(1:10),
  #                    ID1 = c(1,1,1,1,2,2,2,3,3,5),
  #                    ID2 = c(2,3,5,6,3,5,6,5,6,6))
  
  # Choose the starting index. This is an index for a row in the Y data file. 
  # The value in the indexed row should be greater than 1 to accommodate 
  # calculation of antecedent values.
  # Nstart = 1
  # lagplusone = ntimestep+1
  # for( i in 1:lagplusone){ # if we don't have enough space at the start of YIN, move down
  #   nextrow <- i+1
  #   Nstart = ifelse(YIN$ind[i]<lagplusone, nextrow, Nstart)
  # }
  # if(Nstart>1){
  #   YIN <- YIN[Nstart:nrow(YIN),] # trim YIN to reflect new Nstart
  #   Nstart = 1 # make Nstart 1 again
  # }
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
  Y  = pull(YIN, varname)
  Yday = YIN$ind
  
  # Choose the ending index. 
  Nend   = nrow(YIN)
  
  # Prepare data for JAGS -- covariates are scaled
  data = list(
              Nstart = Nstart,
              Nend = Nend,
              Nlag = 7,
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
              #SWC = as.vector(scale(as.numeric(dataIN$SWC),center=TRUE,scale=TRUE)),
              PAR = as.vector(scale(as.numeric(dataIN$PAR),center=TRUE,scale=TRUE)),
              # covariate timesteps into the past
              # this code is flexible in case we want to combine timesteps,
              # or average over multiple timesteps
              C1 = C1, #stop times for covariates
              C2 = C2, #start times for covariates
              T1 = T1, #stop times for covariates over longer timescales
              T2 = T2 #start times for covariates over longer timescales
              )

  
  # Get initials for WUE precision based on the log of previous model versions
  sig.Y = (sd(Y)**2)
  #temp
  if(varname=="WUE_GPP" & sitename=="crk"){ # The crk site has some weird values, so basing initials off of the standard dev is not good enough
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
  X2  = cbind(X1[,1]*X1[,2], X1[,1]*X1[,4], X1[,2]*X1[,4], 
              X1[,5]*X1[,1], X1[,5]*X1[,2], X1[,5]*X1[,4])
  # Fit simple linear model
  fit <- lm(Y[Nstart:Nend] ~ X1[,1] + X1[,2] + X1[,3] + X1[,4] + X1[,5] + # main effects
              X1a[,1] + X1a[,2] + # squared
              X2[,1] + X2[,2] + X2[,3] + X2[,4] + X2[,5] + X2[,6]) # interactions
  # Extract coefficient estimates:
  beta0  = fit$coefficients[1] # the intercept
  beta1  = fit$coefficients[2:6] # main effects
  beta1a = fit$coefficients[7:8] # squared effects
  beta2 = fit$coefficients[9:14] # interactive effects
  
  # Create initials based on the above estimates:
  inits = list(list(beta0 = beta0, beta1 = beta1, beta1a = beta1a, beta2 = beta2, sig.Y = sig.Y), #pink
               list(beta0 = beta0/10, beta1 = beta1/10, beta1a = beta1a/10, beta2 = beta2/10, sig.Y =  sig.Y/5), #green
               list(beta0 = beta0*10, beta1 = beta1*10, beta1a = beta1a*10, beta2 = beta2*10, sig.Y =  sig.Y*5)) #blue
  
  # Initial values: from saved state (if we already ran the model, we want to start where we left off)
  if(newinits==F){
    
    if(file.exists(initfilename)){
      load(initfilename)
      
      if(lowdev == T){
        saved_state <- lowdevrestart(saved_state, vary_by = 10)
      }
      
      initslist <- saved_state[[2]]
      
      # temp
      if(length(saved_state[[2]][[1]][["beta2"]])>6){
        inits2 <- saved_state[[2]]
        inits2[[1]][["beta2"]] <- inits2[[1]][["beta2"]][1:6]
        inits2[[2]][["beta2"]] <- inits2[[2]][["beta2"]][1:6]
        inits2[[3]][["beta2"]] <- inits2[[3]][["beta2"]][1:6]
        
        initslist <- inits2
      }
      
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
             "beta0_p_temp", "beta1_p_temp", "beta1a_p_temp", "beta2_p_temp", # for p-values
             "dYdX", # net sensitivities
             "sig.Y", # variance of Y
             "wV","wT","wSs","wSd","wPAR", # importance weights
             "R2", "Dsum") # model fit
    
# Run model with jagsui package
start<-proc.time() # keep track of run time
jagsui <- jags(data = data,
                   inits = initslist,
                   model.file = modelname,
                   parameters.to.save = params,
                   n.chains = 3,
                   n.adapt = 500,
                   n.thin = ifelse(test==F,3,1),
                   n.iter = ifelse(test==F,50000,200),
                   parallel = ifelse(test==F,TRUE, FALSE))
end<-proc.time()
elapsed<- (end-start)/60
print("jagsui done running; minutes to completion:")
print(elapsed[3])
    
if(overwrite==T){
  save(jagsui, file = jm_codafilename) # save out whole jagsui object for local
}

if(test==T){
  load(jm_codafilename)
}
  
  jm_coda <- jagsui$samples # convert to coda form to work with postjags functions
  
  df_mod <- dumsum(jagsui, type = "jagsUI") # organize the coda object as a dataframe
            
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
