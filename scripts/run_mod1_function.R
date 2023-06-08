#### Function to run JAGS model depending on site and response variable

run_mod1 <- function(dataIN, varname, sitename, newinits = F, overwrite = F, lowdev = F, post_only = F){
  
  # Temp for testing, comment out before using the function
  # dataIN = gebdat
  # varname = "GPP"
  # sitename = "geb"
  # overwrite = T
  # newinits = T
  # lowdev = F
  
  # Define filenames
  modelname <- "./models/model1/modv1.R"
  initfilename <- paste("./models/model1/", sitename, "/inits/inits_", varname, ".RData", sep = "")
  jm_codafilename <- paste("./models/model1/", sitename, "/coda/jm_coda_", varname,".RData", sep = "")
  mcmcfoldername <- paste("./models/model1/", sitename, "/convergence/", varname, sep = "")
  mcmcfilename <- paste("MCMC_", varname, sep = "")
  dffilename <- paste("./models/model1/", sitename,"/coda/df_mod1_", sitename, "_", varname, ".csv", sep = "")
  
  # Create necessary folders if they do not already exist
  if(!file.exists(paste("models/model1/", sitename, sep = ""))) { dir.create(paste("models/model1/", sitename, sep = ""))}
  if(!file.exists(paste("models/model1/", sitename, "/inits", sep = ""))) { dir.create(paste("models/model1/", sitename, "/inits", sep = ""))}
  if(!file.exists(paste("models/model1/", sitename, "/coda", sep = ""))) { dir.create(paste("models/model1/", sitename, "/coda", sep = ""))}
  if(!file.exists(paste("models/model1/", sitename, "/convergence", sep = ""))) { dir.create(paste("models/model1/", sitename, "/convergence", sep = ""))}
  if(!file.exists(mcmcfoldername)) { dir.create(mcmcfoldername)}
  
  
  #####################################################################
  #Part 1: Model setup
  
  # End for ETpart model (all days)
  N = nrow(dataIN)
  
  # Make the response variable file of interest. This file includes 
  # growing-season response (Y) variables of interest along with indices 
  # that link the Y variables back to covariates of interest.
  # NOTE: Kym did this to only test growing season response variables (4/1 to 10/31 for each year)
  YIN = dataIN %>%
    rowid_to_column("ind") %>% # ind: index to link the growing season Y variables back to the appropriate row in the covariate data set
    filter(!is.na(SIF_O2A_sfm)) %>% # filter from 7:30 am to 4:00 pm everyday, when SIF data was taken
    select(ind, GPP, SIF_O2A_sfm, SIF_O2B_sfm, all_of(varname)) # select variables we will be using as Y variables

  # Change Y to the column for the response variable of interest (Gs, GPP, or SIF) 
  Y  = YIN[,varname]
  Yday = YIN$ind
  
  # jIND file provides indices to calculate interactions between covariates
  # Basically a matrix version of X2, defined in the first initialization later in the script
  # key: 1 VPD # 2 Tair # 3 PAR # 4 Sshall # 5 Sdeep
  # X1a = cbind(X1[,1]^2, X1[,2]^2)
  # X2 = X1[,1]*X1[,2], X1[,1]*X1[,4], X1[,1]*X1[,5], X1[,2]*X1[,4], X1[,2]*X1[,5], X1[,4]*X1[,5]
  jIND <- data.frame(j = c(1:6),
                     ID1 = c(1,1,1,2,2,4),
                     ID2 = c(2,4,5,4,5,5))
  
  # an example from a previous SAM version:
  # X2  = cbind(X1[,1]*X1[,2], X1[,1]*X1[,3], X1[,1]*X1[,5], X1[,1]*X1[,6], X1[,2]*X1[,3], X1[,2]*X1[,5], X1[,2]*X1[,6], 
  # X1[,3]*X1[,5], X1[,3]*X1[,6], X1[,5]*X1[,6])
  # jIND <- data.frame(j = c(1:10),
  #                    ID1 = c(1,1,1,1,2,2,2,3,3,5),
  #                    ID2 = c(2,3,5,6,3,5,6,5,6,6))
  
  # Choose the starting index. This is an index for a row in the Y data file. 
  # The value in the indexed row should be greater than 1 to accommodate 
  # calculation of antecedent values.
  Nstart = YIN$ind[1]
  # Choose the ending index. 
  Nend   = length(Y)
  
  # Prepare data for JAGS -- covariates are scaled
  data = list(
              N = N, # Number of rows for input data
              Nstart = Nstart,
              Nend = Nend,
              Nlag = 6,
              Nparms = 5, # Nparms is the number of driving variables included to calculate main effects
              Y = Y,
              Yday = YIN$ind, # Choose column in YIN that provides indices linking response variables with covariates
              ID1 = jIND[,2], 
              ID2 = jIND[,3],
              jlength = nrow(jIND),
              VPD = as.vector(scale(dataIN$VPD,center=TRUE,scale=TRUE)), # scale function takes vector of values, centers and scales by SD
              Tair = as.vector(scale(dataIN$Tair,center=TRUE,scale=TRUE)),
              PAR = as.vector(scale(dataIN$PAR,center=TRUE,scale=TRUE)),
              Sshall = as.vector(scale(dataIN$SWC_shall,center=TRUE,scale=TRUE)),
              Sdeep = as.vector(scale(dataIN$SWC_deep,center=TRUE,scale=TRUE)),
              # covariate timesteps into the past
              # this code is flexible in case we want to combine timesteps
              C1 = c(0, 1, 2, 3, 4, 5), #stop times for covariates
              C2 = c(0, 1, 2, 3, 4, 5) #start times for covariates
              )

  
  # Get initials for WUE precision based on the log of previous model versions
  sig.Y = (sd(Y)**2)
  
  # Initial values are estimated using a linear model. As in the data list (above), 
  # covariates are centered and standardized. Replace name of covariate in quotes
  # with the appropriate column number in the dataIN file.
  X1  = cbind(data$VPD[Yday[Nstart:Nend]], #1
              data$Tair[Yday[Nstart:Nend]], #2
              data$PAR[Yday[Nstart:Nend]], #3
              data$Sshall[Yday[Nstart:Nend]], #4
              data$Sdeep[Yday[Nstart:Nend]] #5
              ) 
  
  # Notes: The code below is indexed numerically, which you will have to pay attention to as you change covariates of interest
  # Squared terms calculated for VPD and Tair
  X1a = cbind(X1[,1]^2, X1[,2]^2) 
  # Put all covariates together;
  # Interactions incorporated into linear model used to estimate initial values
  X2  = cbind(X1[,1]*X1[,2], X1[,1]*X1[,4], X1[,1]*X1[,5], X1[,2]*X1[,4], X1[,2]*X1[,5], 
              X1[,4]*X1[,5])
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
  if(file.exists(initfilename)){
    load(initfilename)
    initslist <- saved_state[[2]]
  }else if(!file.exists(initfilename)){
    initslist <- inits
  }
  
  #####################################################################
  # Part 2: Run JAGS Model (jagsui initializes and runs in one step)
  
  # parameters to track
  params = c("deviance", # deviance
             "beta0","beta1","beta1a","beta2", # intercept, main effects, squared effects, interactive effects
             "dYdX", # net sensitivites
             "sig.Y", # variance of Y
             "wV","wT", "wP","wSs","wSd", # importance weights
             "R2") # model fit
  
  if(post_only==F){ # if we want to run the model
    # Run model with jagsui package
    start<-proc.time() # keep track of run time
    jagsui <- jags(data = data,
                   inits = initslist,
                   model.file = modelname,
                   parameters.to.save = params,
                   n.chains = 3,
                   n.adapt = 1000,
                   n.thin = 10,
                   n.iter = 30000,
                   parallel = TRUE)
    end<-proc.time()
    elapsed<- (end-start)/60
    print("jagsui done running; minutes to completion:")
    print(elapsed[3])
    
    if(overwrite==T){
      save(jagsui, file = jm_codafilename) # save out whole jagsui object for local
    }
    
  } # end post_only==F
  
  if(post_only==T){ # if we just want to load the model
    load(jm_codafilename) # named jagsui
  }
  
  jm_coda <- jagsui$samples # convert to coda form to work with postjags functions
  
  df_mod1 <- data.frame(ID = get_index(jagsui), # identifiers
                                        #param = names(do.call(c, jagsui$mean)), # parameter names
                                        param = gsub('[[:digit:]]+', '', names(do.call(c, jagsui$mean))), # parameter names, removing all numbers
                                        mean = do.call(c, jagsui$mean), 
                                        pc2.5 = do.call(c, jagsui$q2.5), pc97.5 = do.call(c, jagsui$q97.5),
                                        overlap0 = do.call(c, jagsui$overlap0),
                                        gel = do.call(c, jagsui$Rhat))
  df_mod1$param[1] <- "beta0"
  df_mod1$param[2:(1+ncol(X1))] <- "beta1" # rename some betas with numbers to distinguish effect types
  df_mod1$param[(1+ncol(X1)+1):(1+ncol(X1)+ncol(X1a))] <- "beta1a"
  df_mod1$param[(1+ncol(X1)+ncol(X1a)+1):(1+ncol(X1)+ncol(X1a)+ncol(X2))] <- "beta2"
  
  
  write.csv(df_mod1, file = dffilename)
  
  #####################################################################
  # Part 3: Check diagnostics
  
  mcmcplot(jm_coda, parms = c("deviance", # deviance
                              "beta0","beta1","beta1a","beta2", # intercept, main effects, squared effects, interactive effects
                              "sig.Y", # variance of Y
                              "wV","wT", "wP","wSs","wSd", # importance weights
                              "R2"), # model fit
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
  }
  
  # inits to save
  init_names = names(initslist[[1]])
  
  # use get_remove_index function to find which variables to remove
  remove_vars = get_remove_index(init_names, params, type="jagsUI")
  
  newinits <- initfind(jm_coda, OpenBUGS = FALSE)
  newinits[[1]]
  saved_state <- removevars(initsin = newinits, 
                            variables = remove_vars)
  
  if(lowdev == T){
    # saved chain with lowest deviance, and make remaining chains vary around it
    saved_state[[2]][[1]] = saved_state[[2]][[devin]] # Best (low dev) initials for chain 1
    saved_state[[2]][[2]] = lapply(saved_state[[2]][[devin]],"*",2)
    saved_state[[2]][[3]] = lapply(saved_state[[2]][[devin]],"/",2)
  }
  
  #saved_state[[1]]
  
  if(overwrite==T){
    save(saved_state, file = initfilename) #for local
  }
  
  # If converged, run and save replicated data
  # jm_rep <- update(jagsui, parameters.to.save = "Y.rep",
  #                  n.iter = 15000, n.thin = 5)
  # 
  # if(overwrite==T){
  #   save(jm_rep, file = jm_repfilename) #for local
  # }
  # 
  
  
  
  
} # end of function
