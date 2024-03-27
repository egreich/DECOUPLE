### Model 2
### SAM model for WUE (either GPP/T or SIF/T)

model{
  for(i in Nstart:Nend){
    
    # Likelihood for ET data
    ET[i] ~ dnorm(ET.pred[i], tau.Y) # the 'true' ET should be somewhere around the predicted ET, with some precision
    ET.rep[i] ~ dnorm(ET.pred[i], tau.Y)
    ET.pred[i] <- E.model[i] + T.pred[i] # for evaluating model fit, not connected to ET data
    
    # intermediate calculation trick temp - this might just be needed for lae
    #ET.pred2[i] <- ifelse(ET.pred[i]>ET.max, ET.max, ET.pred[i])
    
    # T/ET ratio
    ET.int[i] <- ifelse(ET.pred[i]==0, 0.0000000000000001, ET.pred[i]) # intermediate calculation trick to ensure the denominator is not 0
    T.ratio[i] <- ifelse(ET.pred[i]==0, 0, T.pred[i]/ET.int[i])
    
    # Predicted transpiration
    # This is WUE = Photosynthesis/Transpiration rewritten to solve for T
    WUE.pred.int[i] <- ifelse(WUE.pred[i]==0, 0.0000000000000001, WUE.pred[i]) # intermediate calculation trick to ensure the denominator is not 0
    T.pred[i] <- (1/WUE.pred.int[i])*Phot[i] # Phot (photosynthesis) is the numerator, either GPP or SIF

    # Compute squared difference for calculating posterior predictive loss
    sqdiff[i] <- pow(ET.rep[i]-ET[i],2)
    
    #################### start soil evaporation process-based model ##########################
    
    # Calculate aerodynamic resistance to heat transfer (rah)- allow k to vary btwn 0.35-0.42
    # 0.001 is Z0m (m), the momentum soil roughness, [Yang et al., 2008; Stefan et al., 2015]
    # Notes: look for redundant code, and simplify to make it faster
    rah0[i] <- (1/((vk.pred**2)*ws[i]))*((log(Z/0.001))**2) # s/m
    #rah_unstable[i] <- rah0[i]/((1 + Ri[i])**0.75) # aerodynamic resistance to heat transfer s/m, estimated as in Choudhury et al. [1986]
    rah_stable[i] <- ((1 - 0.75*Ri[i])*((log(Z/0.001))**2))/((vk.pred**2)*ws[i])
    rah[i] <- ifelse( Ri[i] > 0, rah_unstable[i], rah_stable[i])
    
    # calculate parameters for alternative alpha and Bowen
    psi[i] <- psisat.pred * (S[i]/ssat.pred)**(-bch.pred)
    
    # calculate soil resistance (rss)- allow Arss and Brss to vary
    rss1[i] <- ((fc.pred - sres.pred)/(S[i] - sres.pred)) * rssmin # soil resistance for S > sres
    rss2[i] <- exp(8.206 - 4.255*(S[i]/fc.pred)) # another soil resistance derivation
    rss[i] <- ifelse(S[i] > sres.pred, rss1[i], rss2[i]) # the resistance to the diffusion of vapor in large soil pores
    
    # Calculate alpha and Bowen ratios
    alpha[i] <- exp((psi[i]*g)/(10000 * Rwv * Tsoil[i])) # Kelvin equation
    bowen[i] <- ifelse(S[i] > fc.pred, 1, (0.5 - 0.5*cos((S[i] * pi)/fc.pred))**2) # wetness function
    
    # Soil evaporation from E1 (CLM 4.5)
    LE4.5[i] <- ifelse(Tsoil[i] >= 0, bowen[i]*((rho[i]*Cp)/gamma[i])*((alpha[i]*(e.sat[i] - e.a[i]))/rah[i]), 0)
    Esoil4.5[i] <- conv.fact[i]*conv.fact.time*LE4.5[i]
    # Soil evaporation from E2 (CLM 3.5)
    LE3.5[i] <- ifelse(Tsoil[i] >= 0, ((rho[i]*Cp)/gamma[i])*((alpha[i]*(e.sat[i] - e.a[i]))/(rah[i] + rss[i])), 0)
    Esoil3.5[i] <- conv.fact[i]*conv.fact.time*LE3.5[i]
    
    # Intercepted E
    Eint[i] <- (P[i])*(1 - exp(-k.pred*(LAI[i])))
    
    # e.scalar ~ dunif(0,1) - maybe try this if the fit still isn't that great
    E.model[i] <- p*Esoil4.5[i] + (1-p)*Esoil3.5[i] + Eint[i]
    
    #################### end soil evaporation process-based model ##########################
    
    #################### start SAM model ###################################################
    
    # Regression (mean) model
    #WUE.pred[i] <- log.WUE[i]
    WUE.pred[i] <- exp(log.WUE[i]) # force mean to be positive, transform from log scale
    log.WUE[i] <- beta0 + main.effects[i] + squared.terms[i] + interactions[i] 
    
    # Define parts involving main effects, quadratic effects, and 2-way interactions
    main.effects[i]  <- sum(X.effect[,i])
    interactions[i]  <- sum(sum.XX.int[,i])
    squared.terms[i] <- sum(X2.effect[,i])
    
    # Define components of main.effects, squared.terms, and interactions:
    # Main effect parts:
    for(j in 1:Nparms){
      X.effect[j,i]<- beta1[j]*X[j,i]
    }
    
    # Squared terms:
    for(j in 1:2){
      X2.effect[j,i] <- beta1a[j]*pow(X[j,i],2)
    }
    
    # Two-way interaction terms:
    for(j in 1:jlength){
      XX.int[j,i] <- beta2[j]*X[ID1[j],i]*X[ID2[j],i]
      sum.XX.int[j,i] <- sum(XX.int[j,i])
    }
    
    # Creating antecedent covariates; 
    # This matrix of values will end up being used to calculate 
    # the parts involving main effects, interactions, and squared terms
    # in the regression model:
    X[1,i] <- VPDant[i]       ## Also included as squared term
    X[2,i] <- TAant[i]        ## Also included as squared term
    X[3,i] <- Sshall_ant[i]
    X[4,i] <- Sdeep_ant[i]
    X[5,i] <- PAR_ant[i]
    
    # Computed antecedent values. 
    VPDant[i]     <- sum(VPDtemp[i,]) # summing over all lagged j's
    TAant[i]      <- sum(Tairtemp[i,])
    Sshall_ant[i] <- sum(Sshalltemp[i,])
    Sdeep_ant[i] <- sum(Sdeeptemp[i,])
    PAR_ant[i] <- sum(PARtemp[i,])
    
    for(j in 1:Nlag){ # covariates into the past
      VPDtemp[i,j] <- wV[j]*V_temp[i,j]
      V_temp[i,j] <- mean(VPD[(Yday[i]-C1[j]):(Yday[i]-C2[j])]) # mean VPD during that block period
      
      Tairtemp[i,j] <- wT[j]*T_temp[i,j]
      T_temp[i,j] <- max(Tair[(Yday[i]-T1[j]):(Yday[i]-T2[j])])  # max temperature during that block period
      
      Sshalltemp[i,j] <- wSs[j]*Ss_temp[i,j]
      Ss_temp[i,j] <- mean(Sshall[(Yday[i]-T1[j]):(Yday[i]-T2[j])])
      
      Sdeeptemp[i,j] <- wT[j]*Sd_temp[i,j]
      Sd_temp[i,j] <- mean(Sdeep[(Yday[i]-T1[j]):(Yday[i]-T2[j])])
      
      PARtemp[i,j] <- wT[j]*PAR_temp[i,j]
      PAR_temp[i,j] <- mean(PAR[(Yday[i]-C1[j]):(Yday[i]-C2[j])])
    }
    
    # Calculate net sensitivities (derivative) -- derived quantities
    # key:
    # 1 VPD # 2 Tair # 3 Sshall # 4 Sdeep # 5 PAR
    # X1a = cbind(X1[,1]^2, X1[,2]^2)
    # X2  = cbind(X1[,1]*X1[,2], X1[,1]*X1[,4], X1[,2]*X1[,4], X1[,5]*X1[,1], X1[,5]*X1[,2], X1[,5]*X1[,4])
    dYdVPD[i] <- beta1[1] + 2*beta1a[1]*VPDant[i] + beta2[1]*TAant[i] + beta2[2]*Sdeep_ant[i] + beta2[4]*PAR_ant[i]
    dYdT[i]   <- beta1[2] + 2*beta1a[2]*TAant[i] + beta2[1]*VPDant[i] + beta2[3]*Sdeep_ant[i] + beta2[5]*PAR_ant[i]
    dYdSs[i]  <- beta1[3]
    dYdSd[i]  <- beta1[4] + beta2[2]*VPDant[i] + beta2[3]*TAant[i] + beta2[6]*PAR_ant[i]
    dYdPAR[i]  <- beta1[5] + beta2[4]*VPDant[i] + beta2[5]*TAant[i] + beta2[6]*Sdeep_ant[i]
    
    # Put all net sensitivities into one array, for easy monitoring
    dYdX[i,1] <- dYdVPD[i] * WUE.pred[i] #  * WUE.pred[i] converts the net sensitivities off the log scale
    dYdX[i,2] <- dYdT[i] * WUE.pred[i]
    dYdX[i,3] <- dYdSs[i] * WUE.pred[i]
    dYdX[i,4] <- dYdSd[i] * WUE.pred[i]
    dYdX[i,5] <- dYdPAR[i] * WUE.pred[i]
  }
  
  # Relatively non-informative priors for regression parameters:
  
  # Overall intercept:
  beta0 ~ dnorm(0,0.00001)
  beta0_p_temp[1] <- step(beta0) # Bayesian p-values
  
  # Main effects:
  for(j in 1:Nparms){
    beta1[j] ~ dnorm(0,0.00001)
    beta1_p_temp[j] <- step(beta1[j]) # Bayesian p-values
  }
  
  # Quadratic effects
  for(j in 1:2){
    beta1a[j] ~ dnorm(0,0.00001)
    beta1a_p_temp[j] <- step(beta1a[j]) # Bayesian p-values
  }
  
  # Two-way interaction effects:
  for(j in 1:jlength){
    beta2[j] ~ dnorm(0,0.00001)
    beta2_p_temp[j] <- step(beta2[j]) # Bayesian p-values
  }
  
  # Priors for importance weights for each covariate, "delta" or gamma "trick"
  # for imposing Dirichlet(1) priors for the weights:  
  for(j in 1:Nlag){
    # Priors for unnormalized weights
    dV[j]    ~ dgamma(1,1)
    dT[j]    ~ dgamma(1,1)
    dSs[j]   ~ dgamma(1,1)
    dSd[j]   ~ dgamma(1,1)
    dPAR[j]   ~ dgamma(1,1)
    
    # Compute normalized weights:
    wV[j]    <- dV[j]/sum(dV[])
    wT[j]    <- dT[j]/sum(dT[])
    wSs[j]   <- dSs[j]/sum(dSs[])
    wSd[j]   <- dSd[j]/sum(dSd[])
    wPAR[j]   <- dPAR[j]/sum(dPAR[])
  }
  
  #################### end SAM model ###################################################
  
  #################### start priors ##########################################################
  
  # Priors for ET:
  tau.Y ~ dgamma(.1,.1) # since this is associated with the data model for ET.
  sig.Y <- 1/sqrt(tau.Y)
  
  # Priors for stochastic parameters for evap equations
  k.pred ~ dnorm(0.5, 10)T(0,) # k = 0.5 # decay function k for intercepted E
  vk.pred ~ dunif(0.35, 0.42) # the von Karman constant is usually 0.40
  bch.pred ~ dnorm(bch, 0.40)T(0,) # Clapp and Hornberger parameter
  #k.pred ~ dnorm(0.5, 10)T(0,) # k = 0.5 # decay function k for intercepted E
  fc.pred ~ dnorm(fc, 200)T(S.min,1) # upper limit 1
  sres.pred ~ dnorm(sres, 80000)T(0,S.min)# residual soil moisture
  ssat.pred ~ dnorm(ssat, 50)T(0,) # soil moisture at saturation
  psisat.pred ~ dnorm(psisat, 0.015)T(-1000,0)  # parameterized air entry pressure, in mm of water
  
  # given p (proportion "contribution" pf Esoil4.5 to estimated Esoil) a prior, or set = 0.5 in data list
  p ~ dunif(0,1)
  
  #################### end priors ##########################################################
  
  #################### summary stats #######################################################
  
  # sum across days
  Dsum <- sum(sqdiff[Nstart:Nend])
  
  # Compute quantities for calculating Bayesian R2
  var.pred <- pow(sd(ET.pred[Nstart:Nend]),2) # var.pred <- pow(sd(mu[Nstart:Nend]),2)
  var.resid <- 1/tau.Y
  R2 <- var.pred/(var.pred + var.resid)

}