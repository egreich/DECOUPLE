### DEPART model
### Model for partitioning ET into T and E and for calculating WUE

model{
  for(i in Nstart:Nend){
    
    # Likelihood for ET data
    ET[i] ~ dnorm(ET.pred[i], tau.Y) # the 'true' ET should be somewhere around the predicted ET, with some precision
    ET.rep[i] ~ dnorm(ET.pred[i], tau.Y)
    ET.pred[i] <- E.model[i] + T.pred[i] # for evaluating model fit, not connected to ET data
    
    # T/ET ratio
    ET.int[i] <- ifelse(ET.pred[i]==0, 0.0000000000000001, ET.pred[i]) # intermediate calculation trick to ensure the denominator is not 0
    T.ratio[i] <- ifelse(ET.pred[i]==0, 0, T.pred[i]/ET.int[i])
    
    # Predicted transpiration
    # This is WUE = Photosynthesis/Transpiration rewritten to solve for T
    T.pred[i] <- (1/WUE.pred[i])*Phot[i] # Phot (photosynthesis) is the numerator, either GPP or SIF
    
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
    Esoil4.5[i] <- conv.fact[i]*LE4.5[i]
    # Soil evaporation from E2 (CLM 3.5)
    LE3.5[i] <- ifelse(Tsoil[i] >= 0, ((rho[i]*Cp)/gamma[i])*((alpha[i]*(e.sat[i] - e.a[i]))/(rah[i] + rss[i])), 0)
    Esoil3.5[i] <- conv.fact[i]*LE3.5[i]
    
    # Intercepted E
    Eint[i] <- (P[i])*(1 - exp(-k.pred*(LAI[i])))
    
    # e.scalar ~ dunif(0,1) - maybe try this if the fit still isn't that great
    E.model[i] <- p*Esoil4.5[i] + (1-p)*Esoil3.5[i] + Eint[i]
    
  }
    
    #################### end soil evaporation process-based model ##########################
    
    #################### start SAM model ###################################################
    
    # Regression (mean) model
    WUE.pred[i] <- log.WUE[i]
    #WUE.pred[i] <- exp(log.WUE) # force mean to be positive, transform from log scale
    log.WUE[i] <- beta0 + main.effects[i] + squared.terms[i] + interactions[i] 

  
  #################### end SAM model ###################################################
  
  #################### start priors ##########################################################
  
  # Priors for ET:
  tau.Y ~ dgamma(0.1,0.1) # since this is associated with the data model for ET.
  sig.Y <- 1/sqrt(tau.Y)
  
  # Priors for WUE:
  sig.WUE ~ dunif(0,100)
  tau.WUE <- pow(sig.WUE,-2)
  
  # Priors for stochastic parameters for evap equations
  k.pred ~ dnorm(0.5, 10)T(0,) # k = 0.5 # decay function k for intercepted E
  vk.pred ~ dunif(0.35, 0.42) # the von Karman constant is usually 0.40
  bch.pred ~ dnorm(bch, 0.40)T(0,) # Clapp and Hornberger parameter
  #k.pred ~ dnorm(0.5, 10)T(0,) # k = 0.5 # decay function k for intercepted E
  fc.pred ~ dnorm(fc, 200)T(S.min,1) # upper limit 1
  # could pick a precision that's less precise
  # Look at clay fraction across all sites, base precision off that, make more flexible:
  # decrease precision, increase variance (but check with sites)
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