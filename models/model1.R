


model{
  for(i in Nstart:Nend){
    
    # Likelihood for observed T, GPP, or SIF data
    Y[i] ~ dnorm(mu[i], tau.Y) # the 'true' Y should be somewhere around the mean Y, with some precision
    # Replicated data
    Y.rep[i] ~ dnorm(mu[i], tau.Y)

    # Compute log predictive density for WAIC
    #ldx[i] <- logdensity.norm(Y[i], Y.rep[i], tau.Y)
    # Compute predictive density for WAIC:
    #dx[i] <- exp(ldx[i])

    # Compute squared difference for calculating posterior predictive loss
    sqdiff[i] <- pow(Y.rep[i]-Y[i],2)
    
    # Regression (mean) model
    mu[i] <- beta0 + main.effects[i] + squared.terms[i] + interactions[i] 
    
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
    
    for(j in 1:Nlag){ # covariates weeks, months into the past
      VPDtemp[i,j] <- wV[j]*V_temp[i,j]
      V_temp[i,j] <- mean(VPD[(Yday[i]-C1[j]):(Yday[i]-C2[j])]) # mean VPD during that block period
      
      Tairtemp[i,j] <- wT[j]*T_temp[i,j]
      T_temp[i,j] <- mean(Tair[(Yday[i]-C1[j]):(Yday[i]-C2[j])])
      
      Sshalltemp[i,j] <- wSs[j]*Ss_temp[i,j]
      Ss_temp[i,j] <- mean(Sshall[(Yday[i]-C1[j]):(Yday[i]-C2[j])])
      
      Sdeeptemp[i,j] <- wT[j]*Sd_temp[i,j]
      Sd_temp[i,j] <- mean(Sdeep[(Yday[i]-C1[j]):(Yday[i]-C2[j])])
      
      PARtemp[i,j] <- wT[j]*PAR_temp[i,j]
      PAR_temp[i,j] <- mean(PAR[(Yday[i]-C1[j]):(Yday[i]-C2[j])])
    }
    
    # Calculate net sensitivities (derivative) -- derived quantities
    # key:
    # 1 VPD # 2 Tair # 3 Sshall # 4 Sdeep # 5 PAR
    # X1a = cbind(X1[,1]^2, X1[,2]^2)
    # X2 = X1[,1]*X1[,2], X1[,1]*X1[,3], X1[,1]*X1[,4], X1[,2]*X1[,3], X1[,2]*X1[,4], X1[,3]*X1[,4], X1[,5]*X1[,1], X1[,5]*X1[,2], X1[,5]*X1[,3], X1[,5]*X1[,4]
    dYdVPD[i] <- beta1[1] + 2*beta1a[1]*VPDant[i] + beta2[1]*TAant[i] + beta2[2]*Sshall_ant[i] + beta2[3]*Sdeep_ant[i] + beta2[7]*PAR_ant[i]
    dYdT[i]   <- beta1[2] + 2*beta1a[2]*TAant[i] + beta2[1]*VPDant[i] + beta2[4]*Sshall_ant[i] + beta2[5]*Sdeep_ant[i] + beta2[8]*PAR_ant[i]
    dYdSs[i]  <- beta1[3] + beta2[2]*VPDant[i] + beta2[4]*TAant[i] + beta2[6]*Sdeep_ant[i] + beta2[9]*PAR_ant[i]
    dYdSd[i]  <- beta1[4] + beta2[3]*VPDant[i] + beta2[5]*TAant[i] + beta2[6]*Sshall_ant[i] + beta2[10]*PAR_ant[i]
    dYdPAR[i]  <- beta1[5] + beta2[7]*VPDant[i] + beta2[8]*TAant[i] + beta2[9]*Sshall_ant[i] + beta2[10]*Sshall_ant[i]
    
    # Put all net sensitivities into one array, for easy monitoring
    dYdX[i,1] <- dYdVPD[i]
    dYdX[i,2] <- dYdT[i]
    dYdX[i,3] <- dYdSs[i]
    dYdX[i,4] <- dYdSd[i]
    dYdX[i,5] <- dYdPAR[i]
  }
  
  # Relatively non-informative priors for regression parameters:
  
  # Overall intercept:
  beta0 ~ dnorm(0,0.00001)
  
  # Main effects:
  for(j in 1:Nparms){
    beta1[j] ~ dnorm(0,0.00001)
  }
  
  # Quadratic effects
  for(j in 1:2){
    beta1a[j] ~ dnorm(0,0.00001)
  }
  
  # Two-way interaction effects:
  for(j in 1:jlength){
    beta2[j] ~ dnorm(0,0.00001)
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
  

  
  #Prior for standard deviation in data likelihood 
  tau.Y <- pow(sig.Y,-2)
  sig.Y ~ dunif(0,1000)
  
  # sum across days
  Dsum <- sum(sqdiff[Nstart:Nend])
  
  # Compute quantities for calculating Bayesian R2
  var.pred <- pow(sd(mu[Nstart:Nend]),2)
  var.resid <- 1/tau.Y
  R2 <- var.pred/(var.pred + var.resid)
  
  # Priors for resonse variable:
  #tau.Y ~ dgamma(0.1,0.1) # since this is associated with the data model for Y
  #sig.Y <- 1/sqrt(tau.Y)

}