import numpy as np
import xarray as xr
from matplotlib import pyplot as pl
from scipy.optimize import fmin

class QuantReg:
    
    def __init__(self, PolyDeg=1, rho=0.95):
        '''QuantReg(PolyDeg=1,rho=0.95)

        Quantile regression

        Fits a polynomial function (of degree PolyDeg) using quantile regression based on a percentile (rho).
        Based on script by Dr. Phillip M. Feldman, and based on method by Koenker, Roger, and
        Gilbert Bassett Jr. “Regression Quantiles.” Econometrica: Journal of
        the Econometric Society, 1978, 33–50.


        Parameters
        ----------
        PolyDeg : int
            The degree of the polynomial function
        rho : float between 0-1
            The percentile to fit to, must be between 0-1
        '''        

        self.PolyDeg = PolyDeg
        self.rho     = rho
        self.N_coefficients = PolyDeg+1
        

    def model(self, x, beta):
       """
       This example defines the model as a polynomial, where the coefficients of the
       polynomial are passed via `beta`.
       """
       if self.PolyDeg == 0:
           return x*beta
       else:
           return np.polyval(beta, x)

    
    def tilted_abs(self, x, weights):
       """
       OVERVIEW

       The tilted absolute value function is used in quantile regression.


       INPUTS

       x: This parameter represents a value of the independent variable, and in
       general takes any real value (float) or NumPy array of floats.
       """

       return weights * x * (self.rho - (x < 0))

    def objective(self, beta, x, y, weights):
       """
       The objective function to be minimized is the sum of the tilted absolute
       values of the differences between the observations and the model.
       """
       return self.tilted_abs(y - self.model(x, beta), weights).sum()

    def fit(self, x, y, weights=None):
        '''fit(x, y, weights=None)

        fit

        Fits a polynomial function (of degree PolyDeg) using quantile regression based on a percentile (rho).
        Based on script by Dr. Phillip M. Feldman, and based on method by Koenker, Roger, and
        Gilbert Bassett Jr. “Regression Quantiles.” Econometrica: Journal of
        the Econometric Society, 1978, 33–50.


        Parameters
        ----------
        x : list or list like
            independent variable
        y : list or list like
            dependent variable
        weights : list or list like
            Vector to weight each point, must be same size as x

         Returns
        -------
        list
            The resulting parameters in order of degree from low to high
        '''
        
        # Build weights if they don't exits:
        if weights is None:
            weights=np.ones(x.shape)

        # Define starting point for optimization:
        beta_0 = np.zeros(self.N_coefficients)
        if self.N_coefficients >= 2:
           beta_0[1]= 1.0

        # `beta_hat[i]` will store the parameter estimates for the quantile
        # corresponding to `fractions[i]`:
        beta_hat= []

        #for i, fraction in enumerate(fractions):
        beta_hat = fmin(self.objective, x0=beta_0, args=(x, y, weights), xtol=1e-8,
          disp=False, maxiter=3000) 
        return beta_hat
        
        
## Load Data        
#site = "DE-Hai"
#ref       = xval.sel(site=site).training_mask.mean(dim='hour', skipna=False)
#sw_in_pot = precomp_loc.sel(site=site).SW_IN_POT.mean(dim='hour', skipna=False).where(ref>=0.7)
#sw_in     = precomp_loc.sel(site=site).SW_IN.mean(dim='hour', skipna=False).where(ref>=0.7)
import pandas as pd 

data = pd.read_csv("./data_formatted/CH-Lae/CH-Lae_dat.csv")
sw_in_pot = data["SW_IN"]
sw_in = data["SW_IN"]

#### Version with an intercept:

mask = sw_in_pot.notnull() & sw_in.notnull()

## PolyDeg is the degree of polynomial (0 => y=x*m, 1 => y=x*m+b)
## rho is the quantile to fit (0.95 => 95th percentile)
qreg = QuantReg(PolyDeg=1, rho=0.95)
params = qreg.fit(sw_in_pot[mask], sw_in[mask])

pl.figure(figsize=(8,8))
pl.title("With intercept")
xlim = (0, sw_in_pot.max().values)
pl.scatter(sw_in_pot[mask], sw_in[mask], marker='.')
xlin = np.linspace(*xlim)
pl.plot(xlin, qreg.model(xlin, params), 'r', lw=4)
pl.xlim(*xlim)
pl.ylim(*xlim)


#### Version with no intercept (slope only)

qreg = QuantReg(PolyDeg=0, rho=0.95)
params = qreg.fit(sw_in_pot[mask], sw_in[mask])

pl.figure(figsize=(8,8))
pl.title("No intercept")
xlim = (0, sw_in_pot.max().values)
pl.scatter(sw_in_pot[mask], sw_in[mask], marker='.')
xlin = np.linspace(*xlim)
pl.plot(xlin, qreg.model(xlin, params), 'r', lw=4)
pl.xlim(*xlim)
pl.ylim(*xlim)


