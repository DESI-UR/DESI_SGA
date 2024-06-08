'''
Define a class MultLinFit that takes in a list of datasets and covariances, each 
in the same format used by Hyperfit's LinFit.  I.e.,
  - data are 2xN arrays, and datasets is a list of data
  - cov are 2x2xN arrays, and covs is a list of cov
Note that the data sets can be different sizes.

While the LinFit class allows for higher-dimensional fits -- planes and 
hyperplanes in addition to lines -- this class only fits lines, as needed for 
the Tully-Fisher relation.  It is assumed that all data sets share a common 
slope but have different intercepts and scatter parameters.  The fit parameters 
are of the form

\vec{\theta} = (a, b1, b2,...,b_m, sigma1, sigma2,...,sigma_m),

where a is the common slope, b1,...,b_m are the intercepts for the m datasets, 
and sigma1,...,sigma_m are the scatter parameters for each dataset.
'''


################################################################################
# Import modules
#-------------------------------------------------------------------------------
import numpy as np

from scipy.special import loggamma
from scipy.optimize import differential_evolution

import emcee
################################################################################




class MultiLinFit:
    """
    Class to implement linear fits to multiple datasets, assuming a common slope 
    but different intercepts across each set.
    
    Based on the hyperfit algorithm of Robotham and Obreschkow (PASP 2015) and the 
    Python LinFit implementation of Howlett and Gordon
    (https://hyperfit.readthedocs.io/en/latest/).
    
    Attributes
    ----------
    nsets : int
        Number of data sets and covariances entered by user.
    ndims : int
        Dimensionality of the data (expect 2).
    ndata : ndarray
        Array giving the length of every input data vector.
    params : ndarray
        Best-fit linear parameters for the data sets.
    params_scatter : ndarray
        Best-fit scatters along the y-axis for each data set.
    
    Parameters
    ----------
    datasets : list or ndarray
        An Mx2xN array of 2xN data vectors.
    covs : list or ndarray
        An Mx2x2xN array of 2x2xN covariance matrices.
    weights : ndarray
        Array of weights for each data set. Unit weights if not specified.
    vertaxis : float
        Specify which coordinate axis in data is the 'vertical' one. Defaults to last axis (-1).
    """
    
    def __init__(self, datasets, covs, weights=None, vertaxis=-1):
        
        self.nsets = len(datasets)
        self.ndims = np.shape(datasets[0])[0]
        self.ndata = np.array([np.shape(data)[1] for data in datasets])
        self.datasets = datasets
        self.covs = covs
        self.data = None
        self.cov = None
        
        if self.ndata[-1] < 3:
            # The last data set has too few points to have a quantifiable scatter
            self.npars = 1 + self.nsets - 1
            self.params_scatter = np.zeros(self.nsets - 1)
        else:
            self.npars = 1 + self.nsets # slope + intercepts + sigmas
            self.params_scatter = np.zeros(self.nsets)
        self.params = np.zeros(self.npars)
        
        self.weights = [np.ones(n) for n in self.ndata] if weights is None else weights
        self.vertaxis = vertaxis
        
        self.param_bounds = None      # parameter fit bounds for all data sets
        
    # Log posterior function.
    def _lnpost(self, params):
        lnpost = 0.

        for i in range(self.nsets):
            # Loop over individual data sets. 
            self.data = self.datasets[i]
            self.cov  = self.covs[i]
            
            # Set up parameter and bounds arrays for each data set.
            if self.ndata[-1] > 2:
                pars_i = np.array([params[0]] + [params[1+i]] + [params[self.nsets+1+i]])
                bounds_i = [self.param_bounds[0]] + \
                           [self.param_bounds[1+i]] + \
                           [self.param_bounds[self.nsets+1+i]]
            else:
                pars_i = np.array([params[0]] + [params[1+i]])
                bounds_i = [self.param_bounds[0]] + \
                           [self.param_bounds[1+i]]

            # Set up weights for each data set.
            weights = self.weights[i]
            
            # Sum over all data sets.
            lnprior = self._lnprior(pars_i, bounds_i)
            lnlike = self._lnlike(pars_i)                
            lnpost += np.sum(weights * lnlike) + lnprior
        
        return lnpost
            
    # Log prior function.
    def _lnprior(self, params, bounds):
        lnprior = 0.
        for i, (param, bound) in enumerate(zip(params.T, bounds)):
            lnprior += np.where(np.logical_or(param < bound[0], param > bound[1]), -np.inf, 0.0)

        return lnprior
    
    # Log likelihood function.
    def _lnlike(self, params):
        
        if len(params) > 2:
            a, b, sigma = params
        else:
            a, b = params
            sigma = 0

        x, dx2 = self.data[0], self.cov[0,0]
        y, dy2 = self.data[1], self.cov[1,1]
        dxy = self.cov[0,1]
        sy2 = sigma**2 + a**2*dx2 + dy2 - 2*dxy*a
        lnlike = 0.5*np.sum(np.log((a**2 + 1)/sy2) - (a*x - y + b)**2/sy2)

        return lnlike
    
    def bessel_cochran(self, sigma):
        """
        Bessel-Cochran correction of sample scatter to population scatter.
        
        Parameters
        ----------
        sigma : ndarray
            1xM array of scatters for the M input datasets.
        
        Returns
        -------
        sigma_corr : ndarray
            1xM array of corrected scatter parameters.
        """
        
        if self.ndata[-1] < 3:
            ndata = self.ndata[:-1]
        else:
            ndata = self.ndata
            
        sigma_corr = (np.sqrt(0.5 * ndata) * np.exp(loggamma(0.5 * (ndata - self.ndims)) - loggamma(0.5 * (ndata - self.ndims + 1.0)))) * sigma

        return sigma_corr
    
    def optimize(self, bounds, tol=1e-6, verbose=False):
        """
        Find the best-fit line for multiple datasets, assuming a common slope 
        across all sets but independent intercepts and scatters.
        
        Parameters
        ----------
        bounds : sequence
            Bounds for variables [a, b1, ..., bm, sig1, ..., sigm].
        tol : float
            Optimization tolerance.
        verbose : bool
            Print fit result.
            
        Returns
        -------
        params : ndarray
            Array of best-fit slope and intercepts [a, b1, b2, ..., bm]
        params_scatter : ndarray
            Array of vertical axis scatter parameters [sig1, sig2, ... sigm]
        log_posterior : float
            Value of ln(posterior) at the best fit point.
        """
        self.param_bounds = bounds
        res = differential_evolution(lambda *args: -self._lnpost(*args), self.param_bounds, tol=tol)

        if verbose:
            print(res)
            
        if self.ndata[-1] < 3:
            param_lim = self.nsets-1
        else:
            param_lim = self.nsets
            
        self.params = res.x[:-param_lim]
        self.params_scatter = np.fabs(res.x[-param_lim:])
        self.params_scatter = self.bessel_cochran(self.params_scatter)
        return self.params, self.params_scatter, -res.fun
    
    def emcee(self, 
              bounds, 
              max_iter=100000, 
              batchsize=1000, 
              ntau=50.0, 
              tautol=0.05, 
              verbose=False):
        """
        Run MCMC using the emcee EnsembleSampler.
        
        The MCMC is seeded using a randomization of the best-fit values of the 
        common slope, intercepts, and vertical scatters 
        [a, b1, ..., bm, sig1, ..., sigm].
        
        Parameters
        ----------
        bounds : sequence
            Bounds for variables [a, b1, ..., bm, sig1, ..., sigm].
        max_iter : int
            Maximum number of MCMC iterations.
        batchsize : int
            Size of each batch. Convergence checked after each batch.
        ntau : float
            Minimum autocorrelation length to consider for convergence.
        tautol : float
            Maximum fractional deviation between successive autocorrelation 
            lengths for convergence.
        verbose : bool
            Print out convergence statistics and progress bars if True.
            
        Returns
        -------
        mcmc_samples : ndarray
            Array of flattened and burned-in MCMC samples.
        mcmc_lnlike : ndarray
            Log-likelihood values of every MCMC sample.
        """

        # Set up emcee. Start the walkers in a small 1 percent ball around the 
        # best fit.  The best fit will set self.params and self.params_scatter.
        self.optimize(bounds, verbose=verbose)
        ndim = len(self.params) + len(self.params_scatter)
        nwalker = 4 * ndim
        seeds = np.asarray([[(0.01 * np.random.rand() + 0.995) * j for j in np.concatenate([self.params, self.params_scatter])] for _ in range(nwalker)])

        sampler = emcee.EnsembleSampler(nwalker, ndim, self._lnpost)

        old_tau = np.inf
        niter = 0
        converged = 0
        while ~converged:
            sampler.run_mcmc(seeds, nsteps=batchsize, progress=verbose)
            tau = sampler.get_autocorr_time(discard=int(0.5 * niter), tol=0)
            converged = np.all(ntau * tau < niter)
            converged &= np.all(np.abs(old_tau - tau) / tau < tautol)
            old_tau = tau
            begin = None
            niter += 1000
            if verbose:
                print("Niterations/Max Iterations: ", niter, "/", max_iter)
                print("Integrated ACT/Min Convergence Iterations: ", tau, "/", np.amax(ntau * tau))
            if niter >= max_iter:
                break

        # Remove burn-in and and save the samples
        tau = sampler.get_autocorr_time(discard=int(0.5 * niter), tol=0)
        burnin = int(2 * np.max(tau))
        samples = sampler.get_chain(discard=burnin, flat=True).T
        mcmc_samples = samples
        mcmc_lnlike = sampler.get_log_prob(discard=burnin, flat=True)

        return mcmc_samples, mcmc_lnlike