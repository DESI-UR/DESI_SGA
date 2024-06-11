'''
This is a collection of the various methods of fitting a line to data.
'''


################################################################################
# Import modules
#-------------------------------------------------------------------------------
import numpy as np

from hyperfit.linfit import LinFit
from hyperfit_v2 import MultiLinFit

from scipy.optimize import minimize
################################################################################




################################################################################
# Invert the slope and y-intercept
#-------------------------------------------------------------------------------
def param_invert(w0, w1, cov):
    '''
    Convert the slope and y-intercept from the fit
                                x = w0*y + w1
    to the fit
                                y = a*x + b


    PARAMETERS
    ==========

    w0 : float
        slope of the fit x = w0*y + w1
        
    w1 : ndarray of shape (m,)
        y-intercepts of the fit x = w0*y + w1

    cov : ndarray of shape (1+m,1+m)
        Covariance matrix of w0 and w1


    RETURNS
    =======

    a : float
        slope of the fit y = a*x + b
        
    b : ndarray of shape (m,)
        y-intercepts of the fit y = a*x + b

    cov_ab : ndarray of shape (1+m,1+m)
        Covariance matrix of a and b
    '''


    ############################################################################
    # Convert w0, w1 to a, b
    #---------------------------------------------------------------------------
    a = 1./w0
    b = -w1/w0
    ############################################################################


    ############################################################################
    # Convert the covariance matrix from (w0, w1) to (a, b)
    #---------------------------------------------------------------------------
    # Generate a bunch of realizations of w0 and w1
    rng = np.random.default_rng()
    if isinstance(w1, np.ndarray):
        w_rng = rng.multivariate_normal([w0, *w1], cov, 5000).T
        w1_rng = w_rng[1:]
    else:
        w_rng = rng.multivariate_normal([w0, w1], cov, 5000).T
        w1_rng = w_rng[1]
    w0_rng = w_rng[0]

    # Transform them from (w0, w1) space to (a, b) space
    a_rng = 1./w0_rng
    b_rng = -w1_rng/w0_rng

    # Calculate the covariance matrix for (a, b)
    if isinstance(w1, np.ndarray):
        ab_rng = np.stack((a_rng, *b_rng))
    else:
        ab_rng = np.stack((a_rng, b_rng))
    cov_ab = np.cov(ab_rng)
    ############################################################################

    return a, b, cov_ab
################################################################################
################################################################################




################################################################################
# Fitting using hyperfit
#-------------------------------------------------------------------------------
def hyperfit_line(x, y, dx, dy, bounds):
    '''
    Fit a line, y = ax + b, to the data, using the hyperfit package.  This 
    accepts uncertainties in both x and y.


    PARAMETERS
    ==========
    
    x, y : ndarrays of shape (N,)
        (x, y) values to fit the line to.

    dx, dy : ndarrays of shape (N,)
        Uncertainties in x and y.

    bounds : length-3 tuple of (min, max)
        The lower and upper bounds for each of the three fit parameters:
          1. (slope minimum, slope maximum)
          2. (y-intercept minimum, y-intercept maximum)
          3. (scatter minimum, scatter maximum)


    RETURNS
    =======

    a, b, sig : floats
        slope, y-intercept, and scatter best-fit values

    cov : ndarray of shape (3,3)
        Covariance matrix of the fit

    mcmc_samples : ndarray of shape (3,M)
        Flattened values from the MCMC chain for each of the three parameters

    hf : hyperfit object
    '''

    N = len(x) # Number of data points

    ############################################################################
    # Create a covariance matrix
    # 
    # This assumes that x and y are independent.
    #---------------------------------------------------------------------------
    cov = np.empty((2, 2, N))

    for i in range(N):
        cov[:,:,i] = np.array([[dx[i]**2, 0.], [0., dy[i]**2]])
    ############################################################################


    ############################################################################
    # Create the hacked hyperfit object
    #---------------------------------------------------------------------------
    hf = MultiLinFit([x, y], cov)
    ############################################################################


    ############################################################################
    # Run MCMC to fit line
    #---------------------------------------------------------------------------
    mcmc_samples, mcmc_lnlike = hf.emcee(bounds, verbose=True)
    ############################################################################


    ############################################################################
    # Calculate parameter values and covariance matrix
    #---------------------------------------------------------------------------
    a, b, sig = np.median(mcmc_samples, axis=1)

    cov = np.cov(mcmc_samples)
    ############################################################################

    return a, b, sig, cov, mcmc_samples, hf
################################################################################
################################################################################




################################################################################
# Fitting using the hacked version of hyperfit (for multiple data sets)
#-------------------------------------------------------------------------------
def hyperfit_line_multi(x, y, dx, dy, bounds):
    '''
    Fit a line, y = ax + b, to the data, using a hacked version of the hyperfit 
    package.  This accepts uncertainties in both x and y.


    PARAMETERS
    ==========
    
    x, y : lists with length M of ndarrays of shape (N,)
        (x, y) values to fit the line to.
        
        Note that each (xi, yi) pair must have the same length, but all xi do 
        not need to have the same length.

    dx, dy : lists with length M of ndarrays of shape (N,)
        Uncertainties in x and y.
        
        Note that each (dxi, dyi) pair must have the same length (and match the 
        length of the corresponding (xi, yi) pair).

    bounds : length-(2M + 1) tuple of (min, max)
        The lower and upper bounds for each fit parameter:
          1. (slope minimum, slope maximum)
          2. (y-intercept minimum, y-intercept maximum) x M
          3. (scatter minimum, scatter maximum) x M


    RETURNS
    =======

    a : float
        slope best-fit value
        
    b, sig : ndarrays of shape (M,)
        y-intercept and scatter best-fit values for each of the M data sets

    cov : ndarray of shape (2M + 1, 2M + 1)
        Covariance matrix of the fit

    mcmc_samples : ndarray of shape (2M + 1, P)
        Flattened values from the MCMC chain for each of the 2M + 1 parameters

    hf : hyperfit object
    '''
    
    M = len(x) # Number of data sets

    ############################################################################
    # Pack data for MultiLinFit class
    #---------------------------------------------------------------------------
    datasets = [] # list of (2xN) arrays of the data
    covs = []     # list of covariance matrices for each data set
    
    for j in range(M):
        
        #-----------------------------------------------------------------------
        # Create covariance matrix (This assumes that x and y are independent.)
        #-----------------------------------------------------------------------
        N = len(x[j]) # Number of data points in the jth data set
        
        cov = np.empty((2, 2, N))

        for i in range(N):
            cov[:,:,i] = np.array([[dx[j][i]**2, 0.], [0., dy[j][i]**2]])
        #-----------------------------------------------------------------------
        covs.append(cov)
        #-----------------------------------------------------------------------
        # Create data object
        #-----------------------------------------------------------------------
        data = np.empty((2, N))
        
        data[0] = x[j]
        data[1] = y[j]
        #-----------------------------------------------------------------------
        datasets.append(data)
    ############################################################################


    ############################################################################
    # Create the hyperfit object
    #---------------------------------------------------------------------------
    hf = MultiLinFit(datasets, covs)
    ############################################################################


    ############################################################################
    # Run MCMC to fit line
    #---------------------------------------------------------------------------
    mcmc_samples, mcmc_lnlike = hf.emcee(bounds, verbose=True)
    ############################################################################


    ############################################################################
    # Calculate parameter values and covariance matrix
    #---------------------------------------------------------------------------
    best_fits = np.median(mcmc_samples, axis=1)
    a = best_fits[0]
    b = best_fits[1:M+1]
    sig = best_fits[M+1:]

    cov = np.cov(mcmc_samples)
    ############################################################################

    return a, b, sig, cov, mcmc_samples, hf
################################################################################
################################################################################




################################################################################
# Fit by minimizing chi2
#-------------------------------------------------------------------------------
def calculate_chi2(params, x, y, dx, dy):
    '''
    Calculate the chi2 value of the current line.


    PARAMETERS
    ==========

    params : list of length 2
        [slope, y-intercept]

    x, y : ndarrays of shape (N,)
        (x, y) values to fit the line to.

    dx, dy : ndarrays of shape (N,)
        Uncertainties in x and y.


    RETURNS
    =======

    chi2 : float
        Chi2 value of the current best fit
    '''


    a, b = params


    ############################################################################
    # Calculate the values of the current best-fit line
    #---------------------------------------------------------------------------
    y_fit = a*x + b
    ############################################################################


    ############################################################################
    # Calculate chi2 of the current fit
    #---------------------------------------------------------------------------
    chi2 = np.sum((y - y_fit)**2/(dy**2 + (a*dx)**2))
    ############################################################################

    return chi2



def chi2_line(x, y, dx, dy, p0, limits):
    '''
    Fit a line, y = mx + b, to the data by minimizing the chi2.  This accepts 
    uncertainties in both x and y.


    PARAMETERS
    ==========

    x, y : ndarrays of shape (N,)
        (x, y) values to fit the line to.

    dx, dy : ndarrays of shape (N,)
        Uncertainties in x and y.

    p0 : list of length 2
        Initial guesses for the slope and y-intercept

    limits : list of length 2
        Lower and upper limits for the slope and y-intercept:
            [(slope minimum, slope maximum), (y-intercept minimum, y-intercept maximum)]


    RETURNS
    =======

    a, b : floats
        slope and y-intercept of the best-fit

    cov : ndarray of shape (2,2)
        Covariance matrix of the fit parameters

    result : 
        Output from scipy.optimize.minimize
    '''


    ############################################################################
    # Fit using scipy.optimize.minimize
    #---------------------------------------------------------------------------
    result = minimize(calculate_chi2, 
                      p0, 
                      args=(x, y, dx, dy), 
                      bounds=limits)
    ############################################################################


    ############################################################################
    # Extract the best-fit parameters and their uncertainties
    #---------------------------------------------------------------------------
    a = result.x[0]
    b = result.x[1]

    '''
    uncertainties = np.sqrt(np.diag(result.hess_inv.to_dense()))
    da = uncertainties[0]
    db = uncertainties[1]
    '''
    cov = result.hess_inv.todense()
    ############################################################################

    return a, b, cov, result
################################################################################
################################################################################







################################################################################
# Fit using linear least squares
#-------------------------------------------------------------------------------
def lls_line(x, y, dy):
    '''
    Fit a line, y = mx + b, to the data by least squares.  This accepts 
    uncertainties only in y (the dependent variable).


    PARAMETERS
    ==========

    x, y : ndarrays of shape (N,)
        (x, y) values to fit the line to.

    dy : ndarray of shape (N,)
        Uncertainties in y.


    RETURNS
    =======

    a, b : floats
        slope and y-intercept best-fit values

    cov : ndarray of shape (2,2)
        Covariance matrix between a and b
    '''

    A = np.vander(x, 2)

    ATA = np.dot(A.T, A/(dy**2)[:, None])

    w = np.linalg.solve(ATA, np.dot(A.T, y/dy**2))

    # Calculate the covariance matrix (for uncertainties)
    cov = np.linalg.inv(ATA)


    ############################################################################
    # Extract the best-fit parameters and their uncertainties
    #---------------------------------------------------------------------------
    a = w[0]
    b = w[1]
    ############################################################################

    return a, b, cov
################################################################################
################################################################################





################################################################################
# Custom linear least squares with scatter in x
#-------------------------------------------------------------------------------
def chi2_sigmax(params, x, y, dx, dy):
    '''
    Calculate the chi2 for the line, y = ax + b, assuming that there is some 
    additional scatter in the x-direction.


    PARAMETERS
    ==========

    params : list of length 3
        best-fit parameters: [a, b, sigma_x]

    x, y : ndarrays of shape (N,)
        (x, y) data points to which we are fitting the line

    dx, dy : ndarrays of shape (N,)
        Uncertainties in x and y


    RETURNS
    =======

    chi2 : float
        chi2 for the current values in params
    '''


    a, b, sigma_x = params

    Delta_y = y - (a*x + b)

    delta2 = dy**2 + a**2 * (dx**2 + sigma_x**2)

    chi2 = np.sum(Delta_y**2 / delta2)

    return chi2



def chi2_line_sigmax(x, y, dx, dy, p0, limits):
    '''
    Fit a line, y = ax + b, assuming that there is some additional scatter in 
    the x-direction.  This fit accepts uncertainties in both x and y.


    PARAMETERS
    ==========

    x, y : ndarrays of shape (N,)
        (x, y) data points to which we are fitting the line

    dx, dy : ndarrays of shape (N,)
        Uncertainties in x and y

    p0 : list of length 3
        Initial guesses for the slope, y-intercept, and x-scatter

    limits : list of length 3
        Lower and upper limits for the slope, y-intercept, and x-scatter:
            [(slope minimum, slope maximum), 
             (y-intercept minimum, y-intercept maximum), 
             (x-scatter minimum, x-scatter maximum)]


    RETURNS
    =======

    a : float
        best-fit value for the slope

    b : float
        best-fit value for the y-intercept

    sigma_x : float
        best-fit value for the additional scatter in the x-direction

    cov : ndarray of shape (3,3)
        Covariance matrix between a, b, and sigma_x

    result : 
        Output from scipy.optimize.minimize
    '''


    ############################################################################
    # Fit using scipy.optimize.minimize
    #---------------------------------------------------------------------------
    result = minimize(chi2_sigmax, 
                      p0, 
                      args=(x, y, dx, dy), 
                      bounds=limits)
    ############################################################################


    ############################################################################
    # Extract the best-fit parameters
    #---------------------------------------------------------------------------
    a = result.x[0]
    b = result.x[1]
    sigma_x = result.x[2]

    '''
    uncertainties = np.sqrt(np.diag(result.hess_inv.to_dense()))
    da = uncertainties[0]
    db = uncertainties[1]
    '''
    cov = result.hess_inv.todense()
    ############################################################################

    return a, b, sigma_x, cov, result
################################################################################
################################################################################
















