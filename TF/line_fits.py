'''
This is a collection of the various methods of fitting a line to data.
'''


################################################################################
# Import modules
#-------------------------------------------------------------------------------
import numpy as np

from hyperfit.linfit import LinFit

from scipy.optimize import minimize
################################################################################




def param_invert(w0, w1, cov):
    '''
    Convert the slope and y-intercept from the fit
                                x = w0*y + w1
    to the fit
                                y = a*x + b


    PARAMETERS
    ==========

    w0, w1 : floats
        slope (w0) and y-intercept (w1) of the fit x = w0*y + w1

    cov : ndarray of shape (2,2)
        Covariance matrix of w0 and w1


    RETURNS
    =======

    a, b : floats
        slope (a) and y-intercept (b) of the fit y = a*x + b

    cov_ab : ndarray of shape (2,2)
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
    w0_rng, w1_rng = rng.multivariate_normal([w0, w1], cov, 5000).T

    # Transform them from (w0, w1) space to (a, b) space
    a_rng = 1./w0_rng
    b_rng = -w1_rng/w0_rng

    # Calculate the covariance matrix for (a, b)
    ab_rng = np.stack((a_rng, b_rng))
    cov_ab = np.cov(ab_rng)
    ############################################################################

    return a, b, cov_ab




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
    # Create the hyperfit object
    #---------------------------------------------------------------------------
    hf = LinFit([x, y], cov)
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

    return a, b, sig, cov, mcmc_samples
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
    cov = result.hess_inv.to_dense()
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


















