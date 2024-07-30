'''
Run a "fit" of the TFR on simulated data (SGA_TFR_simtest_0x.fits from Segev) to 
test how well hyperfit recovers the slope.
'''

################################################################################
# Import modules
#-------------------------------------------------------------------------------
import numpy as np

import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = False
from matplotlib import cm, colors
from matplotlib.patches import Ellipse
import matplotlib as mpl

from astropy.table import Table

import corner

from line_fits import param_invert, hyperfit_line

import sys
sys.path.insert(1, '/global/u1/k/kadglass/DESI_SGA/TF/SV/')
# sys.path.insert(1, '/Users/kdouglass/Documents/Research/DESI/PV_survey/code/TF/SV')
from help_functions import adjust_lightness
################################################################################




################################################################################
# Constants
#-------------------------------------------------------------------------------
h = 1
H0 = 100*h

c = 3e5

q0 = 0.2

V0 = 2.5 # Set 0-pt of the TFR
################################################################################




################################################################################
# Run through the various files and fit the TFR
#-------------------------------------------------------------------------------
simtest_directory = '/global/cfs/cdirs/desi/science/td/pv/tfgalaxies/SV/sim/'

bounds_itfr = ((-1000.0, 1000.0), (-10.0, 10.0), (1.0e-5, 500.0))

a_itfr = np.zeros(10)
a_err = np.zeros(10)
b_itfr = np.zeros(10)
b_err = np.zeros(10)
s_itfr = np.zeros(10)

for i in range(10):
    
    simtest_num = i + 1
    
    ############################################################################
    # Read in simulated data
    #---------------------------------------------------------------------------
    if simtest_num < 10:
        tdata = Table.read(simtest_directory + 'SGA_TFR_simtest_00' + str(simtest_num) + '.fits')
    else:
        tdata = Table.read(simtest_directory + 'SGA_TFR_simtest_0' + str(simtest_num) + '.fits')
    ############################################################################
    
    
    ############################################################################
    # Fit the ITFR using hyperfit
    #---------------------------------------------------------------------------
    w0, w1, sig_w, cov_w, itfr_mcmc_samples, hf_itfr = hyperfit_line(tdata['R_MAG_SB26'], 
                                                                     np.log10(tdata['V_0p33R26']) - V0, 
                                                                     tdata['R_MAG_SB26_ERR'], 
                                                                     0.434*tdata['V_0p33R26_err']/tdata['V_0p33R26'], 
                                                                     bounds_itfr)
    ############################################################################
    
    
    ############################################################################
    # Calculate the best-fit values
    #---------------------------------------------------------------------------
    a_itfr[i], b_itfr[i], cov_itfr = param_invert(w0, w1, cov_w[:2,:2])
    
    a_err[i] = np.sqrt(cov_itfr[0,0])
    
    b_err[i] = np.sqrt(cov_itfr[1,1])
    
    s_itfr[i] = np.mean(itfr_mcmc_samples[2])
    ############################################################################
################################################################################




################################################################################
# Save the best-fit values
#-------------------------------------------------------------------------------
best_fits = Table()

best_fits['test'] = np.arange(10) + 1
best_fits['slope'] = a_itfr
best_fits['slope_err'] = a_err
best_fits['y_int'] = b_itfr
best_fits['y_int_err'] = b_err
best_fits['scatter'] = s_itfr

best_fits.write('SGA_TFR_simtest_hyperfit_results.txt', 
                format='ascii.commented_header', 
                overwrite=True)
################################################################################