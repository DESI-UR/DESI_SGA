'''
Plot of the log distance ratio v. z for the calibrated TFR galaxies.
'''

################################################################################
# Import modules
#-------------------------------------------------------------------------------
from astropy.table import Table

import numpy as np

import matplotlib.pyplot as plt

# Custom functions
import sys
sys.path.insert(1, '/Users/kdouglass/Documents/Research/DESI/PV_survey/code/TF/')
from help_functions import profile_histogram
################################################################################



################################################################################
# Import data
#-------------------------------------------------------------------------------
# data_directory = '/global/cfs/cdirs/desi/science/td/pv/tfgalaxies/Y1/'
data_directory = '/Users/kdouglass/Documents/Research/data/DESI/Y1/'

SGA_TF = Table.read(data_directory + 'DESI-DR1_TF_pv_cat_v14.fits')

# Plot those in the main cosmology sample differently
sample1 = SGA_TF['MAIN']
################################################################################



################################################################################
# Plot
#-------------------------------------------------------------------------------
_zbins = np.arange(0, 0.105, 0.005)
dz = 0.5*np.diff(_zbins)
zc = 0.5*(_zbins[1:] + _zbins[:-1])

fig = plt.figure()

gs = fig.add_gridspec(2, 2, width_ratios=(4,1), height_ratios=(1,4), 
	                  left=0.1, right=0.9, bottom=0.1, top=0.9, 
	                  wspace=0.05, hspace=0.05)
#-------------------------------------------------------------------------------
# Main plot
#-------------------------------------------------------------------------------
ax = fig.add_subplot(gs[1,0])

ax.errorbar(SGA_TF['Z_DESI_CMB'][~sample1], 
             SGA_TF['LOGDIST'][~sample1], 
             xerr=SGA_TF['ZERR_DESI'][~sample1], 
             yerr=SGA_TF['LOGDIST_ERR'][~sample1],
             fmt='.', 
             color='gray',
             alpha=0.1, 
             ecolor='gray')
ax.errorbar(SGA_TF['Z_DESI_CMB'][sample1], 
             SGA_TF['LOGDIST'][sample1], 
             xerr=SGA_TF['ZERR_DESI'][sample1], 
             yerr=SGA_TF['LOGDIST_ERR'][sample1],
             fmt='.', 
             markersize=4, 
             alpha=0.3, 
             ecolor='gray')

# Plot the weighted mean
N, y_avg, y_std = profile_histogram(SGA_TF['Z_DESI_CMB'][sample1], 
                                    SGA_TF['LOGDIST'][sample1], 
                                    _zbins, 
                                    weights=SGA_TF['LOGDIST_ERR'][sample1]**-2, 
                                    weighted=True)
ax.errorbar(zc, y_avg, xerr=dz, yerr=y_std, fmt='x', 
            color='darkblue', 
            # color='m',
            label='weighted mean', 
            zorder=10)

# Line at eta = 0
ax.hlines(0, 0, 0.2, linestyles='dashed', colors='k', zorder=5)

ax.grid(ls=':')

plt.tick_params(axis='both', which='major', labelsize=12)

ax.set_xlabel(r'$z_{\text{CMB}}$', fontsize=14)
ax.set_ylabel(r'$\eta = \log \left( \frac{D_z}{D_{\text{TFR}}} \right)$', fontsize=14)

ax.set_ylim((-0.5, 0.5))
ax.set_xlim((0, 0.19))

ax.grid(ls=':')
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# redshift histogram
#-------------------------------------------------------------------------------
ax_histx = fig.add_subplot(gs[0,0], sharex=ax)

ax_histx.hist(SGA_TF['Z_DESI_CMB'][sample1], 
              bins=np.arange(0, 0.175, 0.005))
ax_histx.hist(SGA_TF['Z_DESI_CMB'][~sample1], 
              bins=np.arange(0, 0.175, 0.005), 
              color='darkgray')
ax_histx.hist(SGA_TF['Z_DESI_CMB'][sample1], 
              bins=np.arange(0, 0.175, 0.005), 
              color='tab:blue', 
              histtype='step')

ax_histx.tick_params(axis='x', labelbottom=False)
ax_histx.set_ylabel('count')
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# log distance histogram
#-------------------------------------------------------------------------------
ax_histy = fig.add_subplot(gs[1,1], sharey=ax)

ax_histy.hist(SGA_TF['LOGDIST'][sample1], 
              bins=np.arange(-2, 2, 0.05), 
              orientation='horizontal')
ax_histy.hist(SGA_TF['LOGDIST'][~sample1], 
              bins=np.arange(-2, 2, 0.05), 
              color='darkgray', 
              orientation='horizontal')
ax_histy.hist(SGA_TF['LOGDIST'][sample1], 
              bins=np.arange(-2, 2, 0.05), 
              orientation='horizontal', 
              histtype='step',
              color='tab:blue')

ax_histy.tick_params(axis='y', labelleft=False)
ax_histy.set_xlabel('count');
#-------------------------------------------------------------------------------

# plt.show()

plt.savefig('../../../figures/Y1_papers/iron_logdist-v-z_cutsAlex_dz0p005_weightsVmax-1_20250926.png', 
            dpi=150, 
            facecolor='none', 
            bbox_inches='tight');
################################################################################