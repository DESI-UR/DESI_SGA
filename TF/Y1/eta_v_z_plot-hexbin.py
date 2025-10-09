'''
Plot of the log distance ratio v. z for the calibrated TFR galaxies, with the 
galaxies shown as hexbins for their density (instead of scattered points).
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

gs = fig.add_gridspec(2, 3, width_ratios=(4,1,0.3), height_ratios=(1,4), 
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
             ecolor='gray', 
             zorder=1)

hx = ax.hexbin(SGA_TF['Z_DESI_CMB'][sample1], 
          SGA_TF['LOGDIST'][sample1], 
          cmap='plasma', 
          mincnt=1, 
          vmin=1, 
          vmax=70, 
          gridsize=(30,25), 
          extent=(0, 0.2, -0.5, 0.5), 
          zorder=2)

# Plot the weighted mean
N, y_avg, y_std = profile_histogram(SGA_TF['Z_DESI_CMB'][sample1], 
                                    SGA_TF['LOGDIST'][sample1], 
                                    _zbins, 
                                    weights=SGA_TF['LOGDIST_ERR'][sample1]**-2, 
                                    weighted=True)
ax.errorbar(zc, y_avg, xerr=dz, yerr=y_std, fmt='x', 
            color='w', 
            # color='m',
            label='weighted mean', 
            zorder=3)

# Line at eta = 0
ax.hlines(0, 0, 0.2, linestyles='dashed', colors='lightgray', zorder=4)

ax.grid(ls=':')

plt.tick_params(axis='both', which='major', labelsize=12)

ax.set_xlabel(r'$z_{\text{CMB}}$', fontsize=14)
ax.set_ylabel(r'$\eta = \log \left( \frac{D_z}{D_{\text{TFR}}} \right)$', fontsize=14)

ax.set_ylim((-0.21, 0.21))
ax.set_xlim((0, 0.13))

ax.grid(ls=':')
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# colorbar
#-------------------------------------------------------------------------------
ax_cb = fig.add_subplot(gs[:,2])

plt.colorbar(mappable=hx, cax=ax_cb, label='Number of galaxies')
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# redshift histogram
#-------------------------------------------------------------------------------
ax_histx = fig.add_subplot(gs[0,0], sharex=ax)

ax_histx.hist(SGA_TF['Z_DESI_CMB'][sample1], 
              bins=np.arange(0, 0.175, 0.005), 
              color='darkblue')
ax_histx.hist(SGA_TF['Z_DESI_CMB'][~sample1], 
              bins=np.arange(0, 0.175, 0.005), 
              color='darkgray')
ax_histx.hist(SGA_TF['Z_DESI_CMB'][sample1], 
              bins=np.arange(0, 0.175, 0.005), 
              color='darkblue', 
              histtype='step')

ax_histx.tick_params(axis='x', labelbottom=False)
ax_histx.set_ylabel('count')
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# log distance histogram
#-------------------------------------------------------------------------------
ax_histy = fig.add_subplot(gs[1,1], sharey=ax)

ax_histy.hist(SGA_TF['LOGDIST'][sample1], 
              bins=np.arange(-2, 2, 0.02), 
              orientation='horizontal', 
              color='darkblue')
ax_histy.hist(SGA_TF['LOGDIST'][~sample1], 
              bins=np.arange(-2, 2, 0.02), 
              color='darkgray', 
              orientation='horizontal')
ax_histy.hist(SGA_TF['LOGDIST'][sample1], 
              bins=np.arange(-2, 2, 0.02), 
              orientation='horizontal', 
              histtype='step',
              color='darkblue')

ax_histy.tick_params(axis='y', labelleft=False)
ax_histy.set_xlabel('count');
#-------------------------------------------------------------------------------

# plt.show()

plt.savefig('../../../figures/Y1_papers/iron_logdist-v-z_cutsAlex_dz0p005_weightsVmax-1_20250926-hexbin.png', 
            dpi=150, 
            facecolor='none', 
            bbox_inches='tight');
################################################################################