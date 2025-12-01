'''
Plot histograms of the log-distance ratio and uncertainty for the TF sample.
'''

################################################################################
# Import modules
#-------------------------------------------------------------------------------
from astropy.table import Table

import numpy as np

import matplotlib.pyplot as plt
################################################################################



################################################################################
# Import data
#-------------------------------------------------------------------------------
# data_directory = '/global/cfs/cdirs/desi/science/td/pv/tfgalaxies/Y1/'
data_directory = '/Users/kdouglass/Documents/Research/data/DESI/Y1/'

SGA_TF = Table.read(data_directory + 'DESI-DR1_TF_pv_cat_v13.fits')

# Plot those in the main cosmology sample differently
sample1 = SGA_TF['MAIN']
################################################################################



################################################################################
# Plot
#-------------------------------------------------------------------------------
plt.figure(tight_layout=True)

#-------------------------------------------------------------------------------
# Histogram of the log-distance ratio
#-------------------------------------------------------------------------------
ax = plt.subplot(121)

eta_bins = np.arange(-1, 1, 0.05)

plt.hist(SGA_TF['LOGDIST'][sample1], 
         bins=eta_bins, 
         color='darkblue')
plt.hist(SGA_TF['LOGDIST'][~sample1], 
         bins=eta_bins, 
         color='darkgray')
plt.hist(SGA_TF['LOGDIST'][sample1], 
         bins=eta_bins, 
         histtype='step',
         color='darkblue')

plt.ylabel('number of galaxies', fontsize=16)
plt.xlabel(r'$\eta = \log \left( \frac{D_z}{D_{\text{TFR}}} \right)$', 
	       fontsize=16)

plt.tick_params(axis='both', which='major', labelsize=14)
#-------------------------------------------------------------------------------
# Histogram of the uncertainties in the log-distance ratio
#-------------------------------------------------------------------------------
ax2 = plt.subplot(122, sharey=ax)

sigma_bins = np.arange(0.09, 0.2, 0.003)

plt.hist(SGA_TF['LOGDIST_ERR'][sample1], 
         bins=sigma_bins, 
         color='darkblue', 
         label='Main Sample')
plt.hist(SGA_TF['LOGDIST_ERR'][~sample1], 
         bins=sigma_bins, 
         color='darkgray', 
         label='Dwarfs')
plt.hist(SGA_TF['LOGDIST_ERR'][sample1], 
         bins=sigma_bins, 
         histtype='step',
         color='darkblue')

plt.legend()

plt.xlabel(r'$\sigma_{\eta}$', fontsize=16)

ax2.tick_params(axis='y', labelleft=False)
plt.tick_params(axis='x', which='major', labelsize=14)
#-------------------------------------------------------------------------------


# plt.show()

plt.savefig('../../../figures/Y1_papers/iron_logdist_hists_v13.png', 
            dpi=150, 
            facecolor='none');
################################################################################