'''
Plot the N(z) distribution of the TF sample.
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

z_bins = np.arange(0, 0.175, 0.005)

plt.hist(SGA_TF['Z_DESI_CMB'][sample1], 
         bins=z_bins, 
         color='darkblue', 
         label='Main Sample')
plt.hist(SGA_TF['Z_DESI_CMB'][~sample1], 
         bins=z_bins, 
         color='darkgray', 
         label='Dwarfs')
plt.hist(SGA_TF['Z_DESI_CMB'][sample1], 
         bins=z_bins, 
         color='darkblue', 
         histtype='step')

plt.legend()

plt.xlabel(r'$z_{\text{CMB}}$', fontsize=16)
plt.ylabel('number of galaxies', fontsize=16)

plt.tick_params(axis='both', which='major', labelsize=14)

# plt.show()

plt.savefig('../../../figures/Y1_papers/iron_Nz-distribution_v13.png', 
            dpi=150, 
            facecolor='none');
################################################################################