'''
Plot the N(z) distribution of the TF sample.
'''


################################################################################
# Import modules
#-------------------------------------------------------------------------------
from astropy.table import Table
from astropy.io import fits

import numpy as np

import matplotlib.pyplot as plt
################################################################################



################################################################################
# Import data
#-------------------------------------------------------------------------------
# data_directory = '/global/cfs/cdirs/desi/science/td/pv/tfgalaxies/Y1/'
# data_directory = '/Users/kdouglass/Documents/Research/data/DESI/Y1/'

# v13 (cosmology catalog)
# SGA_TF = Table.read(data_directory + 'DESI-DR1_TF_pv_cat_v13.fits')

SGA_TF = Table.read('SGA_iron_jointTFR_moduli-v18_20260708.fits')

# Plot those in the main cosmology sample differently
sample1 = SGA_TF['MAIN']
################################################################################



################################################################################
# Plot
#-------------------------------------------------------------------------------
plt.figure(tight_layout=True)

z_bins = np.arange(0, 0.175, 0.005)

N, edges,_ = plt.hist(SGA_TF['Z_DESI_CMB'][sample1], 
                      bins=z_bins, 
                      color='darkblue', 
                      label='Main Sample')
N_dwarf,_,_ = plt.hist(SGA_TF['Z_DESI_CMB'][~sample1], 
                       bins=z_bins, 
                       color='darkgray', 
                       label='Dwarfs & Outliers')
plt.hist(SGA_TF['Z_DESI_CMB'][sample1], 
         bins=z_bins, 
         color='darkblue', 
         histtype='step')

plt.legend()

plt.xlabel(r'$z_{\text{CMB}}$', fontsize=16)
plt.ylabel('number of galaxies', fontsize=16)

plt.tick_params(axis='both', which='major', labelsize=14)

plt.show()

# plt.savefig('../../../figures/Y1_papers/iron_Nz-distribution_v18.png', 
#             dpi=150, 
#             facecolor='none');
################################################################################


################################################################################
# Save figure data
#-------------------------------------------------------------------------------
# Build header
#-------------------------------------------------------------------------------
hdr = fits.Header()

hdr['DESI_DR'] = 'DR1'
hdr['FIGURE'] = 12

empty_primary = fits.PrimaryHDU(header=hdr)
#-------------------------------------------------------------------------------
N_hdu = fits.BinTableHDU(data=Table([edges[:-1], N, N_dwarf], 
                         names=['BIN_EDGE', 'N_MAIN', 'N_DWARF']))
#-------------------------------------------------------------------------------
hdul = fits.HDUList([empty_primary, N_hdu])

hdul.writeto('paper_figures/Fig12/fig12_data.fits', overwrite=True)
################################################################################