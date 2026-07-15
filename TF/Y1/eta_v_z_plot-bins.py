'''
Plot of the log distance ratio v. z for the calibrated TFR galaxies, but only 
plot the binned values.
'''

################################################################################
# Import modules
#-------------------------------------------------------------------------------
from astropy.io import fits
from astropy.table import Table
from astropy.cosmology import FlatLambdaCDM
import astropy.constants as const

import numpy as np

import matplotlib.pyplot as plt

# Custom functions
import sys
sys.path.insert(1, '/Users/kdouglass/Documents/Research/DESI/PV_survey/DESI_SGA/TF/')
from help_functions import profile_histogram
################################################################################



################################################################################
# Import data
#-------------------------------------------------------------------------------
# data_directory = '/global/cfs/cdirs/desi/science/td/pv/tfgalaxies/Y1/'
data_directory = '/Users/kdouglass/Documents/Research/data/DESI/Y1/'

SGA_TF = Table.read(data_directory + 'DESI-DR1_TF_pv_cat_v13.fits') #<-- v13 - cosmology catalog
# SGA_TF = Table.read('SGA_iron_jointTFR_moduli-v18_20260708.fits')

SGA_TF_latest = Table.read('SGA_iron_jointTFR_moduli-v18_20260708.fits')

# Plot those in the main cosmology sample differently
sample1 = SGA_TF['MAIN']
sample1_latest = SGA_TF_latest['MAIN']
################################################################################



################################################################################
# Cosmology
#-------------------------------------------------------------------------------
LightSpeed = const.c.to('km/s').value

h = 1
H0 = 100*h

cosmo = FlatLambdaCDM(H0=H0, Om0=0.3151)
################################################################################



################################################################################
# PV lines (code taken from Cullen)
#-------------------------------------------------------------------------------
velarray = np.arange(-400, 401, 100)
zarray = np.linspace(0.001, 0.11, 200)

dzarray = cosmo.comoving_distance(zarray).value
dharray = cosmo.comoving_distance(np.outer(1.0/(1.0 + velarray/LightSpeed), 
                                  (1.0 + zarray)) - 1.0).value
deltamarray = np.log10(dzarray/dharray)

rotation = [20.0, 25.0, 30.0, 45.0, 0.0, -45.0, -30.0, -25.0, -20.0]
labels = ["-400", "-300", "-200", "-100", "0", "100", "200", "300", "400"]
xcoord = np.array([26000.0, 20000.0, 15000.0, 8000.0, -1000.0, 8000.0, 15000.0, 20000.0, 26000.0])
coord = np.searchsorted(zarray, xcoord/LightSpeed)
ycoord = np.array([deltamarray[i,j] for i, j in enumerate(coord)])
################################################################################




################################################################################
# Bin galaxies
#-------------------------------------------------------------------------------
_zbins = np.arange(0, 0.105, 0.005)
dz = 0.5*np.diff(_zbins)
zc = 0.5*(_zbins[1:] + _zbins[:-1])

N, y_avg, y_std = profile_histogram(SGA_TF['Z_DESI_CMB'][sample1], 
                                    SGA_TF['LOGDIST'][sample1], 
                                    _zbins, 
                                    weights=SGA_TF['LOGDIST_ERR'][sample1]**-2, 
                                    weighted=True)

N_alt, y_avg_alt, y_std_alt = profile_histogram(SGA_TF_latest['Z_DESI_CMB'][sample1_latest], 
                                                SGA_TF_latest['LOGDIST'][sample1_latest], 
                                                _zbins, 
                                                weights=SGA_TF_latest['LOGDIST_ERR'][sample1_latest]**-2, 
                                                weighted=True)

################################################################################



"""
################################################################################
# Save data for figure
#-------------------------------------------------------------------------------
# Build the header
#-------------------------------------------------------------------------------
hdr = fits.Header()

hdr['DESI_DR'] = 'DR1'
# hdr['FIGURE'] = 9
hdr['FIGURE'] = 14

empty_primary = fits.PrimaryHDU(header=hdr)
#-------------------------------------------------------------------------------
# Binned galaxy data
#-------------------------------------------------------------------------------
data_table_hdu = fits.BinTableHDU(data=Table([zc, dz, y_avg, y_std], 
                                             names=['Z', 'DZ', 'ETA', 'ETA_ERR']))
'''
pv_lines_hdu = fits.BinTableHDU(data=Table([zarray, deltamarray.T], 
                                           names=['Z', 'ETA']))

pv_labels_hdu = fits.BinTableHDU(data=Table([velarray, xcoord, ycoord, rotation, labels], 
                                            names=['V', 'X', 'Y', 'ROT', 'LABEL']))

hdul = fits.HDUList([empty_primary, data_table_hdu, pv_lines_hdu, pv_labels_hdu])

hdul.writeto('paper_figures/Fig9/fig9_data.fits', overwrite=True)
'''
hdul = fits.HDUList([empty_primary, data_table_hdu])

hdul.writeto('paper_figures/Fig14/fig14_data.fits', overwrite=True)
################################################################################
exit()
"""


################################################################################
# Plot
#-------------------------------------------------------------------------------
fig = plt.figure(tight_layout=True)

#-------------------------------------------------------------------------------
# Plot the weighted mean
#-------------------------------------------------------------------------------
# Cosmology sample
plt.errorbar(zc, y_avg, xerr=dz, yerr=y_std, fmt='o', 
             color='turquoise', 
             label='Original')

# Latest catalog
plt.errorbar(zc, y_avg_alt, xerr=dz, yerr=y_std_alt, fmt='x', 
             color='darkblue', 
             label='Updated')

plt.legend()

#-------------------------------------------------------------------------------
# PV lines (code taken from Cullen)
#-------------------------------------------------------------------------------
colors = 0.8*np.fabs(velarray)/np.amax(np.fabs(velarray)) + 0.2

for v in range(len(velarray)):
    c = plt.cm.YlOrRd(colors[v])

    plt.plot(zarray, deltamarray[v,:], 
             color=c, 
             linestyle='-', 
             alpha=0.7, 
             zorder=0)

    if (v != 4):
        plt.text(xcoord[v]/LightSpeed, ycoord[v], 
                 labels[v], 
                 color=c, 
                 fontsize=12, 
                 rotation=rotation[v], 
                 ha="center", 
                 va="center", 
                 bbox=dict(boxstyle="square", ec="w", fc="w"), 
                 zorder=1)

# Line at eta = 0
plt.hlines(0, 0, 0.2, linestyles='dashed', colors='k', zorder=5)
#-------------------------------------------------------------------------------

plt.grid(ls=':')

plt.tick_params(axis='both', which='major', labelsize=14)

plt.xlabel(r'$z_{\text{CMB}}$', fontsize=16)
plt.ylabel(r'$\eta = \log \left( \frac{D_z}{D_{\text{TFR}}} \right)$', 
           fontsize=16)

plt.ylim((-0.05, 0.05))
# plt.ylim((-0.13, 0.13))
plt.xlim((0, 0.105))

# plt.show()

plt.savefig('../../../figures/Y1_papers/iron_logdist-v-z_bins_v13-v18.png', 
            dpi=150, 
            facecolor='none')
################################################################################