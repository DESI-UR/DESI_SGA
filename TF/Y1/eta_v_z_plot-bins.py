'''
Plot of the log distance ratio v. z for the calibrated TFR galaxies, but only 
plot the binned values.
'''

################################################################################
# Import modules
#-------------------------------------------------------------------------------
from astropy.table import Table
from astropy.cosmology import FlatLambdaCDM
import astropy.constants as const

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
# Cosmology
#-------------------------------------------------------------------------------
LightSpeed = const.c.to('km/s').value

h = 1
H0 = 100*h

cosmo = FlatLambdaCDM(H0=H0, Om0=0.3151)
################################################################################



################################################################################
# Plot
#-------------------------------------------------------------------------------
fig = plt.figure(tight_layout=True)

#-------------------------------------------------------------------------------
# Plot the weighted mean
#-------------------------------------------------------------------------------
_zbins = np.arange(0, 0.105, 0.005)
dz = 0.5*np.diff(_zbins)
zc = 0.5*(_zbins[1:] + _zbins[:-1])

# N, y_avg, y_std = profile_histogram(SGA_TF['Z_DESI_CMB'], 
#                                     SGA_TF['LOGDIST'], 
#                                     _zbins, 
#                                     weights=SGA_TF['LOGDIST_ERR']**-2, 
#                                     weighted=True)
# plt.errorbar(zc, y_avg, xerr=dz, yerr=y_std, fmt='x', 
#              color='darkgray', 
#              label='weighted mean')

N, y_avg, y_std = profile_histogram(SGA_TF['Z_DESI_CMB'][sample1], 
                                    SGA_TF['LOGDIST'][sample1], 
                                    _zbins, 
                                    weights=SGA_TF['LOGDIST_ERR'][sample1]**-2, 
                                    weighted=True)
plt.errorbar(zc, y_avg, xerr=dz, yerr=y_std, fmt='x', 
             color='darkblue', 
             label='weighted mean')
#-------------------------------------------------------------------------------



#-------------------------------------------------------------------------------
# PV lines (code taken from Cullen)
#-------------------------------------------------------------------------------
velarray = np.arange(-400, 401, 100)
zarray = np.linspace(0.005, 0.11, 200)

dzarray = cosmo.comoving_distance(zarray).value
dharray = cosmo.comoving_distance(np.outer(1.0/(1.0 + velarray/LightSpeed), 
								  (1.0 + zarray)) - 1.0).value
deltamarray = np.log10(dzarray/dharray)

rotation = [20.0, 25.0, 30.0, 45.0, 0.0, -45.0, -30.0, -25.0, -20.0]
labels = ["-400", "-300", "-200", "-100", "0", "100", "200", "300", "400"]
xcoord = np.array([26000.0, 20000.0, 15000.0, 8000.0, -1000.0, 8000.0, 15000.0, 20000.0, 26000.0])
coord = np.searchsorted(zarray, xcoord/LightSpeed)
ycoord = np.array([deltamarray[i,j] for i, j in enumerate(coord)])

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

plt.show()

# plt.savefig('../../../figures/Y1_papers/iron_logdist-v-z_bins_v14.png', 
#             dpi=150, 
#             facecolor='none')
################################################################################