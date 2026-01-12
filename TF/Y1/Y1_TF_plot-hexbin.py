'''
This file generates the Y1 TF plot with hexbin coloring instead of scattered 
points.
'''

################################################################################
# Import modules
#-------------------------------------------------------------------------------
import numpy as np

from astropy.cosmology import Planck18, LambdaCDM
from astropy.table import Table
import astropy.units as u

import pickle

import matplotlib.pyplot as plt
################################################################################



################################################################################
# Read in best-fit pickle file
#-------------------------------------------------------------------------------
temp_infile = open('cov_ab_iron_jointTFR_v14.pickle', 'rb')
cov_ab, tfr_samples, logV0, zmin, zmax, dz, zbins = pickle.load(temp_infile)
temp_infile.close()
################################################################################



################################################################################
# Read in galaxies
#-------------------------------------------------------------------------------
# data_directory = '/global/cfs/cdirs/desi/science/td/pv/tfgalaxies/Y1/'
data_directory = '/Users/kdouglass/Documents/Research/data/DESI/Y1/'

SGA_TF = Table.read(data_directory + 'DESI-DR1_TF_pv_cat_v14.fits')

# Plot those in the main cosmology sample differently
main = SGA_TF['MAIN']
################################################################################



################################################################################
# Define the cosmology
#-------------------------------------------------------------------------------
h = 1
H0 = 100*h

cosmo = LambdaCDM(H0=H0, 
                  Om0=Planck18.Om0, 
                  Tcmb0=Planck18.Tcmb0, 
                  Neff=Planck18.Neff, 
                  m_nu=Planck18.m_nu, 
                  Ob0=Planck18.Ob0, 
                  Ode0=Planck18.Ode0)
################################################################################



################################################################################
# Convert the y-intercepts to absolute magnitudes
#-------------------------------------------------------------------------------
# Center redshift values of each bin
# NOTE: zc should really use zbins[:-1], but we are dropping the last bin 
# (0.1-0.105) because we mistakenly used it while calibrating
zc = 0.5*dz + zbins[:-2]

# Distance modulus for each redshift bin center
mu_zc = cosmo.distmod(zc)

# Extract slope
slope = np.median(tfr_samples[0])
slope_err = np.sqrt(cov_ab[0,0])

# Each redshift bin has its own 0pt
# To put it in absolute-magnitude space, we'll convert it to an absolute 
# magnitude using the middle of the redshift bin
# NOTE: Again, ZP should really be median(tfr_samples[1:-1], axis=1), but we are 
# dropping the 0.1-0.105 redshift bin that we mistakenly used while calibrating
ZP = np.median(tfr_samples[1:-2], axis=1) - mu_zc.value
ZP_err = np.sqrt(np.diagonal(cov_ab[1:-2,1:-2])) # Should include z-bin width to this uncertainty

sig = np.median(tfr_samples[-1])

logv = np.linspace(-1*np.ones(len(zbins)-2), 3.5*np.ones(len(zbins)-2), 100)
absmag = slope*(logv - logV0) + ZP
################################################################################



################################################################################
# Plot
#-------------------------------------------------------------------------------
plt.figure(figsize=(5,7), tight_layout=True)

plt.grid(ls=':')

# plt.fill_between(logv, line_err[0], line_err[1], color='lightgray')

plt.errorbar(np.log10(SGA_TF['V_0p4R26'][~main]), 
             SGA_TF['R_ABSMAG_SB26'][~main], 
             xerr=0.434*SGA_TF['V_0p4R26_ERR'][~main]/SGA_TF['V_0p4R26'][~main],
             yerr=SGA_TF['R_ABSMAG_SB26_ERR'][~main], 
             fmt='.',
             color='gray',
             alpha=0.1, 
             ecolor='gray', 
             zorder=1)

plt.hexbin(np.log10(SGA_TF['V_0p4R26'][main]), 
           SGA_TF['R_ABSMAG_SB26'][main], 
           cmap='plasma', 
           mincnt=1, 
           vmin=1, 
           vmax=70, 
           gridsize=(70,80), 
           extent=(-0.1, 3.1, -25, -12.25), 
           zorder=2)

plt.colorbar(label='Number of galaxies')

plt.plot(logv, absmag, 'k', zorder=3)
plt.plot(logv, absmag + sig, 'k:', zorder=4)
plt.plot(logv, absmag - sig, 'k:', zorder=5)

plt.xlim([0.5, 3.1])
plt.ylim([-12.25, -24.5])

plt.xlabel('log($V(0.4R_{26})$ [km/s])', fontsize=14)
plt.ylabel('$M_r^{0.1} (26) - 5\log h$', fontsize=14);

ax = plt.gca()
ax.tick_params(axis='both', which='major', labelsize=12);

# plt.show()

plt.savefig('../../../figures/Y1_papers/iron_TFR_dz0p005_weightsVmax-1_cutsAlex_20250926-hexbin.png', 
            dpi=150, 
            facecolor='none')
################################################################################