'''
This file generates the TF calibration plot with each redshift bin in its own 
subplot.  It does not contain any of the fitting procedures.
'''

################################################################################
# Import modules
#-------------------------------------------------------------------------------
import numpy as np

from astropy.io import fits
from astropy.table import Table
import astropy.units as u

import pickle

import matplotlib.pyplot as plt
################################################################################



################################################################################
# Read in best-fit pickle file
#-------------------------------------------------------------------------------
temp_infile = open('cov_ab_iron_jointTFR_varyV0-dwarfsAlex_z0p1_zbins0p005_weightsVmax-1_dVsys_KAD-20250813.pickle', 
                   'rb')
cov_tfr, tfr_mcmc_samples, logV0, zmin, zmax, dz, zbins = pickle.load(temp_infile)
temp_infile.close()
################################################################################



################################################################################
# Read in galaxies
#-------------------------------------------------------------------------------
# data_directory = '/global/cfs/cdirs/desi/science/td/pv/tfgalaxies/Y1/'
data_directory = '/Users/kdouglass/Documents/Research/data/DESI/Y1/'

TF_Y1 = Table.read(data_directory + 'DESI-DR1_TF_pv_cat_v13.fits')
################################################################################



################################################################################
# Define the calibration sample
#-------------------------------------------------------------------------------
# Only keep those in the main sample
# This takes care of Alex's velocity and magnitude cut
main = TF_Y1['MAIN']

# Inclination > 45
q0 = 0.2
cosi2 = (TF_Y1['BA']**2 - q0**2) / (1 - q0**2)
i_min = 45. * u.degree
cosi2_max = np.cos(i_min)**2
is_good_incl = cosi2 < cosi2_max

# SSL spirals only
is_good_morph_ML = np.zeros_like(is_good_incl, dtype=bool)
for i in range(len(TF_Y1)):
    if TF_Y1['MORPHTYPE_AI'][i] == 'Spiral':
        is_good_morph_ML[i] = True

# John's VI
is_good_John = TF_Y1['JOHN_VI'].mask

good_gals = TF_Y1[main & is_good_incl & is_good_morph_ML & is_good_John]
################################################################################



################################################################################
# Assign galaxies to redshift bins
#-------------------------------------------------------------------------------
good_gals['Z_BIN_IDX'] = np.digitize(good_gals['Z_DESI_CMB'], zbins, right=True)

# Remove galaxies that fall outside the calibration range
no_use = (good_gals['Z_BIN_IDX'] == 0) | (good_gals['Z_BIN_IDX'] == len(zbins))
gals = good_gals[~no_use]

# Find list of all zbin indices
_zbin_ids = np.unique(gals['Z_BIN_IDX'])
m = len(_zbin_ids)

# Pack the galaxies into arrays
logV, logV_err = [], []
mag, mag_err = [], []

# Loop over the redshift bins
for k, _zbin_id in enumerate(_zbin_ids):

    select_zbin = gals['Z_BIN_IDX'] == _zbin_id

    logV.append(np.log10(gals['V_0p4R26'][select_zbin]) - logV0)
    logV_err.append(0.434*gals['V_0p4R26_ERR'][select_zbin] / gals['V_0p4R26'][select_zbin])

    mag.append(gals['R_MAG_SB26_CORR'][select_zbin])
    mag_err.append(gals['R_MAG_SB26_ERR_CORR'][select_zbin])
################################################################################


'''
for i in range(len(zbins) + 1):
    if i == 0:
        print(f'{i:2d}  z <= {zbins[i]:0.3f}  {np.sum(zbin_indices == i):3d} galaxies')
    elif i == len(zbins):
        print(f'{i:2d}  z > {zbins[i-1]:0.3f}  {np.sum(zbin_indices == i):3d} galaxies')
    else:
        print(f'{i:2d}  {zbins[i-1]:0.3f} < z <= {zbins[i]:0.3f}  {np.sum(zbin_indices == i):3d} galaxies')
'''

'''
################################################################################
# Save figure data for paper
#-------------------------------------------------------------------------------
# Header
#-------------------------------------------------------------------------------
hdr = fits.Header()

hdr['DESI_DR'] = 'DR1'
hdr['FIGURE'] = 6

empty_primary = fits.PrimaryHDU(header=hdr)
#-------------------------------------------------------------------------------
# Data
#-------------------------------------------------------------------------------
hdulist = [empty_primary]

for i in range(m):
    gal_hdu = fits.BinTableHDU(data=Table([logV[i], logV_err[i], mag[i], mag_err[i]], 
                                          names=('LOGV', 'LOGV_ERR', 'R_MAG', 'R_MAG_ERR')))
    hdulist.append(gal_hdu)

hdul = fits.HDUList(hdulist)

hdul.writeto('paper_figures/Fig6/fig6_data.fits', overwrite=True)
#-------------------------------------------------------------------------------
# Header
#-------------------------------------------------------------------------------
hdr2 = fits.Header()

hdr2['DESI_DR'] = 'DR1'
hdr2['FIGURE'] = 7
hdr2['LOG_V0'] = float(f'{logV0:.3f}')
#-------------------------------------------------------------------------------
# Data
#-------------------------------------------------------------------------------
mcmcfit_hdul = fits.PrimaryHDU(header=hdr2, data=tfr_mcmc_samples)

mcmcfit_hdul.writeto('paper_figures/Fig7/fig7_data.fits', overwrite=True)
################################################################################
exit()
'''


################################################################################
# Pull out best-fit MCMC samples
#-------------------------------------------------------------------------------
a_mcmc = np.percentile(tfr_mcmc_samples[0], [16., 50., 84])

b_mcmc = []
for k in range(1, m+1):
    b_mcmc.append(np.percentile(tfr_mcmc_samples[k], [16., 50., 84.]))
b_mcmc = np.asarray(b_mcmc)
################################################################################



################################################################################
# Plot
#-------------------------------------------------------------------------------
fig, axs = plt.subplots(nrows=2, ncols=7, 
                        sharex=True, sharey=True, 
                        figsize=(14,6), 
                        tight_layout=True)

_logv = np.arange(0, 3, 0.1) - logV0

a_   = a_mcmc[1]
b_   = b_mcmc[:,1]

# color = iter(plt.cm.cool(np.linspace(0,1,m)))
color = iter(plt.cm.winter(np.linspace(0,1,m)))
for i in range(m):
    c = next(color)
    
    row = int(i/7)
    col = i%7
    
    eb = axs[row,col].errorbar(logV0 + logV[i],
                               mag[i],
                               xerr=logV_err[i],
                               yerr=mag_err[i],
                               fmt='.', color=c, alpha=0.1)
    axs[row,col].plot(_logv + logV0, a_*_logv + b_[i], color=c)
    
    axs[row,col].set(xlim=[1.7, 2.5], 
                     ylim=[17.2, 12.5], 
                     title=f'{zbins[i]:.3f} < z $\leq$ {zbins[i+1]:.3f}')

# Delete extra axes
# fig.delaxes(axs[-1,-1])

fig.supxlabel(r'$\log{(V(0.4R_{26})~[\mathrm{km/s}]}$)', fontsize=16)
fig.supylabel(r'$m_r^{0.1} (26)$', fontsize=16);

# plt.show()

plt.savefig('../../../figures/Y1_papers/TF_Y1_zbin_calibration-array_v13-winter.png', 
            dpi=150, 
            facecolor='none')
################################################################################