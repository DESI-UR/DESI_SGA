'''
Script to create the figure comparing delta V (0.33R26) v. delta phi (MaNGA - 
SGA) for SDSS DR17.  Based off of the notebook `PV_MaNGA_comp.ipynb`.

This figure is published in the DESI SV TFR paper.
'''

################################################################################
# Import modules
#-------------------------------------------------------------------------------
import numpy as np
import numpy.ma as ma

from astropy.io import fits
from astropy.table import Table, join
from astropy.coordinates import SkyCoord
import astropy.constants as const
import astropy.units as u
from astropy import wcs

import scipy.stats as stats
from scipy.optimize import minimize, curve_fit

import numdifftools as ndt

from hyperfit.linfit import LinFit

import corner

import os

# Path to RotationCurve libraries. Update as needed.
#rotcurvepath = os.path.join(os.environ['HOME'], 'desi/RotationCurves/spirals')
rotcurvepath = os.path.join(os.environ['HOME'], 'RotationCurves/spirals')
if not os.path.exists(rotcurvepath):
    raise FileNotFoundError(f'{rotcurvepath} does not exist.')

import sys
sys.path.insert(1, rotcurvepath)
from dark_matter_mass_v1 import rot_fit_BB

from tqdm import tqdm

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib.lines as mlines
################################################################################



################################################################################
# Plotting settings
#-------------------------------------------------------------------------------
mpl.rc('font', size=12)
################################################################################



################################################################################
# Constants
#-------------------------------------------------------------------------------
h = 1
H0 = 100*h*u.km/u.s/u.Mpc

c = const.c.to('km/s')

q0 = 0.2

MANGA_SPAXEL_SIZE = 0.5*u.arcsec
################################################################################



################################################################################
# Import data
#-------------------------------------------------------------------------------
# DESI
#-------------------------------------------------------------------------------
tf_targets = Table.read('../../../TF/SV/SGA-2020_fuji_Vrot_VI_dVsys_photsys_corr.fits')

SGA_dict = {}
for i in range(len(tf_targets)):
    SGA_dict[tf_targets['SGA_ID'][i]] = i
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# SDSS MaNGA cross-match file
#-------------------------------------------------------------------------------
SGA_MaNGA = Table.read('../../MaNGA_SGA_crossmatch_2022-06-28.txt', 
                       format='ascii.commented_header')
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# SDSS MaNGA best fits
#
# Nitya's fitting from 2023 on DR17
#-------------------------------------------------------------------------------
MaNGA_fits = Table.read('../Elliptical_sphdisk_refitspirals_BPT_illustris_v11.fits')
#-------------------------------------------------------------------------------
################################################################################



################################################################################
# Combine tables
#-------------------------------------------------------------------------------
# Add the rotational velocities
#-------------------------------------------------------------------------------
SGA_MaNGA_pv = join(SGA_MaNGA, 
                    tf_targets['SGA_ID', 'Z_DESI', 'ZERR_DESI', 'D26', 'BA', 'PA', 'V_0p33R26', 'V_0p33R26_ERR'], 
                    keys='SGA_ID')

SGA_MaNGA_pv['R26'] = 0.5*SGA_MaNGA_pv['D26']

SGA_MaNGA_pv['SKY_FIBER_DIST'] = (0.33*SGA_MaNGA_pv['R26']*u.arcmin).to('radian')
#-------------------------------------------------------------------------------
# Add the MaNGA best-fit values to the table
#-------------------------------------------------------------------------------
SGA_MaNGA_pv_fits = join(SGA_MaNGA_pv, 
                         MaNGA_fits['plateifu', 'v_max', 'v_max_err', 'r_turn', 'alpha', 'ba', 'ba_err', 'nsa_elpetro_ba', 'phi', 'phi_err', 'nsa_elpetro_phi', 'nsa_z'], 
                         keys='plateifu')
        
# Flip all -99 values to NaN
for col_name in SGA_MaNGA_pv_fits.colnames:
    bad_values = SGA_MaNGA_pv_fits[col_name] == -999
    if np.any(bad_values):
        SGA_MaNGA_pv_fits[col_name][bad_values] = np.nan

good_V = np.isfinite(SGA_MaNGA_pv_fits['v_max']) & np.isfinite(SGA_MaNGA_pv_fits['V_0p33R26'])
#-------------------------------------------------------------------------------
# Remove those galaxies from the final sample that have V(R26) < 0.9Vmax, and 
# have Vmax > 1000 km/s
#-------------------------------------------------------------------------------
# 1 - Convert R26 to kpc for each galaxy
dist_to_galaxy = SGA_MaNGA_pv_fits['Z_DESI']*c/H0
R26_kpc = dist_to_galaxy.to('kpc')*np.tan((SGA_MaNGA_pv_fits['R26']*u.arcmin).to(u.rad))

# 2 - Compute V(R26)
SGA_MaNGA_pv_fits['Vfit_R26'] = rot_fit_BB(R26_kpc.data, 
                                           [SGA_MaNGA_pv_fits['v_max'], SGA_MaNGA_pv_fits['r_turn'], SGA_MaNGA_pv_fits['alpha']])

# 3 - Filter out those with V(R26) < 0.9Vmax
goodVmax = SGA_MaNGA_pv_fits['Vfit_R26'] >= 0.9*SGA_MaNGA_pv_fits['v_max']

# 4 - Filter out those with Vmax > 1000 km/s
lowVmax = SGA_MaNGA_pv_fits['v_max'] < 1000.

final_sample = good_V & goodVmax & lowVmax
#-------------------------------------------------------------------------------
################################################################################



################################################################################
# Compute V(0.33R26) for the MaNGA fit, correcting for (b/a) differences
#-------------------------------------------------------------------------------
# Convert R26 to kpc for each galaxy
dist_to_galaxy_NSA = SGA_MaNGA_pv_fits['nsa_z']*c/H0
R26_kpc_NSA = dist_to_galaxy_NSA.to('kpc')*np.tan((SGA_MaNGA_pv_fits['R26']*u.arcmin).to(u.rad))

SGA_MaNGA_pv_fits['Vfit_0p33R26'] = rot_fit_BB((0.33*R26_kpc_NSA).data, 
                                               [SGA_MaNGA_pv_fits['v_max'], 
                                                SGA_MaNGA_pv_fits['r_turn'], 
                                                SGA_MaNGA_pv_fits['alpha']])

delta_pa = np.abs((SGA_MaNGA_pv_fits['phi']%180.) - SGA_MaNGA_pv_fits['PA'])
delta_pa[delta_pa > 90.] = 180. - delta_pa[delta_pa > 90.]
SGA_MaNGA_pv_fits['delta_phi'] = delta_pa

cosi2_sga = (SGA_MaNGA_pv_fits['BA']**2 - q0**2)/(1 - q0**2)
cosi2_manga = (SGA_MaNGA_pv_fits['ba']**2 - q0**2)/(1 - q0**2)

cosi2_sga[cosi2_sga < 0] = 0.
cosi2_manga[cosi2_manga < 0] = 0.

SGA_MaNGA_pv_fits['Vfit_corr_0p33R26'] = (np.sin(np.arccos(np.sqrt(cosi2_manga)))/np.sin(np.arccos(np.sqrt(cosi2_sga))))*SGA_MaNGA_pv_fits['Vfit_0p33R26']

#-------------------------------------------------------------------------------
# Compute the uncertainty in the velocity at V(0.33 R26) for the MaNGA fit
#-------------------------------------------------------------------------------
SGA_MaNGA_pv_fits['Vfit_corr_0p33R26_err'] = np.nan

# Convert b/a to i (to match what is in the Hessians)
SGA_MaNGA_pv_fits['i_map'] = np.arccos(np.sqrt((SGA_MaNGA_pv_fits['ba']**2 - q0**2)/(1 - q0**2)))

for i in range(len(SGA_MaNGA_pv_fits)):
    
    if np.isfinite(SGA_MaNGA_pv_fits['v_max'][i]) and np.isfinite(SGA_MaNGA_pv_fits['V_0p33R26'][i]):
        
        # gal_ID = str(SGA_MaNGA_pv_fits['plate'][i]) + '-' + str(SGA_MaNGA_pv_fits['ifudsgn'][i])
        gal_ID = SGA_MaNGA_pv_fits['plateifu'][i]

        try:
            Hessian = np.load('/global/u1/k/kadglass/RotationCurves/spirals/DRP_map_Hessians/dr17/' + gal_ID + '_Hessian.npy')
            hess_inv_all = 2*np.linalg.inv(Hessian)
            
            # Reconstruct the inverse Hessian to contain just the parameters that we need
            hess_inv = np.zeros((5,5))
            hess_inv[-4:,-4:] = hess_inv_all[-4:,-4:] # copies phi, Vmax, Rturn, and alpha
            hess_inv[0,0] = hess_inv_all[1,1] # copies i
            hess_inv[0,-4:] = hess_inv_all[1,-4:] # copies off-diagonal elements for i
            hess_inv[-4:,0] = hess_inv_all[-4:,1] # copies off-diagonal elements for i

            N_samples = 10000

            random_sample = np.random.multivariate_normal(mean=[SGA_MaNGA_pv_fits['i_map'][i],
                                                                SGA_MaNGA_pv_fits['phi'][i], 
                                                                SGA_MaNGA_pv_fits['v_max'][i], 
                                                                SGA_MaNGA_pv_fits['r_turn'][i], 
                                                                SGA_MaNGA_pv_fits['alpha'][i]], 
                                                          cov=hess_inv, 
                                                          size=N_samples)

            # Remove bad samples (those with negative values for any of the parameters)
            is_good_random = np.all(random_sample > 0, axis=1)
            good_randoms = random_sample[is_good_random, :]

            # Calculate values of curve at this location
            y_sample = rot_fit_BB(R26_kpc[i].value, [good_randoms[:,-3], 
                                                     good_randoms[:,-2], 
                                                     good_randoms[:,-1]])
            
            # Adjust for differences in i and phi
            delta_pa_sample = np.abs((good_randoms[:,1]%180.) - SGA_MaNGA_pv_fits['PA'][i])
            delta_pa_sample[delta_pa_sample > 90.] = 180. - delta_pa_sample[delta_pa_sample > 90.]
            
            y_sample = np.cos(delta_pa_sample*np.pi/180)*(np.sin(good_randoms[:,0])/np.sin(np.arccos(np.sqrt(cosi2_sga[i]))))*y_sample

            SGA_MaNGA_pv_fits['Vfit_corr_0p33R26_err'][i] = np.std(y_sample, axis=0)
            
        except (FileNotFoundError, np.linalg.LinAlgError) as error:
            print(gal_ID, error)
#-------------------------------------------------------------------------------
################################################################################



'''
################################################################################
# Plot Delta V (0.33 R26) v. Delta phi (MaNGA - SGA)
#-------------------------------------------------------------------------------
plt.figure(tight_layout=True)

plt.errorbar(np.abs(SGA_MaNGA['delta_phi'][final_sample]), 
             SGA_MaNGA['V_0p33R26'][final_sample] - SGA_MaNGA['Vfit_corr_0p33R26'][final_sample], 
             xerr=SGA_MaNGA['phi_err_map'][final_sample], 
             yerr=np.sqrt(SGA_MaNGA['V_0p33R26_ERR'][final_sample]**2 + SGA_MaNGA['Vfit_corr_0p33R26_err'][final_sample]**2), 
             fmt='o', 
             # alpha=0.5,
             ecolor='lightskyblue')
plt.hlines(0, 0., 90., linestyles='dotted', colors='k')

plt.xlim(0, 70)
#plt.ylim(-325, 75)

plt.tick_params(axis='both', which='major', labelsize=14)

plt.xlabel('$\Delta \phi$ (MaNGA - SGA) [deg]', fontsize=16)
plt.ylabel('$\Delta V(0.33R_{26})$ [km/s]', fontsize=16)

plt.savefig('../../Figures/MaNGA_dr17/MaNGA_fuji_deltaV0p33R26-vs-deltaPhi_20240424.png', 
            dpi=150)
################################################################################
'''

################################################################################
# Plot Delta V (0.33 R26) v. Delta phi (MaNGA - SGA) with a scatter plot for 
# Delta V
#-------------------------------------------------------------------------------
fig = plt.figure()

gs = fig.add_gridspec(1, 2, width_ratios=(4, 1), 
                      left=0.1, right=0.9, bottom=0.1, top=0.9, 
                      wspace=0.05)

ax = fig.add_subplot(gs[0, 0])
ax.errorbar(np.abs(SGA_MaNGA_pv_fits['delta_phi'][final_sample]), 
            SGA_MaNGA_pv_fits['V_0p33R26'][final_sample] - SGA_MaNGA_pv_fits['Vfit_corr_0p33R26'][final_sample],
            xerr=SGA_MaNGA_pv_fits['phi_err'][final_sample], 
            yerr=np.sqrt(SGA_MaNGA_pv_fits['V_0p33R26_ERR'][final_sample]**2 + SGA_MaNGA_pv_fits['Vfit_corr_0p33R26_err'][final_sample]**2), 
            fmt='o', 
            # alpha=0.5,
            ecolor='lightskyblue')
ax.hlines(0, 0., 90., linestyles='dotted', colors='k')

ax.set_xlim(0, 70)
#plt.ylim(-325, 75)

plt.tick_params(axis='both', which='major', labelsize=14)

ax.set_xlabel('$\Delta \phi$ (MaNGA - SGA) [deg]', fontsize=16)
ax.set_ylabel('$\Delta V(0.33R_{26})$ [km/s]', fontsize=16)


ax_histy = fig.add_subplot(gs[0, 1], sharey=ax)
ax_histy.hist(SGA_MaNGA_pv_fits['V_0p33R26'][final_sample] - SGA_MaNGA_pv_fits['Vfit_corr_0p33R26'][final_sample], 
              bins=np.arange(-100, 150, 10), 
              orientation='horizontal')

ax_histy.tick_params(axis='y', labelleft=False)

plt.savefig('../../../Figures/MaNGA_dr17/MaNGA_fuji_deltaV0p33R26-vs-deltaPhi_wHist_20260113.png', 
            dpi=150, facecolor='none', bbox_inches='tight')
# bbox_inches='tight' is being used to keep savefig() from cutting off the axis 
# labels; tight_layout wasn't working with these subplots.
################################################################################