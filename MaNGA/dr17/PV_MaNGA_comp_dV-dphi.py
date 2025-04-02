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
from astropy.table import Table
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
tf_targets = Table.read('/global/cfs/cdirs/desi/science/td/pv/desi_pv_tf_fuji_healpix.fits')
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# SGA
#-------------------------------------------------------------------------------
SGA = Table.read('/global/cfs/cdirs/cosmo/data/sga/2020/SGA-2020.fits', 
                 'ELLIPSE')

SGA_dict = {}
for i in range(len(SGA)):
    SGA_dict[SGA['SGA_ID'][i]] = i
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# SDSS MaNGA cross-match file
#-------------------------------------------------------------------------------
SGA_MaNGA = Table.read('../MaNGA_SGA_crossmatch_2022-06-28.txt', 
                       format='ascii.commented_header')
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# SDSS MaNGA best fits
#
# Nitya's fitting from 2023 on DR17
#-------------------------------------------------------------------------------
MaNGA_fits = Table.read('master_table_Halpha_BB_HI_H2_MxCG_R90_CMD.txt', 
                        format='ascii.commented_header')
#-------------------------------------------------------------------------------
################################################################################



################################################################################
# Calculate the rotational velocities for the DESI galaxies
#-------------------------------------------------------------------------------
# Separate the DESI data into center and off-center observations
#-------------------------------------------------------------------------------
tf_targets['SKY_FIBER_DIST'] = 0.
tf_targets['SKY_FIBER_DIST_R26'] = 0.

# For each SGA galaxy that has at least one observation, calculate the distance 
# for all of that galaxy's targets
for sga_id in np.unique(tf_targets['SGA_ID']):
    
    # Identify all galaxy targets on this galaxy
    obs_idx = tf_targets['SGA_ID'] == sga_id
    
    # Find galaxy index in SGA catalog
    sga_idx = SGA_dict[sga_id]
    
    #---------------------------------------------------------------------------
    # Calculate distance between each observation and the center of the SGA 
    # galaxy
    #---------------------------------------------------------------------------
    SGA_coords = SkyCoord(ra=SGA['RA'][sga_idx], 
                          dec=SGA['DEC'][sga_idx], 
                          unit=u.degree)
    target_coords = SkyCoord(ra=tf_targets['RA'][obs_idx], 
                             dec=tf_targets['DEC'][obs_idx], 
                             unit=u.degree)
    
    sep2d = target_coords.separation(SGA_coords)
    
    tf_targets['SKY_FIBER_DIST'][obs_idx] = sep2d
    tf_targets['SKY_FIBER_DIST_R26'][obs_idx] = 2*sep2d.to('arcmin')/(SGA['D26'][sga_idx]*u.arcmin)
    #---------------------------------------------------------------------------
    
centers_boolean = tf_targets['SKY_FIBER_DIST_R26'] < 0.1

centers = tf_targets[centers_boolean]
axis = tf_targets[~centers_boolean]
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Clean the DESI center observations
#
# Only keep those with
#  - DELTACHI2 > 25
#  - ZWARN == 0
#-------------------------------------------------------------------------------
good_centers = centers[(centers['DELTACHI2'] > 25) & (centers['ZWARN'] == 0)]
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# If an SGA galaxy has at least one observation at its center, set the redshift 
# of that galaxy
#-------------------------------------------------------------------------------
SGA['Z_DESI'] = np.nan
SGA['ZERR_DESI'] = np.nan

weights = 1./(good_centers['ZERR']**2)

for sga_id in np.unique(good_centers['SGA_ID']):
    
    # Find all the center observations of this galaxy
    obs_idx = good_centers['SGA_ID'] == sga_id
    
    # Find the row in SGA for this galaxy
    SGA_idx = SGA_dict[sga_id]
    
    # Set the redshift of this galaxy to be the weighted average redshift of all 
    # good center observations
    SGA['Z_DESI'][SGA_idx] = np.average(good_centers['Z'][obs_idx], 
                                        weights=weights[obs_idx])
    SGA['ZERR_DESI'][SGA_idx] = np.sqrt(1./np.sum(weights[obs_idx]))


SGA_MaNGA['Z_DESI'] = np.nan
SGA_MaNGA['ZERR_DESI'] = np.nan
SGA_MaNGA['R26'] = np.nan
SGA_MaNGA['BA'] = np.nan
SGA_MaNGA['PA'] = np.nan

for i in range(len(SGA_MaNGA)):
    
    # Find the row in SGA for this galaxy
    SGA_idx = SGA_dict[SGA_MaNGA['SGA_ID'][i]]
    
    # Set the redshift of this galaxy
    SGA_MaNGA['Z_DESI'][i] = SGA['Z_DESI'][SGA_idx]
    SGA_MaNGA['ZERR_DESI'][i] = SGA['ZERR_DESI'][SGA_idx]
    
    # Transfer R26 over to the SGA_MaNGA table
    SGA_MaNGA['R26'][i] = 0.5*SGA['D26'][SGA_idx]
    
    # Transfer b/a over to the SGA_MaNGA table
    SGA_MaNGA['BA'][i] = SGA['BA'][SGA_idx]
    
    # Transfer phi over to the SGA_MaNGA table
    SGA_MaNGA['PA'][i] = SGA['PA'][SGA_idx]
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Calculate the rotational velocity
#-------------------------------------------------------------------------------
axis['SKY_FIBER_DIST'] = 0.
axis['SKY_FIBER_DIST_R26'] = 0.
axis['V_ROT'] = np.nan
axis['V_ROT_ERR'] = np.nan


# For each SGA galaxy that has at least one center observation, calculate the 
# distance for all of that galaxy's targets
for sga_gal in np.unique(centers['SGA_ID']):
    
    # Identify all galaxy targets on this galaxy
    obs_idx = axis['SGA_ID'] == sga_gal
    
    # Find galaxy index in SGA catalog
    sga_idx = SGA_dict[sga_gal]
    
    #---------------------------------------------------------------------------
    # Calculate distance between each observation and the center
    #---------------------------------------------------------------------------
    center_coords = SkyCoord(ra=SGA['RA'][sga_idx], 
                             dec=SGA['DEC'][sga_idx], 
                             unit=u.degree)
    target_coords = SkyCoord(ra=axis['RA'][obs_idx], 
                             dec=axis['DEC'][obs_idx], 
                             unit=u.degree)
    
    sep2d = target_coords.separation(center_coords)
    
    axis['SKY_FIBER_DIST'][obs_idx] = sep2d.to('radian')
    axis['SKY_FIBER_DIST_R26'][obs_idx] = 2*sep2d.to('arcmin')/(SGA['D26'][sga_idx]*u.arcmin)
    #---------------------------------------------------------------------------
    
    
    #---------------------------------------------------------------------------
    # Calculate rotational velocity
    #---------------------------------------------------------------------------
    # Use the average redshift of all center observations for the systemic velocity
    z_center = np.mean(SGA['Z_DESI'][sga_idx])
    z_err_center2 = SGA['ZERR_DESI'][sga_idx]**2

    # Calculate rotational velocity for all observations of the galaxy
    axis['V_ROT'][obs_idx] = c*(axis['Z'][obs_idx] - z_center)
    axis['V_ROT_ERR'][obs_idx] = c*np.sqrt(axis['ZERR'][obs_idx]**2 + z_err_center2)
    #---------------------------------------------------------------------------
    
    
    #---------------------------------------------------------------------------
    # Correct rotational velocities for inclination angle
    #---------------------------------------------------------------------------
    cosi2 = (SGA['BA'][sga_idx]**2 - q0**2)/(1 - q0**2)
    
    # Galaxies with b/a < q0
    if cosi2 < 0:
        cosi2 = 0
    
    axis['V_ROT'][obs_idx] /= np.sin(np.arccos(np.sqrt(cosi2)))
    #---------------------------------------------------------------------------
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Just keep those velocities measured at 0.33 R26
#-------------------------------------------------------------------------------
axis_0p3 = axis[(axis['SKY_FIBER_DIST_R26'] > 0.3) & (axis['SKY_FIBER_DIST_R26'] < 0.4)]
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Remove "bad" galaxies
#
# Those with
#  - 10 > V > 1000 km/s
#  - Delta V / Vmin < 5
#-------------------------------------------------------------------------------
axis_0p3_goodV = axis_0p3[(np.abs(axis_0p3['V_ROT']) < 1000) & (np.abs(axis_0p3['V_ROT']) > 10)]

good_deltaV = np.ones(len(axis_0p3_goodV), dtype=bool)

for sga_id in np.unique(axis_0p3_goodV['SGA_ID']):
    
    # Identify all galaxy targets on this galaxy
    obs_idx = axis_0p3_goodV['SGA_ID'] == sga_id
    
    n_obs = np.sum(obs_idx)
    
    if n_obs > 1:
        
        Vmin = np.min(np.abs(axis_0p3_goodV['V_ROT'][obs_idx]))
        Vmax = np.max(np.abs(axis_0p3_goodV['V_ROT'][obs_idx]))
        
        v_norm_min = np.abs(axis_0p3_goodV['V_ROT'][obs_idx])/Vmin
        v_norm_max = np.abs(axis_0p3_goodV['V_ROT'][obs_idx])/Vmax
        
        diff_matrix = np.abs(axis_0p3_goodV['V_ROT'][obs_idx]).reshape(n_obs, 1) - np.abs(axis_0p3_goodV['V_ROT'][obs_idx]).reshape(1, n_obs)
        
        diff_matrix_norm = diff_matrix/Vmin
        
        if np.any(np.abs(diff_matrix_norm) > 5.):
            
            # Remove all observations with DELTACHI2 < 25
            # Note: This also typically removes observations with ZWARN != 0
            deltachi2_idx = axis_0p3_goodV['DELTACHI2'] >= 25
            
            good_deltaV[obs_idx & ~deltachi2_idx] = False
            
            good_obs_idx = obs_idx & deltachi2_idx
            
            n_obs_good = np.sum(good_obs_idx)
            
            # Check to make sure that, if there are still multiple observations, they all satisfy our relative velocity criteria
            if n_obs_good > 1:
                
                Vmin = np.min(np.abs(axis_0p3_goodV['V_ROT'][good_obs_idx]))
                
                diff_matrix = np.abs(axis_0p3_goodV['V_ROT'][good_obs_idx]).reshape(n_obs_good, 1) - np.abs(axis_0p3_goodV['V_ROT'][good_obs_idx]).reshape(1, n_obs_good)
                
                diff_matrix_norm = diff_matrix/Vmin
                
                if np.any(np.abs(diff_matrix_norm) > 5.):

                    # Set all of these so that we don't look at this galaxy
                    good_deltaV[good_obs_idx] = False
                    
axis_0p3_good = axis_0p3_goodV[good_deltaV]
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Remove those that have been visually inspected to be suspicious
#-------------------------------------------------------------------------------
VI_remove = Table.read('../../TF/SV/fuji_VI.txt', format='ascii.commented_header')

remove_targets = np.zeros(len(axis_0p3_good), dtype=bool)

for targetid in VI_remove['TARGETID']:
    
    remove_targets = remove_targets & (axis_0p3_good['TARGETID'] == targetid)
    
VI_axis_0p3_good = axis_0p3_good[~remove_targets]
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Add the rotational velocities to the SDSS MaNGA - SGA cross-match file
#-------------------------------------------------------------------------------
SGA_MaNGA['V_0p33R26'] = np.nan
SGA_MaNGA['V_0p33R26_ERR'] = np.nan

SGA_MaNGA['SKY_FIBER_DIST'] = np.nan

weights = 1./(VI_axis_0p3_good['V_ROT_ERR']**2)

for i in range(len(SGA_MaNGA)):
    
    # Does this galaxy have any observations?
    i_obs = VI_axis_0p3_good['SGA_ID'] == SGA_MaNGA['SGA_ID'][i]
    
    if np.sum(i_obs) > 0:
        
        # Average all velocities at this radius
        SGA_MaNGA['V_0p33R26'][i] = np.average(np.abs(VI_axis_0p3_good['V_ROT'][i_obs]), 
                                               weights=weights[i_obs])
        SGA_MaNGA['V_0p33R26_ERR'][i] = np.sqrt(1./np.sum(weights[i_obs]))
        
        # Copy over the distance from the center for this observation
        SGA_MaNGA['SKY_FIBER_DIST'][i] = np.average(VI_axis_0p3_good['SKY_FIBER_DIST'][i_obs])
#-------------------------------------------------------------------------------
################################################################################



################################################################################
# Add the MaNGA best-fit values to the table
#-------------------------------------------------------------------------------
SGA_MaNGA['Vmax_map'] = np.nan
SGA_MaNGA['Vmax_err_map'] = np.nan
SGA_MaNGA['Rturn_map'] = np.nan
SGA_MaNGA['alpha_map'] = np.nan

SGA_MaNGA['ba_map'] = np.nan
SGA_MaNGA['ba_err_map'] = np.nan
SGA_MaNGA['ba_NSA'] = np.nan

SGA_MaNGA['phi_map'] = np.nan
SGA_MaNGA['phi_err_map'] = np.nan
SGA_MaNGA['phi_NSA'] = np.nan

SGA_MaNGA['Z_NSA'] = np.nan

for i in range(len(SGA_MaNGA)):
    
    gal_id = SGA_MaNGA['plateifu'][i]
    
    # Find galaxy row in MaNGA fits table
    plate_bool = MaNGA_fits['MaNGA_plate'] == SGA_MaNGA['plate'][i]
    ifu_bool = MaNGA_fits['MaNGA_IFU'] == SGA_MaNGA['ifudsgn'][i]
    
    i_fit = plate_bool & ifu_bool
    
    # Copy best-fit parameter values from fit table to galaxy table
    if (np.sum(i_fit) > 0): #and (gal_id not in []):
        SGA_MaNGA['Vmax_map'][i] = MaNGA_fits['Vmax_map'][i_fit]
        SGA_MaNGA['Vmax_err_map'][i] = MaNGA_fits['Vmax_err_map'][i_fit]
        SGA_MaNGA['Rturn_map'][i] = MaNGA_fits['Rturn_map'][i_fit]
        SGA_MaNGA['alpha_map'][i] = MaNGA_fits['alpha_map'][i_fit]
        
        SGA_MaNGA['ba_map'][i] = MaNGA_fits['ba_map'][i_fit]
        SGA_MaNGA['ba_err_map'][i] = MaNGA_fits['ba_err_map'][i_fit]
        SGA_MaNGA['ba_NSA'][i] = MaNGA_fits['NSA_ba'][i_fit]
        
        SGA_MaNGA['phi_map'][i] = MaNGA_fits['phi_map'][i_fit]
        SGA_MaNGA['phi_err_map'][i] = MaNGA_fits['phi_err_map'][i_fit]
        SGA_MaNGA['phi_NSA'][i] = MaNGA_fits['NSA_phi'][i_fit]
        
        SGA_MaNGA['Z_NSA'][i] = MaNGA_fits['NSA_redshift'][i_fit]
        
# Flip all -99 values to NaN
SGA_MaNGA['Vmax_map'][SGA_MaNGA['Vmax_map'] == -999] = np.nan

good_V = np.isfinite(SGA_MaNGA['Vmax_map']) & np.isfinite(SGA_MaNGA['V_0p33R26'])


#-------------------------------------------------------------------------------
# Remove those galaxies from the final sample that have V(R26) < 0.9Vmax, and 
# have Vmax > 1000 km/s
#-------------------------------------------------------------------------------
# 1 - Convert R26 to kpc for each galaxy
dist_to_galaxy = SGA_MaNGA['Z_DESI']*c/H0
R26_kpc = dist_to_galaxy.to('kpc')*np.tan((SGA_MaNGA['R26']*u.arcmin).to(u.rad))

# 2 - Compute V(R26)
SGA_MaNGA['Vfit_R26'] = rot_fit_BB(R26_kpc.data, [SGA_MaNGA['Vmax_map'], SGA_MaNGA['Rturn_map'], SGA_MaNGA['alpha_map']])

# 3 - Filter out those with V(R26) < 0.9Vmax
goodVmax = SGA_MaNGA['Vfit_R26'] >= 0.9*SGA_MaNGA['Vmax_map']

# 4 - Filter out those with Vmax > 1000 km/s
lowVmax = SGA_MaNGA['Vmax_map'] < 1000.
'''
# 5 - Filter out those with alpha > 99
good_alpha = SGA_MaNGA['alpha_map'] < 99.

# 6 - Filter out those with large uncertainties in Vmax
goodVmax2 = SGA_MaNGA['Vmax_err_map']/SGA_MaNGA['Vmax_map'] <= 2
''';
final_sample = good_V & goodVmax & lowVmax #& good_alpha & goodVmax2
#-------------------------------------------------------------------------------
################################################################################



################################################################################
# Compute V(0.33R26) for the MaNGA fit, correcting for (b/a) differences
#-------------------------------------------------------------------------------
# Convert R26 to kpc for each galaxy
dist_to_galaxy = SGA_MaNGA['Z_DESI']*c/H0
p33R26_kpc = dist_to_galaxy.to('kpc')*np.tan(SGA_MaNGA['SKY_FIBER_DIST'])

SGA_MaNGA['Vfit_0p33R26'] = rot_fit_BB(p33R26_kpc.data, 
                                       [SGA_MaNGA['Vmax_map'], 
                                        SGA_MaNGA['Rturn_map'], 
                                        SGA_MaNGA['alpha_map']])

delta_pa = np.abs((SGA_MaNGA['phi_map']%180.) - SGA_MaNGA['PA'])
delta_pa[delta_pa > 90.] = 180. - delta_pa[delta_pa > 90.]
SGA_MaNGA['delta_phi'] = delta_pa

cosi2_sga = (SGA_MaNGA['BA']**2 - q0**2)/(1 - q0**2)
cosi2_manga = (SGA_MaNGA['ba_map']**2 - q0**2)/(1 - q0**2)

cosi2_sga[cosi2_sga < 0] = 0.
cosi2_manga[cosi2_manga < 0] = 0.

# DESI_corrected = (1./np.cos(delta_pa*np.pi/180.))*(np.sin(np.arccos(np.sqrt(cosi2_sga)))/np.sin(np.arccos(np.sqrt(cosi2_manga))))*SGA_MaNGA['V_0p33R26']

SGA_MaNGA['Vfit_corr_0p33R26'] = (np.sin(np.arccos(np.sqrt(cosi2_manga)))/np.sin(np.arccos(np.sqrt(cosi2_sga))))*SGA_MaNGA['Vfit_0p33R26']

#-------------------------------------------------------------------------------
# Compute the uncertainty in the velocity at V(0.33 R26) for the MaNGA fit
#-------------------------------------------------------------------------------
SGA_MaNGA['Vfit_corr_0p33R26_err'] = np.nan

# Convert b/a to i (to match what is in the Hessians)
SGA_MaNGA['i_map'] = np.arccos(np.sqrt((SGA_MaNGA['ba_map']**2 - q0**2)/(1 - q0**2)))

for i in range(len(SGA_MaNGA)):
    
    if np.isfinite(SGA_MaNGA['Vmax_map'][i]) and np.isfinite(SGA_MaNGA['V_0p33R26'][i]):
        
        gal_ID = str(SGA_MaNGA['plate'][i]) + '-' + str(SGA_MaNGA['ifudsgn'][i])

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

            random_sample = np.random.multivariate_normal(mean=[SGA_MaNGA['i_map'][i],
                                                                SGA_MaNGA['phi_map'][i], 
                                                                SGA_MaNGA['Vmax_map'][i], 
                                                                SGA_MaNGA['Rturn_map'][i], 
                                                                SGA_MaNGA['alpha_map'][i]], 
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
            delta_pa_sample = np.abs((good_randoms[:,1]%180.) - SGA_MaNGA['PA'][i])
            delta_pa_sample[delta_pa_sample > 90.] = 180. - delta_pa_sample[delta_pa_sample > 90.]
            
            y_sample = np.cos(delta_pa_sample*np.pi/180)*(np.sin(good_randoms[:,0])/np.sin(np.arccos(np.sqrt(cosi2_sga[i]))))*y_sample

            SGA_MaNGA['Vfit_corr_0p33R26_err'][i] = np.std(y_sample, axis=0)
            
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

gs = fig.add_gridspec(1, 2, width_ratios=(4, 1), left=0.1, right=0.9, bottom=0.1, top=0.9, wspace=0.05)

ax = fig.add_subplot(gs[0, 0])
ax.errorbar(np.abs(SGA_MaNGA['delta_phi'][final_sample]), 
            SGA_MaNGA['V_0p33R26'][final_sample] - SGA_MaNGA['Vfit_corr_0p33R26'][final_sample],
            xerr=SGA_MaNGA['phi_err_map'][final_sample], 
            yerr=np.sqrt(SGA_MaNGA['V_0p33R26_ERR'][final_sample]**2 + SGA_MaNGA['Vfit_corr_0p33R26_err'][final_sample]**2), 
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
ax_histy.hist(SGA_MaNGA['V_0p33R26'][final_sample] - SGA_MaNGA['Vfit_corr_0p33R26'][final_sample], 
              bins=np.arange(-100, 150, 10), 
              orientation='horizontal')

ax_histy.tick_params(axis='y', labelleft=False)

plt.savefig('../../Figures/MaNGA_dr17/MaNGA_fuji_deltaV0p33R26-vs-deltaPhi_wHist_20240424.png', 
            dpi=150, bbox_inches='tight')
# bbox_inches='tight' is being used to keep savefig() from cutting off the axis 
# labels; tight_layout wasn't working with these subplots.
################################################################################