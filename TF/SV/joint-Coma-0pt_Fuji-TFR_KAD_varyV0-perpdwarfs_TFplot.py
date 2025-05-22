'''
This file is the same as joint-Coma-0pt_Fuji-TFR_KAD_varyV0-perpdwarfs.ipynb, 
but should only be used to generate the final TFR plot; it does not contain any 
of the fitting procedures.
'''


################################################################################
# Load modules
#-------------------------------------------------------------------------------
import numpy as np

import matplotlib.pyplot as plt
plt.rcParams["text.usetex"] = False
from matplotlib import cm, colors
from matplotlib.patches import Ellipse
# import matplotlib as mpl

from astropy.io import fits
from astropy.table import Table, join
from astropy.coordinates import SkyCoord
from astropy import units as u
import astropy.constants as const
# from astropy.visualization.wcsaxes import SphericalCircle

import corner

# import os

# import requests

import pickle

import sys
sys.path.insert(1, '/global/u1/k/kadglass/DESI_SGA/TF/')
# #sys.path.insert(1, '/Users/kellydouglass/Documents/Research/DESI/Targets/code/TF/')
from help_functions import adjust_lightness
# from line_fits import param_invert, hyperfit_line
from TF_photoCorrect import BASS_corr, MW_dust, k_corr, internal_dust
from z_CMB_convert import convert_z_frame
################################################################################




################################################################################
# Constants
#-------------------------------------------------------------------------------
h = 1
H0 = 100*h

c = const.c.to('km/s')

q0 = 0.2
################################################################################




################################################################################
# Import data
#-------------------------------------------------------------------------------
# fuji
#-------------------------------------------------------------------------------
tfuji = Table.read('/global/cfs/projectdirs/desi/science/td/pv/tfgalaxies/desi_pv_tf_fuji_healpix.fits')
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# SGA
#-------------------------------------------------------------------------------
SGA = Table.read('../SGA_distances.fits')

SGA_dict = {}

for i in range(len(SGA)):
    
    SGA_dict[SGA['SGA_ID'][i]] = i
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Best-fit
#-------------------------------------------------------------------------------
temp_infile = open('cov_ab_fuji_joint_TFR_varyV0-perpdwarfs0_AnthonyUpdates_weightsVmax-1_KAD.pickle', 
                   'rb')
cov_ab, tfr_samples, V0 = pickle.load(temp_infile)
temp_infile.close()

m_fit = np.median(tfr_samples[0])
b_fit = np.median(tfr_samples[1:3], axis=1)
sig_fit = np.median(tfr_samples[3:], axis=1)

w0_fit = 1/m_fit
#-------------------------------------------------------------------------------
################################################################################




################################################################################
# Separate the fuji data into center and off-center observations
#-------------------------------------------------------------------------------
tfuji['SKY_FIBER_DIST'] = 0.
tfuji['SKY_FIBER_DIST_R26'] = 0.

# For each SGA galaxy that has at least one observation, calculate the distance 
# for all of that galaxy's targets
for sga_id in np.unique(tfuji['SGA_ID']):
    
    # Identify all galaxy targets on this galaxy
    obs_idx = tfuji['SGA_ID'] == sga_id
    
    # Find galaxy index in SGA catalog
    sga_idx = SGA_dict[sga_id]
    
    #---------------------------------------------------------------------------
    # Calculate distance between each observation and the center of the SGA 
    # galaxy
    #---------------------------------------------------------------------------
    SGA_coords = SkyCoord(ra=SGA['RA'][sga_idx], 
                          dec=SGA['DEC'][sga_idx], 
                          unit=u.degree)
    target_coords = SkyCoord(ra=tfuji['RA'][obs_idx], 
                             dec=tfuji['DEC'][obs_idx], 
                             unit=u.degree)
    
    sep2d = target_coords.separation(SGA_coords)
    
    tfuji['SKY_FIBER_DIST'][obs_idx] = sep2d
    tfuji['SKY_FIBER_DIST_R26'][obs_idx] = 2*sep2d.to('arcmin')/(SGA['D26'][sga_idx]*u.arcmin)
    #---------------------------------------------------------------------------


centers_boolean = tfuji['SKY_FIBER_DIST_R26'] < 0.1

fuji_centers = tfuji[centers_boolean]
fuji_axis = tfuji[~centers_boolean]
################################################################################




################################################################################
# Clean the fuji center observations
#
# Only keep those with
#  - DELTACHI2 > 25
#  - ZWARN == 0
#-------------------------------------------------------------------------------
good_centers = fuji_centers[(fuji_centers['DELTACHI2'] > 25) & (fuji_centers['ZWARN'] == 0)]
################################################################################




################################################################################
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
################################################################################




################################################################################
# Coma cluster membership
#-------------------------------------------------------------------------------
# Data table #3 from Tully (2015)
#-------------------------------------------------------------------------------
hdu = fits.open('../Tully15-Table3.fits')
table3 = Table(hdu[1].data)
hdu.close()
#-------------------------------------------------------------------------------


Coma_nest = 100001

Coma_row_t3 = table3['Nest'] == Coma_nest

R2t_Coma = table3['R2t'][Coma_row_t3][0]
sigma_Coma = table3['sigV'][Coma_row_t3][0]
mu_Coma = table3['DM'][Coma_row_t3][0]


Coma_coords = SkyCoord(table3['SGLON'][Coma_row_t3]*u.degree, 
                       table3['SGLAT'][Coma_row_t3]*u.degree, 
                       frame='supergalactic').icrs

zHelio_Coma = convert_z_frame(table3['<Vcmba>'][Coma_row_t3][0]/c.value, 
                              Coma_coords.ra.deg, 
                              Coma_coords.dec.deg, 
                              corrtype='-full')[0]

V_Coma = 100 * 10**(0.2*(mu_Coma - 25)) / (1 + zHelio_Coma)


#-------------------------------------------------------------------------------
# Calculate the projected distance between the Coma cluster and each SGA galaxy
#-------------------------------------------------------------------------------
# First, we need to convert R2t from Mpc to an angle, using the group's velocity
# Note that we are NOT assuming that the size of the cluster is a small angle!!
R2t_Coma_angle_1p5 = np.arctan(1.5*R2t_Coma/(V_Coma/H0))*u.radian
R2t_Coma_angle_3 = np.arctan(3*R2t_Coma/(V_Coma/H0))*u.radian

SGA_coords = SkyCoord(SGA['RA'], SGA['DEC'], unit='deg')

sep = Coma_coords.separation(SGA_coords)

SGA_in_Coma1 = (sep < R2t_Coma_angle_1p5) & (SGA['Z_DESI']*c.value > V_Coma - 3*sigma_Coma) & (SGA['Z_DESI']*c.value < V_Coma + 3*sigma_Coma)

SGA_in_Coma2 = (sep >= R2t_Coma_angle_1p5) & (sep < R2t_Coma_angle_3) & (SGA['Z_DESI']*c.value > V_Coma - 2*sigma_Coma) & (SGA['Z_DESI']*c.value < V_Coma + 2*sigma_Coma)

SGA_in_Coma = SGA_in_Coma1 | SGA_in_Coma2
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Keep all observations of each galaxy that are within the Coma cluster
#-------------------------------------------------------------------------------
SGA_ID_in_Coma = SGA['SGA_ID'][SGA_in_Coma]

centers_inComa = good_centers[np.in1d(good_centers['SGA_ID'], SGA_ID_in_Coma)]

axis_inComa = fuji_axis[np.in1d(fuji_axis['SGA_ID'], SGA_ID_in_Coma)]
#-------------------------------------------------------------------------------
################################################################################



################################################################################
# 0-pt calibrators
#
# Which objects with center observations also have independent distances?
#-------------------------------------------------------------------------------
distances = (SGA['DM_Stellar'] != -1) | (SGA['DM1_SN'] != -1)
centers = np.isfinite(SGA['Z_DESI'])

#-------------------------------------------------------------------------------
# Keep all observations of each galaxy that have independent distances
#-------------------------------------------------------------------------------
SGA_ID_dist = SGA['SGA_ID'][distances & centers]

centers_dist = good_centers[np.in1d(good_centers['SGA_ID'], SGA_ID_dist)]

axis_dist = fuji_axis[np.in1d(fuji_axis['SGA_ID'], SGA_ID_dist)]
#-------------------------------------------------------------------------------
################################################################################




################################################################################
# Calculate the rotational velocity
#-------------------------------------------------------------------------------
# for Coma galaxies
#-------------------------------------------------------------------------------
axis_inComa['SKY_FIBER_DIST'] = 0.
axis_inComa['SKY_FIBER_DIST_R26'] = 0.
axis_inComa['V_ROT'] = np.nan
axis_inComa['V_ROT_ERR'] = np.nan


# For each SGA galaxy that has at least one center observation, calculate the 
# distance for all of that galaxy's targets
for sga_gal in np.unique(centers_inComa['SGA_ID']):
    
    # Identify all galaxy targets on this galaxy
    obs_idx = axis_inComa['SGA_ID'] == sga_gal
    
    # Find galaxy index in SGA catalog
    sga_idx = SGA_dict[sga_gal]
    
    #---------------------------------------------------------------------------
    # Calculate distance between each observation and the center
    #---------------------------------------------------------------------------
    center_coords = SkyCoord(ra=SGA['RA'][sga_idx], 
                             dec=SGA['DEC'][sga_idx], 
                             unit=u.degree)
    target_coords = SkyCoord(ra=axis_inComa['RA'][obs_idx], 
                             dec=axis_inComa['DEC'][obs_idx], 
                             unit=u.degree)
    
    sep2d = target_coords.separation(center_coords)
    
    axis_inComa['SKY_FIBER_DIST'][obs_idx] = sep2d
    axis_inComa['SKY_FIBER_DIST_R26'][obs_idx] = 2*sep2d.to('arcmin')/(SGA['D26'][sga_idx]*u.arcmin)
    #---------------------------------------------------------------------------
    
    
    #---------------------------------------------------------------------------
    # Calculate rotational velocity
    #---------------------------------------------------------------------------
    # Use the average redshift of all center observations for the systemic velocity
    z_center = np.mean(SGA['Z_DESI'][sga_idx])
    z_err_center2 = SGA['ZERR_DESI'][sga_idx]**2

    # Calculate rotational velocity for all observations of the galaxy
    z_rot = (1 + axis_inComa['Z'][obs_idx])/(1 + z_center) - 1
    axis_inComa['V_ROT'][obs_idx] = c*z_rot
    axis_inComa['V_ROT_ERR'][obs_idx] = c*np.sqrt((axis_inComa['ZERR'][obs_idx]/(1 + z_center))**2 + z_err_center2*((1 + axis_inComa['Z'][obs_idx])/(1 + z_center)**2))
    #---------------------------------------------------------------------------
    
    
    #---------------------------------------------------------------------------
    # Correct rotational velocities for inclination angle
    #---------------------------------------------------------------------------
    cosi2 = (SGA['BA'][sga_idx]**2 - q0**2)/(1 - q0**2)
    
    # Galaxies with b/a < q0
    if cosi2 < 0:
        cosi2 = 0
    
    axis_inComa['V_ROT'][obs_idx] /= np.sin(np.arccos(np.sqrt(cosi2)))
    #---------------------------------------------------------------------------
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# For 0-pt calibrators
#-------------------------------------------------------------------------------
axis_dist['SKY_FIBER_DIST'] = 0.
axis_dist['SKY_FIBER_DIST_R26'] = 0.
axis_dist['V_ROT'] = np.nan
axis_dist['V_ROT_ERR'] = np.nan


# For each SGA galaxy that has at least one center observation, calculate the 
# distance for all of that galaxy's targets
for sga_gal in np.unique(centers_dist['SGA_ID']):
    
    # Identify all galaxy targets on this galaxy
    obs_idx = axis_dist['SGA_ID'] == sga_gal
    
    # Find galaxy index in SGA catalog
    sga_idx = SGA_dict[sga_gal]
    
    #---------------------------------------------------------------------------
    # Calculate distance between each observation and the center
    #---------------------------------------------------------------------------
    center_coords = SkyCoord(ra=SGA['RA'][sga_idx], 
                             dec=SGA['DEC'][sga_idx], 
                             unit=u.degree)
    target_coords = SkyCoord(ra=axis_dist['RA'][obs_idx], 
                             dec=axis_dist['DEC'][obs_idx], 
                             unit=u.degree)
    
    sep2d = target_coords.separation(center_coords)
    
    axis_dist['SKY_FIBER_DIST'][obs_idx] = sep2d
    axis_dist['SKY_FIBER_DIST_R26'][obs_idx] = 2*sep2d.to('arcmin')/(SGA['D26'][sga_idx]*u.arcmin)
    #---------------------------------------------------------------------------
    
    
    #---------------------------------------------------------------------------
    # Calculate rotational velocity
    #---------------------------------------------------------------------------
    # Use the average redshift of all center observations for the systemic velocity
    z_center = np.mean(SGA['Z_DESI'][sga_idx])
    z_err_center2 = SGA['ZERR_DESI'][sga_idx]**2

    # Calculate rotational velocity for all observations of the galaxy
    # axis_dist['V_ROT'][obs_idx] = c*(axis_dist['Z'][obs_idx] - z_center)
    # axis_dist['V_ROT_ERR'][obs_idx] = c*np.sqrt(axis_dist['ZERR'][obs_idx]**2 + z_err_center2)
    z_rot = (1 + axis_dist['Z'][obs_idx])/(1 + z_center) - 1
    axis_dist['V_ROT'][obs_idx] = c*z_rot
    axis_dist['V_ROT_ERR'][obs_idx] = c*np.sqrt((axis_dist['ZERR'][obs_idx]/(1 + z_center))**2 + z_err_center2*((1 + axis_dist['Z'][obs_idx])/(1 + z_center)**2))
    #---------------------------------------------------------------------------
    
    
    #---------------------------------------------------------------------------
    # Correct rotational velocities for inclination angle
    #---------------------------------------------------------------------------
    cosi2 = (SGA['BA'][sga_idx]**2 - q0**2)/(1 - q0**2)
    
    # Galaxies with b/a < q0
    if cosi2 < 0:
        cosi2 = 0
    
    axis_dist['V_ROT'][obs_idx] /= np.sin(np.arccos(np.sqrt(cosi2)))
    #---------------------------------------------------------------------------
#-------------------------------------------------------------------------------
################################################################################




################################################################################
# Cut for galaxies suitable for calibrating the TFR
#
# Requirements:
#  - 10 < Vrot < 1000 km/s at 0.33R26
#  - Delta V / Vmin <= 5
#  - i > 45 degrees
#  - spiral-type morphology
#  - passes visual inspection
#-------------------------------------------------------------------------------
# Velocity cut
#-------------------------------------------------------------------------------
# Coma
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
r0p3 = (axis_inComa['SKY_FIBER_DIST_R26'] > 0.3) & (axis_inComa['SKY_FIBER_DIST_R26'] < 0.4)

Vgood = (np.abs(axis_inComa['V_ROT']) < 1000) & (np.abs(axis_inComa['V_ROT']) > 10)

good_axis_inComa = axis_inComa[r0p3 & Vgood]
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# 0-pt calibrators
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
r0p3_0pt = (axis_dist['SKY_FIBER_DIST_R26'] > 0.3) & (axis_dist['SKY_FIBER_DIST_R26'] < 0.4)

Vgood_0pt = (np.abs(axis_dist['V_ROT']) < 1000) & (np.abs(axis_dist['V_ROT']) > 10)

good_axis_dist = axis_dist[r0p3_0pt & Vgood_0pt]
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Relative velocity cut
#-------------------------------------------------------------------------------
# Coma
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
good_deltaV = np.ones(len(good_axis_inComa), dtype=bool)

for sga_id in np.unique(good_axis_inComa['SGA_ID']):
    
    # Identify all galaxy targets on this galaxy
    obs_idx = good_axis_inComa['SGA_ID'] == sga_id
    
    n_obs = np.sum(obs_idx)
    
    if n_obs > 1:
        
        Vmin = np.min(np.abs(good_axis_inComa['V_ROT'][obs_idx]))
        Vmax = np.max(np.abs(good_axis_inComa['V_ROT'][obs_idx]))
        
        v_norm_min = np.abs(good_axis_inComa['V_ROT'][obs_idx])/Vmin
        v_norm_max = np.abs(good_axis_inComa['V_ROT'][obs_idx])/Vmax
        
        diff_matrix = np.abs(good_axis_inComa['V_ROT'][obs_idx]).reshape(n_obs, 1) - np.abs(good_axis_inComa['V_ROT'][obs_idx]).reshape(1, n_obs)
        
        diff_matrix_norm = diff_matrix/Vmin
        
        if np.any(np.abs(diff_matrix_norm) > 5.):
            
            # Remove all observations with DELTACHI2 < 25
            # Note: This also typically removes observations with ZWARN != 0
            deltachi2_idx = good_axis_inComa['DELTACHI2'] >= 25
            
            good_deltaV[obs_idx & ~deltachi2_idx] = False
            
            good_obs_idx = obs_idx & deltachi2_idx
            
            n_obs_good = np.sum(good_obs_idx)
            
            # Check to make sure that, if there are still multiple observations, 
            # they all satisfy our relative velocity criteria
            if n_obs_good > 1:
                
                Vmin = np.min(np.abs(good_axis_inComa['V_ROT'][good_obs_idx]))
                
                diff_matrix = np.abs(good_axis_inComa['V_ROT'][good_obs_idx]).reshape(n_obs_good, 1) - np.abs(good_axis_inComa['V_ROT'][good_obs_idx]).reshape(1, n_obs_good)
                
                diff_matrix_norm = diff_matrix/Vmin
                
                if np.any(np.abs(diff_matrix_norm) > 5.):

                    # Set all of these so that we don't look at this galaxy
                    good_deltaV[good_obs_idx] = False
                    
good_deltaV_axis_inComa = good_axis_inComa[good_deltaV]
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# 0-pt calibrators
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
good_deltaV_0pt = np.ones(len(good_axis_dist), dtype=bool)

for sga_id in np.unique(good_axis_dist['SGA_ID']):
    
    # Identify all galaxy targets on this galaxy
    obs_idx = good_axis_dist['SGA_ID'] == sga_id
    
    n_obs = np.sum(obs_idx)
    
    if n_obs > 1:
        
        Vmin = np.min(np.abs(good_axis_dist['V_ROT'][obs_idx]))
        Vmax = np.max(np.abs(good_axis_dist['V_ROT'][obs_idx]))
        
        v_norm_min = np.abs(good_axis_dist['V_ROT'][obs_idx])/Vmin
        v_norm_max = np.abs(good_axis_dist['V_ROT'][obs_idx])/Vmax
        
        diff_matrix = np.abs(good_axis_dist['V_ROT'][obs_idx]).reshape(n_obs, 1) - np.abs(good_axis_dist['V_ROT'][obs_idx]).reshape(1, n_obs)
        
        diff_matrix_norm = diff_matrix/Vmin
        
        if np.any(np.abs(diff_matrix_norm) > 5.):
            
            # Remove all observations with DELTACHI2 < 25
            # Note: This also typically removes observations with ZWARN != 0
            deltachi2_idx = good_axis_dist['DELTACHI2'] >= 25
            
            good_deltaV_0pt[obs_idx & ~deltachi2_idx] = False
            
            good_obs_idx = obs_idx & deltachi2_idx
            
            n_obs_good = np.sum(good_obs_idx)
            
            # Check to make sure that, if there are still multiple observations, they all satisfy our relative velocity criteria
            if n_obs_good > 1:
                
                Vmin = np.min(np.abs(good_axis_dist['V_ROT'][good_obs_idx]))
                
                diff_matrix = np.abs(good_axis_dist['V_ROT'][good_obs_idx]).reshape(n_obs_good, 1) - np.abs(good_axis_dist['V_ROT'][good_obs_idx]).reshape(1, n_obs_good)
                
                diff_matrix_norm = diff_matrix/Vmin
                
                if np.any(np.abs(diff_matrix_norm) > 5.):
                    
                    # Set all of these so that we don't look at this galaxy
                    good_deltaV_0pt[good_obs_idx] = False

good_deltaV_axis_dist = good_axis_dist[good_deltaV_0pt]
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Inclination angle cut
#-------------------------------------------------------------------------------
i_min = 45. # degrees
cosi2_max = np.cos(i_min*np.pi/180.)**2

SGA['cosi2'] = (SGA['BA']**2 - q0**2)/(1 - q0**2)
SGA['cosi2'][SGA['cosi2'] < 0] = 0
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Coma
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
good_deltaV_axis_inComa['iSGA'] = -1

for i in range(len(good_deltaV_axis_inComa)):
    
    # Find galaxy in SGA
    sga_idx = SGA_dict[good_deltaV_axis_inComa['SGA_ID'][i]]
    
    good_deltaV_axis_inComa['iSGA'][i] = sga_idx
    
good_deltaV_axis_inComa['cosi2'] = SGA['cosi2'][good_deltaV_axis_inComa['iSGA']]

edge = good_deltaV_axis_inComa['cosi2'] <= cosi2_max

good_edge_axis_inComa = good_deltaV_axis_inComa[edge]
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# 0-pt calibrators
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
good_deltaV_axis_dist['iSGA'] = -1

for i in range(len(good_deltaV_axis_dist)):
    
    # Find galaxy in SGA
    sga_idx = SGA_dict[good_deltaV_axis_dist['SGA_ID'][i]]
    
    good_deltaV_axis_dist['iSGA'][i] = sga_idx
    
good_deltaV_axis_dist['cosi2'] = SGA['cosi2'][good_deltaV_axis_dist['iSGA']]

edge_0pt = good_deltaV_axis_dist['cosi2'] <= cosi2_max

good_edge_axis_dist = good_deltaV_axis_dist[edge_0pt]
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Morphology cut
#-------------------------------------------------------------------------------
# Coma
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
good_edge_axis_inComa['MORPHTYPE'] = SGA['MORPHTYPE'][good_edge_axis_inComa['iSGA']]

spirals = np.zeros(len(good_edge_axis_inComa), dtype=bool)

for i in range(len(good_edge_axis_inComa)):
    
    try:    
        if (good_edge_axis_inComa['MORPHTYPE'][i][0] == 'S') and (good_edge_axis_inComa['MORPHTYPE'][i][:2] != 'S0'):
            spirals[i] = True
    except IndexError:
        print(good_edge_axis_inComa['MORPHTYPE'][i])

good_edge_spirals_axis_inComa = good_edge_axis_inComa[spirals]
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# 0-pt calibrators
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
good_edge_axis_dist['MORPHTYPE'] = SGA['MORPHTYPE'][good_edge_axis_dist['iSGA']]

spirals_0pt = np.zeros(len(good_edge_axis_dist), dtype=bool)

for i in range(len(good_edge_axis_dist)):
    
    try:    
        if (good_edge_axis_dist['MORPHTYPE'][i][0] == 'S') and (good_edge_axis_dist['MORPHTYPE'][i][:2] != 'S0'):
            spirals_0pt[i] = True
        elif good_edge_axis_dist['MORPHTYPE'][i] == 'N/A':
            spirals_0pt[i] = True
    except IndexError:
        print(good_edge_axis_dist['MORPHTYPE'][i])

good_edge_spirals_axis_dist = good_edge_axis_dist[spirals_0pt]
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Visual inspection cut
#-------------------------------------------------------------------------------
VI_remove = Table.read('fuji_VI.txt', format='ascii.commented_header')

remove_targets = np.zeros(len(good_edge_spirals_axis_inComa), dtype=bool)

for targetid in VI_remove['TARGETID']:
    
    remove_targets = remove_targets & (good_edge_spirals_axis_inComa['TARGETID'] == targetid)
    
VI_good_edge_spirals_axis_inComa = good_edge_spirals_axis_inComa[~remove_targets]

# Both of the 0-pt calibrators pass VI.
#-------------------------------------------------------------------------------
################################################################################




################################################################################
# Compute the weighted average velocity for those galaxies with more than one 
# observation at 0.33R26
#-------------------------------------------------------------------------------
# Coma
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
SGA['V_0p33R26'] = np.nan
SGA['V_0p33R26_err'] = np.nan

weights = 1./(VI_good_edge_spirals_axis_inComa['V_ROT_ERR']**2)

for sga_id in np.unique(VI_good_edge_spirals_axis_inComa['SGA_ID']):
    
    # Identify all galaxy targets on this galaxy
    obs_idx = VI_good_edge_spirals_axis_inComa['SGA_ID'] == sga_id
    
    SGA['V_0p33R26'][SGA_dict[sga_id]] = np.average(np.abs(VI_good_edge_spirals_axis_inComa['V_ROT'][obs_idx]), 
                                                    weights=weights[obs_idx])

    SGA['V_0p33R26_err'][SGA_dict[sga_id]] = np.sqrt(1./np.sum(weights[obs_idx]))

# Make a catalog of just Coma galaxies with velocities
SGA_TF = SGA[np.isfinite(SGA['V_0p33R26']) & (SGA['R_MAG_SB26'] > 0)]
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# 0-pt calibrators
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
SGA['V_0p33R26'] = np.nan
SGA['V_0p33R26_err'] = np.nan

good_edge_spirals_axis_dist['R_MAG_SB26'] = np.nan
good_edge_spirals_axis_dist['R_MAG_SB26_ERR'] = np.nan

weights = 1./(good_edge_spirals_axis_dist['V_ROT_ERR']**2)

for sga_id in np.unique(good_edge_spirals_axis_dist['SGA_ID']):
    
    # Identify all galaxy targets on this galaxy
    obs_idx = good_edge_spirals_axis_dist['SGA_ID'] == sga_id
    
    SGA['V_0p33R26'][SGA_dict[sga_id]] = np.average(np.abs(good_edge_spirals_axis_dist['V_ROT'][obs_idx]), 
                                                    weights=weights[obs_idx])

    SGA['V_0p33R26_err'][SGA_dict[sga_id]] = np.sqrt(1./np.sum(weights[obs_idx]))

# Make a catalog of just 0-pt galaxies with velocities
SGA_0pt = SGA[np.isfinite(SGA['V_0p33R26']) & (SGA['R_MAG_SB26'] > 0)]
################################################################################




################################################################################
# Photometric corrections
#-------------------------------------------------------------------------------
# Survey offsets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
gals_directory = '/global/cfs/cdirs/desi/science/td/pv/tfgalaxies/SV/'
# gals_directory = '/Users/kdouglass/Documents/Research/data/DESI/SV/'
gals_filename = 'SGA-2020_fuji_Vrot_photsys.fits'
gals = Table.read(gals_directory + gals_filename)
gal_photsys = gals['SGA_ID', 'PHOTSYS']

# Add PHOTOSYS column to target tables
SGA_TF = join(SGA_TF, gal_photsys, join_type='left', keys='SGA_ID')
SGA_0pt = join(SGA_0pt, gal_photsys, join_type='left', keys='SGA_ID')

# Calculate corrections
SGA_sys_corr, SGA_sys_corr_err = BASS_corr(SGA_TF['PHOTSYS'])
zpt_sys_corr, zpt_sys_corr_err = BASS_corr(SGA_0pt['PHOTSYS'])
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# MW dust corrections
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Import E(B-V) dust map
ebv_directory = '/global/cfs/cdirs/desicollab/users/rongpu/dust/desi_ebv/public_data/maps/'
# ebv_directory = '/Users/kdouglass/Documents/Research/data/DESI/'
ebv_filename = 'desi_dust_gr_512.fits'
ebv_map = Table.read(ebv_directory + ebv_filename)

# Calculate corrections
SGA_dust_corr, SGA_dust_corr_err = MW_dust(SGA_TF['RA'], SGA_TF['DEC'], ebv_map)
zpt_dust_corr, zpt_dust_corr_err = MW_dust(SGA_0pt['RA'], SGA_0pt['DEC'], ebv_map)

# Flip NaN values to 0
SGA_dust_corr_err[np.isnan(SGA_dust_corr_err)] = 0
zpt_dust_corr_err[np.isnan(zpt_dust_corr_err)] = 0
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# K-corrections
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
SGA_kcorr = k_corr(SGA_TF['Z_DESI'], 
                   [SGA_TF['G_MAG_SB26'], SGA_TF['R_MAG_SB26'], SGA_TF['Z_MAG_SB26']], 
                   [SGA_TF['G_MAG_SB26_ERR'], SGA_TF['R_MAG_SB26_ERR'], SGA_TF['Z_MAG_SB26_ERR']])
zpt_kcorr = k_corr(SGA_0pt['Z_DESI'], 
                   [SGA_0pt['G_MAG_SB26'], SGA_0pt['R_MAG_SB26'], SGA_0pt['Z_MAG_SB26']], 
                   [SGA_0pt['G_MAG_SB26_ERR'], SGA_0pt['R_MAG_SB26_ERR'], SGA_0pt['Z_MAG_SB26_ERR']])
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Internal dust extinction
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
temp_infile = open('fuji_internalDust_mcmc-20250205.pickle', 'rb')
dust_mcmc_samples,_ = pickle.load(temp_infile)
temp_infile.close()

internalDust_coeffs = np.median(dust_mcmc_samples, axis=1)
internalDust_coeffs_err = np.zeros(len(internalDust_coeffs))
internalDust_coeffs_err[0] = np.std(dust_mcmc_samples[0][(-1.5 < dust_mcmc_samples[0]) & (dust_mcmc_samples[0] < 0)])
internalDust_coeffs_err[1] = np.std(dust_mcmc_samples[1][(0 < dust_mcmc_samples[1]) & (dust_mcmc_samples[1] < 1)])

SGA_internalDust_corr, SGA_internalDust_corr_err = internal_dust(SGA_TF['BA'], 
                                                                 internalDust_coeffs, 
                                                                 internalDust_coeffs_err)
zpt_internalDust_corr, zpt_internalDust_corr_err = internal_dust(SGA_0pt['BA'], 
                                                                 internalDust_coeffs, 
                                                                 internalDust_coeffs_err)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Apply corrections
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
SGA_TF['R_MAG_SB26_CORR'] = SGA_TF['R_MAG_SB26'] - SGA_dust_corr[1] + SGA_sys_corr + SGA_kcorr[:,1] - SGA_internalDust_corr
SGA_0pt['R_MAG_SB26_CORR'] = SGA_0pt['R_MAG_SB26'] - zpt_dust_corr[1] + zpt_sys_corr + zpt_kcorr[:,1] - zpt_internalDust_corr

SGA_TF['R_MAG_SB26_ERR_CORR'] = np.sqrt(SGA_TF['R_MAG_SB26_ERR']**2 + SGA_dust_corr_err[1]**2 + SGA_sys_corr_err**2 + SGA_internalDust_corr_err**2)
SGA_0pt['R_MAG_SB26_ERR_CORR'] = np.sqrt(SGA_0pt['R_MAG_SB26_ERR']**2 + zpt_dust_corr_err[1]**2 + zpt_sys_corr_err**2 + zpt_internalDust_corr_err**2)
################################################################################




################################################################################
# Compute the absolute magnitudes for the 0-pt calibrators based on their 
# distance moduli
#-------------------------------------------------------------------------------
SGA_0pt['R_ABSMAG_SB26'] = SGA_0pt['R_MAG_SB26_CORR'] - SGA_0pt['DM1_SN']
SGA_0pt['R_ABSMAG_SB26_err'] = np.sqrt(SGA_0pt['R_MAG_SB26_ERR_CORR']**2 + SGA_0pt['e_DM1_SN']**2)
################################################################################




# THERE ARE NO DWARF GALAXIES IN THIS CALIBRATION
################################################################################
# Identify the dwarf galaxies
#-------------------------------------------------------------------------------
# First, define the line perpendicular to the calibration
#-------------------------------------------------------------------------------
logV_n17 = (-17 - b_fit[1])/m_fit + V0
b_perp = -17 + w0_fit*(logV_n17 - V0)
#-------------------------------------------------------------------------------
'''
#-------------------------------------------------------------------------------
# Identify the dwarf galaxies
#-------------------------------------------------------------------------------
dwarfs = (SGA_TF['R_MAG_SB26_CORR'] - SGA_TF['R_MAG_SB26_ERR_CORR']) > (-w0_fit*(np.log10(SGA_TF['V_0p33R26']) - V0) + b_perp + b_fit[0] - b_fit[1])

SGA_TF_bright = SGA_TF[~dwarfs]
#-------------------------------------------------------------------------------
'''
################################################################################




################################################################################
# Plot the calibrated TFR for Coma
#-------------------------------------------------------------------------------
xvals = np.linspace(0.5, 3., 1000)
yvals = np.zeros((len(b_fit), len(xvals)))

for i in range(len(b_fit)):
    yvals[i] = m_fit * (xvals - V0) + b_fit[i]
    
yvals_perp = -w0_fit*(xvals - V0) + (b_perp + b_fit[0] - b_fit[1])


#-------------------------------------------------------------------------------
# Get the MCMC 1-sigma quantiles to plot with the fit
#-------------------------------------------------------------------------------
y_chain1 = np.outer(xvals - V0, tfr_samples[0]) + tfr_samples[1]
y_chain1_quantiles = np.quantile(y_chain1, [0.1587, 0.8414], axis=1)

y_chain2 = np.outer(xvals - V0, tfr_samples[0]) + tfr_samples[2]
y_chain2_quantiles = np.quantile(y_chain2, [0.1587, 0.8414], axis=1)
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Pack info into data
#-------------------------------------------------------------------------------
'''
data1 = [np.log10(SGA_TF_bright['V_0p33R26']), SGA_TF_bright['R_MAG_SB26_CORR']]
x1_err = 0.434*SGA_TF_bright['V_0p33R26_err']/SGA_TF_bright['V_0p33R26']
y1_err = SGA_TF_bright['R_MAG_SB26_ERR_CORR']
'''
data1 = [np.log10(SGA_TF['V_0p33R26']), SGA_TF['R_MAG_SB26_CORR']]
x1_err = 0.434*SGA_TF['V_0p33R26_err']/SGA_TF['V_0p33R26']
y1_err = SGA_TF['R_MAG_SB26_ERR_CORR']
corr1_xy = np.zeros_like(x1_err)

data2 = [np.log10(SGA_0pt['V_0p33R26']), SGA_0pt['R_ABSMAG_SB26']]
x2_err = 0.434*SGA_0pt['V_0p33R26_err']/SGA_0pt['V_0p33R26']
y2_err = SGA_0pt['R_ABSMAG_SB26_err']
corr2_xy = np.zeros_like(x2_err)
#-------------------------------------------------------------------------------

'''
#-------------------------------------------------------------------------------
# Dwarf galaxies
#-------------------------------------------------------------------------------
data_dwarfs = [np.log10(SGA_TF['V_0p33R26'][dwarfs]), 
               SGA_TF['R_MAG_SB26'][dwarfs]]
x_err_dwarfs = 0.434*SGA_TF['V_0p33R26_err'][dwarfs]/SGA_TF['V_0p33R26'][dwarfs]
y_err_dwarfs = SGA_TF['R_MAG_SB26_ERR'][dwarfs]
corr_xy_dwarfs = np.zeros_like(x_err_dwarfs)
#-------------------------------------------------------------------------------
'''

#-------------------------------------------------------------------------------
# Generate ellipses
#-------------------------------------------------------------------------------
ells1 = [
    Ellipse(
        xy=[data1[0][i], data1[1][i]],
        width=2*y1_err[i],
        height=2*x1_err[i],
        angle=np.rad2deg(np.arccos(corr1_xy[i])),
    )
    for i in range(len(data1[0]))
]

ells2 = [
    Ellipse(
        xy=[data2[0][i], data2[1][i]],
        width=2*y2_err[i],
        height=2*x2_err[i],
        angle=np.rad2deg(np.arccos(corr2_xy[i])),
    )
    for i in range(len(data2[0]))
]
'''
ells_dwarfs = [
    Ellipse(
        xy=[data_dwarfs[0][i], data_dwarfs[1][i]],
        width=2*y_err_dwarfs[i],
        height=2*x_err_dwarfs[i],
        angle=np.rad2deg(np.arccos(corr_xy_dwarfs[i])),
    )
    for i in range(len(data_dwarfs[0]))
]
'''
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Make the plot
#-------------------------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(2, 
                               figsize=(4,4.8), 
                               tight_layout=True, 
                               sharex=True, 
                               height_ratios=[3, 5])
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Coma
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Dwarf region
ax2.fill_between(xvals, 
                 yvals_perp, 
                 18, 
                 color='gainsboro')

# Perpendicular line (used to define dwarfs)
# plt.plot(xvals, yvals_perp, ':', c='darkgray', lw=1.3)

# Calibrated TFR uncertainty
ax2.fill_between(xvals, 
                 y_chain1_quantiles[0], 
                 y_chain1_quantiles[1], 
                 color="darkgray")

# Calibrated TFR w/ Coma intercept
ax2.plot(xvals, yvals[0], "k", lw=1.3)
ax2.plot(xvals, yvals[0] - sig_fit[0], "k--", lw=1.3)
ax2.plot(xvals, yvals[0] + sig_fit[0], "k--", lw=1.3)

# Uncertainty ellipses for Coma galaxies
for i, e in enumerate(ells1):
    ax2.add_artist(e)
    e.set_color(adjust_lightness('tab:blue', amount=1.75))
'''
# Uncertainty ellipses for dwarf Coma galaxies
for i, e in enumerate(ells_dwarfs):
    ax2.add_artist(e)
    e.set_color(adjust_lightness('gray', amount=1.75))
'''
# Coma galaxies
ax2.plot(data1[0], data1[1], 'x', label='Coma ({} galaxies)'.format(len(data1[0])))
'''
# Dwarf galaxies in Coma
ax2.plot(data_dwarfs[0], data_dwarfs[1], 'x', c='gray')
'''

ax2.set_xlabel(r"$\log{(V(0.33R_{26}) [\mathrm{km/s}])}$", fontsize=14)
ax2.set_ylabel(r"$m_r^{0.05}(26)$", fontsize=14)

ax2.legend(loc='upper left')

ax2.set_xlim(0.5, 3)
ax2.set_ylim(17.5, 12)
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 0-pt calibrators
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Calibrated TFR uncertainty
ax1.fill_between(xvals, 
                 y_chain2_quantiles[0], 
                 y_chain2_quantiles[1], 
                 color='darkgray')

# Calibrated TFR
ax1.plot(xvals, yvals[1], 'k', lw=1.3)

for i, e in enumerate(ells2):
    ax1.add_artist(e)
    e.set_facecolor(adjust_lightness('tab:orange', amount=1.25))
    
ax1.plot(data2[0], data2[1], 'x', c='tab:orange', 
         label='0-pt ({} galaxies)'.format(len(data2[0])))

ax1.set_ylabel(r'$M_r^{0.05} (26) - 5\log h$', fontsize=14)

ax1.legend(loc='upper left')

ax1.set_ylim(-19.5, -23)
# ax1.set_aspect('equal', adjustable='box', share=True)
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Save the figure
#-------------------------------------------------------------------------------
plt.savefig('../../Figures/SV/fuji_joint-Coma-0pt_TFR_varyV0-perpdwarfs_AnthonyUpdates_weightsVmax-1_20250522.png', 
            dpi=150, 
            facecolor='none')
#-------------------------------------------------------------------------------
################################################################################