'''
This file is the same as Coma-Fuji-ITFR-nodwarfs_KAD.ipynb, but should only be 
used to generate the final plot; it does not contain any of the fitting 
procedures.
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
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u
# from astropy.visualization.wcsaxes import SphericalCircle

# import corner

# import os

# import requests

import pickle

from help_functions import adjust_lightness

# import sys
# sys.path.insert(1, '/global/u1/k/kadglass/DESI_SGA/TF/')
# #sys.path.insert(1, '/Users/kellydouglass/Documents/Research/DESI/Targets/code/TF/')
# from line_fits import param_invert, hyperfit_line
################################################################################




################################################################################
# Constants
#-------------------------------------------------------------------------------
h = 1
H0 = 100*h

c = 3e5

q0 = 0.2

V0 = 2.5
################################################################################




################################################################################
# Import data
#-------------------------------------------------------------------------------
# fuji
#-------------------------------------------------------------------------------
tfuji = Table.read('/global/cfs/projectdirs/desi/science/td/pv/desi_pv_tf_fuji_healpix.fits')
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
sigma_Coma = table3['sigP'][Coma_row_t3][0]
mu_Coma = table3['DM'][Coma_row_t3][0]


Coma_coords = SkyCoord(table3['SGLON'][Coma_row_t3]*u.degree, 
                       table3['SGLAT'][Coma_row_t3]*u.degree, 
                       frame='supergalactic')

d_Coma = 10*10**(0.2*mu_Coma) # pc
V_Coma = 100*(d_Coma*1e-6)    # km/s


#-------------------------------------------------------------------------------
# Calculate the projected distance between the Coma cluster and each SGA galaxy
#-------------------------------------------------------------------------------
# First, we need to convert R2t from Mpc to an angle, using the group's velocity
# Note that we are NOT assuming that the size of the cluster is a small angle!!
R2t_Coma_angle_1p5 = np.arctan(1.5*R2t_Coma/(d_Coma*1e-6))*u.radian
R2t_Coma_angle_3 = np.arctan(3*R2t_Coma/(d_Coma*1e-6))*u.radian

SGA_coords = SkyCoord(SGA['RA'], SGA['DEC'], unit='deg')

sep = Coma_coords.separation(SGA_coords)

SGA_in_Coma1 = (sep < R2t_Coma_angle_1p5) & (SGA['Z_DESI']*c > V_Coma - 3*sigma_Coma) & (SGA['Z_DESI']*c < V_Coma + 3*sigma_Coma)

SGA_in_Coma2 = (sep >= R2t_Coma_angle_1p5) & (sep < R2t_Coma_angle_3) & (SGA['Z_DESI']*c > V_Coma - 2*sigma_Coma) & (SGA['Z_DESI']*c < V_Coma + 2*sigma_Coma)

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
# Calculate the rotational velocity
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
    axis_inComa['V_ROT'][obs_idx] = c*(axis_inComa['Z'][obs_idx] - z_center)
    axis_inComa['V_ROT_ERR'][obs_idx] = c*np.sqrt(axis_inComa['ZERR'][obs_idx]**2 + z_err_center2)
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
################################################################################




################################################################################
# Cut for Coma galaxies suitable for calibrating the TFR
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
r0p3 = (axis_inComa['SKY_FIBER_DIST_R26'] > 0.3) & (axis_inComa['SKY_FIBER_DIST_R26'] < 0.4)

Vgood = (np.abs(axis_inComa['V_ROT']) < 1000) & (np.abs(axis_inComa['V_ROT']) > 10)

good_axis_inComa = axis_inComa[r0p3 & Vgood]
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Relative velocity cut
#-------------------------------------------------------------------------------
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
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Inclination angle cut
#-------------------------------------------------------------------------------
SGA['cosi2'] = (SGA['BA']**2 - q0**2)/(1 - q0**2)
SGA['cosi2'][SGA['cosi2'] < 0] = 0

good_deltaV_axis_inComa['iSGA'] = -1

for i in range(len(good_deltaV_axis_inComa)):
    
    # Find galaxy in SGA
    sga_idx = SGA_dict[good_deltaV_axis_inComa['SGA_ID'][i]]
    
    good_deltaV_axis_inComa['iSGA'][i] = sga_idx
    
good_deltaV_axis_inComa['cosi2'] = SGA['cosi2'][good_deltaV_axis_inComa['iSGA']]

i_min = 45. # degrees

cosi2_max = np.cos(i_min*np.pi/180.)**2

edge = good_deltaV_axis_inComa['cosi2'] <= cosi2_max

good_edge_axis_inComa = good_deltaV_axis_inComa[edge]
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Morphology cut
#-------------------------------------------------------------------------------
good_edge_axis_inComa['MORPHTYPE'] = SGA['MORPHTYPE'][good_edge_axis_inComa['iSGA']]

spirals = np.zeros(len(good_edge_axis_inComa), dtype=bool)

for i in range(len(good_edge_axis_inComa)):
    
    try:    
        if (good_edge_axis_inComa['MORPHTYPE'][i][0] == 'S') and (good_edge_axis_inComa['MORPHTYPE'][i][:2] != 'S0'):
            spirals[i] = True
    except IndexError:
        print(good_edge_axis_inComa['MORPHTYPE'][i])

good_edge_spirals_axis_inComa = good_edge_axis_inComa[spirals]
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Visual inspection cut
#-------------------------------------------------------------------------------
VI_remove = Table.read('fuji_VI.txt', format='ascii.commented_header')

remove_targets = np.zeros(len(good_edge_spirals_axis_inComa), dtype=bool)

for targetid in VI_remove['TARGETID']:
    
    remove_targets = remove_targets & (good_edge_spirals_axis_inComa['TARGETID'] == targetid)
    
VI_good_edge_spirals_axis_inComa = good_edge_spirals_axis_inComa[~remove_targets]
#-------------------------------------------------------------------------------
################################################################################




################################################################################
# Compute the weighted average velocity for those galaxies with more than one 
# observation at 0.33R26
#-------------------------------------------------------------------------------
SGA['V_0p33R26'] = np.nan
SGA['V_0p33R26_err'] = np.nan

weights = 1./(VI_good_edge_spirals_axis_inComa['V_ROT_ERR']**2)

for sga_id in np.unique(VI_good_edge_spirals_axis_inComa['SGA_ID']):
    
    # Identify all galaxy targets on this galaxy
    obs_idx = VI_good_edge_spirals_axis_inComa['SGA_ID'] == sga_id
    
    SGA['V_0p33R26'][SGA_dict[sga_id]] = np.average(np.abs(VI_good_edge_spirals_axis_inComa['V_ROT'][obs_idx]), 
                                                    weights=weights[obs_idx])

    SGA['V_0p33R26_err'][SGA_dict[sga_id]] = np.sqrt(1./np.sum(weights[obs_idx]))
################################################################################




################################################################################
# Make a catalog of just those galaxies with velocities
#-------------------------------------------------------------------------------
SGA_TF = SGA[np.isfinite(SGA['V_0p33R26']) & (SGA['R_MAG_SB26'] > 0)]
################################################################################




################################################################################
# Identify the dwarf galaxies
#-------------------------------------------------------------------------------
temp_infile = open('cov_ab_fuji_zero-point_ITFR_nodwarf1_KAD.pickle', 'rb')
cov_ab, itfr_samples, tfr_samples = pickle.load(temp_infile)
temp_infile.close()

#-------------------------------------------------------------------------------
# First, calculate the absolute magnitudes for the galaxies based on the TF 
# calibration
#-------------------------------------------------------------------------------
m_fit = np.median(tfr_samples[0])
b_fit = np.median(tfr_samples[1])

SGA_TF['R_ABSMAG_SB26'] = m_fit*(np.log10(SGA_TF['V_0p33R26']) - V0) + b_fit
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Calculate the uncertainty in the absolute magnitudes based on the TF 
# calibration
#-------------------------------------------------------------------------------
rng = np.random.default_rng()

m_random = tfr_samples[0]
b_random = tfr_samples[1]

N_samples = len(m_random)

SGA_TF['R_ABSMAG_SB26_err'] = np.nan

for i in range(len(SGA_TF)):
    
    v_random = rng.normal(SGA_TF['V_0p33R26'][i], 
                          SGA_TF['V_0p33R26_err'][i], 
                          size=N_samples)
    
    Ms = m_random*(np.log10(v_random) - V0) + b_random
    
    SGA_TF['R_ABSMAG_SB26_err'][i] = np.nanstd(Ms)
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Identify the dwarf galaxies
#-------------------------------------------------------------------------------
dwarfs = (SGA_TF['R_ABSMAG_SB26'] - SGA_TF['R_ABSMAG_SB26_err']) > -17

SGA_TF_bright = SGA_TF[~dwarfs]
#-------------------------------------------------------------------------------
################################################################################




################################################################################
# Plot the calibrated TFR for Coma
#-------------------------------------------------------------------------------
temp_infile = open('mcmc_fuji_Coma_ITFR_nodwarfs1_KAD.pickle', 'rb')
(cov_w, cov_itfr, itfr_mcmc_samples, tfr_mcmc_samples) = pickle.load(temp_infile)
temp_infile.close()


b_Coma = np.median(tfr_mcmc_samples[1])


xvals = np.linspace(1., 3., 1000)
yvals = m_fit * (xvals - V0) + b_Coma


#-------------------------------------------------------------------------------
# Get the MCMC 1-sigma quantiles to plot with the fit
#-------------------------------------------------------------------------------
y_chain = np.outer(xvals - V0, tfr_mcmc_samples[0]) + tfr_mcmc_samples[1]
y_chain_quantiles = np.quantile(y_chain, [0.1587, 0.8414], axis=1)
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Pack info into data
#-------------------------------------------------------------------------------
data = [np.log10(SGA_TF_bright['V_0p33R26']), SGA_TF_bright['R_MAG_SB26']]
x_err = 0.434*SGA_TF_bright['V_0p33R26_err']/SGA_TF_bright['V_0p33R26']
y_err = SGA_TF_bright['R_MAG_SB26_ERR']
corr_xy = np.zeros_like(x_err)
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Dwarf galaxies
#-------------------------------------------------------------------------------
data_dwarfs = [np.log10(SGA_TF['V_0p33R26'][dwarfs]), 
               SGA_TF['R_MAG_SB26'][dwarfs]]
x_err_dwarfs = 0.434*SGA_TF['V_0p33R26_err'][dwarfs]/SGA_TF['V_0p33R26'][dwarfs]
y_err_dwarfs = SGA_TF['R_MAG_SB26_ERR'][dwarfs]
corr_xy_dwarfs = np.zeros_like(x_err_dwarfs)
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Generate ellipses
#-------------------------------------------------------------------------------
ells = [
    Ellipse(
        xy=[data[0][i], data[1][i]],
        width=2*y_err[i],
        height=2*x_err[i],
        angle=np.rad2deg(np.arccos(corr_xy[i])),
    )
    for i in range(len(data[0]))
]

ells_dwarfs = [
    Ellipse(
        xy=[data_dwarfs[0][i], data_dwarfs[1][i]],
        width=2*y_err_dwarfs[i],
        height=2*x_err_dwarfs[i],
        angle=np.rad2deg(np.arccos(corr_xy_dwarfs[i])),
    )
    for i in range(len(data_dwarfs[0]))
]
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Make the plot
#-------------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(4,5), tight_layout=True)

ax.fill_between(xvals, 
                y_chain_quantiles[0], 
                y_chain_quantiles[1], 
                color="lightgray")

for i, e in enumerate(ells):
    ax.add_artist(e)
    e.set_color(adjust_lightness('tab:blue', amount=1.75))
    
for i, e in enumerate(ells_dwarfs):
    ax.add_artist(e)
    e.set_color(adjust_lightness('gray', amount=1.75))
    
ax.plot(data[0], data[1], 'x')

ax.plot(data_dwarfs[0], data_dwarfs[1], 'x', c='gray')

ax.plot(xvals, yvals, c="k", marker="None", ls="-", lw=1.3)#, alpha=0.9)

'''
ax.plot(xvals - hf.vert_scat, 
        yvals, 
        c="k", 
        marker="None", 
        ls="--", 
        lw=1.3)#, alpha=0.9)
ax.plot(xvals + hf.vert_scat, 
        yvals, 
        c="k", 
        marker="None", 
        ls="--", 
        lw=1.3)#, alpha=0.9)
'''

ax.set_xlabel(r"$\log{(V_\mathrm{0.33R_{26}} [\mathrm{km/s}])}$", fontsize=14)
ax.set_ylabel(r"$m_r(26)$", fontsize=14)

ax.set_title("Coma Cluster ({} galaxies)".format(len(SGA_TF_bright)), 
             fontsize = 14)

ax.set_xlim(0.5, 3)
ax.set_ylim(18, 13)
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Save the figure
#-------------------------------------------------------------------------------
plt.savefig('../../Figures/SV/fuji_Coma_TFR_nodwarfs1_noquantiles_20230818.png', 
            dpi=150, 
            facecolor='none')
#-------------------------------------------------------------------------------
################################################################################