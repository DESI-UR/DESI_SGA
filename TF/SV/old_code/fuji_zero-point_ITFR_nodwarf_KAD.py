'''
This file is a replica of the notebook with the same name.  However, it does not 
include any of the fitting - this file is to be used when adjusting any of the 
figure parameters / styling.
'''

################################################################################
# Import modules
#-------------------------------------------------------------------------------
import numpy as np

from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u

import corner

# import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

import pickle

from help_functions import adjust_lightness

# import sys
# sys.path.insert(1, '/global/u1/k/kadglass/DESI_SGA/TF/')
# from line_fits import param_invert
################################################################################



'''
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
# Read in our version of the SGA that includes distances from the Extragalactic
# Distance Database. (This file was made with the data_match.ipynb notebook.)
#-------------------------------------------------------------------------------
SGA = Table.read('../SGA_distances.fits')

SGA_dict = {}

for i in range(len(SGA)):
    SGA_dict[SGA['SGA_ID'][i]] = i
#-------------------------------------------------------------------------------
################################################################################




################################################################################
# Separate fuji data into center and off-center observations
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
#  - ZWARN = 0
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
# Filter SGA to keep only those objects with center observations and 
# independent distances
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
# Calculate the rotational velocities
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
    axis_dist['V_ROT'][obs_idx] = c*(axis_dist['Z'][obs_idx] - z_center)
    axis_dist['V_ROT_ERR'][obs_idx] = c*np.sqrt(axis_dist['ZERR'][obs_idx]**2 + z_err_center2)
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
################################################################################




################################################################################
# Cut for galaxies with distances suitable for calibrating the TFR
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
r0p3 = (axis_dist['SKY_FIBER_DIST_R26'] > 0.3) & (axis_dist['SKY_FIBER_DIST_R26'] < 0.4)

Vgood = (np.abs(axis_dist['V_ROT']) < 1000) & (np.abs(axis_dist['V_ROT']) > 10)

good_axis_dist = axis_dist[r0p3 & Vgood]
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Relative velocity cut
#-------------------------------------------------------------------------------
good_deltaV = np.ones(len(good_axis_dist), dtype=bool)

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
            
            good_deltaV[obs_idx & ~deltachi2_idx] = False
            
            good_obs_idx = obs_idx & deltachi2_idx
            
            n_obs_good = np.sum(good_obs_idx)
            
            # Check to make sure that, if there are still multiple observations, 
            # they all satisfy our relative velocity criteria
            if n_obs_good > 1:
                
                Vmin = np.min(np.abs(good_axis_dist['V_ROT'][good_obs_idx]))
                
                diff_matrix = np.abs(good_axis_dist['V_ROT'][good_obs_idx]).reshape(n_obs_good, 1) - np.abs(good_axis_dist['V_ROT'][good_obs_idx]).reshape(1, n_obs_good)
                
                diff_matrix_norm = diff_matrix/Vmin
                
                if np.any(np.abs(diff_matrix_norm) > 5.):
                    
                    # Set all of these so that we don't look at this galaxy
                    good_deltaV[good_obs_idx] = False
                    
good_deltaV_axis_dist = good_axis_dist[good_deltaV]
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Inclination angle cut
#-------------------------------------------------------------------------------
SGA['cosi2'] = (SGA['BA']**2 - q0**2)/(1 - q0**2)
SGA['cosi2'][SGA['cosi2'] < 0] = 0

good_deltaV_axis_dist['iSGA'] = -1

for i in range(len(good_deltaV_axis_dist)):
    
    # Find galaxy in SGA
    sga_idx = SGA_dict[good_deltaV_axis_dist['SGA_ID'][i]]
    
    good_deltaV_axis_dist['iSGA'][i] = sga_idx
    
good_deltaV_axis_dist['cosi2'] = SGA['cosi2'][good_deltaV_axis_dist['iSGA']]

i_min = 45. # degrees

cosi2_max = np.cos(i_min*np.pi/180.)**2

edge = good_deltaV_axis_dist['cosi2'] <= cosi2_max

good_edge_axis_dist = good_deltaV_axis_dist[edge]
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Morphology cut
#-------------------------------------------------------------------------------
good_edge_axis_dist['MORPHTYPE'] = SGA['MORPHTYPE'][good_edge_axis_dist['iSGA']]

spirals = np.zeros(len(good_edge_axis_dist), dtype=bool)

for i in range(len(good_edge_axis_dist)):
    
    try:    
        if (good_edge_axis_dist['MORPHTYPE'][i][0] == 'S') and (good_edge_axis_dist['MORPHTYPE'][i][:2] != 'S0'):
            spirals[i] = True
    except IndexError:
        print(good_edge_axis_dist['MORPHTYPE'][i])

good_edge_spirals_axis_dist = good_edge_axis_dist[spirals]
#-------------------------------------------------------------------------------

SGA_idx = []

for SGA_id in np.unique(good_edge_spirals_axis_dist['SGA_ID']):
    
    SGA_idx.append(SGA_dict[SGA_id])

#-------------------------------------------------------------------------------
# Visual inspection
#
# No targets need to be removed - both objects and their observations pass VI.
#-------------------------------------------------------------------------------
################################################################################




################################################################################
# Compute the weighted average velocity for those galaxies with more than one
# observation at 0.33R26
#-------------------------------------------------------------------------------
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

# Make a catalog of just those galaxies with velocities
SGA_0pt = SGA[np.isfinite(SGA['V_0p33R26']) & (SGA['R_MAG_SB26'] > 0)]
################################################################################




################################################################################
# Compute the absolute magnitudes based on the distance measurements
#
# Mr - 5logh = mr - mu - 5logh
#
# where h is the reduced Hubble constant used to calibrate the distance modulus, 
# mu.
#
# Both of our galaxies have distance moduli from Stahl et al. (2021), which 
# appears to use h = 1.  Therefore, the distance moduli that we have to use are 
# mu_h = mu - 5logh.
#-------------------------------------------------------------------------------
SGA_0pt['R_ABSMAG_SB26'] = SGA_0pt['R_MAG_SB26'] - SGA_0pt['DM1_SN']
SGA_0pt['R_ABSMAG_SB26_err'] = np.sqrt(SGA_0pt['R_MAG_SB26_ERR']**2 + SGA_0pt['e_DM1_SN']**2)
################################################################################




################################################################################
# Read in the calibrated slope and 0-pt
#-------------------------------------------------------------------------------
# Slope info
#-------------------------------------------------------------------------------
temp_infile = open('mcmc_fuji_Coma_ITFR_nodwarfs1_KAD.pickle', 'rb')
(cov_w, cov_itfr, itfr_mcmc_samples, tfr_mcmc_samples) = pickle.load(temp_infile)
temp_infile.close()
#-------------------------------------------------------------------------------
'''

#-------------------------------------------------------------------------------
# 0-pt info
#-------------------------------------------------------------------------------
temp_infile = open('cov_ab_fuji_zero-point_ITFR_nodwarf1_KAD.pickle', 'rb')
cov_ab, itfr_samples, tfr_samples = pickle.load(temp_infile)
temp_infile.close()
#-------------------------------------------------------------------------------
################################################################################



'''
################################################################################
# Calculate the calibrated slope, 0-pt, and uncertainties
#-------------------------------------------------------------------------------
m_fit = np.median(tfr_samples[0])
b_fit = np.median(tfr_samples[1])

w0_random = itfr_samples[0]
w1_fit_array = itfr_samples[1]

xvals = np.linspace(1, 3, 1000)
yvals = m_fit * (xvals - V0) + b_fit

N_samples = len(w0_random)

yvals_random = np.zeros((N_samples, len(xvals)))

for i in range(N_samples):
    m_random_i, b_fit_i,_ = param_invert(w0_random[i], 
                                         w1_fit_array[i], 
                                         np.cov(itfr_samples))
    yvals_random[i] = (xvals - V0) * m_random_i + b_fit_i
    
edges = np.nanpercentile(yvals_random, [16, 84], axis=0)
################################################################################
'''



################################################################################
# Corner plot
#-------------------------------------------------------------------------------
fig = corner.corner(tfr_samples.T, bins=30, smooth=1,
                    range=[[-13, -4], [-23, -21.75]],   # Range for a, b. Adjust as needed.
                    labels=['$a$', '$b$'],
                    levels=(1-np.exp(-0.5), 1-np.exp(-2)),
                    quantiles=[0.16, 0.5, 0.84],
                    color='tab:blue',
                    hist_kwargs={'histtype':'stepfilled', 'alpha':0.3},
                    plot_datapoints=False,
                    fill_contours=True,
                    show_titles=True,
                    title = {'0-pt calibration'},
                    title_kwargs={"fontsize": 14})

plt.savefig('../../Figures/SV/fuji_0pt_corner_nodwarfs1_20240130.png', 
            dpi=150, 
            facecolor='none')
################################################################################



'''
################################################################################
# Scatter plot
#-------------------------------------------------------------------------------
# Pack info into data
#-------------------------------------------------------------------------------
data = [np.log10(SGA_0pt['V_0p33R26']), SGA_0pt['R_ABSMAG_SB26']]
x_err = 0.434*SGA_0pt['V_0p33R26_err']/SGA_0pt['V_0p33R26']
y_err = SGA_0pt['R_ABSMAG_SB26_err']
corr_xy = np.zeros_like(x_err)
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
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Make the plot
#-------------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(4,5), tight_layout=True)

ax.fill_between(xvals, edges[0], edges[1], color="lightgray")

for i, e in enumerate(ells):
    ax.add_artist(e)
    e.set_color(adjust_lightness('tab:blue', amount=1.75))
    
ax.plot(data[0], data[1], 'x')

ax.plot(xvals, yvals, c="k", marker="None", ls="-", lw=1.3)

ax.set_xlabel(r"$\log{(V_\mathrm{0.33R_{26}} [\mathrm{km/s}])}$", fontsize=14)
ax.set_ylabel(r"$M_r(26) - 5$log$h$", fontsize=14)

ax.set_title("0-pt calibration ({} galaxies)".format(len(SGA_0pt)), 
             fontsize = 14)

ax.set_xlim(2, 2.5)
ax.set_ylim(-19, -22)
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Save the figure
#-------------------------------------------------------------------------------
plt.savefig('../../Figures/SV/fuji_0pt_TFR_nodwarf1_20230821.png', 
            dpi=150, 
            facecolor='none')
#-------------------------------------------------------------------------------
################################################################################
'''