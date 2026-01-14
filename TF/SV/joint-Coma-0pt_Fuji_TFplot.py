'''
This file is the same as joint-Coma-0pt_Fuji-TFR_KAD_varyV0-perpdwarfs.ipynb, 
but should only be used to generate the final Coma + 0pt calibrated TFR plot; it 
does not contain any of the fitting procedures.
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
from astropy.table import Table, join, unique
from astropy.coordinates import SkyCoord
from astropy import units as u
import astropy.constants as const
# from astropy.visualization.wcsaxes import SphericalCircle

import corner

# import os

# import requests

import pickle

import sys
# sys.path.insert(1, '/global/u1/k/kadglass/DESI_SGA/TF/')
sys.path.insert(1, '/Users/kdouglass/Documents/Research/DESI/PV_Survey/DESI_SGA/TF/')
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

# For error propagation
rng = np.random.default_rng()
N_samples = 10000
################################################################################




################################################################################
# Import data
#-------------------------------------------------------------------------------
# Fuji
#-------------------------------------------------------------------------------
tfuji = Table.read('SGA-2020_fuji_Vrot_VI_dVsys_photsys_corr.fits')
#-------------------------------------------------------------------------------
# SGA distances
#-------------------------------------------------------------------------------
SGA_dist = Table.read('../SGA_distances.fits')
#-------------------------------------------------------------------------------
# Merge the two tables
#-------------------------------------------------------------------------------
SGA_dist_unique_cols = list(set(SGA_dist.colnames).difference(tfuji.colnames))

SGA = join(tfuji, SGA_dist[['SGA_ID'] + SGA_dist_unique_cols], keys='SGA_ID')

SGA_dict = {}
for i in range(len(SGA)):
    SGA_dict[SGA['SGA_ID'][i]] = i
#-------------------------------------------------------------------------------
# Best-fit
#-------------------------------------------------------------------------------
temp_infile = open('cov_ab_fuji_joint_TFR_varyV0-perpdwarfs0_AnthonyUpdates_weightsVmax-1_dVsys_corr_KAD-20260114.pickle', 
                   'rb')
cov_ab, tfr_samples, V0 = pickle.load(temp_infile)
temp_infile.close()

m_fit = np.median(tfr_samples[0])
b_fit = np.median(tfr_samples[1:3], axis=1)
sig_fit = np.median(tfr_samples[3:], axis=1)
#-------------------------------------------------------------------------------
################################################################################




###############################################################################
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
# Keep all observations of each galaxy that are within the Coma cluster
#-------------------------------------------------------------------------------
SGA_ID_in_Coma = SGA['SGA_ID'][SGA_in_Coma]
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
# Only use Stahl distances
#-------------------------------------------------------------------------------
Stahl = SGA['SN_Catalog'] == 'Stahl-SNIa'
#-------------------------------------------------------------------------------
# Keep all observations of each galaxy that have independent distances
#-------------------------------------------------------------------------------
SGA_ID_dist = SGA['SGA_ID'][distances & centers & Stahl]
#-------------------------------------------------------------------------------
################################################################################




################################################################################
# Cut for galaxies suitable for calibrating the TFR
#
# Requirements:
#  - 10 < Vrot < 1000 km/s at 0.33R26
#  - Delta V / Vmin <= 5
#  - Velocity sign consistency
#  - i > 45 degrees
#  - spiral-type morphology
#  - passes visual inspection
#-------------------------------------------------------------------------------
# Inclination angle cut
#-------------------------------------------------------------------------------
i_min = 45. # degrees
cosi2_max = np.cos(i_min*np.pi/180.)**2

SGA['cosi2'] = (SGA['BA']**2 - q0**2)/(1 - q0**2)
SGA['cosi2'][SGA['cosi2'] < 0] = 0

edge = SGA['cosi2'] <= cosi2_max

# Coma
edge_inComa = SGA[SGA_in_Coma & edge]

# 0-pt calibrators
edge_dist = SGA[distances & centers & Stahl & edge]
#-------------------------------------------------------------------------------
# Morphology cut
#-------------------------------------------------------------------------------
# Coma
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
spirals = np.zeros(len(edge_inComa), dtype=bool)

for i in range(len(edge_inComa)):
    
    try:    
        if (edge_inComa['MORPHTYPE'][i][0] == 'S') and (edge_inComa['MORPHTYPE'][i][:2] != 'S0'):
            spirals[i] = True
    except IndexError:
        print(edge_inComa['MORPHTYPE'][i])

edge_spirals_inComa = edge_inComa[spirals]
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# 0-pt calibrators
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
spirals_0pt = np.zeros(len(edge_dist), dtype=bool)

for i in range(len(edge_dist)):
    
    try:    
        if (edge_dist['MORPHTYPE'][i][0] == 'S') and (edge_dist['MORPHTYPE'][i][:2] != 'S0'):
            spirals_0pt[i] = True
        elif edge_dist['MORPHTYPE'][i] == 'N/A':
            spirals_0pt[i] = True
    except IndexError:
        print(edge_dist['MORPHTYPE'][i])

edge_spirals_dist = edge_dist[spirals_0pt]
#-------------------------------------------------------------------------------

# Make a catalog of just Coma galaxies with velocities
SGA_TF = edge_spirals_inComa[np.isfinite(edge_spirals_inComa['V_0p33R26']) & (edge_spirals_inComa['R_MAG_SB26'] > 0)]

# Make a catalog of just 0-pt galaxies with velocities
SGA_0pt = edge_spirals_dist[np.isfinite(edge_spirals_dist['V_0p33R26']) & (edge_spirals_dist['R_MAG_SB26'] > 0)]
################################################################################




################################################################################
# Photometric corrections
#-------------------------------------------------------------------------------
# Survey offsets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Calculate corrections
SGA_sys_corr, SGA_sys_corr_err = BASS_corr(SGA_TF['PHOTSYS'])
zpt_sys_corr, zpt_sys_corr_err = BASS_corr(SGA_0pt['PHOTSYS'])
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# MW dust corrections
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Import E(B-V) dust map
# ebv_directory = '/global/cfs/cdirs/desicollab/users/rongpu/dust/desi_ebv/public_data/maps/'
ebv_directory = '/Users/kdouglass/Documents/Research/data/DESI/'
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
SGA_0pt['R_ABSMAG_SB26_ERR'] = np.sqrt(SGA_0pt['R_MAG_SB26_ERR_CORR']**2 + SGA_0pt['e_DM1_SN']**2)
################################################################################




################################################################################
# Save figure data for paper
#-------------------------------------------------------------------------------
# Header
#-------------------------------------------------------------------------------
hdr = fits.Header()

hdr['DESI_DR'] = 'EDR'
hdr['FIGURE'] = 7

empty_primary = fits.PrimaryHDU(header=hdr)
#-------------------------------------------------------------------------------
# Data
#-------------------------------------------------------------------------------
coma_hdu = fits.BinTableHDU(data=SGA_TF['R_MAG_SB26_CORR', 'R_MAG_SB26_ERR_CORR', 'V_0p33R26', 'V_0p33R26_ERR'], 
                            name='CAL')
coma_hdu.columns['R_MAG_SB26_CORR'].name = 'R_MAG'
coma_hdu.columns['R_MAG_SB26_ERR_CORR'].name = 'R_MAG_ERR'
coma_hdu.columns['V_0p33R26'].name = 'VROT'
coma_hdu.columns['V_0p33R26_ERR'].name = 'VROT_ERR'

zpt_hdu = fits.BinTableHDU(data=SGA_0pt['R_ABSMAG_SB26', 'R_ABSMAG_SB26_ERR', 'V_0p33R26', 'V_0p33R26_ERR'],
                           name='0PT')
zpt_hdu.columns['R_ABSMAG_SB26'].name = 'R_ABSMAG'
zpt_hdu.columns['R_ABSMAG_SB26_ERR'].name = 'R_ABSMAG_ERR'
zpt_hdu.columns['V_0p33R26'].name = 'VROT'
zpt_hdu.columns['V_0p33R26_ERR'].name = 'VROT_ERR'

hdul = fits.HDUList([empty_primary, coma_hdu, zpt_hdu])

hdul.writeto('paper_figures/Fig7/fig7_data-rev1.fits', overwrite=True)
#-------------------------------------------------------------------------------
# Header
#-------------------------------------------------------------------------------
hdr2 = fits.Header()

hdr2['DESI_DR'] = 'EDR'
hdr2['FIGURE'] = 8
hdr2['LOG_V0'] = float(f'{V0:.3f}')
#-------------------------------------------------------------------------------
# Data
#-------------------------------------------------------------------------------
mcmcfit_hdu = fits.PrimaryHDU(header=hdr2, data=tfr_samples)

mcmcfit_hdu.writeto('paper_figures/Fig8/fig8_data-rev1.fits', overwrite=True)
################################################################################
# exit()


# THERE ARE NO DWARF GALAXIES IN THIS CALIBRATION
################################################################################
# Identify the dwarf galaxies
#-------------------------------------------------------------------------------
# First, define the line perpendicular to the calibration
#-------------------------------------------------------------------------------
logV_n17 = (-17 - b_fit[1])/m_fit + V0
b_perp = -17 + (logV_n17 - V0)/m_fit
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
    
yvals_perp = -(xvals - V0)/m_fit + (b_perp + b_fit[0] - b_fit[1])


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
x1_err = 0.434*SGA_TF_bright['V_0p33R26_ERR']/SGA_TF_bright['V_0p33R26']
y1_err = SGA_TF_bright['R_MAG_SB26_ERR_CORR']
'''
data1 = [np.log10(SGA_TF['V_0p33R26']), SGA_TF['R_MAG_SB26_CORR']]
x1_err = 0.434*SGA_TF['V_0p33R26_ERR']/SGA_TF['V_0p33R26']
y1_err = SGA_TF['R_MAG_SB26_ERR_CORR']
corr1_xy = np.zeros_like(x1_err)

data2 = [np.log10(SGA_0pt['V_0p33R26']), SGA_0pt['R_ABSMAG_SB26']]
x2_err = 0.434*SGA_0pt['V_0p33R26_ERR']/SGA_0pt['V_0p33R26']
y2_err = SGA_0pt['R_ABSMAG_SB26_ERR']
corr2_xy = np.zeros_like(x2_err)
#-------------------------------------------------------------------------------

'''
#-------------------------------------------------------------------------------
# Dwarf galaxies
#-------------------------------------------------------------------------------
data_dwarfs = [np.log10(SGA_TF['V_0p33R26'][dwarfs]), 
               SGA_TF['R_MAG_SB26'][dwarfs]]
x_err_dwarfs = 0.434*SGA_TF['V_0p33R26_ERR'][dwarfs]/SGA_TF['V_0p33R26'][dwarfs]
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
# plt.savefig('../../Figures/SV/fuji_joint-Coma-0pt_TFR_varyV0-perpdwarfs_AnthonyUpdates_weightsVmax-1_dVsys_20250620.png', 
plt.savefig('../../../figures/SV_papers/fuji_joint-Coma-0pt_TFR_varyV0-perpdwarfs_AnthonyUpdates_weightsVmax-1_dVsys_corr_20260114.png', 
            dpi=150, 
            facecolor='none')
#-------------------------------------------------------------------------------
################################################################################



