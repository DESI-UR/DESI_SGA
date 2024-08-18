'''
This is essentially a copy of the notebook by the same name.  It is intended 
for remaking plots - there is no fitting occurring in this script.
'''

################################################################################
# Import modules
#-------------------------------------------------------------------------------
import numpy as np

import matplotlib.pyplot as plt
from matplotlib import cm, colors
import matplotlib as mpl

# We have more than 10 clusters, so we need to change the default color map so that they are all plotted with different colors
from cycler import cycler
plt.rcParams['axes.prop_cycle'] = cycler('color', plt.get_cmap('tab20').colors)

from astropy.table import Table

import os

from corner import corner

import pickle
################################################################################


################################################################################
# Matplotlib plot settings
#-------------------------------------------------------------------------------
plt.rcParams['text.usetex'] = False
plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
################################################################################


################################################################################
# Import data
#-------------------------------------------------------------------------------
# SGA Y1 (Iron) galaxies w/ rotational velocities that are in clusters
#-------------------------------------------------------------------------------
calpath = '.'
calfits = os.path.join(calpath, 'SGA-2020_iron_Vrot_cluster_calib.fits')
caltab = Table.read(calfits)
#-------------------------------------------------------------------------------
# SGA Y1 (Iron) galaxies w/ rotational velocities and SNIa distances
#-------------------------------------------------------------------------------
zeroptpath = '.'
zeroptab = Table.read(os.path.join(zeroptpath, 
                                   'SGA-2020_iron_Vrot_0pt_calib.fits'))
################################################################################


################################################################################
# Import best-fit TFR
#-------------------------------------------------------------------------------
temp_infile = open('cov_ab_iron_jointTFR_varyV0-perpdwarfs_KAD.pickle', 'rb')
cov_tfr, tfr_mcmc_samples, logV0 = pickle.load(temp_infile)
temp_infile.close()
################################################################################


################################################################################
# Separate the data into clusters, 0pt calibrators
#-------------------------------------------------------------------------------
# Extract successful clusters
_cluster_ids = np.unique(caltab['CLUSTER_ID'])
m = len(_cluster_ids)

# Pack the results into arrays data and cov.
logV, logV_err = [], []
mag, mag_err = [], []

# Start with the zero-point calibrators
logV.append(np.log10(zeroptab['V_0p4R26']) - logV0)
logV_err.append(0.434*zeroptab['V_0p4R26_err'] / zeroptab['V_0p4R26'])
mag.append(zeroptab['R_ABSMAG_SB26'])
mag_err.append(zeroptab['R_ABSMAG_SB26_err'])

# Loop over the SGA clusters
for k, _cluster_id in enumerate(_cluster_ids):
    select_cluster = np.in1d(caltab['CLUSTER_ID'], _cluster_id)
    logV.append(np.log10(caltab['V_0p4R26'][select_cluster]) - logV0)
    logV_err.append(0.434*caltab['V_0p4R26_err'][select_cluster] / caltab['V_0p4R26'][select_cluster])
    mag.append(caltab['R_MAG_SB26'][select_cluster])
    mag_err.append(caltab['R_MAG_SB26_ERR'][select_cluster])
################################################################################


################################################################################
# Corner plot
#-------------------------------------------------------------------------------
labels  = ['$a$']
labels += ['$b_{0\mathrm{pt}}$'] + [f'$b_{{ {k+1} }}$' for k in np.arange(m)]
labels += [r'$\sigma$']

fig = corner(tfr_mcmc_samples.T, bins=25, smooth=1,
             labels=labels,
             label_kwargs={'fontsize':18},
             labelpad=0.1,
             levels=(1-np.exp(-0.5), 1-np.exp(-2)),
             quantiles=[0.16, 0.5, 0.84],
             color='tab:blue',
             hist_kwargs={'histtype':'stepfilled', 'alpha':0.3},
             plot_datapoints=False,
             fill_contours=True,
             show_titles=True,
             title_kwargs={"fontsize": 18, 'loc':'left', 'pad':10})

for ax in fig.get_axes():
    ax.tick_params(axis='both', which='major', labelsize=16)

# fig.savefig('../../Figures/Y1/TF_Y1_cluster_calibration_binaryML_corner0.png', dpi=150, facecolor='none')
plt.close()
################################################################################


################################################################################
# TFR calibration plot
#-------------------------------------------------------------------------------
sigma_mcmc = np.percentile(tfr_mcmc_samples[-1], [16., 50., 84.])
a_mcmc = np.percentile(tfr_mcmc_samples[0], [16., 50., 84])
b_mcmc = []
for k in range(1, m+2):
    b_mcmc.append(np.percentile(tfr_mcmc_samples[k], [16., 50., 84.]))
b_mcmc = np.asarray(b_mcmc)

# fig, axes = plt.subplots(2,2, figsize=(10,7), sharex=True,
                           # gridspec_kw={'height_ratios':[4,1], 'hspace':0.04, 'wspace':0.25})
fig, axes = plt.subplots(1,2, figsize=(10,7), 
                         gridspec_kw={'hspace':0.04, 'wspace':0.25})

a_   = a_mcmc[1]
b0pt = b_mcmc[0,1]
b_   = b_mcmc[:,1][1:]

#-------------------------------------------------------------------------------
# Plot clusters
ax = axes[0]#,0]

_logv = np.arange(0, 3, 0.1) - logV0
for k in range(m):
    eb = ax.errorbar(x=logV[k+1] + logV0, y=mag[k+1], xerr=logV_err[k+1], yerr=mag_err[k+1],
                     fmt='.', label=f'{_cluster_ids[k]}')

    ax.plot(_logv + logV0, a_*_logv + b_[k], color=eb[0].get_color(), ls='--', alpha=0.5, label='fit')

ax.set(xlim=[1.25, 2.75],
       ylim=[18, 13])
ax.set_xlabel(r'$\log{(V_{0.4R_{26}}~[\mathrm{km/s}]}$)', fontsize=14)
ax.set_ylabel(r'$m_r (26)$', fontsize=14)
ax.legend(loc='upper left', fontsize=9, ncol=2);
#-------------------------------------------------------------------------------
# Plot calibrators
ax = axes[1]#[0,1]
eb = ax.errorbar(x=logV[0] + logV0, y=mag[0], xerr=logV_err[0], yerr=mag_err[0],
                 fmt='.', label=f'0-pt')

ax.plot(_logv + logV0, a_*_logv + b0pt, color=eb[0].get_color(), ls='--', label='fit')
ax.set(xlim=[1.25, 2.75],
       ylim=[-16.75, -23])
ax.set_xlabel(r'$\log{(V_{0.4R_{26}}~[\mathrm{km/s}]}$)', fontsize=14)
ax.set_ylabel(r'$M_r(26) = m_r(26) - \mu$', fontsize=14)
ax.legend(loc='upper left', fontsize=9, ncol=2)
'''
#-------------------------------------------------------------------------------
# Plot residuals: clusters
ax = axes[1,0]

for k in range(m):
    logv_obs = logV[k+1]
    m_obs = mag[k+1]
    m_exp = (a_*logv_obs + b_[k])
    eb = ax.errorbar(x=logv_obs + logV0, y=(m_exp-m_obs)/m_exp, xerr=logV_err[k+1], yerr=mag_err[k+1],
                     fmt='.', label=f'{_cluster_ids[k]}')

ax.axhline(0, ls=':', color='k', alpha=0.3)

ax.set(xlim=[1.25, 2.75],
       ylim=[-0.4, 0.4])
ax.set_xlabel(r'$\log{(V_{0.4R_{26}}~[\mathrm{km/s}]}$)', fontsize=14)
ax.set_ylabel(r'$\Delta m_r/m_r$', fontsize=14)

#-------------------------------------------------------------------------------
# Plot residuals: calibrators
ax = axes[1,1]

logv_obs = logV[0]
m_obs = mag[0]
m_exp = (a_*logv_obs + b0pt)

b = ax.errorbar(x=logv_obs + logV0, y=(m_exp-m_obs)/m_exp, xerr=logV_err[0], yerr=mag_err[0],
                     fmt='.', label=f'0-pt')
ax.axhline(0, ls=':', color='k', alpha=0.3)

ax.set(xlim=[1.25, 2.75],
       ylim=[-0.4, 0.4])
ax.set_xlabel(r'$\log{(V_{0.4R_{26}}~[\mathrm{km/s}]}$)', fontsize=14)
ax.set_ylabel(r'$\Delta M_r/M_r$', fontsize=14)
'''
fig.subplots_adjust(left=0.1, bottom=0.1, top=0.9, right=0.9);
# fig.savefig('../../Figures/Y1/TF_Y1_cluster_calibration_0pt_fit0.png', 
fig.savefig('../../Figures/Y1/TF_Y1_cluster_calibration_0pt_fit0_noResiduals.png', 
            dpi=150, 
            facecolor='none')
################################################################################