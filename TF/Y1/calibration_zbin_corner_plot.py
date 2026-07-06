'''
This file generates the TF calibration corner plot.  It does not contain any of 
the fitting procedures.
'''

################################################################################
# Import modules
#-------------------------------------------------------------------------------
import numpy as np

import pickle

from corner import corner

import matplotlib.pyplot as plt
################################################################################



################################################################################
# Read in best-fit pickle file
#-------------------------------------------------------------------------------
temp_infile = open(#'cov_ab_iron_jointTFR_varyV0-dwarfsAlex_z0p1_zbins0p005_weightsVmax-1_dVsys_KAD-20250813.pickle', 
                   'cov_ab_iron_jointTFR_v14.pickle',
                   'rb')
cov_tfr, tfr_mcmc_samples, logV0, zmin, zmax, dz, zbins = pickle.load(temp_infile)
temp_infile.close()

# Number of redshift bins
m = len(zbins) - 1
# m = len(zbins) - 2
# NOTE: m should really be just len(zbins)-1, but we calibrated the v13 catalog 
# with an extra redshift bin (0.01-0.015) that we want to not show in the paper.

# Create a boolean index array to mask the redshift bin we don't want to show
# mask = np.ones(tfr_mcmc_samples.shape[0], dtype=bool)
# mask[-2] = False
################################################################################



################################################################################
# Limits for each parameter in the corner plot
#-------------------------------------------------------------------------------
values = np.median(tfr_mcmc_samples, axis=1)
sigmas = np.sqrt(np.diag(cov_tfr))

multiple = 3.5

limits = []

for i in range(len(values)):
    # if i != len(values) - 2: # Need to skip the last z-bin
    limits.append([values[i] - multiple*sigmas[i], 
                   values[i] + multiple*sigmas[i]])
################################################################################



################################################################################
# Plot
#-------------------------------------------------------------------------------
labels  = ['$a$']
labels += [f'$b_{{ {k+1} }}$' for k in np.arange(m)]
labels += [r'$\sigma$']

fig = plt.figure(figsize=(50,50))

corner(tfr_mcmc_samples.T,#[mask].T, 
       bins=25, smooth=1, 
       range=limits, 
       labels=labels,
       label_kwargs={'fontsize':18},
       levels=(1-np.exp(-0.5), 1-np.exp(-2)),
       quantiles=[0.16, 0.5, 0.84],
       color='tab:blue',
       hist_kwargs={'histtype':'stepfilled', 'alpha':0.3},
       plot_datapoints=False,
       fill_contours=True,
       show_titles=True,
       title_fmt='.3f', 
       title_kwargs={"fontsize": 18, 'loc':'left', 'pad':10}, 
       fig=fig
       );

#-------------------------------------------------------------------------------
# Add a table of the best-fit parameters to the upper right corner
#-------------------------------------------------------------------------------
ax_table = fig.add_axes([0.55, 0.65, 0.35, 0.25])
ax_table.axis('off')

table_text = []
for i in range(len(labels)):
    table_text.append([f'${values[i]:.3f} \pm {sigmas[i]:.3f}$'])

param_table = ax_table.table(cellText=table_text, 
                             rowLabels=labels, 
                             loc='center', 
                             cellLoc='center', 
                             fontsize=48)
param_table.scale(0.35, 6)
#-------------------------------------------------------------------------------

for ax in fig.get_axes():
    ax.tick_params(axis='both', which='major', labelsize=16)

plt.show()

# plt.savefig('../../../figures/Y1_papers/TFcorner_Y1_v13alt_table.png', 
#             dpi=150, 
#             facecolor='none')
################################################################################



################################################################################
# Save figure data for paper
#-------------------------------------------------------------------------------
temp_infile = open('paper_figures/Fig7/fig7_data.pickle', 'wb')

pickle.dump(tfr_mcmc_samples, temp_infile)

temp_infile.close()
################################################################################

