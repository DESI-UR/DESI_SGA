'''
This file is the same as joint-Coma-0pt_Fuji-TFR_KAD_varyV0-perpdwarfs.ipynb, 
but should only be used to generate the final corner plot; it does not contain 
any of the fitting procedures.
'''


################################################################################
# Load modules
#-------------------------------------------------------------------------------
import numpy as np

import matplotlib.pyplot as plt
plt.rcParams["text.usetex"] = False

import corner

import pickle
################################################################################




################################################################################
# Import data
#-------------------------------------------------------------------------------
# Best-fit
#-------------------------------------------------------------------------------
temp_infile = open('cov_ab_fuji_joint_TFR_varyV0-perpdwarfs1_KAD.pickle', 'rb')
cov_ab, tfr_samples, V0 = pickle.load(temp_infile)
temp_infile.close()
#-------------------------------------------------------------------------------
################################################################################




################################################################################
# Plot the corner plot
#-------------------------------------------------------------------------------
fig = corner.corner(tfr_samples.T, bins=30, smooth=1,
                    range=[[-8.2, -6.8], [15.5, 15.75], [-19.2, -18.25], [0.95, 1.2]],
                    labels=['$a$', '$b_{Coma}$', '$b_{0pt}$', r'$\sigma_{Coma}$'],
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
#-------------------------------------------------------------------------------
# Save the figure
#-------------------------------------------------------------------------------
plt.savefig('../../Figures/SV/fuji_joint_Coma_corner_20241009.png', 
            dpi=150, 
            facecolor='none', 
            bbox_inches='tight')
#-------------------------------------------------------------------------------
################################################################################