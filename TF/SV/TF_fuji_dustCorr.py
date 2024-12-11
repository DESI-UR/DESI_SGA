'''
This is a copy of TF_fuji_dustCorr.ipynb, to be used only to generate plots.  
There is no fitting that occurs in this file.
'''

################################################################################
# Import modules
#-------------------------------------------------------------------------------
################################################################################



################################################################################
# Plot mr v. b/a with fit
#-------------------------------------------------------------------------------
plt.figure(tight_layout=True)

plt.plot(gals['BA'], gals['R_MAG_SB26_CORR'], '.', alpha=0.2)
# plt.errorbar(ba, mr_bins, yerr=np.sqrt(cov_bins[1,1,:]), fmt='x', c='tab:blue')
plt.errorbar(ba + 0.05, mr_bins, yerr=np.sqrt(cov_bins[1,1,:]), fmt='x', c='darkblue')

plt.plot(np.array([0,1]), A_ba*np.array([0,1]) + (B + mr_median), 'k')

plt.xlim([0, 1])
plt.ylim([18.5, 11.75])

plt.tick_params('both', which='major', labelsize=14)

plt.xlabel('$b/a$', fontsize=16)
plt.ylabel('$m_r$', fontsize=16)

plt.savefig('../../Figures/SV/fuji_internalDustCorr_20241115.png', dpi=150)
################################################################################