# import packages ------------------------------------------------
from astropy.table import Table
from astropy import units as u
from astropy.constants import c

import matplotlib as mpl
import matplotlib.pyplot as plt

from scipy.optimize import minimize

import numpy as np

import numdifftools as ndt

from curve_fit_fxns import *
#==================================================================

# initialize tables and empty columns ============================================
loa_rotvel = Table.read('/pscratch/sd/d/dbustos/rot_curves/shredded_rkpc.fits')

SGA = Table.read('/global/cfs/cdirs/cosmo/data/sga/2020/SGA-2020.fits', 'ELLIPSE')

SGA_dict = {}
for i in range(len(SGA)):
    SGA_dict[SGA['SGA_ID'][i]] = i

loa_rotvel['chi2_reduced'] = np.nan
loa_rotvel['vmax_fit'] = np.nan
loa_rotvel['rturn_fit'] = np.nan
loa_rotvel['vmax_err'] = np.nan
loa_rotvel['rturn_err'] = np.nan
#===================================================================================

# find curve fit & uncertainty & plot ----------------------------------------------
for sga_id in np.unique(loa_rotvel['SGA_ID']):
    
    targ_id = loa_rotvel[loa_rotvel['SGA_ID'] == sga_id]

    # grab all fibers with velocity < 1000 km/s and passed VI
    # valid_fibers = targ_id[(abs(targ_id['Velocity'])<1000) & (targ_id['bad_map']==0)]
    valid_fibers = targ_id[(abs(targ_id['Velocity'])<1000)]
    
    # grab radius
    r_kpc = valid_fibers['r_kpc']
    
    # make sure there are still 3 points to curve fit
    # if (len(valid_fibers) < 3) or (len(np.unique(r_kpc.round(5))) < 3):
    #     continue

    # absolute velocities
    velocity = abs(valid_fibers['Velocity'])

    v_err = valid_fibers['V_err']

    # normalize r
    idx = np.argmax(velocity)

    r26_ratio = valid_fibers['DIST_R26']
    r26 = r_kpc[idx]/r26_ratio[idx]
    
#----------------------------------------------------------------------------
#create a pseudo-center fiber if there isn't one to assist with curve fitting
#----------------------------------------------------------------------------

    # if np.all(r_kpc) != 0:
    #     z_cen_err = valid_fibers['Z_center'][0]
    #     velocity = np.append(velocity,0)
    #     r_kpc = np.append(r_kpc,0)
    #     v_err = np.append(v_err,c*z_cen_err)


#-----------------------------------------
# curve fitting
#-----------------------------------------
     
    # bounds for v max, r turn
    bounds = [(0,1000),(0.01,np.max(r_kpc))]

    # initial guess for v max, r turn
    initial_guess = [velocity[idx], r_kpc[idx]]

    #-------------
    # curve fit ----------------------------------
    #-------------
    min_fxn = minimize(fun = chi2, 
                 x0 = initial_guess, 
                 args = (velocity, v_err, r_kpc), 
                 bounds = bounds, 
                 method = 'Powell')
    #---------------------------------------------
    
    # vmax_fit, rturn_fit, alpha_fit
    min_fits = min_fxn.x
    
    #--------------------
    # get reduced chi 2 -------------------------------
    #--------------------
    chi2_fit= min_fxn.fun
    
    data_pts = len(velocity)

    prms = 2
    
    reduced_chi2 = chi2_reduced(chi2_fit, data_pts, prms)
    #------------------------------------------------
    
    

#-------------------------------
# uncertainty from curve fit
#-------------------------------
    
    hessian = ndt.Hessian(chi2)
    hess = hessian(min_fits,velocity, v_err, r_kpc)

    # get covariance matrix for errors
    # make empty matrix of nans to prevent from plotting
    try:
    #covariance matrix
        hess_inv = 2*np.linalg.inv(hess)
    # vmax_err, rturn_err, alpha_err
        fit_params_err = np.sqrt(np.diag(np.abs(hess_inv)))
    
    except np.linalg.LinAlgError:
        # Do an alternate to the above
        hess_inv = np.full((2,2),np.nan)
        
#----------------
# plot
#----------------
    fig = plt.figure()
    ax = fig.add_subplot()

    #---------------
    # curve fit -------------------------------------------
    #---------------
    r1 = np.linspace(0,r26,500)

    v_r = v_rot(r1,min_fits[0],min_fits[1])
    #------------------------------------------------------

    # sample size
    size = 1000
    
    #------------------------------
    # uncertainty on curve fit   ---------------------------------------------------------------
    #------------------------------
    if ~np.isnan(hess_inv).any():
        samples=np.random.multivariate_normal(min_fits,hess_inv,size=size)

        # make sure all samples are positive
        good_samples = samples[(samples > 0).all(axis = 1)]
    
        # empty array for sample curve fit to go in
        v_sample = np.zeros((len(good_samples),500))
        
        # get curve fit for each sample distribution
        for i in range(len(good_samples)):
            indx = samples[i]
            v_x = v_rot(r1, indx[0], indx[1])
            v_sample[i] = v_x
        
        #take standard deviation for all v_rot along each r
        std_dev = np.std(v_sample[~np.isnan(v_sample).any(axis=1)], axis = 0)

        # plot uncertainty
        ax.fill_between(r1/r26, v_r - std_dev, v_r + std_dev, alpha = .12, color = 'chartreuse')
        
        # plot curve fit
        ax.plot(r1/r26, v_r, color = 'mediumseagreen', zorder =1 )
        
    else:
        # plot curve fit
        ax.plot(r1/r26, v_r, color = 'mediumseagreen', zorder = 1)
        
        # note the covariance matrix returned NAN
        ax.annotate('invalid hessian', xy = (10,10), xycoords = 'figure pixels')
    #------------------------------------------------------------------------------------------------

    
    

    #------------------------------
    # v vs r points and errorbars -------------------------------------------------------------------------------------
    #------------------------------
    ax.errorbar(r_kpc/r26, velocity, yerr = v_err, ls = 'none', ecolor = 'xkcd:lighter purple', capsize = 5, zorder=2, alpha = .6)
    ax.scatter(r_kpc/r26, velocity, c = 'xkcd:vivid purple', zorder = 3)
    # ------------------------------------------------------------------------------------------------------------------
    
    ax.grid(ls=':')
    ax.set_title('SGA ID: {}'.format(sga_id))
    ax.set_xlim(None, 1.01)
    ax.set_ylim(-1, None)
    ax.set(xlabel='$r/R_{26}$', ylabel='$v_{rot}$ $($ $km/s$ $)$')

    
    
#-----------------------------------------  
# save everything to table or $PSCRATCH
#----------------------------------------- 
    
    fig.savefig('/pscratch/sd/d/dbustos/curve_fit/' + 'sga_{}_curve_fit.png'.format(sga_id), dpi=120)
    
    fig.clear()
    plt.close(fig)
    
    np.save('/pscratch/sd/d/dbustos/hessian/' + str(sga_id) + '_hessian.npy', hess)
    
    table_idx = (loa_rotvel['SGA_ID'] == sga_id) & (abs(loa_rotvel['Velocity']) < 1000)
    loa_rotvel['chi2_reduced'][table_idx] = reduced_chi2

    loa_rotvel['vmax_fit'][table_idx] = min_fits[0]
    loa_rotvel['rturn_fit'][table_idx] = min_fits[1]

    loa_rotvel['vmax_err'] = fit_params_err[0]
    loa_rotvel['rturn_err'] = fit_params_err[1]
    

    # sanity check ---------------------------------
    # print('sga_id:',str(sga_id))
    # print('minimize fits (v, r): ', min_fits)
    # print('fit_params_err (v, r): ', fit_params_err)
    # print('reduced chi2: ', reduced_chi2)

# save new table
loa_rotvel.write('/pscratch/sd/d/dbustos/rot_curves/loa_rotvel_curvefit.fits',format='fits',overwrite=True)