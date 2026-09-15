# Y3 Tully-Fisher Analysis

This directory contains the code used to calibrate the DESI DR2 TFR with observations from the DESI Peculiar Velocity Survey, a secondary targeting program in DESI.


1. `loa_rot_vel.ipynb` - This notebook computes the rotational velocity for as many galaxies within the Loa (DR2) sample as possible.
    * Inputs:
        * `desi_pv_tf_loa_healpix.fits` (produced from a SQL query of the DESI database)
    * Output: `SGA-2020_loa_Vrot_v*.fits`
        * Center observations have been cleaned (`DELTACHI2` > 25, `ZWARN` = 0)
        * Rotational velocities at 0.4$R_{26}$ satisfy 10 < $V$ < 1000 km/s and $\Delta V/V_{min} \leq 5$
        * Rotational velocities at 0.4$R_{26}$ have the same sign on the same side of the galaxy, and opposite signs on opposite sides
        * Any galaxy which has been removed due to VI is also not included (VI done with `TF_Y3_VI.ipynb`)
        * 7 km/s statistical uncertainty added to all reported Redrock uncertainties
        * 0.06 systematic uncertainty in b/a added to all galaxies
        * Velocities rescaled to a redshift of z=0.05 to account for cosmological surface brightness dimming

2. `TF_loa_internal-dustCorr.ipynb` - This notebook fits for the correlation between the observed apparent magnitude and the axis ratio of the galaxies in the Loa sample to correct for internal dust extinction.
    * Input: `SGA-2020_loa_Vrot_v*.fits` (produced from `loa_rot_vel.ipynb`)
    * Output: `loa_internalDust_nokcorr.pickle` (contains MCMC samples and median $m_r$ from linear fit to $m_r$ v. $b/a$)

3. `TF_Y3_zbin_calibration_v5_irrcal.ipynb` - This notebook calibrates the Tully Fisher relation using redshift bins 
    * Inputs:
        * `SGA-2020_loa_Vrot_v5.fits` (produced from `loa_rot_vel.ipynb`)
        * `loa_internalDust_nokcorr.pickle` (contains MCMC samples and median $m_r$ for internal dust correction)
        * `TFY3_Classification.csv` (contains morphology classifications from SSL binary classifiers)
    * Output:
        * `cov_ab_loa_jointTFR_ellipse_v5a_spirals.pickle` (contains covariance matrix, MCMC samples, and log $V_0$ value from calibration)
        * `cov_ab_loa_jointTFR_ellipse_v5a_irregulars.pickle` (contains covariance matrix, MCMC samples, and log $V_0$ value from calibration)
        * `SGA_loa_jointTFR_v5a.fits` (Main catalog)
