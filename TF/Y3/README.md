# Y3 Tully-Fisher Analysis

This directory contains the code used to calibrate the DESI DR2 TFR with observations from the DESI Peculiar Velocity Survey, a secondary targeting program in DESI.

# Primary Directory Outline 

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

3. `TF-Y3_calibration_v7a.ipynb` - This notebook calibrates the Tully Fisher relation using redshift bins 
    * Inputs:
        * `SGA-2020_loa_Vrot_v5.fits` (produced from `loa_rot_vel.ipynb`)
        * `loa_internalDust_nokcorr.pickle` (contains MCMC samples and median $m_r$ for internal dust correction)
        * `TFY3_Classification.csv` (contains morphology classifications from SSL binary classifiers)
    * Output:
        * `cov_ab_loa_jointTFR_ellipse_v*.pickle` (contains covariance matrix, MCMC samples, and log $V_0$ value from calibration)
        * `TF_Y3_TFR_fit_params_v*.fits` (table with best fit parameters: slope, per-bin intercepts, and uncertainties)
        * `DESI-DR2_TF_pv_cat_v*.fits` (Main catalog)
     
# Additional Notebooks

## Mocks

* `mocks/TFR_DR2_mock_gen.ipynb` - applies an equivalent calibration on Abacus mock data for use in downstream cosmology

## Other 
* `DR1_comparison` - Compares galaxy properties and fit between DESI DR1 and DR2 TF samples
* `Loa_Vcomp.ipynb` - Looks at how the velocities compare when we have observations made on both sides of the galaxy's center
* `mock_scale_tests.ipynb` - uses the calibration on mock data to determine the optimal size of the ellipse used for identifying the main galaxy population
* `Morphology_VI.ipynb` - Visual inspection of galaxy morphologies
* `TF_CF4_hyperfit.ipynb` - Apply our calibration method to CF4 data
* `TF_Y3_scatter.ipynb` - Looks at how inclination angle uncertainties would impact the TFR scatter
* `relative_PV_check.ipynb` - Looks at the relative contribution of peculiar velocities as a function of redshift
* `rot_curve_corrections` - Corrects the rotational velocity measurements, accounting for cosmological SB dimming in fiber placement (incorporated into `loa_rot_vel.ipynb`)

