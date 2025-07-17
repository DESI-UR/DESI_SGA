This directory contains the code used to calibrate the DESI Y1 TFR with observations from the DESI Peculiar Velocity Survey, a secondary targeting program in DESI.

# Primary directory outline

1. `iron_rot_vel.ipynb` - This notebook computes the rotational velocity for as many galaxies within the Iron (Y1) sample as possible.
    * Inputs:
        * `desi_pv_tf_iron_healpix.fits` (produced from a SQL query of the DESI database)
        * Additional file with `PHOTSYS` column (and `SGA_ID`)
    * Output: `SGA-2020_iron_Vrot_dVsys_VI.fits`
        * Center observations have been cleaned (`DELTACHI2` > 25, `ZWARN` = 0)
        * Rotational velocities at 0.4$R_{26}$ satisfy 10 < $V$ < 1000 km/s and $\Delta V/V_{min} \leq 5$
        * Any galaxy which has been removed due to VI is also not included
        * 7 km/s statistical uncertainty added to all reported Redrock uncertainties

2. `TF_iron_internal-dustCorr.ipynb` - This notebook fits for the correlation between the observed apparent magnitude and the axis ratio of the galaxies in the Iron sample to correct for internal dust extinction.
    * Input: `SGA-2020_iron_Vrot_dVsys_VI_photsys.fits` (produced from `iron_rot_vel.ipynb`)
    * Output: `iron_internalDust_z0p1_mcmc.pickle` (contains MCMC samples and median $m_r$ from linear fit to $m_r$ v. $b/a$)

3. `TF_Y1_cluster_calibration_AnthonyUpdates2_weightsVmax-1_KAD.ipynb` - This notebook matches Iron galaxies to known galaxy clusters and calibrates the TFR based on those clusters and any galaxies with known independent distance measurements.
    * Inputs:
        * `SGA-2020_iron_Vrot_dVsys_VI_photsys.fits` (produced from `iron_rot_vel.ipynb`)
        * `iron_internalDust_z0p1_mcmc.pickle` (contains MCMC samples and median $m_r$ for internal dust correction)
        * `TFY1_Classification.csv` (contains morphology classifications from SSL binary classifiers)
        * `PVsample_group_lookup+SNe.csv` (contains distance moduli for all Y1 TF+FP galaxies)
        * `new_Pplus_group_match+SGA_avg.csv` (contains distance moduli of Y1 TF+FP galaxies that host a SNIa)
    * Output: `cov_ab_iron_jointTFR_varyV0-perpdwarfs0_*.pickle` (contains covariance matrix, MCMC samples, and log $V_0$ value from calibration)

4. `TF_iron-jointTFR-varyV0-perpdwarf-AnthonyUpdates-weightsVmax-1.ipynb` - This notebook calculates the distance moduli for the Iron galaxies using the calibrated TFR.
    * Inputs:
        * `SGA-2020_iron_Vrot_VI_photsys.fits` (produced from `iron_rot_vel.ipynb`)
        * `iron_internalDust_z0p1_mcmc.pickle` (contains MCMC samples and median $m_r$ for internal dust correction)
        * `cov_ab_iron_jointTFR_varyV0-perpdwarfs0_*.pickle` (contains covariance matrix, MCMC samples, and log $V_0$ value from calibration)
    * Output: `SGA_iron_jointTFR-varyV0-perpdwarf_*_moduli.fits` (header contains calibrated TFR values)

5. Peculiar velocities are calculated using the `pv_calc.py` script, found in the main TF directory (one level up from this one).

# Additional notebooks

`Iron_EDD.ipynb` - This notebook matches the Iron galaxies with rotational velocities with the Extragalactic Distance Database to identify potential 0-pt calibrators.

`iron_PVs.ipynb` - This notebook looks at the poperties of the peculiar velocities measured for the Iron galaxies using the calibrated TFR.

`TF_Y1_cluster_calibration_KAD-Stahl010pt.ipynb` - This notebook is a copy of `TF_Y1_cluster_calibration_KAD.ipynb` described above (step #3), but calibrates the 0-pt with those Iron galaxies with distance moduli from Stahl01's SNIa calibration.

`TF_iron-jointTFR-varyV0-perpdwarf_H0.ipynb` - This notebook fits for $H_0$ using the Iron galaxies and the calibrated TFR.

`TF_Y1_targets.ipynb` - This notebook is a preliminary glance at the statistics of the Iron sample.