This directory contains the code used to calibrate the DESI Y1 TFR with observations from the DESI Peculiar Velocity Survey, a secondary targeting program in DESI.

# Primary directory outline

1. `iron_rot_vel.ipynb` - This notebook computes the rotational velocity for as many galaxies within the Iron (Y1) sample as possible.
    * Inputs:
        * `desi_pv_tf_iron_healpix.fits` (produced from a SQL query of the DESI database)
        * Additional file with `PHOTSYS` column (and `SGA_ID`)
    * Output: `SGA-2020_iron_Vrot_dVsys_VI_photsys_v2.fits`
        * Center observations have been cleaned (`DELTACHI2` > 25, `ZWARN` = 0)
        * Rotational velocities at 0.4$R_{26}$ satisfy 10 < $V$ < 1000 km/s and $\Delta V/V_{min} \leq 5$
        * Rotational velocities at 0.4$R_{26}$ have the same sign on the same side of the galaxy, and opposite signs on opposite sides
        * Any galaxy which has been removed due to VI is also not included (VI done with `TF_Y1_VI.ipynb`)
        * 7 km/s statistical uncertainty added to all reported Redrock uncertainties

2. `TF_iron_internal-dustCorr.ipynb` - This notebook fits for the correlation between the observed apparent magnitude and the axis ratio of the galaxies in the Iron sample to correct for internal dust extinction.
    * Input: `SGA-2020_iron_Vrot_dVsys_VI_photsys_v2.fits` (produced from `iron_rot_vel.ipynb`)
    * Output: `iron_internalDust_z0p1_mcmc-20260426.pickle` (contains MCMC samples and median $m_r$ from linear fit to $m_r$ v. $b/a$)

3. `TF_Y1_zbin_calibration_weightsVmax-1_cutsAlex_KAD.ipynb` - This notebook calibrates the Tully Fisher relation using redshift bins 
    * Inputs:
        * `SGA-2020_iron_Vrot_dVsys_VI_photsys_v2.fits` (produced from `iron_rot_vel.ipynb`)
        * `iron_internalDust_z0p1_mcmc-20260424.pickle` (contains MCMC samples and median $m_r$ for internal dust correction)
        * `TFY1_Classification.csv` (contains morphology classifications from SSL binary classifiers)
    * Output: `cov_ab_iron_jointTFR_varyV0-dwarfsAlex_z0p1_zbins0p005_weightsVmax-1_dVsys_*.pickle` (contains covariance matrix, MCMC samples, and log $V_0$ value from calibration)

4. `TF_iron-jointTFR-dwarfAlex-zbin0p005-weightsVmax-1.ipynb` - This notebook calculates the distance moduli for the Iron galaxies using the calibrated TFR.
    * Inputs:
        * `SGA-2020_iron_Vrot_dVsys_VI_photsys_v2.fits` (produced from `iron_rot_vel.ipynb`)
        * `iron_internalDust_z0p1_mcmc.pickle` (contains MCMC samples and median $m_r$ for internal dust correction)
        * `cov_ab_iron_jointTFR_varyV0-dwarfsAlex_z0p1_zbins0p005_weightsVmax-1_dVsys_*.pickle` (contains covariance matrix, MCMC samples, and log $V_0$ value from calibration)
    * Output: `SGA_iron_jointTFR-varyV0-dwarfAlex_*_moduli.fits` (header contains calibrated TFR values)

5. Peculiar velocities are calculated using the `pv_calc.py` script, found in the main TF directory (one level up from this one).

6. The $\eta$ covariance matrix for the main cosmological sample is calculated in `TF_Y1_cov.ipynb`

# Additional notebooks

## Cosmology

`combine_fp_tf_copy.ipynb` - This notebook is a copy of Cullan's notebook that combines the FP and TF catalogs, renormalized the $\eta$ values from both to have an average of 0.

`Cullan_debug.ipynb`

`iron_PV_H0_study.ipynb` - This notebook looks at the influence of $H_0$ on $\eta$ v. $z$.

`TF_iron-jointTFR-varyV0-perpdwarf_H0.ipynb` - This notebook fits for $H_0$ using the Iron galaxies and the calibrated TFR.

`Y1_y-intercepts.ipynb` - Do the $y$-intercepts of the $z$-bin calibration follow Hubble's Law?

## Zero-pointing

`Iron_EDD.ipynb` - This notebook matches the Iron galaxies with rotational velocities with the Extragalactic Distance Database to identify potential 0-pt calibrators.

## Catalog properties

`iron_PVs.ipynb` - This notebook looks at the poperties of the peculiar velocities measured for the Iron galaxies using the calibrated TFR.

`iron_Vcomp.ipynb` - This notebook looks at how the velocities compare when we have observations made on both sides of the galaxy's center.

`TF_iron-BayesTFR-dwarfAlex-field.ipynb` - This notebook looks at $\eta$ v. $z$ for Alex's Bayesian TFR calibration.

`TF_iron_onsky.ipynb` - This notebook plots the on-sky distribution of galaxies in the catalog.

`TF_PV_corr_checks.ipynb` - A companion notebook to `TF_iron_internal-dustCorr.ipynb`, this notebook checks to confirm that all correlations have been zeroed out by the internal dust correction.

`TF_redshift_trends.ipynb` - This notebook looks at $eta$ v. $z$ for different subsets of the TF catalog.

`TF_Y1_morph_stats.ipynb` - What are the properties of the TF catalog in each of the different SSL morphological classes?

`TF_Y1_targets.ipynb` - This notebook is a preliminary glance at the statistics of the Iron sample.

`TF_Y1_v15_stats.ipynb` - General properties of the TF catalog.

`TF_Y1_VI-John_compare.ipynb` - How do John's VI results compare with the SSL morphology classification?

## Alternate calibration notebooks

`generate_ellipse_contours.ipynb` - This notebook fits the $z$-bin TFRs with ellipses as an alternate definition of the main cosmology sample.

`iron_rot_vel_alex.ipynb` - This notebook replicates Alex's method of calculating rotational velocities (assuming all observed redshifts are drawn from a PDF)

`TF_Y1_calibration_v15-KAD.ipynb` - This notebook is an alternate of `TF_Y1_zbin_calibration_weightsVmax-1_cutsAlex_KAD.ipynb`, corresponding to the v15-catalog calibration.

`TF_Y1_cluster_calibration_KAD-Stahl010pt.ipynb` - This notebook is a copy of `TF_Y1_cluster_calibration_KAD.ipynb` described above (step #3), but calibrates the 0-pt with those Iron galaxies with distance moduli from Stahl01's SNIa calibration.

`TF_Y1_clustering.ipynb` - This notebook looks at using a k-means clustering algorithm to define the main TF sample.

## Comparison with other catalogs

`Cosmicflows4_comparison_KAD.ipynb` and `Cosmicflows4_comparison.ipynb` - These notebooks compare the absolute magnitudes of our TFR calibration with those galaxies common to both our work and CF4.

`DESI-DR1-CF4.ipynb` - This notebook compares the distance modulus values reported in the CF4 catalog with our TFR calibration.

`catalog_trends.ipynb` and `Y1_calibration_comparison.ipynb` - Compares different catalogs.

`FitResults-Copy1.ipynb` - This notebook compares Alex's TF calibration with ours.
