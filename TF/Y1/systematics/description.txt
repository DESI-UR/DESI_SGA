File descriptions: 

1) DESI-DR1_TF_pv_cat_v14.fits
This is the original catalog.

2) DESI-DR1_TF_pv_cat_v14mod_Dcut.fits
This catalog includes all of the same cuts as before, but also applies a sixe cut of D_26 > 0.5' to the calibration and the full samples.

3) DESI-DR1_TF_pv_cat_v14mod_deltaV.fits
This catalog applies a cut of deltaV/Vmin < 2.5 rather than 5, which is propagated through the calibration and full sample.

4) DESI-DR1_TF_pv_cat_v14mod_nodeltaV.fits
This catalog does not apply a deltaV/Vmin cut which is propagated through the calibration and full sample.

5) DESI-DR1_TF_pv_cat_v14mod_hdbscan.fits
This catalog replaces Alex's velocity and magnitude cuts and replaces it by using HDBScan to pull out the central galaxy population.

6) DESI-DR1_TF_pv_cat_v14mod_JohnVIcal.fits
This catalog applies John's VI identification of good spiral galaxies to the calibration sample instead of the SSL morphological classifications used for the calibration sample in the original catalog.

7) DESI-DR1_TF_pv_cat_v14mod_JohnVIfull.fits
This catalog applies John's VI identification of good spiral galaxies to the calibration and full samples instead of the SSL morphological classifications used for the calibration sample in the original catalog.

8) DESI-DR1_TF_pv_cat_v14mod_spirals.fits
This catalog applies the SSL morphological classification cut to the calibration and full samples as opposed to just the calibration sample like the original catalog.
