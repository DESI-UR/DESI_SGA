# Directory Outline

1. `rot_curve_targets_loa.ipynb` - Gets the distances between center observation and all other fibers. It files for galaxies that have at least three fibers, and a center observation or two symmetric points. Or galaxies that have >10 unique points
  * Inputs
    * `desi_pv_loa_healpix.fits`
    * `desi_pv_fuji_healpix.fits`
    * `shredded_tables.txt`
    * `SGA-2020.fits` (`ELLIPSE`)
  * Output: `loa_targs_v2.fits`
    * Galaxies that pass criteria selection are marked with a 1 in `Selection` column
    * `DIST` - distance from the fiber to the center of the galaxy (in arcmin)
    * `DIST_R26` - the distance from the fiber to the center of the galaxy divided by R26 distance (in arcmin)
    * `PA` - the position angle of the galaxy
    * `C_TO_F_ANGLE` - the angle of separation from the center of the galaxy to the fiber
    * `ANGLE_OFF_AXIS` - the angle of the fiber from the major axis 

2. `projected_velocity.ipynb` - find projected velocities of shredded galaxies to reduce number of galaxies before VI
  * Inputs
    * `loa_targs_v2.fits`
    * `SGA-2020.fits` (`ELLIPSE`)
  * Output: `loa_targs_v2.fits`
    * `ZERR_MOD` - redshift error accounting for redrock 7km/s uncertainty
    * `unique_obs` - the number of unique observations in the galaxy
    * `Velocity` - projected velocity
    * `V_err` - velocity error
    * `Z_center` - center redshift

3. `cutouts_loa.ipynb` - This notebook generates cutouts for the galaxies and places the fibers on the cutout.
  * Inputs
    * `loa_targs.fits` (produced from `rot_curve_targets_loa.ipynb`)
    * `SGA-2020.fits` (`ELLIPSE`)
    * `get_cutouts.py`
      * This file generates the cutout and saves it to specified directory
  * Outputs
    * Cutouts - save to `$PSCRATCH` directory
    * Cutouts with fibers - save to `$PSCRATCH` directory

4. `shredded_VI.ipynb` - This notebook shows the cutouts, to make visual inspection faster.
  * Inputs
    * `loa_targs.fits`
    * `SGA-2020.fits` (`ELLIPSE`)
    * `loa_targs_missing.fits` (during VI, I realized the shredded galaxies did not contain all the fibers, so I did a crossmatch, redid steps 1 & 2 under `finding_mising_fibers.ipynb`, and added it to VI)
  * Outputs
    * `sga_ids.npy`
    * `redo_cutout.npy`
    * `target_ids.npy`
    * `matched_ids.npy`
    * `manual_VI_ids.npy`

5. `Filing_bad_fibers.ipynb` - After VI, all bad fibers are marked in this notebook with a 1 under `bad_fiber`
   * Input:
     * `loa_targs.fits`
     * `loa_targs_missing.fits`
   * Output: `shredded_VI.fits`
     * `deproj_dist` - deprojected distance in degrees
     * `deproj_r26` - deprojected distance in r26

6. `velocity_maps.ipynb` - This notebook calculates the rotational velocity all points in a galaxy and makes a velocity map
   * Inputs
     * `shredded_VI.fits`
     * `SGA-2020.fits` (`ELLIPSE`)
     * `get_cutouts.py` (Note: use the same directory that you have previously generated the cutouts from in `cutouts_loa.ipynb`)
     * `galaxy_selection.py`
     * `velocity_map_fxns.py`
  * Outputs
    * `loa_rot_velocity.fits`
      * `ZERR_MOD` corrects for redrock 7km/s systematic uncertainty
      * `Velocity` is the rotational velocity of each point (in km/s)
      * `V_err` is the error for rotational velocity (in km/s)
      * `Z_center` is the redshift of the center of the galaxy
      * `c_or_s` marks if the velocity was found using the center (marked as 0), or symmetric points (marked as 1)
    * velocity maps (saved to `$PSCRATCH`)

7. `VI_velmaps.ipynb` - similar to the first VI, this notebook helps VI go faster. All fibers and velocity maps that don't look right are removed
  * Inputs
    * `loa_rot_velocity.fits`
    * `SGA-2020.fits` (`ELLIPSE`)
  * Outputs
    * `shredded_vel_VI.fits`
   
8. `curve_fit.ipynb` - this produces the rotation curve fit for remaining galaxies
  * Inputs
    * `shredded_vel_VI.fits`
    * `SGA-2020.fits` (`ELLIPSE`)
  * Outputs
    * `shredded_rkpc.fits`
      * `r_kpc` - radius in units kpc
    * `loa_rotvel_curvefit.fits`
      * `chi2_reduced` - reduced chi^2
      * `vmax_fit` - fit param for maximum velocity
      * `rturn_fit` - fit param for radius at which curve goes from increasing to flat
      * `vmax_err` - error for vmax
      * `rturn_err` - error for rturn

9. `crossmatch.ipynb` - catalog crossmatching with SGA to get stellar and HI masses
  * Inputs
    * `SGA-2020.fits` (`ELLIPSE`)
    * `dr2_galaxy_sedfitting_v1.0.fits`
    * `SGA-2020_ALFALFA.fits`
    * `a100.code12.table2.190808.csv`
  * Outputs
    * `sga_xmatch.fits`
    * `sga_xmatch_HI.fits`

10. `total_mass.ipynb` - constructs total mass from curve fits
  * Inputs
    * `loa_rotvel_curvefit.fits`
  * Outputs
    * `sga_total_mass.fits`
   
11. `mass_plots` - constructs BTFR and visible-total mass
  * Inputs
    * `loa_rotvel_curvefit.fits`
    * `sga_total_mass.fits`
  * Outputs
    