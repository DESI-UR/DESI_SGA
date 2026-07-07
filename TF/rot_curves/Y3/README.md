# Directory Outline

1. `rot_curve_targets_loa.ipynb` - Gets the distances between center observation and all other fibers. It files for galaxies that have at least three fibers, and a center observation or two symmetric points.
  * Inputs
    * `desi_pv_loa_healpix.fits`
    * `desi_pv_fuji_healpix.fits`
    * `SGA-2020.fits` (`ELLIPSE`)
  * Output: `loa_targs.fits`
    * Galaxies that pass criteria selection are marked with a 1 in `Selection` column
    * `DIST` is distance from the fiber to the center of the galaxy (in arcmin)
    * `DIST_R26` is the distance from the fiber to the center of the galaxy divided by R26 distance (in arcmin)

2. `cutouts_loa.ipynb` - This notebook generates cutouts for the galaxies and places the fibers on the cutout.
  * Inputs
    * `loa_targs.fits` (produced from `rot_curve_targets_loa.ipynb`)
    * `SGA-2020.fits` (`ELLIPSE`)
    * `get_cutouts.py`
      * This file generates the cutout and saves it to specified directory
  * Outputs
    * Cutouts - save to `$PSCRATCH` directory
    * Cutouts with fibers - save to `$PSCRATCH` directory

3. `Filing_bad_fibers.ipynb` - After VI, all bad fibers are marked in this notebook with a 1 under `bad_fiber`
   * See also Google Spreadsheet in drive called 'fiber jury' for a comprehensive list of which fibers I removed and why
   * Input: `loa_targs.fits`
   * Output: `loa_targs_edited.fits`

4. `velocity_maps.ipynb` - This notebook calculates the rotational velocity all points in a galaxy and makes a velocity map
   * Inputs
     * `loa_targs_edited.fits` (produced from `Filing_bad_fibers.ipynb`)
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

5. `rot_vel_filing.ipynb` - after VI, all fibers and velocity maps that don't look right are removed
   * See also Google Spreadsheet in drive called 'Velocity Map jury' for a comprehensive list of which fibers I removed and why
   * Input: `loa_rot_velocity.fits`
   * Output: `loa_rotvel.fits`
     * `bad_map` indicates if the fibers are good (marked as 0) or not (marked as 1)
