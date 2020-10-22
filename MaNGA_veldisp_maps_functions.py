

import numpy as np
import numpy.ma as ma

from astropy.io import fits
from astropy.table import Table

import matplotlib.pyplot as plt

from MaNGA_veldisp_maps_plotting_functions import plot_diagnostic_panel, plot_veldisp_r

import sys
sys.path.insert(1, '/Users/kellydouglass/Documents/Research/Rotation_curves/RotationCurves/spirals/')
from DRP_rotation_curve_plottingFunctions import plot_rband_image


MANGA_SPAXEL_SIZE = 0.5 # spaxel size (0.5")

c = 3e5 # km/s
H0 = 100 # km/s/Mpc


def extract_veldisp(gal_ID, z, ba, phi, r50, maps_filename,
                    plot_diagnostics=True, IMAGE_DIR=None, IMAGE_FORMAT='eps'):
    '''
    Extract the velocity dispersions as a function of distance from the center 
    of the galaxy.


    PARAMETERS
    ==========

    gal_ID : string
        <plate> - <IFU> of galaxy

    z : float
        galaxy redshift

    ba : float
        galaxy semimajor to semiminor axis ratio

    phi : float
        Galaxy rotation angle E of N [degrees]

    r50 : float
        50% light radius of galaxy (from elliptical petrosian fit) [arcsec]

    maps_filename : string
        Location of data map fits file

    plot_diagnostics : boolean
        Flag to determine whether or not to plot the various plot diagnostics.  
        Default is True (plot all figures).

    IMAGE_DIR : string
        File path to which diagnostic images are saved.  
        Default value is None (do not save images).

    IMAGE_FORMAT : string
        Saved image file format.  Default format is eps.
    '''

    ############################################################################
    # Read in stellar velocity dispersion map
    #---------------------------------------------------------------------------
    maps = fits.open(maps_filename)

    # Extract average r-band image
    r_band = maps['SPX_MFLUX'].data

    # Extract (and correct for instrumental resolution effects) stellar velocity 
    # dispersion map
    star_sigma = np.sqrt(maps['STELLAR_SIGMA'].data**2 - maps['STELLAR_SIGMACORR'].data**2)
    # See https://www.sdss.org/dr16/manga/manga-data/working-with-manga-data/ for 
    # correction details.

    # Stellar velocity dispersion mask extension
    star_sigma_mask_extension = maps['STELLAR_SIGMA'].header['QUALDATA']
    #---------------------------------------------------------------------------
    # Mask arrays
    #---------------------------------------------------------------------------
    mStar_sigma = ma.array(star_sigma, mask=maps[star_sigma_mask_extension].data > 0)
    mr_band = ma.masked_equal(r_band, 0)
    ############################################################################




    ############################################################################
    # Check if all of the array is masked.
    #---------------------------------------------------------------------------
    num_masked_spaxels = np.sum(mStar_sigma.mask) - np.sum(mr_band.mask)
    frac_masked_spaxels = num_masked_spaxels/np.sum(~mr_band.mask)

    unmasked_data = False

    if frac_masked_spaxels < 1:
        unmasked_data = True
    ############################################################################




    ############################################################################
    # Find coordinates of center spaxel
    #
    # The center of the galaxy is defined to be the same as the brightest spaxel 
    # in the galaxy.
    #---------------------------------------------------------------------------
    if unmasked_data:
        optical_center = np.asarray( mr_band.max() == mr_band).nonzero()

        x_center = optical_center[1][0]
        y_center = optical_center[0][0]
    else:
        rows, cols = mr_band.shape

        x_center = int(0.5*rows)
        y_center = int(0.5*cols)

        optical_center = ([y_center], [x_center])
    ############################################################################




    ############################################################################
    # Normalize velocity dispersion values by central value
    #---------------------------------------------------------------------------
    mStar_sigma_norm = mStar_sigma/mStar_sigma[optical_center]
    ############################################################################




    ############################################################################
    # Calculate distance from each spaxel to center spaxel, in spaxel units
    #---------------------------------------------------------------------------
    y, x = np.indices(mStar_sigma.shape)

    distance_spaxels = np.hypot(x - x_center, y - y_center)

    deproj_distance_spaxels = deproject_r(x - x_center, y - y_center, distance_spaxels, ba, phi)
    ############################################################################




    ############################################################################
    # Convert distances to arcseconds, kpc
    #---------------------------------------------------------------------------
    distance_arcsec = distance_spaxels*MANGA_SPAXEL_SIZE
    deproj_distance_arcsec = deproj_distance_spaxels*MANGA_SPAXEL_SIZE
    deproj_distance_arcsec_norm = deproj_distance_arcsec/r50


    dist_to_galaxy_kpc = ( z * c / (H0/1000))
    spaxel_scale_factor = dist_to_galaxy_kpc * np.tan( MANGA_SPAXEL_SIZE*(1/60)*(1/60)*(np.pi/180))
    distance_kpc = spaxel_scale_factor*distance_spaxels
    deproj_distance_kpc = spaxel_scale_factor*deproj_distance_spaxels
    ############################################################################




    ############################################################################
    # Create table of distances (spaxels, arcseconds, and kpc) and stellar 
    # velocity dispersion value (raw, normalized to center spaxel value)
    #---------------------------------------------------------------------------
    veldisp_table = Table()
    veldisp_table['vel_disp'] = mStar_sigma.flatten()
    veldisp_table['vel_disp_norm'] = mStar_sigma_norm.flatten()
    veldisp_table['r_spaxels'] = distance_spaxels.flatten()
    veldisp_table['r_deproj_spaxels'] = deproj_distance_spaxels.flatten()
    veldisp_table['r_arcsec'] = distance_arcsec.flatten()
    veldisp_table['r_deproj_arcsec'] = deproj_distance_arcsec.flatten()
    veldisp_table['r_deproj_arcsec_norm'] = deproj_distance_arcsec_norm.flatten()
    veldisp_table['r_kpc'] = distance_kpc.flatten()
    veldisp_table['r_deproj_kpc'] = deproj_distance_kpc.flatten()
    ############################################################################


    
    if plot_diagnostics:
        
        plot_diagnostic_panel(gal_ID, mr_band, mStar_sigma, mStar_sigma_norm, 
                              veldisp_table, deproj_distance_arcsec_norm, IMAGE_DIR=IMAGE_DIR)
        
        
        plot_veldisp_r( veldisp_table['r_deproj_arcsec_norm'], 
                        veldisp_table['vel_disp_norm'], 
                        gal_ID, 
                        IMAGE_DIR=IMAGE_DIR)
        
        if IMAGE_DIR is None:
            plt.show()


    return veldisp_table




################################################################################
################################################################################
################################################################################


def deproject_r(x_coords, y_coords, distance_spaxels, axis_ratio, phi):
    '''
    Deproject the radial distances from the center of the galaxy.


    PARAMETERS
    ==========

    x_coords, y_coords : ndarray of shape (n,n)
        Index coordinate values for each spaxel.

    distance_spaxels : ndarray of shape (n,n)
        Radial distance from the galaxy's center spaxel, in spaxel units

    axis_ratio : float
        Galaxy ratio of semimajor to semiminor axis

    phi : float
        Axis of rotation of the semimajor axis E of N [degrees]


    RETURNS
    =======

    deproj_dist : ndarray of shape (n,n)
        Deprojected radial distance from galaxy's center spaxel
    '''


    tanTheta_prime = x_coords/y_coords

    Theta = np.pi - phi - np.arctan(tanTheta_prime)

    cosTheta = np.cos(Theta)
    sinTheta = np.sin(Theta)

    deproj_dist = distance_spaxels*np.sqrt(cosTheta*cosTheta \
                  + sinTheta*sinTheta/(axis_ratio*axis_ratio))

    return deproj_dist





