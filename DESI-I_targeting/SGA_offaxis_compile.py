
'''
The goal is to make a compiled list of all the off-axis targets for the largest 
SGA galaxies.
'''

import json

import os

import numpy as np

from astropy.table import Table
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy import wcs



# Pixel resolution
pix_scale = 0.25 # arcsec/pixel
pix_scale_arcmin = pix_scale/60
pix_scale_degrees = pix_scale/3600 # deg/pixel



################################################################################
# Files
#-------------------------------------------------------------------------------
galaxy_filenames = []
target_filenames = []
clean_target_filenames = []

for i in range(13):
    galaxy_filenames.append('SGA_large_galaxies_' + str(i) + '.fits')
    clean_target_filenames.append('../target_files/SGA_off-axis_targets_' + str(i) + '_cleaned.txt')
    target_filenames.append('../target_files/SGA_off-axis_targets_' + str(i) + '.txt')
################################################################################




def xra(ra0, dec0, x, y):

    sep = np.sqrt(x**2 + y**2)
    c = np.arctan(sep)

    return ra0 + np.arctan2(x*np.sin(c), 
                            sep*np.cos(dec0)*np.cos(c) - y*np.sin(dec0)*np.sin(c))


def ydec(dec0, x, y):

    sep = np.sqrt(x**2 + y**2)
    c = np.arctan(sep)

    return np.arcsin(np.cos(c)*np.sin(dec0) + (y*np.sin(c)*np.cos(dec0))/sep)


def xy_to_radec(x, y, ra0=0, dec0=0, x0=0, y0=0, xscale=1, yscale=1):
    '''
    Converts x/y (pixels) to RA/DEC (degrees) position using the TAN Gnomonic 
    projection system.

    From https://rdrr.io/github/AngusWright/LAMBDAR/src/R/xy.to.radec.R


    PARAMETERS
    ==========

    x, y : float
        Point of interest

    ra0, dec0 : float
        anchor point

    x0, y0 : float
        anchor point

    xscale, yscale : float
        scale (degrees per pixel)
    '''


    ra0 = ra0*np.pi/180.
    dec0 = dec0*np.pi/180.

    xscale = xscale*np.pi/180.
    yscale = yscale*np.pi/180.

    x = (x0 - x)*np.tan(xscale)
    y = (y0 - y)*np.tan(yscale)

    RA = xra(ra0, dec0, x, y)*180./np.pi
    DEC = ydec(dec0, x, y)*180./np.pi

    return RA, DEC




################################################################################
# Build table of targets
#-------------------------------------------------------------------------------
target_table = Table()

GALAXY = []
RA = []
DEC = []

for i in range(len(galaxy_filenames)):

    ############################################################################
    # Open files
    #---------------------------------------------------------------------------
    galaxy_table = Table.read(galaxy_filenames[i], format='fits')

    if os.path.isfile(clean_target_filenames[i]):
        infile = open(clean_target_filenames[i], 'r')
    else:
        infile = open(target_filenames[i], 'r')

    targets = json.load(infile)
    infile.close()
    ############################################################################


    ############################################################################
    # Build look-up dictionary for galaxy table
    #---------------------------------------------------------------------------
    galaxy_index = {}

    for i in range(len(galaxy_table)):
        galaxy_index[galaxy_table['GALAXY'][i]] = i
    ############################################################################


    ############################################################################
    # Put targets into lists
    #---------------------------------------------------------------------------
    for galaxy in targets:

        ########################################################################
        # Find galaxy in galaxy table
        #-----------------------------------------------------------------------
        i_gal = galaxy_index[galaxy]
        ########################################################################


        ########################################################################
        # Get center coordinates of galaxy (image)
        #-----------------------------------------------------------------------
        ra_center = galaxy_table['SGA_RA'][i_gal]
        dec_center = galaxy_table['SGA_DEC'][i_gal]
        ########################################################################
        

        ########################################################################
        # Determine size of image needed
        #-----------------------------------------------------------------------
        major_axis = galaxy_table['DIAM'][i_gal]

        major_axis_pixels = major_axis/pix_scale_arcmin

        img_size = int(major_axis_pixels + 100)
        ########################################################################


        ########################################################################
        # Download image file and create WCS object from header
        #-----------------------------------------------------------------------
        img_url = 'https://www.legacysurvey.org/viewer/cutout.fits?ra={}&dec={}&%22/pix={}&layer=ls-dr9&size={}'.format(ra_center, dec_center, pix_scale, img_size)

        try:
            hdu = fits.open(img_url, format='fits')
        except:
            if os.path.isfile('../large_gal_images/' + galaxy + '.fits'):
                hdu = fits.open('../large_gal_images/' + galaxy + '.fits', 
                                format='fits')
            else:
                if galaxy == 'NGC0205':
                    print(galaxy, 'has no image!')
                    continue
                    
                print(galaxy, img_url)
                raise

        gal_header = hdu[0].header
        w = WCS(gal_header)
        ########################################################################


        ########################################################################
        # Transform pixel coordinates to ra, dec
        #-----------------------------------------------------------------------
        galaxy_target_pixels = np.array(targets[galaxy])

        if galaxy_target_pixels.shape[0] == 0:
            print(galaxy, 'has no off-axis targets!')
            continue

        galaxy_target_pixels = galaxy_target_pixels[:,::-1]

        galaxy_target_pixels[:,1] = np.abs(galaxy_target_pixels[:,1] - gal_header['IMAGEH'])

        zeros = np.zeros((galaxy_target_pixels.shape[0], 1))

        galaxy_target_pixels = np.concatenate((galaxy_target_pixels, zeros), axis=1)

        coords = w.wcs_pix2world(galaxy_target_pixels, 0)

        GALAXY.extend([galaxy]*coords.shape[0])
        RA.append(coords[:,0])
        DEC.append(coords[:,1])
        ########################################################################
    ############################################################################


#-------------------------------------------------------------------------------
# Add lists to table
#-------------------------------------------------------------------------------
target_table['GALAXY'] = GALAXY
#target_table['ROW'] = ROW
#target_table['COL'] = COL
target_table['RA'] = np.concatenate(RA)
target_table['DEC'] = np.concatenate(DEC)
#-------------------------------------------------------------------------------


fits_table = Table()
fits_table['RA'] = np.concatenate(RA)
fits_table['DEC'] = np.concatenate(DEC)
################################################################################




################################################################################
# Save results
#-------------------------------------------------------------------------------
target_table.write('../target_files/SGA_off-axis_targets.txt', 
                   format='ascii.commented_header', 
                   overwrite=True)

fits_table.write('../target_files/SGA_off-axis_targets.fits', 
                   format='fits', 
                   overwrite=True)
################################################################################



