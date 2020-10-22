

from astropy.table import Table

import numpy as np

from MaNGA_veldisp_maps_functions import extract_veldisp

import os


maps_directory = '/Users/kellydouglass/Documents/Research/data/SDSS/MaNGA/'



################################################################################
# Galaxies to analyze
#-------------------------------------------------------------------------------
run_all_galaxies = True

gal_IDs = ['7443-1902']
################################################################################



################################################################################
# Read in general data
#-------------------------------------------------------------------------------
drpall_filename = maps_directory + 'drpall-v2_4_3.fits'

drpall = Table.read(drpall_filename, format='fits')

if run_all_galaxies:
    gal_IDs = drpall['plateifu']
    redshifts = drpall['z']
    axis_ratios = drpall['nsa_elpetro_ba']
    rotation_angles = drpall['nsa_elpetro_phi']
    radii_50 = drpall['nsa_elpetro_th50_r']
else:
    redshifts = np.zeros(len(gal_IDs))
    axis_ratios = np.zeros(len(gal_IDs))
    rotation_angles = np.zeros(len(gal_IDs))
    radii_50 = np.zeros(len(gal_IDs))

    for i in range(len(gal_IDs)):

        gal_idx = np.asarray(drpall['plateifu'] == gal_IDs[i]).nonzero()

        redshifts[i] = drpall['z'][gal_idx]
        axis_ratios[i] = drpall['nsa_elpetro_ba'][gal_idx]
        rotation_angles[i] = drpall['nsa_elpetro_phi'][gal_idx]
        radii_50[i] = drpall['nsa_elpetro_th50_r'][gal_idx]
################################################################################



for i in range(len(gal_IDs)):

    gal_ID = gal_IDs[i]
    z = redshifts[i]
    ba = axis_ratios[i]
    phi = rotation_angles[i]
    r50 = radii_50[i]

    plate, fiber = gal_ID.split('-')

    maps_filename = maps_directory + 'MaNGA_DR16/HYB10-GAU-MILESHC/' + plate \
                    + '/manga-' + gal_ID + '-MAPS-HYB10-GAU-MILESHC.fits.gz'

    if os.path.exists(maps_filename):

        print('Processing galaxy', gal_ID, '(', i, 'of', len(gal_IDs), ')')

        veldisp_table = extract_veldisp(gal_ID, z, ba, phi, r50, maps_filename, 
                                        IMAGE_DIR='../figures/')

        ############################################################################
        # Save table of distances (spaxels, arcseconds, and kpc) and stellar 
        # velocity dispersion value (raw, normalized to center spaxel value)
        #---------------------------------------------------------------------------
        veldisp_filename = '../veldisp_data_files/' + gal_ID + '_veldisp.txt'

        veldisp_table.write(veldisp_filename, format='ascii.commented_header', 
                            overwrite=True)
        ############################################################################