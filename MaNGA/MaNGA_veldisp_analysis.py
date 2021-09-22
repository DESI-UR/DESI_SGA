

from astropy.table import QTable,Table

import numpy as np

from MaNGA_veldisp_analysis_functions import find_radii

import os


data_table_directory = '../veldisp_data_files/'



################################################################################
# Galaxies to analyze
#-------------------------------------------------------------------------------
run_all_galaxies = True

gal_IDs = ['8713-6101']
################################################################################



################################################################################
# Read in general data
#-------------------------------------------------------------------------------
drpall_master_filename = '/Users/kellydouglass/Documents/Research/Rotation_curves/RotationCurves/spirals/DRPall-master_file.txt'

drpall_master = QTable.read(drpall_master_filename, format='ascii.ecsv')
################################################################################



################################################################################
# Maximum variances to find
#-------------------------------------------------------------------------------
max_variances = [0.01, 0.05, 0.1]

if run_all_galaxies:

    for value in max_variances:
        drpall_master['veldisp_' + str(value)] = np.zeros(len(drpall_master))

    # Build list of gal_IDs for calculations
    gal_IDs = []

    for i in range(len(drpall_master)):
        gal_IDs.append(str(drpall_master['MaNGA_plate'][i]) + '-' + str(drpall_master['MaNGA_IFU'][i]))
################################################################################



################################################################################
#-------------------------------------------------------------------------------
for i in range(len(gal_IDs)):

    gal_ID = gal_IDs[i]

    data_table_filename = data_table_directory + gal_ID + '_veldisp.txt'

    if os.path.exists(data_table_filename):

        print('Processing galaxy', gal_ID, '(', i+1, 'of', len(gal_IDs), ')')

        max_radii = find_radii(data_table_filename, max_variances)

        if run_all_galaxies:
            for j in range(len(max_variances)):
                drpall_master['veldisp_' + str(max_variances[j])][i] = max_radii[j]

        else:

            gal_table = Table()
            gal_table['veldisp_var'] = max_variances
            gal_table['max_radii'] = max_radii
            gal_table.pprint()
################################################################################



################################################################################
# Save table
#-------------------------------------------------------------------------------
if run_all_galaxies:
    drpall_filename = '../DRPall-master_file_veldisp.txt'

    drpall_master.write(drpall_filename, format='ascii.ecsv', overwrite=True)
################################################################################




