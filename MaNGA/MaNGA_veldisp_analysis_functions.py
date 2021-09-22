


from astropy.table import Table

import numpy as np



def find_radii(data_table_filename, max_variations):
    '''
    Calculate the normalized radius at which the distribution widens to more than 
    some fraction of the central velocity dispersion. 


    PARAMETERS
    ==========

    data_table_filename : string
        Location of galaxy's velocity dispersion data file (output from 
        MaNGA_veldisp_maps.py)

    max_variations : list or ndarray
        Maximum variation allowed in the velocity dispersion


    RETURNS
    =======

    max_radii : ndarray with the same length as max_variations
        First radius where the 
    '''

    ############################################################################
    # Read in galaxy's velocity dispersion data table
    #---------------------------------------------------------------------------
    data_table = Table.read(data_table_filename, format='ascii.commented_header')
    ############################################################################



    ############################################################################
    # In steps of 0.1 arcsec, calculate the width of the distribution of the 
    # normalized velocity dispersion
    #---------------------------------------------------------------------------
    bin_width = 0.1 # arcsec
    n_bins = int(np.nanmax(data_table['r_deproj_arcsec_norm'])/bin_width)

    max_radii = np.zeros(len(max_variations))

    if n_bins > 0:

        radii_bins = np.linspace(0, np.nanmax(data_table['r_deproj_arcsec_norm']), n_bins)

        bin_indices = np.digitize(data_table['r_deproj_arcsec_norm'], radii_bins)

        averages = np.zeros(n_bins)
        std_dev = np.zeros(n_bins)

        for i in range(n_bins):
            bin_boolean = bin_indices == i

            vel_disp_bin_values = data_table['vel_disp_norm'][bin_boolean]

            if np.sum(vel_disp_bin_values.mask) < len(vel_disp_bin_values):
                avg = np.nanmean(vel_disp_bin_values)
                stdev = np.nanstd(vel_disp_bin_values)
            else:
                avg = 1
                stdev = -1

            averages[i] = avg
            std_dev[i] = stdev
        ########################################################################
        # Search through the results and find the first time that the standard  
        # deviation is greater than some value
        #-----------------------------------------------------------------------
        for i in range(len(max_variations)):

            points_up = (averages + std_dev) > 1 + max_variations[i]
            points_down = (averages - std_dev) > 1 - max_variations[i]

            j = np.argmax(points) - 1

            max_radii[i] = radii_bins[j]
        ########################################################################
    ############################################################################



    return max_radii