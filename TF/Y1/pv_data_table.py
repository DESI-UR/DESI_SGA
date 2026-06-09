'''
Create data table of galaxies w/ PVs with the data to be published with the 
paper.
'''


################################################################################
# Import modules
#-------------------------------------------------------------------------------
import astropy.units as u
from astropy.table import Table
from astropy.coordinates import Angle, SkyCoord
from astropy.io import fits

import numpy as np
################################################################################



################################################################################
# User input
#-------------------------------------------------------------------------------
# Galaxy data file directory
data_directory = '/Users/kdouglass/Documents/Research/data/DESI/Y1/'

# Galaxy data file name
data_filename = 'DESI-DR1_TF_pv_cat_v14.fits'

# Output data file name
out_filename = 'DESI-DR1_TF_pv_cat_v14-pub.fits'

# Columns to include in data table
col_names = ['SGA_ID', 
             'RA', 
             'DEC', 
             'Z_DESI', 
             'D26', 
             'R_MAG_SB26', 
             'V_0p4R26', 
             'MU_TF', 
             'LOGDIST',
             'MAIN']
err_dict = {'Z_DESI':'ZERR_DESI', 
            'R_MAG_SB26':'R_MAG_SB26_ERR', 
            'V_0p4R26':'V_0p4R26_ERR', 
            'MU_TF':'MU_TF_ERR', 
            'LOGDIST':'LOGDIST_ERR'}

# Units for each column with units
# Note that the same units will be used for the uncertainty column
unit_dict = {'RA':u.degree, 
             'DEC':u.degree, 
             'D26':u.arcmin, 
             'R_MAG_SB26':u.mag, 
             'V_0p4R26':u.km/u.s, 
             'MU_TF':u.mag}
################################################################################




################################################################################
# Read in galaxy data
#-------------------------------------------------------------------------------
hdul = fits.open(data_directory + data_filename)

hdr = hdul[0].header

data_table = Table(hdul[1].data)

hdul.close()
################################################################################



################################################################################
# Build output table
#
# Note that all output columns will have ALL CAP column names
#-------------------------------------------------------------------------------
out_table = Table()

for name in col_names:
    out_table[name.upper()] = data_table[name]
    
    if name in err_dict.keys():
        out_table[err_dict[name].upper()] = data_table[err_dict[name]]

    # Add units to columns
    if (out_table[name.upper()].unit is None) and (name in unit_dict.keys()):
        out_table[name.upper()].unit = unit_dict[name]

        if name in err_dict.keys():
            out_table[err_dict[name].upper()].unit = unit_dict[name]
################################################################################



################################################################################
# Write table to file
#-------------------------------------------------------------------------------
empty_primary = fits.PrimaryHDU(header=hdr)
table_hdu = fits.BinTableHDU(data=out_table)
hdul = fits.HDUList([empty_primary, table_hdu])

hdul.writeto(data_directory + out_filename, overwrite=True)
################################################################################