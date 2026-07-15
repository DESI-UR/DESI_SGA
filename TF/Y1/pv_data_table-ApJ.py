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
# data_directory = '/Users/kdouglass/Documents/Research/data/DESI/Y1/'
data_directory = ''

# Galaxy data file name
# data_filename = 'DESI-DR1_TF_pv_cat_v13.fits'
data_filename = 'SGA_iron_jointTFR_moduli-v18_20260708.fits'

# Output data file name
out_filename = 'DESI-DR1_TF_pv_cat_v18-ApJ.fits'

# Columns to include in data table
# Also include the number of significant digits to be used for the column
col_names = {'SGA_ID':0, 
             'RA':4, 
             'DEC':4, 
             'Z_DESI':5, 
             'D26':2, 
             'R_MAG_SB26':2, 
             'V_0p4R26':0, 
             'MU_TF':2, 
             'LOGDIST':2,
             'MAIN':0}
err_dict = {'Z_DESI':'ZERR_DESI', 
            'R_MAG_SB26':'R_MAG_SB26_ERR', 
            'V_0p4R26':'V_0p4R26_ERR', 
            'MU_TF':'MU_TF_ERR', 
            'LOGDIST':'LOGDIST_ERR'}
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

for name in col_names.keys():
    out_table[name] = data_table[name]
    
    if name in err_dict.keys():
        out_table[err_dict[name]] = data_table[err_dict[name]]

        # Round column to same decimal used for parent column
        out_table.round({err_dict[name]:col_names[name]})

# Round table values to specified sig. figs, as requested by the ApJ data editor
out_table.round(col_names)
################################################################################



################################################################################
# Write table to file
#-------------------------------------------------------------------------------
empty_primary = fits.PrimaryHDU(header=hdr)
table_hdu = fits.BinTableHDU(data=out_table)
hdul = fits.HDUList([empty_primary, table_hdu])

hdul.writeto(data_directory + out_filename, overwrite=True)
################################################################################