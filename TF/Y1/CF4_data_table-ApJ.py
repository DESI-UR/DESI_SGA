'''
Create data table of galaxies matched with CF4 with the data to be published in 
ApJ.
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
# data_filename = 'CF4_DESI-v13_overlap.fits'
data_filename = 'CF4_DESI-v18_overlap.fits'

# Output data file name
out_filename = 'CF4_DESI-v18_overlap-ApJ.fits'

# Columns to include in data table
# Also include the number of significant digits to be used for the column
col_names = {'SGA_ID':0, 
             'RA':4, 
             'DEC':4, 
             'Z_DESI':5, 
             'D26':2, 
             'R_MAG_SB26':2, 
             'CF_rmag':2,
             'V_0p4R26':0, 
             'CF_Wmx':0,
             'MU_TF':2, 
             'CF_DM-r':2,
             'MAIN':0}
err_dict = {'Z_DESI':'ZERR_DESI', 
            'R_MAG_SB26':'R_MAG_SB26_ERR', 
            'V_0p4R26':'V_0p4R26_ERR', 
            'CF_Wmx':'CF_Wmx_err',
            'MU_TF':'MU_TF_ERR', 
            'CF_DM-r':'CF_DM-r_err'}
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

# Rename columns
out_table['CF_rmag'].name = 'R_MAG_CF4'
out_table['CF_Wmx'].name = 'WMX_CF4'
out_table['CF_Wmx_err'].name = 'WMX_CF4_ERR'
out_table['MU_TF'].name = 'MU_DESI'
out_table['MU_TF_ERR'].name = 'MU_DESI_ERR'
out_table['CF_DM-r'].name = 'MU_CF4'
out_table['CF_DM-r_err'].name = 'MU_CF4_ERR'
################################################################################



################################################################################
# Write table to file
#-------------------------------------------------------------------------------
empty_primary = fits.PrimaryHDU(header=hdr)
table_hdu = fits.BinTableHDU(data=out_table)
hdul = fits.HDUList([empty_primary, table_hdu])

hdul.writeto(data_directory + out_filename, overwrite=True)
################################################################################