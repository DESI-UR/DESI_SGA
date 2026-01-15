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
# Galaxy data file name
# data_filename = 'SGA_fuji_jointTFR-varyV0-perpdwarf_moduli_pec-Watkins15.fits'
# data_filename = 'SGA_fuji_jointTFR-varyV0-perpdwarf-zCMB_dVsys_moduli_pec-Watkins15.fits'
data_filename = 'SGA_fuji_jointTFR-varyV0-perpdwarf-zCMB_dVsys_corr_moduli-20260114_pec-Watkins15.fits'

# Output data file name
out_filename = 'fuji_TF_pv.fits'

# Columns to include in data table
col_names = ['SGA_ID', 
             'RA', 
             'DEC', 
             'Z_DESI', 
             'D26', 
             'R_MAG_SB26', 
             'V_0p33R26', 
             'MU_TFbright', 
             'V_PEC']
err_dict = {'Z_DESI':'ZERR_DESI', 
            'R_MAG_SB26':'R_MAG_SB26_ERR', 
            'V_0p33R26':'V_0p33R26_ERR', 
            'MU_TFbright':'MU_TFbright_ERR', 
            'V_PEC':'VERR_PEC'}
################################################################################




################################################################################
# Read in galaxy data
#-------------------------------------------------------------------------------
hdul = fits.open(data_filename)

hdr = hdul[0].header

data_table = Table(hdul[1].data)

hdul.close()
################################################################################



################################################################################
# Build output table
#-------------------------------------------------------------------------------
out_table = Table()

for name in col_names:
    out_table[name] = data_table[name]
    
    if name in err_dict.keys():
        out_table[err_dict[name]] = data_table[err_dict[name]]

# Rename columns
out_table['MU_TFbright'].name = 'MU_TF'
out_table['MU_TFbright_ERR'].name = 'MU_TF_ERR'
################################################################################



################################################################################
# Write table to file
#-------------------------------------------------------------------------------
empty_primary = fits.PrimaryHDU(header=hdr)
table_hdu = fits.BinTableHDU(data=out_table)
hdul = fits.HDUList([empty_primary, table_hdu])

hdul.writeto(out_filename, overwrite=True)
################################################################################