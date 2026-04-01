'''
Create Fig. 11 data file from EDR TF PV catalog.
'''

################################################################################
# Import modules
#-------------------------------------------------------------------------------
from astropy.io import fits
from astropy.table import Table
################################################################################



################################################################################
# Read in full catalog
#-------------------------------------------------------------------------------
# hdul = fits.open('../../SGA_fuji_jointTFR-varyV0-perpdwarf-zCMB_dVsys_moduli.fits')
hdul = fits.open('../../SGA_fuji_jointTFR-varyV0-perpdwarf-zCMB_dVsys_corr_moduli-20260114.fits')

gals = Table(hdul[1].data)

hdul.close()
################################################################################



################################################################################
# Save just the velocity and absolute magnitude info
#-------------------------------------------------------------------------------
# Build the header
#-------------------------------------------------------------------------------
hdr = fits.Header()

hdr['DESI_DR'] = 'EDR'
hdr['FIGURE'] = 11

empty_primary = fits.PrimaryHDU(header=hdr)
#-------------------------------------------------------------------------------
table_hdu = fits.BinTableHDU(data=gals['Z_DESI_CMB', 'ZERR_DESI', 'LOGDIST', 'LOGDIST_ERR', 'DWARF'])

table_hdu.columns['Z_DESI_CMB'].name = 'Z'
table_hdu.columns['ZERR_DESI'].name = 'Z_ERR'

hdul = fits.HDUList([empty_primary, table_hdu])

hdul.writeto('fig11_data.fits')
################################################################################