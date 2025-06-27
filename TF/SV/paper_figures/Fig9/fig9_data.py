'''
Create Fig. 9 data file from EDR TF PV catalog.
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
hdul = fits.open('../../SGA_fuji_jointTFR-varyV0-perpdwarf-zCMB_dVsys_moduli.fits')

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
hdr['FIGURE'] = 9

empty_primary = fits.PrimaryHDU(header=hdr)
#-------------------------------------------------------------------------------
table_hdu = fits.BinTableHDU(data=gals['V_0p33R26', 'V_0p33R26_ERR', 'R_ABSMAG_SB26_CORR', 'R_ABSMAG_SB26_CORR_ERR', 'DWARF'])

table_hdu.columns['V_0p33R26'].name = 'VROT'
table_hdu.columns['V_0p33R26_ERR'].name = 'VROT_ERR'
table_hdu.columns['R_ABSMAG_SB26_CORR'].name = 'R_ABSMAG'
table_hdu.columns['R_ABSMAG_SB26_CORR_ERR'].name = 'R_ABSMAG_ERR'

hdul = fits.HDUList([empty_primary, table_hdu])

hdul.writeto('fig9_data.fits', overwrite=True)
################################################################################