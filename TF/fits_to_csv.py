'''
Convert a .fits file into a .csv file.
'''

################################################################################
# Import modules
#-------------------------------------------------------------------------------
from astropy.table import Table
################################################################################


################################################################################
# Data file to be converted
#-------------------------------------------------------------------------------
data_directory = '/Users/kdouglass/Documents/Research/data/DESI/Y1/'

fits_filename = 'DESI-DR1_TF_pv_cat_v14.fits'

columns_to_save = ['SGA_ID', 
				   'RA', 
				   'DEC', 
				   'Z_DESI', 'ZERR_DESI', 
				   'D26', 
				   'R_MAG_SB26', 'R_MAG_SB26_ERR', 
				   'V_0p4R26', 'V_0p4R26_ERR', 
				   'MU_TF', 'MU_TF_ERR', 
				   'LOGDIST', 'LOGDIST_ERR', 
				   'MAIN']
################################################################################


################################################################################
# Construct filename of output file
#-------------------------------------------------------------------------------
filename, ext = fits_filename.split('.')

csv_filename = filename + '.csv'
################################################################################


################################################################################
# Read in original file and save it as a csv
#-------------------------------------------------------------------------------
data = Table.read(data_directory + fits_filename)

data[columns_to_save].write(data_directory + csv_filename, 
							format='ascii.csv', 
							overwrite=True)
################################################################################