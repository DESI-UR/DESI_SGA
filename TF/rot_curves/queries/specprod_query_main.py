from specprod_query import get_rot_curve_targets
from astropy.table import Table, vstack, join, unique


# set data release

redux = 'iron'

# load sga params table

sga_params_table = Table.read('/pscratch/sd/n/nravi/pv_rot_curves/sga_params.fits')
sga_params = [tuple(x) for x in sga_params_table['SGA_ID','RA', 'DEC', 'MAJOR', 'RATIO', 'PA'][:1297]]

# output
out_fn = '/pscratch/sd/n/nravi/pv_rot_curves/pv_tf_rotcurve_targets_' + redux + '_1.fits'

get_rot_curve_targets(redux, sga_params, table_out=out_fn, verbose=False)