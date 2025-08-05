'''
Create LaTeX data table of the cluster galaxies with the data formatted for the 
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
data_filename = 'SGA-2020_iron_Vrot_cluster_calib_z0p1_Anthony2_dVsys.fits'

# Output LaTeX file name
latex_filename = 'iron_cluster_cal_galaxies_20250725.tex'

# Columns to include in LaTeX table
col_names = ['SGA_ID', 
             'RA', 
             'DEC', 
             'Z_DESI', 
             'D26', 
             'R_MAG_SB26', 
             'V_0p4R26']
err_dict = {'Z_DESI':'ZERR_DESI', 
            'R_MAG_SB26':'R_MAG_SB26_ERR', 
            'V_0p4R26':'V_0p4R26_ERR'}

# Column header for table
colhead = '\\tablehead{\\colhead{SGA-2020 ID} & \\colhead{R.A.} & \\colhead{Decl.} & \\colhead{Redshift} & \\colhead{$D(26)$} & \\colhead{$m_r(26)$} & \\colhead{$V(0.4R_{26})$} \\\[-0.5em] & [deg] & [deg] & & [arcmin] & [\\text{AB mag}] & [\\text{km/s}]}'

# Table foot (caption, footnotes)
tabfoot = '\\tablecomments{{List} of the \\NclusterGals galaxies in the \\Nclusters clusters used for calibrating the slope of the TFR.  Sky positions and diameters of the 26 mag arcsec$^{-2}$ isophote in the $r$-band are from the SGA-2020 \\citep{SGA}.  Redshifts are measured from the DESI DR1 spectra, and rotational velocities at $0.4R_{26}$ are computed as described in Sec.~\\ref{sec:measure_rot_vel}.}'

# Table name
tab_name = 'Cluster galaxies used for TFR slope calibration'

# Table label
tab_label = 'tab:cluster_cal'
################################################################################




################################################################################
# Read in galaxy data
#-------------------------------------------------------------------------------
data_table = Table.read(data_filename)
#-------------------------------------------------------------------------------
# Build dictionary to index into data_table
#-------------------------------------------------------------------------------
data_table_dict = {}

for i in range(len(data_table)):
    data_table_dict[data_table['SGA_ID'][i]] = i
################################################################################



################################################################################
# Build output table
#-------------------------------------------------------------------------------
out_table = Table()

for name in col_names:
    out_table[name] = data_table[name]
    
    if name in err_dict.keys():
        out_table[err_dict[name]] = data_table[err_dict[name]]
################################################################################



################################################################################
# Table meta-data
#-------------------------------------------------------------------------------
# Format functions and dictionary
#-------------------------------------------------------------------------------
def latex_1err(error):
    err = '{:.1f}'.format(error)
    return '\\pm{0}'.format(err)

def latex_3err(error):
    err = '{:.3f}'.format(error)
    return '\\pm{0}'.format(err)

def latex_zerr(error):
    # err = '{:.2f}'.format(1e6*error)
    # return '$\\pm$({0}'.format(err) + '$\\times 10^{-6})$'
    err = '{:.0f}'.format(1e5*error)
    return '({0})'.format(err)


format_dict = {'SGA_ID':'%7d',
               'RA':'{:.6f}', 
               'DEC':'{:.6f}', 
               'Z_DESI':'{:.5f}', 
               'ZERR_DESI':latex_zerr,
               'D26':'{:.2f}', 
               'R_MAG_SB26':'{:.2f}', 
               'R_MAG_SB26_ERR':latex_3err,
               'V_0p4R26':'{:.1f}',
               'V_0p4R26_ERR':latex_1err}
#-------------------------------------------------------------------------------
# Format table
#-------------------------------------------------------------------------------
for name in col_names:
    out_table[name].format = format_dict[name]
    
    if name in err_dict.keys():
        out_table[err_dict[name]].format = format_dict[err_dict[name]]
################################################################################



################################################################################
# Write table to file
#-------------------------------------------------------------------------------
out_table.write(latex_filename, 
                format='aastex', 
                caption=(tab_name + '\\label{' + tab_label + '}'), 
                latexdict={'preamble':'\\tablewidth{0pt}\n' + colhead, 
                           'tablefoot':tabfoot}, 
                overwrite=True)
################################################################################