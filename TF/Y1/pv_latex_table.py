'''
Create LaTeX data table of galaxies w/ PVs with the data formatted for the DR1 
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
# Galaxy data directory
data_directory = '/Users/kdouglass/Documents/Research/data/DESI/Y1/'

# Galaxy data file name
data_filename = 'DESI-DR1_TF_pv_cat_v14.fits'

# Output LaTeX file name
latex_filename = 'iron_TF_pv_short.tex'

# Columns to include in LaTeX table
col_names = ['SGA_ID', 
             'RA', 
             'DEC', 
             'Z_DESI', 
             'D26', 
             'R_MAG_SB26', 
             'V_0p4R26',
             'MU_TF',
             'LOGDIST', 
             'MAIN'
             ]
err_dict = {'Z_DESI':'ZERR_DESI', 
            'R_MAG_SB26':'R_MAG_SB26_ERR', 
            'V_0p4R26':'V_0p4R26_ERR', 
            'MU_TF':'MU_TF_ERR', 
            'LOGDIST':'LOGDIST_ERR'
            }

# Column header for table
colhead = '\\tablehead{\\colhead{SGA-2020} & \\colhead{R.A.} & \\colhead{Decl.} & \\colhead{Redshift} & \\colhead{$D(26)$} & \\colhead{$m_r(26)$} & \\colhead{$V(0.4R_{26})$} & \\colhead{$\mu$} & \\colhead{$\eta$} & \\colhead{Main} \\\[-0.5em] \\colhead{ID} & \\colhead{[deg]} & \\colhead{[deg]} & & \colhead{[arcmin]} & \colhead{[AB mag]} & \\colhead{[km/s]} & \\colhead{[AB mag]} & & \\colhead{sample}}'

# Table foot (caption, footnotes)
tabfoot = '\\tablecomments{{Ten} of the \\Ntot galaxies in the DESI DR1 TF catalog.  Sky positions and diameters of the 26 mag arcsec$^{-2}$ $r$-band isophote are from the SGA-2020 \\citep{SGA}.  Redshifts are measured from the DESI DR1 spectra, and rotational velocities at $0.4R_{26}$ are computed as described in Sec.~\\ref{sec:measure_rot_vel}.  Distance moduli are calculated from the calibrated TFR, and the log distance ratios are based on the difference between the observed and predicted distance moduli.  Table~\\ref{tab:pv} is published in its entirety online in a machine-readable format.  A portion is shown here for guidance regarding its form and content.}'

# Table name
tab_name = 'DESI DR1 TF catalog'

# Table label
tab_label = 'tab:pv'
################################################################################




################################################################################
# Read in galaxy data
#-------------------------------------------------------------------------------
data_table = Table.read(data_directory + data_filename)
################################################################################



################################################################################
# Build output table
#-------------------------------------------------------------------------------
out_table = Table()

for name in col_names:
    out_table[name] = data_table[name][:10]
    
    if name in err_dict.keys():
        out_table[err_dict[name]] = data_table[err_dict[name]][:10]
################################################################################



################################################################################
# Table meta-data
#-------------------------------------------------------------------------------
# Format functions and dictionary
#-------------------------------------------------------------------------------
def latex_1err(error):
    err = '{:.1f}'.format(error)
    return '$\\pm${0}'.format(err)
    
def latex_2err(error):
    err = '{:.2f}'.format(error)
    return '$\\pm${0}'.format(err)

def latex_3err(error):
    err = '{:.3f}'.format(error)
    return '$\\pm${0}'.format(err)

def latex_zerr(error):
    # err = '{:.2f}'.format(1e6*error)
    # return '$\\pm$({0}'.format(err) + '$\\times 10^{-6})$'
    err = '{:.0f}'.format(1e5*error)
    return '({0})'.format(err)
    

format_dict = {'SGA_ID':'%7d',
               'RA':'{:.4f}', 
               'DEC':'{:.4f}', 
               'Z_DESI':'{:.5f}', 
               'ZERR_DESI':latex_zerr,
               'D26':'{:.2f}', 
               'R_MAG_SB26':'{:.2f}', 
               'R_MAG_SB26_ERR':latex_2err,
               'V_0p4R26':'{:.1f}',
               'V_0p4R26_ERR':latex_1err, 
               'MU_TF':'{:.2f}', 
               'MU_TF_ERR':latex_2err, 
               'LOGDIST':'{:.2f}', 
               'LOGDIST_ERR':latex_2err, 
               'MAIN':'%s'}
#-------------------------------------------------------------------------------
# Format table
#-------------------------------------------------------------------------------
for name in col_names:
    out_table[name].format = format_dict[name]
    
    if name in err_dict.keys():
        out_table[err_dict[name]].format = format_dict[err_dict[name]]
#-------------------------------------------------------------------------------
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