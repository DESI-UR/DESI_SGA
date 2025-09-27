'''
Create LaTeX data table of the SNe galaxies with the data formatted for the 
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
data_filename = 'SGA-2020_iron_Vrot_VI_0pt-Anthony2_calib_z0p1.fits'

# Output LaTeX file name
latex_filename = 'iron_SNe_cal_galaxies.tex'

# Columns to include in LaTeX table
col_names = ['SGA_ID', 
             'RA', 
             'DEC', 
             'Z_DESI', 
             'D26', 
             'R_MAG_SB26', 
             'V_0p4R26', 
             'MU_SN', 
             'CID']
err_dict = {'Z_DESI':'ZERR_DESI', 
            'R_MAG_SB26':'R_MAG_SB26_ERR', 
            'V_0p4R26':'V_0p4R26_ERR', 
            'MU_SN':'MU_SN_ERR'}

# Column header for table
colhead = '\\tablehead{\\colhead{SGA-2020 ID} & \\colhead{R.A.} & \\colhead{Decl.} & \\colhead{$z_{\\rm CMB}$} & \\colhead{$D(26)$} & \\colhead{$m_r(26)$} & \\colhead{$V(0.4R_{26})$} & \\colhead{$\mu$} & \\colhead{SN ID} \\\[-0.5em] & [\text{deg}] & [\text{deg}] &  & \\colhead{[arcmin]} & [\\text{AB mag}] & \\colhead{[\\text{km/s}]} & [\\text{mag}] & }'

# Table foot (caption, footnotes)
tabfoot = '\\tablecomments{{List} of the 24 galaxies used for calibrating the zero-point of the TFR.  Sky positions and diameters of the 26 mag arcsec$^{-2}$ isophote in the $r$-band are from the SGA-2020 \\citep{SGA}.  Redshifts are measured from the DESI DR1 spectra, and rotational velocities at $0.4R_{26}$ are computed as described in Sec.~\\ref{sec:measure_rot_vel}.  Distance moduli are from \\cite{Pantheon}. \\tablenotetext{a}{SNe with multiple observations and/or multiple SNe in the same galaxy that were averaged accounting for covariance.} \\tablenotetext{b}{Primary source, i.e.~the SN occurred in the TF galaxy.}}'

# Table name
tab_name = 'Galaxies used for TFR zero-point calibration'

# Table label
tab_label = 'tab:SNe_cal'
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

def latex_2err(error):
    err = '{:.2f}'.format(error)
    return '\\pm{0}'.format(err)

def latex_3err(error):
    err = '{:.3f}'.format(error)
    return '\\pm{0}'.format(err)

def latex_zerr(error):
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
               'V_0p4R26_ERR':latex_1err, 
               'MU_SN':'{:.2f}', 
               'MU_SN_ERR':latex_2err, 
               'CID':'%s'}
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