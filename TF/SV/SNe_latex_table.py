'''
Create LaTeX data table of the SNe galaxies with the data formatted for the 
table.
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
data_filename = 'SGA_distances_0pt_fuji.fits'

# Output LaTeX file name
latex_filename = 'fuji_SNe_cal_galaxies.tex'

# Columns to include in LaTeX table
col_names = ['SGA_ID', 
             'RA', 
             'DEC', 
             'Z_DESI', 
             'D26', 
             'R_MAG_SB26', 
             'V_0p33R26', 
             'DM1_SN', 
             'SN']
err_dict = {'Z_DESI':'ZERR_DESI', 
            'R_MAG_SB26':'R_MAG_SB26_ERR', 
            'V_0p33R26':'V_0p33R26_err', 
            'DM1_SN':'e_DM1_SN'}

# Column header for table
colhead = '\\tablehead{\\colhead{SGA ID} & \\colhead{R.A.} & \\colhead{Decl.} & \\multicolumn{2}{c}{Redshift} & \\colhead{$D(26)$ [arcmin]} & \\multicolumn{2}{c}{$m_r(26)$} & \\multicolumn{2}{c}{$V(0.33R_{26})$ [km/s]} & \\multicolumn{2}{c}{$\\mu$} & \\colhead{SNIa}}'

# Table foot (caption, footnotes)
tabfoot = '\\tablecomments{{List} of the 2 galaxies used for calibrating the zero-point of the TFR.  Sky positions and diameters of the 26 mag arcsec$^{-2}$ isophote in the $r$ band are from the SGA.  Redshifts are measured from the DESI EDR spectra, and rotational velocities at $0.33R_{26}$ are computed as described in Sec.~\\ref{sec:measure_rot_vel}.  Distance moduli are from \\cite{Stahl2021}.}'

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
# Galaxies to include
#-------------------------------------------------------------------------------
SGAIDs_cal = Table()
SGAIDs_cal['SGAID'] = [294387, 464075]

cal_boolean = np.zeros(len(data_table), dtype=bool)

for i in range(len(SGAIDs_cal)):
    
    if SGAIDs_cal['SGAID'][i] in data_table_dict.keys():
        
        cal_boolean[data_table_dict[SGAIDs_cal['SGAID'][i]]] = True
################################################################################



################################################################################
# Build output table
#-------------------------------------------------------------------------------
out_table = Table()

for name in col_names:
    out_table[name] = data_table[name][cal_boolean]
    
    if name in err_dict.keys():
        out_table[err_dict[name]] = data_table[err_dict[name]][cal_boolean]
################################################################################



################################################################################
# Table meta-data
#-------------------------------------------------------------------------------
# Format functions and dictionary
#-------------------------------------------------------------------------------
def latex_ra(angle):
    h,m,s = Angle(angle, unit=u.deg).hms
    return '\\RA{{{0}}}{{{1}}}{{{2}}}{{{3}}}'.format('{:02d}'.format(int(abs(h))), '{:02d}'.format(int(abs(m))), '{:02d}'.format(int(abs(s))), '{:.2f}'.format(abs(s) - int(abs(s)))[1:])

def latex_dec(angle):
    if angle >= 0:
        sign = '+'
    else:
        sign = '-'
    d,m,s = Angle(angle, unit=u.deg).dms
    return sign + '\\dec{{{0}}}{{{1}}}{{{2}}}{{{3}}}'.format('{:02d}'.format(int(abs(d))), '{:02d}'.format(int(abs(m))), '{:02d}'.format(int(abs(s))), '{:.2f}'.format(abs(s) - int(abs(s)))[1:])

def latex_2err(error):
    err = '{:.2f}'.format(error)
    return '$\\pm${0}'.format(err)

def latex_3err(error):
    err = '{:.3f}'.format(error)
    return '$\\pm${0}'.format(err)

def latex_zerr(error):
    err = '{:.2f}'.format(1e6*error)
    return '$\\pm$({0}'.format(err) + '$\\times 10^{-6})$'

format_dict = {'SGA_ID':'%7d',
               'RA':latex_ra, 
               'DEC':latex_dec,
               'Z_DESI':'{:.6f}', 
               'ZERR_DESI':latex_zerr,
               'D26':'{:.2f}', 
               'R_MAG_SB26':'{:.2f}', 
               'R_MAG_SB26_ERR':latex_3err,
               'V_0p33R26':'{:.2f}',
               'V_0p33R26_err':latex_2err, 
               'DM1_SN':'{:.2f}', 
               'e_DM1_SN':latex_2err, 
               'SN':'%s'}
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