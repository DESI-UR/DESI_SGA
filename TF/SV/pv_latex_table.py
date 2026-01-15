'''
Create LaTeX data table of galaxies w/ PVs with the data formatted for the SV 
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

# Output LaTeX file name
latex_filename = 'fuji_TF_pv_short-20250115.tex'

# Columns to include in LaTeX table
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

# Column header for table
colhead = '\\tablehead{\\colhead{SGA-2020} & \\colhead{R.A.} & \\colhead{Decl.} & \\colhead{Redshift} & \\colhead{$D(26)$} & \\colhead{$m_r(26)$} & \\colhead{$V(0.33R_{26})$} & \\colhead{$\mu$} & \\colhead{$V_{\\rm pec}$} \\\[-0.5em] \\colhead{ID} & \\colhead{[deg]} & \\colhead{[deg]} & & \\colhead{[arcmin]} & \\colhead{[AB mag]} & \\colhead{[km/s]} & \\colhead{[AB mag]} & \\colhead{[km/s]}}'

# Table foot (caption, footnotes)
tabfoot = '\\tablecomments{{Five} of the \\Nbright galaxies in DESI EDR with peculiar velocities measured using the calibrated TFR.  Sky positions and diameters of the 26 mag arcsec$^{-2}$ isophote in the $r$-band are from the SGA-2020 \\citep{SGA}.  Redshifts are measured from the DESI EDR spectra, and rotational velocities at $0.33R_{26}$ are computed as described in Sec.~\\ref{sec:measure_rot_vel}.  Distance moduli are calculated from the calibrated TFR, and peculiar velocities are based on the difference between the observed redshift and that inferred from the distance moduli following \\cite{Watkins2015}.  Table~\\ref{tab:pv} is published in its entirety online in a machine-readable format.  A portion is shown here for guidance regarding its form and content.}'

# Table name
tab_name = 'DESI EDR Peculiar Velocities from the TFR'

# Table label
tab_label = 'tab:pv'
################################################################################




################################################################################
# Read in galaxy data
#-------------------------------------------------------------------------------
data_table = Table.read(data_filename)
################################################################################



################################################################################
# Build output table
#-------------------------------------------------------------------------------
out_table = Table()

for name in col_names:
    out_table[name] = data_table[name][:5]
    
    if name in err_dict.keys():
        out_table[err_dict[name]] = data_table[err_dict[name]][:5]
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
    # err = '{:.2f}'.format(1e6*error)
    # return '$\\pm$({0}'.format(err) + '$\\times 10^{-6})$'
    err = '{:.0f}'.format(1e5*error)
    return '({0})'.format(err)

def latex_verr(error):
    err = '{:.0f}'.format(error)
    return '\\pm{0}'.format(err)
    

format_dict = {'SGA_ID':'%7d',
               'RA':'{:.6f}', #latex_ra, 
               'DEC':'{:.6f}', #latex_dec, 
               'Z_DESI':'{:.5f}', 
               'ZERR_DESI':latex_zerr,
               'D26':'{:.2f}', 
               'R_MAG_SB26':'{:.2f}', 
               'R_MAG_SB26_ERR':latex_3err,
               'V_0p33R26':'{:.1f}',
               'V_0p33R26_ERR':latex_1err, 
               'MU_TFbright':'{:.2f}', 
               'MU_TFbright_ERR':latex_2err, 
               'V_PEC':'{:.0f}', 
               'VERR_PEC':latex_verr}
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