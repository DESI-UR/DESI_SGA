'''
Functional versions of the photometric corrections seen in the 
TF_photoCorrect.ipynb notebook.
'''

################################################################################
# Load modules
#-------------------------------------------------------------------------------
import numpy as np

from astropy.table import Table
from astropy_healpix import HEALPix
from astropy.coordinates import SkyCoord, ICRS
import astropy.units as u

import kcorrect.kcorrect
################################################################################




################################################################################
# Imaging survey systematics
#
# Khaled found that there is a systematic offset in the magnitudes between BASS 
# and DECaLS:
#                      m_BASS - m_DECaLS = 0.0234
# The RMS deviation is 0.02 mag.
#
# This function adjusts the northern (i.e., BASS) photometry to match DECaLS.
#-------------------------------------------------------------------------------
def BASS_corr(gals):
    '''
    Khaled found that there is a systematic offset in the magnitudes between 
    BASS and DECaLS:
                        m_BASS - m_DECaLS = 0.0234
    The RMS deviation is 0.02 mag.
    
    This function adjusts the northern (i.e., BASS) photometry to match DECaLS.
    '''
    # Check if the PHOTOSYS column is in the provided astropy table
    if 'PHOTSYS' not in gals.colnames:
        print('Please add the PHOTSYS column to identify the photometric survey for each object.')
        return None
    
    corr = np.zeros(len(gals))
    corr_err = np.zeros(len(gals))
    
    bass_gal = gals['PHOTSYS'] == 'N'
    
    corr[bass_gal] += 0.0234
    corr_err[bass_gal] += 0.02
    
    return corr, corr_err
################################################################################




################################################################################
# MW dust extinction
#-------------------------------------------------------------------------------
def MW_dust(gals, ebv_map):
    '''
    Correct for the dust extinction due to the MW's dust.  This uses Rongpu's 
    DESI dust map.
    '''
    
    # Create a dictionary of healpixel coordinates for the map
    ebv_map_dict = {}
    for i in range(len(ebv_map)):
        ebv_map_dict[ebv_map['HPXPIXEL'][i]] = i
        
    # Ratio of total to selective extinction (rederived from Schalfly11 and 
    # listed in 
    # https://www.legacysurvey.org/dr10/catalogs/#galactic-extinction-coefficients
    Rr = 2.165
    
    # nside and order values are specific to the E(B-V) map used
    hp = HEALPix(nside=512, order='ring', frame=ICRS())
    
    # Compute HEALPix index of each galaxy
    gal_skyCoords = SkyCoord(ra=gals['RA']*u.deg, dec=gals['DEC']*u.deg)
    gal_hpCoords = hp.skycoord_to_healpix(gal_skyCoords)
    
    # Extract E(B-V) values for each galaxy
    EBV = np.zeros(len(gals))
    EBV_err = np.zeros(len(gals))
    for i in range(len(gals)):
        if gal_hpCoords[i] in ebv_map_dict.keys():
            i_ebv = ebv_map_dict[gal_hpCoords[i]]
            EBV[i] = ebv_map['EBV_GR'][i_ebv]
            EBV_err[i] = ebv_map['EBV_GR_ERR'][i_ebv]
            
    # Compute dust extinction correction
    # A_dust = Rr E(B-V)
    Adust = Rr*EBV
    Adust_err = Rr*EBV_err
    
    return Adust, Adust_err
################################################################################




################################################################################
# K-correct
#-------------------------------------------------------------------------------
def mag2maggies(m):
    '''
    Convert magnitudes to AB maggies
    '''
    return np.power(10, -0.4*m)


def err2ivar(error, maggie):
    '''
    Convert magnitude uncertainties to maggie inverse variances
    '''
    mag_ivar = 1/error**2
    return mag_ivar/np.square(0.4*np.log(10)*maggie)


def k_corr(gals, z_corr=0.05):
    '''
    Compute K-corrections
    '''
    filters = ['decam_g', 'decam_r', 'decam_z']
    kc = kcorrect.kcorrect.Kcorrect(responses=filters)
    
    magnitude = [gals['G_MAG_SB26'], gals['R_MAG_SB26'], gals['Z_MAG_SB26']]
    maggies = mag2maggies(np.array(magnitude).T)
    
    mag_errors = [gals['G_MAG_SB26_ERR'], gals['R_MAG_SB26_ERR'], gals['Z_MAG_SB26_ERR']]
    ivar = err2ivar(np.array(mag_errors).T, maggies)
    
    coeffs = kc.fit_coeffs(redshift=gals['Z_DESI'], maggies=maggies, ivar=ivar)
    
    k = kc.kcorrect(redshift=gals['Z_DESI'], coeffs=coeffs, band_shift=z_corr)
    
    return k
################################################################################