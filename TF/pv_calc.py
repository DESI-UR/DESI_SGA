'''
Calculate the peculiar velocities of galaxies based on some independent distance 
estimator (e.g., Tully-Fisher relation).
'''

################################################################################
# Import modules
#-------------------------------------------------------------------------------
import numpy as np

from astropy.table import Table
import astropy.constants as const
################################################################################



################################################################################
# Constants
#-------------------------------------------------------------------------------
c = const.c.to('km/s')

# Number of random samples to generate for uncertainties
N_samples = 10000

rng = np.random.default_rng()
################################################################################



################################################################################
# Import catalog of galaxies for which to calculate the peculiar velocities
#-------------------------------------------------------------------------------
data_directory = 'SV/'

filename = 'SGA_fuji_ITFR_moduli.fits'

galaxies = Table.read(data_directory + filename)
################################################################################



################################################################################
# Compute cosmological redshift
#-------------------------------------------------------------------------------
def z_cosmo(z_obs, M_z, M_ind):
    '''
    Compute the cosmological redshift of an object.
    
    PARAMETERS
    ==========
    z_obs : float or ndarray of shape (n,)
        observed redshift of object
    M_z : float or ndarray of shape (n,)
        absolute magnitude of object based on observed redshift
    M_ind : float or ndarray of shape (n,)
        absolute magnitude of object based on independent distance calibrator
    
    RETURNS
    =======
    z_cos : float or ndarray of shape (n,)
        redshift due to Hubble expansion (no peculiar motion included)
    '''
    z_cos = z_obs * 10**(0.2*(M_z - M_ind))
    return z_cos


galaxies['Z_COSMO'] = z_cosmo(galaxies['Z_DESI'], 
                              galaxies['R_ABSMAG_SB26'], 
                              galaxies['R_ABSMAG_SB26_TFbright'])
#-------------------------------------------------------------------------------
# Estimate uncertainty
#-------------------------------------------------------------------------------
galaxies['ZERR_COSMO'] = np.nan

for i in range(len(galaxies)):
    
    z_desi_random = rng.normal(galaxies['Z_DESI'][i], 
                               galaxies['ZERR_DESI'][i], 
                               size=N_samples)
    Mz_random = rng.normal(galaxies['R_ABSMAG_SB26'][i], 
                           galaxies['R_ABSMAG_SB26_err'][i], 
                           size=N_samples)
    Mtfr_random = rng.normal(galaxies['R_ABSMAG_SB26_TFbright'][i], 
                             galaxies['R_ABSMAG_SB26_TFbright_err'][i], 
                             size=N_samples)
    
    z_cosmo_random = z_cosmo(z_desi_random, Mz_random, Mtfr_random)
    
    galaxies['ZERR_COSMO'][i] = np.nanstd(z_cosmo_random)
#-------------------------------------------------------------------------------
################################################################################



################################################################################
# Compute peculiar redshift and velocity
#-------------------------------------------------------------------------------
def z_peculiar(z_obs, z_cos):
    '''
    Compute the peculiar redshift.
    
    PARAMETERS
    ==========
    z_obs : float or ndarray of shape (n,)
        observed redshift
    z_cos : float or ndarray of shape (n,)
        cosmological redshift
    
    RETURNS
    =======
    z_pec : float or ndarray of shape (n,)
        peculiar redshift
    '''
    z_pec = (z_obs - z_cos)/(1 + z_cos)
    return z_pec

galaxies['Z_PEC'] = z_peculiar(galaxies['Z_DESI'], galaxies['Z_COSMO'])

galaxies['V_PEC'] = c*galaxies['Z_PEC']
#-------------------------------------------------------------------------------
# Estimate uncertainties
#-------------------------------------------------------------------------------
galaxies['ZERR_PEC'] = np.nan

for i in range(len(galaxies)):
    
    z_desi_random = rng.normal(galaxies['Z_DESI'][i], 
                               galaxies['ZERR_DESI'][i], 
                               size=N_samples)
    z_cosmo_random = rng.normal(galaxies['Z_COSMO'][i], 
                                galaxies['ZERR_COSMO'][i], 
                                size=N_samples)
    
    z_pec_random = z_peculiar(z_desi_random, z_cosmo_random)
    
    galaxies['ZERR_PEC'][i] = np.nanstd(z_pec_random)


galaxies['VERR_PEC'] = c*galaxies['ZERR_PEC']
#-------------------------------------------------------------------------------
################################################################################



################################################################################
# Save file
#-------------------------------------------------------------------------------
updated_filename = filename[:-5] + '_pec.fits'

galaxies.write(data_directory + updated_filename, format='fits', overwrite=True)
################################################################################