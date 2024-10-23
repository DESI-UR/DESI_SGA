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
import astropy.units as u
################################################################################



################################################################################
# Constants
#-------------------------------------------------------------------------------
# h = 1 # EDR
h = 0.7742 # DR1
H0 = 100*h*u.km/u.s/u.Mpc

c = const.c.to('km/s')

# Number of random samples to generate for uncertainties
N_samples = 10000

rng = np.random.default_rng()
################################################################################



################################################################################
# Import catalog of galaxies for which to calculate the peculiar velocities
#-------------------------------------------------------------------------------
# data_directory = 'SV/'
data_directory = 'Y1/'

# filename = 'SGA_fuji_ITFR_moduli.fits'
# filename = 'SGA_fuji_jointTFR-varyV0-perpdwarf_moduli.fits'
filename = 'SGA_iron_jointTFR-varyV0-perpdwarf-fitH0_moduli.fits'

galaxies = Table.read(data_directory + filename)

if data_directory == 'SV/':
    mu_colname = 'mu_TFbright'
elif data_directory == 'Y1/':
    mu_colname = 'mu_TF'
################################################################################



################################################################################
# Compute peculiar velocities
#-------------------------------------------------------------------------------
def dist(mu):
    '''
    Compute the distance to an object based on a distance modulus.
    
    
    PARAMETERS
    ==========
    mu : float or ndarray of shape (n,)
        Distance moduli
    
    RETURNS
    =======
    d : float or ndarray of shape (n,)
        Distance in units of Mpc (or Mpc/h, depending on the units of the 
        distance modulus)
    '''
    d = 10 * 10**(0.2*mu) # pc or pc/h
    return d/1e6 # Mpc or Mpc/h


gal_d = dist(galaxies[mu_colname]) # Mpc/h or Mpc

galaxies['V_PEC'] = c*galaxies['Z_DESI'] - H0*(gal_d*u.Mpc) # km/s
#-------------------------------------------------------------------------------
# Estimate uncertainty
#-------------------------------------------------------------------------------
galaxies['VERR_PEC'] = np.nan

for i in range(len(galaxies)):
    
    z_desi_random = rng.normal(galaxies['Z_DESI'][i], 
                               galaxies['ZERR_DESI'][i], 
                               size=N_samples)
    mu_random = rng.normal(galaxies[mu_colname][i], 
                           galaxies[mu_colname + '_err'][i], 
                           size=N_samples)
    
    d_random = dist(mu_random)
    
    Vpec_random = c*z_desi_random - H0*(d_random*u.Mpc)
    
    galaxies['VERR_PEC'][i] = np.nanstd(Vpec_random.data)
#-------------------------------------------------------------------------------
################################################################################


"""
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
"""


################################################################################
# Save file
#-------------------------------------------------------------------------------
updated_filename = filename[:-5] + '_pec.fits'

galaxies.write(data_directory + updated_filename, format='fits', overwrite=True)
################################################################################