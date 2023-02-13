'''
Helper functions for retrieving galaxy cutouts, spectra for visual inspection
'''

import os

import requests

from astropy.wcs import WCS
from astropy.table import Table

import numpy as np

from desispec.io import read_spectra
from desispec.coaddition import coadd_cameras
from desispec.spectra import stack as specstack



def get_cutout(targetid, ra, dec, size, verbose=False):
    """
    Grab and cache legacy survey cutouts.
    
    Parameters
    ----------
    targetid : int
        DESI target ID.
    ra : float
        Right ascension (degrees).
    dec : float
        Declination (degrees).
    verbose : bool
        Add some status messages if true.
        
    Returns
    -------
    img_name : str
        Name of JPG cutout file written after query.
    w : astropy.wcs.WCS
        World coordinate system for the image.
    """
    # Either load an existing image or download a cutout.
    img_name = 'cache/coma_{}.jpg'.format(targetid)
    
    if os.path.exists(img_name):
        if verbose:
            print('{} exists.'.format(img_name))
    else:
        img_url = 'https://www.legacysurvey.org/viewer/cutout.jpg?ra={}&dec={}&zoom=14&layer=ls-dr9&size={}&sga'.format(ra, dec, size)
        if verbose:
            print('Get {}'.format(img_url))
            
        with open(img_name, 'wb') as handle: 
            response = requests.get(img_url, stream=True) 
            if not response.ok: 
                print(response) 
            for block in response.iter_content(1024): 
                if not block: 
                    break 
                handle.write(block)
                
    # Set up the WCS.
    wcs_input_dict = {
        'CTYPE1': 'RA---TAN',
        'CUNIT1': 'deg',
        'CDELT1': -0.262/3600,
        'CRPIX1': size/2 + 0.5,
        'CRVAL1': ra,
        'NAXIS1': size,
        'CTYPE2': 'DEC--TAN',
        'CUNIT2': 'deg',
        'CDELT2': 0.262/3600,
        'CRPIX2': size/2 + 0.5,
        'CRVAL2': dec,
        'NAXIS2': size
    }
    w = WCS(wcs_input_dict)
    
    return img_name, w





def get_spectra_for_sga(targtab, sga_id):
    """
    Access spectra corresponding to all observations of a given SGA galaxy.
    
    
    Parameters
    ----------
    targtab : astropy.table.Table
        List of SGA + TFT observations.
        
    sga_id : int
        SGA ID of the central galaxy we're querying.
    """
    prefix = os.environ['DESI_SPECTRO_REDUX']
    
    targets = targtab[targtab['SGA_ID'] == sga_id]
    
    targetids = targets['TARGETID']
    
    coadd_spectra = None
    
    healpixfile = prefix + '/' + targets['FILENAME'][0]
    # Replace "redrock" with "spectra" in filename
    healpixfile_parts = healpixfile.split('redrock')
    filename = 'spectra'.join(healpixfile_parts)
    
    for targetid in targetids:
        
        # Grab the spectra corresponding to the TARGETID we queried.
        fmap = Table.read(filename, 'FIBERMAP')
        in1d = np.in1d(fmap['TARGETID'], targetid)

        # If we have something, read in the spectra and coadd across b,r,z cameras.
        if np.any(in1d):
            spectra = read_spectra(filename)[in1d]
            coadds = coadd_cameras(spectra)
            coadds.scores = None

            # Store the coadded spectra in a data structure.
            if coadd_spectra is None:
                coadd_spectra = coadds
            else:
                coadd_spectra = specstack([coadd_spectra, coadds])
                    
    return coadd_spectra