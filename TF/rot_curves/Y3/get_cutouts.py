#Generate cutouts from legacy survey
import requests

from astropy.wcs import WCS

import os

pix_scale = 0.25 # arcsec/pixel

def get_cutout(sgaid, ra, dec, size, dr = 11, zoom = 14, dir = '', verbose=False):
    """Grab and save legacy survey cutouts.
    
    Parameters
    ----------
    sgaid : int
        SGA galaxy ID.
    ra : float
        Right ascension (degrees).
    dec : float
        Declination (degrees).
    size: int
        Number of pixels to make cutout
    dr: int, optional
        data release to get cutout from. Options are 9, 10, and 11. Defaults to 11
    zoom: int, optional
        How zoomed in you want the cutout. Defaults to 14 
    dir : 'str', optional
        Enter the directory you want the cutout to save in. Defaults to cache
    verbose : bool, optional
        Add some status messages if true
        
    Returns
    -------
    img_name : str
        Name of JPG cutout file written after query.
    w : astropy.wcs.WCS
        World coordinate system for the image.
    """
    # Either load an existing image or download a cutout.
    img_name = (dir) + 'sga_{}.jpg'.format(sgaid)

    if dr == 11:
        layer = 'ls-dr11-early-v2'
    if dr == 10:
        layer = 'ls-dr10'
    if dr == 9:
        layer = 'ls-dr9'
    
    if os.path.exists(img_name):
        if verbose:
            print('{} exists.'.format(img_name))
    else:
        img_url = 'https://www.legacysurvey.org/viewer/cutout.jpg?ra={}&dec={}&%22/pix={}&layer={}&size={}&zoom={}&sga'.format(ra, dec, pix_scale, layer, size, zoom)
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