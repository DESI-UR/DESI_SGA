import numpy as np

import os

import requests

from astropy.wcs import WCS

import matplotlib as mpl
import matplotlib.pyplot as plt


################################################################################
# Plot formatting
#-------------------------------------------------------------------------------
mpl.rc('font', size=12)
#mpl.rc('axes', titlesize='small')
mpl.rc('figure', max_open_warning=0)
################################################################################



def plot_radec_DESI(table):
    '''
    Mollweide projection plot adapted to astro coordinates.
    
    
    PARAMETERS
    ==========
    
    table : astropy.table.Table
        Data table with secondary target info
        
    
    RETURNS
    =======
    
    fig : matplotlib.Figure
        Figure object to let user apply further plot manipulation
    '''
    
    fig, ax = plt.subplots(1,1, 
                           figsize=(8,4), 
                           subplot_kw={'projection':'mollweide'})
    
    ############################################################################
    # Loop through unique classes
    #---------------------------------------------------------------------------
    class_names = np.unique(table['PVTYPE'])
    
    for class_name in class_names:
        select = table['PVTYPE'] == class_name
        
        ########################################################################
        # Convert ra,dec to radians
        # Rotate the ra so that the plot goes 360->0 left to right
        #-----------------------------------------------------------------------
        _ra = np.radians(180. - table['TARGET_RA'][select])
        _dec = np.radians(table['TARGET_DEC'][select])
        ########################################################################
        
        ax.scatter(_ra, _dec, alpha=0.5, s=5, label=class_name)
    ############################################################################
    
    
    ############################################################################
    # Clean up the plot
    #---------------------------------------------------------------------------
    ax.set(xticks=np.radians([-150, -120, -90, -60, -30, 0, 30, 60, 90, 120, 150]), 
           xticklabels=['22h', '20h', '18h', '16h', '14h', '12h', '10h', '8h', '6h', '4h', '2h'])
    ax.tick_params(axis='x', labelsize=10)
    ax.tick_params(axis='y', labelsize=10)
    ax.grid(ls=':')
    
    ax.legend(fontsize=8, loc='lower right')
    fig.tight_layout()
    ############################################################################
    
    return fig



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