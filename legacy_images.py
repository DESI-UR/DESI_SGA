import os
import requests

def get_cutout(targetid, ra, dec, verbose=False):
    '''
    Extracts postage-stamp image from the Legacy Surveys viewer


    PARAMETERS
    ==========

    targetid : 

    ra : float
        Right ascension of target center [degrees]

    dec : float
        Declination of target center [degrees]

    verbose : boolean
        Determines whether or not to print status statements.  Default is False 
        (do not print).


    RETURNS
    =======

    img_name : string
        File name of postage stamp image extracted from Legacy Surveys imaging
    '''

    ############################################################################
    # Initialize returned image file name
    #---------------------------------------------------------------------------
    img_name = 'large_gals/{}.jpg'.format(targetid)
    ############################################################################


    ############################################################################
    # Extract image from Legacy Surveys viewer
    #---------------------------------------------------------------------------
    if os.path.exists(img_name):
        if verbose:
            print('{} exists.'.format(img_name))

    else:
        if verbose:
            print('Accessing {}'.format(img_name))

        img_url = 'https://www.legacysurvey.org/viewer/cutout.jpg?ra={}&dec={}&%22/pix=0.25&layer=ls-dr9&size=1000'.format(ra, dec)
        print(img_url)

        with open(img_name, 'wb') as handle:

            response = requests.get(img_url, stream=True)

            if not response.ok: 
                print(response)

            for block in response.iter_content(1024):

                if not block: 
                    break

                handle.write(block)
    ############################################################################

    return img_name






