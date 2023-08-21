import numpy as np

import matplotlib.colors as mc

import colorsys



################################################################################
#-------------------------------------------------------------------------------
def firstdigit(n):
    """
    Return the first digit of a number.
    
    Parameters
    ----------
    n : int, float, or ndarray
        Number or list of numbers.
    
    Returns
    -------
    digit : int
        First digit of the number.
    """
    
    digit = np.trunc(n * 10**(-np.trunc(np.log10(n)))).astype(int)
    
    return digit
################################################################################
################################################################################




################################################################################
#-------------------------------------------------------------------------------
def adjust_lightness(color, amount=0.5):
    '''
    Change the shade of the color.
    
    
    PARAMETERS
    ==========
    
    color : string or length-3 tuple
        Matplotlib color name OR rgb tuple
    
    amount : float
        Amount by which to lighten or darken the color.  Values < 1 will 
        return a darker shade, while values > 1 will return a lighter shade.  
        Default is 0.5.
    '''
    
    try:
        c = mc.cnames[color]
    except:
        c = color
    
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    
    return colorsys.hls_to_rgb(c[0], max(0, min(1, amount*c[1])), c[2])
################################################################################
################################################################################