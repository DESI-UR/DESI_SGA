import numpy as np

from scipy.stats import binned_statistic

import matplotlib.colors as mc

import colorsys


################################################################################
#-------------------------------------------------------------------------------
def profile_histogram(x, y, xbins, 
                      yerr=None, 
                      weights=None, 
                      median=False, 
                      weighted=False):
    """
    Compute a profile histogram from scattered data.
    
    Parameters
    ----------
    x : list or ndarray
        Ordinates (independent variable).
    y : list or ndarray
        Coordinates (dependent variable).
    xbins : list or ndarray
        Bin edges for the independent variable.
    yerr : list or ndarray
        Uncertainties on the dependent variable. Assumed independent.
    weights : list or ndarray
        If not None (and weighted=True), will use this instead of yerr to weight 
        the summary statistics.
    median : bool
        If true, compute median as central value; else, the (weighted) mean.
    weighted : bool
        Weight the summary statistics, either by the uncertainty in y or the 
        provided weights.
        
    Returns
    -------
    N : ndarray
        Unweighted counts per bin.
    h : ndarray
        Summary statistic (mean or median) of independent variable per bin.
    e : ndarray
        Uncertainty on the summary statistic per bin.
    """
    
    N = binned_statistic(x, y, bins=xbins, statistic='count').statistic

    if weighted:
        if (yerr is None) and (weights is None):
            raise ValueError('need to define either yerr or weights if using weighted fit.')

        if weights is None:
            # weight based on yerr
            w = 1/yerr**2
        else:
            w = weights
        W, H, E = binned_statistic(x, [w, w*y, w*y**2], 
                                   bins=xbins, statistic='sum').statistic
        h = H/W
        e = 1/np.sqrt(W)
    else:
        mean, mean2 = binned_statistic(x, [y, y**2], 
                                       bins=xbins, statistic='mean').statistic
        h = mean
        e = np.sqrt((mean2 - mean**2) / (N - 1))

    if median:
        h = binned_statistic(x, y, bins=xbins, statistic='median').statistic
    
    return N, h, e
################################################################################
################################################################################



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