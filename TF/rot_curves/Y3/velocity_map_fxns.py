import numpy as np

def sin_i(ba, q0):
    """get inclination angle:
    Parameters
    -----------
    ba: axis ratio

    Returns
    -------
    sin of inclination angle
    """
    
    cos2 = (ba**2-q0**2)/(1-q0**2)
    if cos2 < 0:
        cos2 = 0
    return(np.sqrt(1-cos2))

#quality fit, input is table to get from, and index into table
def criteria(table, idx, zwarn, deltachi2):
    """Quality cut:
    Parameters
    -----------
    Table: table to get from

    Index: index into table

    Zwarn: value zwarn should equal for cut

    deltachi2: minimum value for cut
    """
    return (table['ZWARN'][idx] == zwarn) & (table['DELTACHI2'][idx] > deltachi2)

def criteria_sym(table, zwarn, deltachi2)
    """Quality cut w/o index:
    Parameters
    -----------
    Table: table to get from

    Index: index into table

    Zwarn: value zwarn should equal for cut

    deltachi2: minimum value for cut
    """
    return (table['ZWARN'] == zwarn) & (table['DELTACHI2'] > deltachi2)