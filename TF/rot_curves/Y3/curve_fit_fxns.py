from scipy.optimize import minimize

import numpy as np

def v_rot(r, v_max, r_turn):
    """ 
    Parameters
    ----------
    r: Radius
    
    v_max: Maximum velocity

    r_turn: radius at which velocity plateaus

    a: sharpness of curve
    
    """
    return (v_max) * np.tanh(r/r_turn)


def chi2(params, v_data, v_err, r):
    """ Calculate chi2 for rotation curve fitting
    Parameters
    ----------
    params: v_max, r_turn, a

    v_data: velocity of data value

    v_err: Error of velocity values

    r: radius
    """
    v_max, r_turn = params
    v_model = (v_max) * np.tanh(r/r_turn)
    return np.sum(((v_data-v_model)/v_err)**2)

def chi2_reduced(chi2, num_data, num_params):
    """ Calculate reduced chi2
    Parameters
    ----------
    chi 2: chi 2 value

    num_data: number of data points

    num_params: number of parameters

    Output
    -----------
    Reduced chi2 value (chi 2 / (number of data points - number of free parameters))
    """
    dof = num_data - num_params
    return chi2/dof

def curve_fit(func, x0, args, bounds):
    """ Minimize function using scipy minimize

    Parameters
    -----------
    func: function to minimize

    x0: initial guess

    args: parameters

    bounds: bounds
    """
    return minimize(fun = func, 
             x0 = x0, 
             args = args, 
             bounds = bounds, 
             method = 'Powell')