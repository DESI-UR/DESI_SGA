

import matplotlib.pyplot as plt

#from marvin.tools.image import Image

import os
import gc

import sys
sys.path.insert(1, '/Users/kellydouglass/Documents/Research/Rotation_curves/RotationCurves/spirals/')
from DRP_rotation_curve_plottingFunctions import plot_rband_image



def plot_veldisp(veldisp, gal_ID, IMAGE_DIR=None, IMAGE_FORMAT='eps', ax=None):
    '''
    Creates a plot of the stellar velocity dispersion map


    Parameters:
    ===========

    veldisp : numpy array of shape (n,n)
        stellar velocity dispersion map

    gal_ID : string
        [MaNGA plate] - [MaNGA IFU]

    IMAGE_DIR : string
        Path of directory to store images

    IMAGE_FORMAT : string
        Format of saved image.  Default is eps

    ax : matplotlib.pyplot axis object
        Axes handle on which to create plot
    '''


    ###########################################################################
    if ax is None:
        fig, ax = plt.subplots()

    veldisp_im = ax.imshow( veldisp, origin='lower')

    cbar = plt.colorbar( veldisp_im, ax=ax)
    cbar.ax.set_ylabel('$\sigma_*$ [km/s]')

    ax.set_title( gal_ID + ' $\sigma_*$')
    ax.set_xlabel('spaxel')
    ax.set_ylabel('spaxel')
    ###########################################################################


    
    if IMAGE_DIR is not None:
        #######################################################################
        # Create output directory if it does not already exist
        #----------------------------------------------------------------------
        if not os.path.isdir( IMAGE_DIR + '/veldisp'):
            os.makedirs( IMAGE_DIR + '/veldisp')
        #######################################################################

        #######################################################################
        # Save figure
        #----------------------------------------------------------------------
        plt.savefig( IMAGE_DIR + '/veldisp/' + gal_ID + '_veldisp.' + IMAGE_FORMAT, 
                     format=IMAGE_FORMAT)
        #######################################################################

        #######################################################################
        # Figure cleanup
        #----------------------------------------------------------------------
        plt.cla()
        plt.clf()
        plt.close()
        del cbar, veldisp_im
        gc.collect()
        #######################################################################


###############################################################################
###############################################################################
###############################################################################


def plot_veldisp_norm(veldisp_norm, gal_ID, IMAGE_DIR=None, IMAGE_FORMAT='eps', ax=None):
    '''
    Creates a plot of the normalized stellar velocity dispersion map


    Parameters:
    ===========

    veldisp_norm : numpy array of shape (n,n)
        normalized stellar velocity dispersion map

    gal_ID : string
        [MaNGA plate] - [MaNGA IFU]

    IMAGE_DIR : string
        Path of directory to store images

    IMAGE_FORMAT : string
        Format of saved image.  Default is eps

    ax : matplotlib.pyplot axis object
        Axes handle on which to create plot
    '''


    ###########################################################################
    if ax is None:
        fig, ax = plt.subplots()

    veldisp_im = ax.imshow( veldisp_norm, origin='lower')

    cbar = plt.colorbar( veldisp_im, ax=ax)
    cbar.ax.set_ylabel('normalized $\sigma_*$')

    ax.set_title( gal_ID + ' $\sigma_*$')
    ax.set_xlabel('spaxel')
    ax.set_ylabel('spaxel')
    ###########################################################################


    
    if IMAGE_DIR is not None:
        #######################################################################
        # Create output directory if it does not already exist
        #----------------------------------------------------------------------
        if not os.path.isdir( IMAGE_DIR + '/veldisp_norm'):
            os.makedirs( IMAGE_DIR + '/veldisp_norm')
        #######################################################################

        #######################################################################
        # Save figure
        #----------------------------------------------------------------------
        plt.savefig( IMAGE_DIR + '/veldisp_norm/' + gal_ID + '_veldisp_norm.' + IMAGE_FORMAT, 
                     format=IMAGE_FORMAT)
        #######################################################################

        #######################################################################
        # Figure cleanup
        #----------------------------------------------------------------------
        plt.cla()
        plt.clf()
        plt.close()
        del cbar, veldisp_im
        gc.collect()
        #######################################################################


################################################################################
################################################################################
################################################################################



def plot_veldisp_r(x, v, gal_ID, IMAGE_DIR=None, IMAGE_FORMAT='eps', ax=None):
    '''
    Plot velocity dispersion as a function of distance from the center spaxel


    Parameters:
    ===========

    x : ndarray of shape (n,)
        distance from the center spaxel

    v : ndarray of shape (n,)
        (normalized) velocity dispersion at the given distance from the center 
        spaxel

    gal_ID : string
        MaNGA plate number - MaNGA fiberID number

    IMAGE_DIR : string
        Path of directory to store images

    IMAGE_FORMAT : string
        Format of saved image

    ax : matplotlib.pyplot figure axis object
        Axes handle on which to create plot

    '''


    if ax is None:
        fig, ax = plt.subplots( figsize=(5, 5))


    ############################################################################
    ax.set_title(gal_ID)
    ax.plot( x, v, '.')

    ax.tick_params( axis='both', direction='in')
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    ax.set_xlabel('Deprojected distance from the center spaxel [arcseconds]')
    ax.set_ylabel('Normalized $\sigma_*$')
    ############################################################################


    if IMAGE_DIR is not None:
        ########################################################################
        # Create output directory if it does not already exist
        #-----------------------------------------------------------------------
        if not os.path.isdir( IMAGE_DIR + '/veldisp_r'):
            os.makedirs( IMAGE_DIR + '/veldisp_r')
        ########################################################################

        ########################################################################
        # Save figure
        #-----------------------------------------------------------------------
        plt.savefig( IMAGE_DIR + "/veldisp_r/" + gal_ID + "_veldisp_r." + IMAGE_FORMAT,
                     format=IMAGE_FORMAT)
        ########################################################################

        ########################################################################
        # Figure cleanup
        #-----------------------------------------------------------------------
        plt.cla()
        plt.clf()
        plt.close()
        gc.collect()
        ########################################################################


###############################################################################
###############################################################################
###############################################################################


def plot_diagnostic_panel( gal_ID, 
                           r_band, 
                           masked_veldisp, 
                           masked_veldisp_norm, 
                           data_table, 
                           distance_map,
                           IMAGE_DIR=None, IMAGE_FORMAT='eps'):
    '''
    Plot a one by two paneled image containing the entire r-band array and the 
    masked stellar velocity dispersion array.


    Parameters:
    ===========

    gal_ID : string
        MaNGA plate number - MaNGA fiberID number

    r_band : numpy array of shape (n,n)
        r_band flux map

    masked_veldisp : numpy array of shape (n,n)
        Masked stellar velocity dispersion map

    masked_veldisp_norm : numpy array of shape (n,n)
        Masked stellar velocity dispersion map normalized by the central value

    data_table : astropy table of length n*n
        Data table with columns of stellar velocity dispersion and distance to 
        the center spaxel 

    IMAGE_DIR : string
        Path of directory to store images.  Default is None (does not save 
        figure)

    IMAGE_FORMAT : string
        Format of saved image.  Default is 'eps'
    '''


    panel_fig, ((image_panel, r_band_panel), 
                (mveldisp_panel, veldisp_r_panel)) = plt.subplots(2, 2)
    panel_fig.set_figheight( 10)
    panel_fig.set_figwidth( 10)
    plt.suptitle( gal_ID + " Diagnostic Panel", y=1.05, fontsize=16)

    '''
    image = Image(plateifu=gal_ID)
    image.plot(fig=panel_fig, ax=image_panel)
    '''
    '''
    dist_im = image_panel.imshow( distance_map, origin='lower')
    cbar = plt.colorbar( dist_im, ax=image_panel)
    cbar.ax.set_ylabel('$r_{depro}$ [arcsec]')
    image_panel.set_title( gal_ID + ' $r_{depro}$')
    image_panel.set_xlabel('spaxel')
    image_panel.set_ylabel('spaxel')
    '''
    plot_rband_image( r_band, gal_ID, ax=r_band_panel)

    plot_veldisp( masked_veldisp, gal_ID, ax=mveldisp_panel)

    plot_veldisp_r( data_table['r_deproj_arcsec'], data_table['vel_disp_norm'], gal_ID, 
                    ax=veldisp_r_panel)

    panel_fig.tight_layout()



    if IMAGE_DIR is not None:
        ########################################################################
        # Create output directory if it does not already exist
        #-----------------------------------------------------------------------
        if not os.path.isdir( IMAGE_DIR + '/diagnostic_panels'):
            os.makedirs( IMAGE_DIR + '/diagnostic_panels')
        ########################################################################

        ########################################################################
        # Save figure
        #-----------------------------------------------------------------------
        plt.savefig( IMAGE_DIR + "/diagnostic_panels/" + gal_ID + "_diagnostic_panel." + IMAGE_FORMAT,
                    format=IMAGE_FORMAT)
        ########################################################################

        ########################################################################
        # Figure cleanup
        #-----------------------------------------------------------------------
        plt.cla()
        plt.clf()
        plt.close( panel_fig)
        del panel_fig, image_panel, r_band_panel, mveldisp_panel, veldisp_r_panel
        gc.collect()
        ########################################################################