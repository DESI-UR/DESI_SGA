

import matplotlib.pyplot as plt

#from marvin.tools.image import Image

import os
import gc

import sys
sys.path.insert(1, '/Users/kellydouglass/Documents/Research/Rotation_curves/RotationCurves/spirals/')
from DRP_rotation_curve_plottingFunctions import plot_rband_image



def plot_Ha(Ha_flux, gal_ID, IMAGE_DIR=None, IMAGE_FORMAT='eps', ax=None):
    '''
    Creates a plot of the H-alpha flux map


    Parameters:
    ===========

    Ha_flux : numpy array of shape (n,n)
        H-alpha flux map

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

    Ha_im = ax.imshow( Ha_flux, origin='lower')

    cbar = plt.colorbar( Ha_im, ax=ax)
    cbar.ax.set_ylabel(r'H$\alpha$ flux [10$^{-17}$ erg/s/cm$^2$')

    ax.set_title( gal_ID + r' H$\alpha$ flux')
    ax.set_xlabel('spaxel')
    ax.set_ylabel('spaxel')
    ###########################################################################


    
    if IMAGE_DIR is not None:
        #######################################################################
        # Create output directory if it does not already exist
        #----------------------------------------------------------------------
        if not os.path.isdir( IMAGE_DIR + '/Ha_flux'):
            os.makedirs( IMAGE_DIR + '/Ha_flux')
        #######################################################################

        #######################################################################
        # Save figure
        #----------------------------------------------------------------------
        plt.savefig( IMAGE_DIR + '/Ha_flux/' + gal_ID + '_HaFlux.' + IMAGE_FORMAT, 
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


def plot_Ha_norm(Ha_flux_norm, gal_ID, IMAGE_DIR=None, IMAGE_FORMAT='eps', ax=None):
    '''
    Creates a plot of the normalized H-alpha flux map


    Parameters:
    ===========

    Ha_flux_norm : numpy array of shape (n,n)
        normalized H-alpha flux map

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

    Ha_im = ax.imshow( Ha_flux_norm, origin='lower')

    cbar = plt.colorbar( Ha_im, ax=ax)
    cbar.ax.set_ylabel(r'normalized H$\alpha$ flux')

    ax.set_title( gal_ID + r' H$\alpha$ flux')
    ax.set_xlabel('spaxel')
    ax.set_ylabel('spaxel')
    ###########################################################################


    
    if IMAGE_DIR is not None:
        #######################################################################
        # Create output directory if it does not already exist
        #----------------------------------------------------------------------
        if not os.path.isdir( IMAGE_DIR + '/HaFlux_norm'):
            os.makedirs( IMAGE_DIR + '/HaFlux_norm')
        #######################################################################

        #######################################################################
        # Save figure
        #----------------------------------------------------------------------
        plt.savefig( IMAGE_DIR + '/HaFlux_norm/' + gal_ID + '_HaFlux_norm.' + IMAGE_FORMAT, 
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



def plot_Ha_r(x, f, gal_ID, IMAGE_DIR=None, IMAGE_FORMAT='eps', ax=None):
    '''
    Plot H-alpha flux as a function of distance from the center spaxel


    Parameters:
    ===========

    x : ndarray of shape (n,)
        distance from the center spaxel

    f : ndarray of shape (n,)
        (normalized) H-alpha flux at the given distance from the center spaxel

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
    ax.plot( x, f, '.')

    ax.tick_params( axis='both', direction='in')
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    ax.set_xlabel('Deprojected distance from the center spaxel [arcseconds]')
    #ax.set_ylabel(r'Normalized H$\alpha$ flux')
    ax.set_ylabel(r'H$\alpha$ flux [10$^{-17}$ erg/s/cm$^2$]')
    ############################################################################


    if IMAGE_DIR is not None:
        ########################################################################
        # Create output directory if it does not already exist
        #-----------------------------------------------------------------------
        if not os.path.isdir( IMAGE_DIR + '/HaFlux_r'):
            os.makedirs( IMAGE_DIR + '/HaFlux_r')
        ########################################################################

        ########################################################################
        # Save figure
        #-----------------------------------------------------------------------
        plt.savefig( IMAGE_DIR + "/HaFlux_r/" + gal_ID + "_HaFlux_r." + IMAGE_FORMAT,
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
                           masked_Ha_flux, 
                           masked_Ha_flux_norm, 
                           data_table, 
                           distance_map,
                           IMAGE_DIR=None, IMAGE_FORMAT='eps'):
    '''
    Plot a one by two paneled image containing the entire r-band array and the 
    masked H-alpha flux array.


    Parameters:
    ===========

    gal_ID : string
        MaNGA plate number - MaNGA fiberID number

    r_band : numpy array of shape (n,n)
        r_band flux map

    masked_Ha_flux : numpy array of shape (n,n)
        Masked H-alpha flux map

    masked_Ha_flux_norm : numpy array of shape (n,n)
        Masked H-alpha flux map normalized by the central value

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
                (mHa_panel, Ha_r_panel)) = plt.subplots(2, 2)
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

    plot_Ha( masked_Ha_flux, gal_ID, ax=mHa_panel)

    plot_Ha_r( data_table['r_deproj_arcsec'], 
               data_table['Ha_flux'],
               #data_table['Ha_flux_norm'], 
               gal_ID, 
               ax=Ha_r_panel)

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
        plt.savefig( IMAGE_DIR + "/diagnostic_panels/" + gal_ID + "_Ha_diagnostic_panel." + IMAGE_FORMAT,
                    format=IMAGE_FORMAT)
        ########################################################################

        ########################################################################
        # Figure cleanup
        #-----------------------------------------------------------------------
        plt.cla()
        plt.clf()
        plt.close( panel_fig)
        del panel_fig, image_panel, r_band_panel, mHa_panel, Ha_r_panel
        gc.collect()
        ########################################################################