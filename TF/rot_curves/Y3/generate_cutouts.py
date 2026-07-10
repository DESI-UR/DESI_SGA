from get_cutouts import get_cutout

from astropy.table import Table
from astropy.visualization.wcsaxes import SphericalCircle
from astropy import units as u

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from astropy.coordinates import SkyCoord, match_coordinates_sky
import matplotlib.patches as patches

from tqdm import tqdm

target_galaxies = Table.read('/pscratch/sd/d/dbustos/rot_curves/loa_targs.fits')

SGA = Table.read('/global/cfs/projectdirs/cosmo/data/sga/2020/SGA-2020.fits', 'ELLIPSE')

SGA_dict = {}
for i in range(len(SGA)):
    SGA_dict[SGA['SGA_ID'][i]] = i

fig1 = plt.figure(figsize=(7,5))
colors = plt.get_cmap('Set2')
fiber_radius = (107./70) * u.arcsec
dr = 10

# for each rotation curve galaxy, grab cut out, draw fibers on image, and save image
for sga_id in tqdm(np.unique(target_galaxies['SGA_ID'])[:10]):
    
    targ_list = target_galaxies[target_galaxies['SGA_ID']==sga_id]
    
    ra, dec = float(SGA['RA'][SGA_dict[sga_id]]), float(SGA['DEC'][SGA_dict[sga_id]])
    
    # D26 in arcmin
    d26 = SGA['D26'][SGA_dict[sga_id]]

    pix = int(2 * d26*60/0.262)

    if (pix < 2500):
        npix = np.minimum(pix,1024)

    elif (pix > 2500):
        npix = np.minimum(pix,3000)

    img_file, wcs = get_cutout(sga_id, ra, dec, dr=dr, dir='/pscratch/sd/d/dbustos/loa_cutouts/cutouts/',size=npix)
    img = mpl.image.imread(img_file)

    ax = fig1.add_subplot(111, projection=wcs)
    ax.imshow(np.flip(img, axis=0))
    ax.set(xlabel='ra', ylabel='dec')
    ax.text(int(0.02*npix), int(0.85*npix), 'SGA_ID: {}'.format(sga_id), fontsize=9, color='yellow')
    overlay = ax.get_coords_overlay('icrs')
    overlay.grid(color='white', ls='dotted');

    # Add the location of the DESI fibers.
    # SDSS fibers are 2" diameter, DESI is 107 um with 70 um/" plate scale.
    r1 = SphericalCircle((ra * u.deg, dec * u.deg), fiber_radius,
                         edgecolor='black', facecolor='none', alpha=0.8, lw=3,
                         transform=ax.get_transform('icrs'))
    r2 = SphericalCircle((ra * u.deg, dec * u.deg), fiber_radius,
                         edgecolor='red', facecolor='none', alpha=0.8, lw=2,
                         transform=ax.get_transform('icrs'))
    ax.add_patch(r1)
    ax.add_patch(r2)

    num_targs = len(targ_list)
    legend_handles = [None] * num_targs

    for i, targ in enumerate(targ_list):
        ra, dec = targ['TARGET_RA'], targ['TARGET_DEC']
        
        edgecolor2 = colors(i)

        # Add the location of the DESI fibers.
        # SDSS fibers are 2" diameter, DESI is 107 um with 70 um/" plate scale.
        r1 = SphericalCircle((ra * u.deg, dec * u.deg), fiber_radius,
                             edgecolor='black', facecolor='none', alpha=1, lw=3,
                             transform=ax.get_transform('icrs'))
        r2 = SphericalCircle((ra * u.deg, dec * u.deg), fiber_radius,
                             edgecolor=edgecolor2, facecolor='none', alpha=1, lw=2,
                             transform=ax.get_transform('icrs'))
        ax.add_patch(r1)
        ax.add_patch(r2)
        
        #ax.text(ra, dec, str(targ['TARGETID']), transform=ax.get_transform('icrs'), color='white', fontsize=6)
        legend_handles[i] = patches.Circle((0,0), facecolor = edgecolor2, ec = edgecolor2, lw = 2, label =  str(targ['TARGETID']))
    
    ax.legend(handles = legend_handles, title = 'SGA ID: {}'.format(sga_id), loc = 'upper right', bbox_to_anchor = (2, 1))
    fig1.subplots_adjust(top=0.85, right=0.85, bottom=0.15, left=0.15)
    
    fig1.savefig('/pscratch/sd/d/dbustos/loa_cutouts/cutouts_fibers/' + 'cutouts_sga_{}.png'.format(sga_id), bbox_inches = 'tight', dpi=120)
    
    fig1.clear()