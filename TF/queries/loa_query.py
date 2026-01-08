import os
import psycopg2

from astropy import units as u
from astropy.table import Table, join, unique, vstack
from astropy.coordinates import SkyCoord, match_coordinates_sky, search_around_sky

from desitarget.targets import decode_targetid, encode_targetid, resolve
from desitarget.io import releasedict, release_to_photsys

from tqdm.notebook import tqdm_notebook

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


from astropy.table import Table, vstack, hstack

import psycopg2

from tqdm.notebook import tqdm_notebook

import time

#Remove the condition that targetid>0 because there are a few entries in the iron table that had a targetid<0.

def get_tf_targets_modified(redux, use_cached=False, verbose=False):
    """Get TF targets from the DESI observations DB for a given spectroscopic reduction.
    
    Parameters
    ----------
    redux : str
        Spectroscopic reduction. E.g., 'everest', 'fuji', ...
    use_cached : bool
        Use cached data rather than re-running the query.
    
    Returns
    -------
    
    tf_targets : Table
        Table of Tully-Fisher observations.
    """
    tf_targets = None

    if os.path.exists(f'tf_targets_{redux}.fits') and use_cached:
        tf_targets = Table.read('cache/tf_targets_{redux}.fits')
    else:
        try:
            db = psycopg2.connect(host='desidb-rr.lbl.gov', database='desidb', user='desi', password='')
            # db = psycopg2.connect(host='decatdb.lbl.gov', database='desidb', user='desi', password='')
            cursor = db.cursor()

            query = f"""SELECT rdx.targetid, rdx.target_ra, rdx.target_dec, zd.z, zd.zerr, zd.spectype, zd.deltachi2, zd.zwarn, pv.pvtype, pv.sga_id
                   FROM {redux}.healpix_fibermap as rdx, static.pv as pv, {redux}.healpix_redshifts as zd
                   WHERE q3c_join(rdx.target_ra, rdx.target_dec, pv.ra, pv.dec, 1./3600.) 
                         AND zd.targetid = rdx.targetid;"""
                         # AND pv.sga_id IS NOT NULL AND (pv.pvtype LIKE 'TFT' or pv.pvtype LIKE 'EXT' or pv.pvtype LIKE 'SGA');"""
            
            if verbose:
                print(query)

            cursor.execute(query)
            rows = cursor.fetchall()
            tf_targets = Table(list(map(list, zip(*rows))),
                               names=['TARGETID', 'TARGET_RA', 'TARGET_DEC', 'Z', 'ZERR', 'SPECTYPE', 'DELTACHI2', 'ZWARN', 'PVTYPE', 'SGA_ID'])

            #- Select only targets with SGA IDs and PV types matching SGA, EXT, and TFT
            select = (tf_targets['SGA_ID'] != None) & \
                     ((tf_targets['PVTYPE'] == 'TFT') | \
                      (tf_targets['PVTYPE'] == 'EXT') | \
                      (tf_targets['PVTYPE'] == 'SGA'))
            tf_targets = tf_targets[select]

            #- Use TARGETID to extract the photometric system used during targeting
            _, _, releases, _, _, _ = decode_targetid(tf_targets['TARGETID'].value)

            photsys = []
            for i, release in enumerate(releases):
                ps = None

                if release in releasedict:
                    ps = release_to_photsys([release])[0].decode('utf-8')
                else:
                    #- Fall-through case: not all SGA center observations are in the main survey.
                    #  In this case, select 'N' or 'S' based on the SGA object's position.
                    ra  = tf_targets['TARGET_RA'][i]
                    dec = tf_targets['TARGET_DEC'][i]
                    c = SkyCoord(ra=ra, dec=dec, unit='degree')

                    #- N: in galactic northern hemisphere and with dec > 32.375. Else, S.
                    isnorth = (c.galactic.b > 0) & (dec > 32.375)
                    ps = 'N' if isnorth else 'S'

                photsys.append(ps)

            #- Complain if the photsys table doesn't match the size of the Vrot table.
            if len(photsys) != len(tf_targets):
                print(f'photsys array of len {len(photsys)} != targets array of len {len(tf_targets)}')

            tf_targets['PHOTSYS'] = photsys

            # tf_targets.write(f'cache/tf_targets_{redux}.fits', overwrite=True)

        except Exception as error:
            print(error)
        finally:
            if db is not None:
                db.close()

    return tf_targets

loa = get_tf_targets_modified('loa', verbose=True)

loa['SGA_ID'] = [str(item) for item in loa['SGA_ID']]
loa.write('../Y3/desi_pv_tf_loa_healpix.fits', overwrite=True)
print(len(loa))