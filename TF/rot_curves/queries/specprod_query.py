import numpy as np
from astropy.table import Table, vstack, join, unique

from astropy.io import fits
from astropy.coordinates import SkyCoord, match_coordinates_sky, search_around_sky
from astropy import units as u
from tqdm import tqdm

import os
import psycopg2
from psycopg2.extras import execute_values
from psycopg2.extensions import register_adapter, AsIs


import datetime

from desitarget.targets import decode_targetid, encode_targetid, resolve
from desitarget.io import releasedict, release_to_photsys

def adapt_numpy_int64(numpy_int64):
    return AsIs(numpy_int64)

def get_rot_curve_targets(redux, sga_params, table_out=None, verbose=False):
    """Get rotation curve targets from the DESI production DB for a given spectroscopic reduction.
        Returns targets in a square around SGA center of size 1.01 D26 x 1.01 D26
    
    Parameters
    ----------
    redux : str
        Spectroscopic reduction. can be 'fuji', 'guadalupe', 'iron', or 'loa'

    sga_params: list of tuples
        list of tuples of the form [(sga_id, ra, dec, major, pa)]

        sga_id : int
            SGA ID from pv table
        
        ra : float
            galaxy RA from SGA [deg]
        dec : float
            galaxy dec from SGA [deg]
        major : float
            semi-major axis (D26/2) from SGA [deg]
        ratio : float
            axis ratio (BA) from SGA
        pa : float
            position angle east-of-north from SGA [deg]
    
    Returns
    -------
    targets : Table
        Table of observations.
    """
    


    targets = None

    try:
        register_adapter(np.int64, adapt_numpy_int64)
        db = psycopg2.connect(host='specprod-db.desi.lbl.gov', database='desi', user='desi', password='')
        cursor = db.cursor()

        # create temporary pv table

        cursor.execute("""CREATE TEMP TABLE pv(
                                sgaid BIGINT,
                                ra FLOAT,
                                dec FLOAT,
                                major FLOAT,
                                ratio FLOAT,
                                pa FLOAT); 
                                """)
        print('temp table created')
        # put sga parameters into the temporary table
        execute_values(cursor, """INSERT INTO pv(sgaid, ra, dec, major, ratio, pa) VALUES %s """, sga_params)

        print('temp table filled')
        # create index for pv table
        
        cursor.execute("""CREATE INDEX ON pv (q3c_ang2ipix(ra, dec));""")
        cursor.execute("""ANALYZE pv;""")

        print('temp table indexed')    
        # get all targets within sga ellipses
        

        # query = f"""SELECT f.targetid, f.tileid, pv.sgaid, f.target_ra, f.target_dec, zp.healpix, zp.survey, zp.program,
        #                 zp.z, zp.zerr, zp.zwarn, zp.chi2, zp.deltachi2, zp.mean_fiber_ra, zp.mean_fiber_dec,
        #                 zp.std_fiber_ra, zp.std_fiber_dec, zp.spectype
        #                 FROM pv, iron.fiberassign AS f
        #                 JOIN iron.zpix AS zp ON f.targetid = zp.targetid
        #                 WHERE q3c_ellipse_join(f.target_ra, f.target_dec, pv.ra, pv.dec, pv.major, pv.ratio, pv.pa);"""

        

        query = f"""SELECT f.targetid, f.tileid, pv.sgaid, f.target_ra, f.target_dec, zp.healpix, zp.survey, zp.program,
                        zp.z, zp.zerr, zp.zwarn, zp.chi2, zp.deltachi2, zp.mean_fiber_ra, zp.mean_fiber_dec,
                        zp.std_fiber_ra, zp.std_fiber_dec, zp.spectype
                        FROM pv, iron.fiberassign AS f
                        JOIN iron.zpix AS zp ON f.targetid = zp.targetid
                        WHERE f.target_ra < (pv.ra + pv.major) AND f.target_ra > (pv.ra - pv.major) AND f.target_dec < (pv.dec + pv.major) AND 
                        f.target_dec > (pv.dec - pv.major);"""
        
        if verbose:
            print(query)

        cursor.execute(query)
        rows = cursor.fetchall()

        if len(rows)>0:
            targets = Table(list(map(list, zip(*rows))), names=['TARGETID', 'TILEID', 'SGA_ID', 'TARGET_RA', 
                                                                'TARGET_DEC', 'HEALPIX', 'SURVEY', 'PROGRAM',
                          'Z', 'ZERR', 'ZWARN', 'CHI2', 'DELTACHI2', 'MEAN_FIBER_RA', 'MEAN_FIBER_DEC',
                         'STD_FIBER_RA', 'STD_FIBER_DEC', 'SPECTYPE'],
            dtype=['int64', 'int64', 'int64', 'float64', 'float64', 'int64', 'S4', 'S6',
                 'float64', 'float64', 'int64', 'float64', 'float64', 'float64', 'float64',
                     'float64', 'float64', 'S6'])


    except Exception as error:
        print(error)
    finally:
        if db is not None:
            db.close()
    
    if table_out is not None:
        targets.write(table_out, format='fits', overwrite=True)