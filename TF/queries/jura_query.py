from astropy.table import Table, vstack, hstack

import psycopg2

from tqdm.notebook import tqdm_notebook

import time


def match_targets(pvtargtab, redux='daily', search='healpix'):
    """Match PV targets against the redshift DB for a particular spectroscopic reduction.
    
    Parameters
    ----------
    pvtargtab : astropy.Table
        Table of PV target info. Specifically need the RA, DEC, PVTYPE, and SGA_ID fields.
    redux : str
        Spectroscopic reduction: e.g., 'daily', 'everest', 'fuji', 'guadalupe', ...
    search : str
        'healpix' to search the HEALPix tables, 'tiles' to search the tiles tables.
        
    Returns
    -------
    desi_targets : astropy.Table
        Joined table of DESI redshifts and PV targets for all matches.
    """
    # Accumulate data in this table.
    desi_targets = None
        
    try:
        db = psycopg2.connect(host='decatdb.lbl.gov', database='desidb', user='desi', password=) # Use the usual DESI password
        cursor = db.cursor()
        # cursor.execute('SET search_path TO da, public;')

        # Loop over all TNS alerts and perform a coordinate match with DESI observations.
        N = len(pvtargtab)
        n = 0
        with tqdm_notebook(total=N) as progress_bar:

            for i, obj in enumerate(pvtargtab):
                ra, dec = obj['RA'], obj['DEC']

                # Enable search in HEALPix tables.
                if search == 'healpix':
                    query = 'SELECT f.targetid,f.target_ra,f.target_dec,h.healpix,h.survey,r.z,r.zerr,r.zwarn,r.deltachi2,h.filename\n' \
                            f'FROM {redux}.healpix_fibermap f\n' \
                            f'INNER JOIN {redux}.healpix h ON f.healpix_id=h.id\n' \
                            f'INNER JOIN {redux}.healpix_redshifts r ON r.healpix_id=h.id AND r.targetid=f.targetid\n' \
                            f'WHERE q3c_radial_query( f.target_ra, f.target_dec, {ra}, {dec}, 1./3600. );'
                    
                    colnames = ['TARGETID', 'TARGET_RA', 'TARGET_DEC', 'HEALPIX', 'SURVEY', 'Z', 'ZERR', 'ZWARN', 'DELTACHI2', 'FILENAME']
                # Enable search in tiles tables.
                elif search == 'tiles':
                    query = 'SELECT f.targetid,f.target_ra,f.target_dec,c.tileid,c.night,r.z,r.zerr,r.zwarn,r.deltachi2,c.filename\n' \
                            f'FROM {redux}.tiles_fibermap f\n' \
                            f'INNER JOIN {redux}.cumulative_tiles c ON f.cumultile_id=c.id\n' \
                            f'INNER JOIN {redux}.tiles_redshifts r ON r.cumultile_id=c.id AND r.targetid=f.targetid\n' \
                            f'WHERE q3c_radial_query( f.target_ra, f.target_dec, {ra}, {dec}, 1./3600. );'
                    colnames = ['TARGETID', 'TARGET_RA', 'TARGET_DEC', 'TILEID', 'NIGHT', 'Z', 'ZERR', 'ZWARN', 'DELTACHI2', 'FILENAME']
                else:
                    raise ValueError(f'Search {search} not recognized; use "healpix" or "tiles."')

                cursor.execute(query)
                rows = cursor.fetchall()

                if rows:
                    # Convert postgresql row output to an astropy Table.
                    data = Table(list(map(list, zip(*rows))),
                                 names=colnames)

                    # hstack the postgresql rows with the PV target info.
                    # The following vstack loop ensures every row gets a match.
                    pv_data = obj
                    if len(data) > 1:
                        for j in range(1, len(data)):
                            pv_data = vstack([pv_data, obj])
                    data = hstack([data, pv_data['PVTYPE', 'SGA_ID', 'RA', 'DEC']])

                    # Accumulate matched targets.
                    if desi_targets is None:
                        desi_targets = data
                    else:
                        desi_targets = vstack([desi_targets, data], join_type='outer')

                if (i+1) % 50 == 0:
                    progress_bar.update(50)
                    n += 50

            if n < N:
                progress_bar.update(N - n)

    except (Exception, psycopg2.Error) as error:
        print(error)
    finally:
        if db is not None:
            db.close()
            
    return desi_targets


if __name__ == "__main__":
    
    start_time = time.time()
    
    # Read in target lists
    pv_ext = Table.read('/global/cfs/projectdirs/desi/science/td/pv/desi_pv/savepath_dr9_corr/pv_ext.fits', hdu=1)
    pv_sga = Table.read('/global/cfs/projectdirs/desi/science/td/pv/desi_pv/savepath_dr9_corr/pv_sga.fits', hdu=1)
    pv_tf = Table.read('/global/cfs/projectdirs/desi/science/td/pv/desi_pv/savepath_dr9_corr/pv_tf.fits', hdu=1)
    print('Read in file')
    
    
    # Run SQL query for each of the target lists
    pv_ext_jura = match_targets(pv_ext, redux='jura')
    pv_sga_jura = match_targets(pv_sga, redux='jura')
    pv_tf_jura = match_targets(pv_tf, redux='jura')
    print('ran queries')
    
    print('Time to execute query:', time.time() - start_time)

    # Combine all target observations into a single table
    pv_jura = vstack([pv_ext_jura, pv_sga_jura, pv_tf_jura])
    
    # Save full observation table to disk
    pv_ext_jura.write('desi_pv_tf_jura_healpix_test2.fits', overwrite=True)

