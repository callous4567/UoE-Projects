# Function to iterate over the provided list
from astroquery.gaia import Gaia
import pickle
import numpy as np
import astropy.units as u

import windows_directories_new

username = "sstrasza"
password = "CENT@rino4657898"
Gaia.login(user=username,
           password=password)

def deltables(_remove_names):
    """
    Remove all these user tables with a try/except.
    :param _remove_names: list of names
    """
    for i in _remove_names:
        try:
            Gaia.delete_user_table(table_name=i, force_removal=True)
        except Exception as e:
            print(e)

def deljobs(_remove_jobs):
    """
    Removes *all* these job ids.
    :param _remove_jobs: array of job integers (unique to our username)
    """
    try:
        Gaia.remove_jobs(list(_remove_jobs))
    except Exception as e:
        print(e)

def clean_jobstables():

    """
    Quick helper to remove all jobs/tables if you're only dealing with a few.
    Split with multiprocessing if anticipating many (>10) tables/jobs for removal.
    :return: True or False if successful.
    """
    try:

        # Get all the jobs and delete them
        remove_jobs = Gaia.list_async_jobs()
        remove_jobs = [d.jobid for d in remove_jobs]
        deljobs(remove_jobs)

        # Get list of all tables on our Gaia, and prune to useful ones (that match our own upload, specifically) for removal
        tables = Gaia.load_tables(only_names=True, include_shared_tables=False)
        remove_names = []
        for table in tables:
            # split by user (to get rid of common tables, i.e. SDSS or DR3/etc)
            namesplit = table.name.split("user_" + username)
            if len(namesplit) >= 2:
                remove_names.append(table.name)
        remove_names = np.array(remove_names, dtype=str)
        deltables(remove_names)


        return True
    except:
        return False

def upload_table(_table, _name):

    """

    Upload an astropy table to Gaia via ADQL and return the fully qualified identifier (user.table)
    :param _table: astropy table
    :param _name: unique name for table
    :return: table_id (identifier for the table = username._name)
    """

    Gaia.upload_table(upload_resource=_table, table_name=_name, table_description="For the thesis! Hallelujah!")
    table_id = "user_" + username + "." + _name

    return table_id

def obtain_dr3ids(_table_id, _match_table_name, dr2_column="dr2_source_id_local", magdif=0.1, angdist=1):

    """

    Will attempt to match DR3 source to the provided dr2_column in provided _table_id (remote.) Returns the table_id
    of the match, alongside deleting the job once complete/etc (will wait on job to avoid spamming external server.)

    **Deletion is hit-or-miss: recommend complete purge after all processes complete. **

    ***Relevant column for match to dr2 is "dr3_source_id"***

    :param _table_id: (stored by server) table id to crossmatch for DR3
    :param _match_table_name: The name (not id) for the job- use an rng or something.
    :param dr2_column: Column to use for DR2 ID. Default dr2_source_id_local
    :param magdif: maximum magnitude difference admissible for match. Default 0.1 mag
    :param angdist: maximum angular difference (mas) admissible for match. Default 1 mas
    :return: table_id of match
    """

    # Set up the job/query
    """
    dr2_source_id_local (equivalent to dr2_source_id)
    dr2_source_id (match found in DR2)
    dr3_source_id (the match to dr2_source_id from gaiadr3.dr2_neighbourhood match)
    """
    job = Gaia.launch_job_async(query=("""
                                          SELECT * 
                                          FROM {0}
                                          JOIN gaiadr3.dr2_neighbourhood 
                                              ON gaiadr3.dr2_neighbourhood.dr2_source_id = {0}.{1}
                                          WHERE ABS(magnitude_difference) < {2} 
                                          AND ABS(angular_distance) < {3}
                                          """).format(_table_id, dr2_column, magdif, angdist),
                                name=_match_table_name)
    Gaia.upload_table_from_job(job)
    Gaia.delete_user_table(table_name=_table_id, force_removal=True)
    Gaia.remove_jobs([job.jobid])
    job_table_id = "user_" + username + ".t" + str(job.jobid)

    return job_table_id

default_columns = ["dr2_source_id", "dr3_source_id", "ra", "dec",
                                              "parallax", "parallax_error", "pmra",
                                              "pmra_error", "pmdec", "pmdec_error", "radial_velocity",
                                              "radial_velocity_error"]
def obtain_download_data(_table_id, _source_join_name, _data_table_name, _savepath,
                         dr3_columns=default_columns, max_parallax=0.2, min_error_ratio=3):

    """
    Will take the provided table ID with dr3_source_id column, and obtain dr3_columns (list- see "default_columns" in
    gaia_utils_standalone, subject to max_parallax and min_error_ratio. Will download table and save it to savepath.
    Will return whether the job was successful (True/False bool.)

    :param _table_id: id of the table from obtain_dr3ids (has dr3 source ids)
    :param _source_join_name: misc name for job (get all dr3 columns)
    :param _data_table_name: misc name for job (get dr3_columns columns)
    :param _savepath: local save path for the obtained table
    :param dr3_columns: list [1,2,3...] of columns wanted, i.e. [pmra, pmdec, ...]
    :param max_parallax: maximum value of parallax to admit in mas. (1/parallax) in arcseconds = distance in parsecs
    :param min_error_ratio: minimum ratio of parallax to error: factor of three is minimum recommended.
    :return: bool on whether process was successful (True/False)
    """

    try:
        # Join our DR3_ID table with the main gaia_source table
        job = Gaia.launch_job_async(query=("""
                                              SELECT * 
                                              FROM {0}
                                              JOIN gaiadr3.gaia_source AS gaiadr3
                                                  ON gaiadr3.source_id = {0}.dr3_source_id
                                                  
                                              WHERE ABS(parallax) < {1}
                                              AND ABS(parallax_over_error) > {2} 
                                              """).format(_table_id, max_parallax, min_error_ratio),
                                    name=_source_join_name)
        Gaia.upload_table_from_job(job)
        Gaia.delete_user_table(table_name=_table_id)
        Gaia.remove_jobs([job.jobid])
        joined_id = "user_" + username + ".t" + job.jobid

        # Set up the string for selection col_1, col_2, col_3... ,col_n
        selection_string = dr3_columns[0]
        for i in range(1,len(dr3_columns)):
            selection_string += (",{}").format(dr3_columns[i])

        # Grab the columns we're interested in
        job = Gaia.launch_job_async(query=("""
                                              SELECT {0}
                                              FROM {1}
                                              """).format(selection_string, joined_id),
                                    name=_data_table_name)
        results = job.get_results()
        Gaia.remove_jobs([job.jobid])
        Gaia.delete_user_table(table_name=joined_id)

        with open(_savepath, "wb") as f:
            pickle.dump(file=f, obj=results)

        return True

    except Exception as e:

        with open(_savepath, "w") as f:
            f.write(str(e))

        return False

def up_obt_download_save(_table, _subtable_index, _final_savepath):

    """

    Upload _table using _subtable_index as a label to keep track of it, saving the final columns (under default_columns)
    to _final_savepath as an astropy table you can pickle_load. True/False return for if it was a success.

    :param _table: the subtable
    :param _subtable_index: the index
    :param _final_savepath: the savepath
    :return: True/False on whether it was a success or a fail.

    """

    # Generate name (using the provided subtable index for this subtable)
    uploaded_id = "local_upload_" + str(_subtable_index)
    id = upload_table(_table, uploaded_id)

    # Generate a random name for the dr3ids-matching job
    dr3id_match_name = "dr3match_" + uploaded_id
    id = obtain_dr3ids(id, dr3id_match_name)

    # Now, two random names for the jobs squirted out for the fullmatch, and the final columns we select
    fullmatch_name = uploaded_id + "_fullmatch"
    selected_columns_name = uploaded_id + "_columns"

    # Return whether job was success (True) or fail (False)
    return obtain_download_data(id, fullmatch_name, selected_columns_name, _savepath=_final_savepath)

def xmatch_gaiadr3(_table,
                   _RA_COLNAME, _DEC_COLNAME,
                   _OBJID_COLNAME, _MAS_TOLERANCE,
                   _COLS_INTEREST):
    """
    Run an xMatch to the Gaia DR3 catalogue as stored by Vizier using the astropy xMatch class.
    The table is clipped into smaller pieces for the sake of avoiding the timeout issue.
    Hit-or-miss as to whether this will work- in general default to xmatch_gaiadr3_GaiaService.
    Not our fault it won't work- CDS xMatch is a public service, and depends on others not abusing it.

    :param _table: --
    :param _RA_COLNAME: the column name for RA
    :param _DEC_COLNAME: the column name for DEC
    :param _OBJID_COLNAME: the column name for objID/index
    :param _MAS_TOLERANCE: tolerance in mas for conesearch
    :param _COLS_INTEREST: array-like with columns of interest as string
    :return: astropy table
    """

    from astroquery.xmatch import XMatch

    if len(_table) >= 2e5:

        # Set up cache directory
        import os, shutil
        cachedir = os.path.join(windows_directories_new.lamodir, "XMATCH_CACHE")
        try: # to make the cache directory
            os.mkdir(cachedir)
        except: # shutil remove the directory + then try again. if this fails, exception is thrown and code breaks.
            shutil.rmtree(cachedir, ignore_errors=True)
            os.mkdir(cachedir)

        # Generate subtables of matches (max 2e6 rows for a vizier query.)
        cachefilespaths = []
        list_of_indices_lists = np.array_split(np.arange(0, len(_table), 1),
                                               int(len(_table)/2e5) + 1)
        try:
            for num,indices in enumerate(list_of_indices_lists):
                print("Currently matching subtable index ", num)
                filesavepath = os.path.join(cachedir, "XMATCHCACHE_FILE_" + str(num) + ".txt")
                cachefilespaths.append(filesavepath)
                table = XMatch.query(cat1=_table[indices],
                                     cat2='vizier:I/355/gaiadr3',
                                     max_distance=_MAS_TOLERANCE*u.mas,
                                     colRA1=_RA_COLNAME,
                                     colDec1=_DEC_COLNAME,
                                     colRA2="ra",
                                     colDec2="dec",
                                     area="allsky",
                                     cache=False,
                                     selection='best')
                print("The length of the original cat for this was ", len(_table[indices]), "with xMatch of length ", len(table))
                #print("Colnames of the obtained Vizier are... ", table.colnames)
                table = table[_COLS_INTEREST]

                with open(filesavepath, "wb") as f:
                    pickle.dump(obj=table, file=f)
        except Exception as e:
            print("An exception occurred when handling the subtables- purging XMatch Cache and stalling,"
                  "See xmatch_gaiadr3 in gaia_utils_standalone.")
            shutil.rmtree(cachedir, ignore_errors=True)
            print(Exception)
            import time
            time.sleep(99999)

        # Load all subfiles/etc
        master_table = []
        for cachefile in cachefilespaths:
            with open(cachefile, "rb") as f:
                master_table.append(pickle.load(f))
        from astropy.table import vstack
        master_table = vstack(master_table)
        shutil.rmtree(cachedir, ignore_errors=True)

    else:

        master_table = XMatch.query(
            cat1=_table,
            cat2='vizier:I/355/gaiadr3',
            max_distance=_MAS_TOLERANCE * u.mas,
            colRA1=_RA_COLNAME,
            colDec1=_DEC_COLNAME,
            colRA2="ra",
            colDec2="dec",
            area="allsky",
            cache=False,
            selection='best'
        )[_COLS_INTEREST]

    return master_table

def _xmatch_gaiadr3_conesearch_gmatch(
        _table,
        _COLS_INTEREST,
        _MAS_TOLERANCE,
        _G_TOLERANCE=None,
        _RA_COLNAME=None,
        _DEC_COLNAME=None,
        _G_COLNAME=None,
        _DO_G_TOLERANCE=False,
        _USE_LITE=True
):
    """

    NOTE: If your table exceeds around 1e6 rows, use xmatch_gaiadr3_conesearch_gmatch.
    NOTE: gmatch is optional. See _DO_G_TOLERANCE=True to enable it.
    NOTE: Ensure your table has valid values for all elements relevant
    NOTE: _USE_LITE enables/disables gaia_source_lite usage, which speeds up the matching process tremendously
    NOTE: If you need to optimize bandwidth/memory, split the query in two (i.e. SELECT from the catalogue externally.)

    Run an xMatch to gaiadr3.gaia_source(_lite) using ADQL queries with a conesearch with _MAS_TOLERANCE. Optionally,
    enable the choice of requiring stars to be within _G_TOLERANCE in Gaia G-band magnitude. Optionally, use gaia_lite.

    TODO: This seems to work as of 18/11/22, i.e. matches do indeed have the correct result.

    :param _table: astropy
    :param _COLS_INTEREST: columns of interest to select
    :param _MAS_TOLERANCE: in mas for cone search
    :param _G_TOLERANCE: in magnitudes
    :param _RA_COLNAME: internally renamed to ra
    :param _DEC_COLNAME: internally renamed to dec
    :param _G_COLNAME: interally renamed to g_obs
    :param _DO_G_TOLERANCE: bool True/False on whether to require g-band magnitudes to match, too
    :param _USE_LITE: bool True/False on whether to use gaia_source_lite
    :return: match table
    """

    clean_jobstables()

    if len(_table) > 1e6 + 1:
        raise RuntimeWarning("TABLE SIZE EXCEEDS 1e6 : Can not recommend running this with such large tables!")

    if _RA_COLNAME != None:
        _table.rename_column(_RA_COLNAME, "ra_obs")
    if _DEC_COLNAME != None:
        _table.rename_column(_DEC_COLNAME, "dec_obs")
    if _G_COLNAME != None:
        _table.rename_column(_G_COLNAME, "g_obs")

    uploaded_id = upload_table(_table, "xmatch_data")

    if _USE_LITE != True:
        ext = "gaiadr3.gaia_source"
    else:
        ext = "gaiadr3.gaia_source_lite"

    if _DO_G_TOLERANCE != True:
        query = ("""
            SELECT *
            FROM {0} as base
            JOIN {2} as ext
            ON 1 = CONTAINS(
                POINT('ICRS', ra_obs, dec_obs),
                CIRCLE('ICRS', ra, dec, {1})
            )
            ORDER BY base.objid ASC 
        """).format(uploaded_id, _MAS_TOLERANCE/(1000*3600), ext)
    else: # query including the requirement that g-band magnitudes are closer than the tolerance
        query = ("""
            SELECT *
            FROM {0} as base
            JOIN {3} as ext
            ON ABS(base.g_obs - ext.phot_g_mean_mag) < {2} 
            AND 1 = CONTAINS(
                POINT('ICRS', base.ra_obs, base.dec_obs),
                CIRCLE('ICRS', ext.ra, ext.dec, {1}))
            ORDER BY base.objid ASC 
        """).format(uploaded_id, _MAS_TOLERANCE/(1000*3600), _G_TOLERANCE, ext)


    job = Gaia.launch_job_async(query=query,
                                name="run_xmatch")
    results = job.get_results()[_COLS_INTEREST]

    return results

def xmatch_gaiadr3_conesearch_gmatch(
        _table,
        _COLS_INTEREST,
        _MAS_TOLERANCE,
        _G_TOLERANCE=None,
        _RA_COLNAME=None,
        _DEC_COLNAME=None,
        _G_COLNAME=None,
        _DO_G_TOLERANCE=False,
        _USE_LITE=True
):

    """
    Run an xMatch to Gaia DR3 using a cone-search- much slower than the Vizier xMatch service- requiring that the match
    has a similar G-band magnitude (Gaia G.) Note that if you fail to propagate the magnitude into Gaia G, you will end
    up with possible fuzzyness- use at your own risk if this is the case.

    Note: Ensure that the table G-band mags are all valid. Not all stars in LAMOST for example have valid values for
    the ugr(=approx,G)iz magnitudes.

    Note: _COLS_INTEREST MUST MATCH GAIA_SOURCE_LITE IN THIS CASE! WE ARE NOT USING THE SAME COLS AS FOR VIZIER!!!!!!!!

    Note: Only 5455339 sources DO NOT seem to have mean G-band magnitudes. Still, modify the query to ensure that the
    magnitude is indeed present.

    :param _table: astropy
    :param _COLS_INTEREST: columns of interest to select
    :param _G_TOLERANCE: in magnitudes
    :param _MAS_TOLERANCE: in mas for cone search
    :param _RA_COLNAME: internally renamed to ra
    :param _DEC_COLNAME: internally renamed to dec
    :param _USE_LITE: bool True/False on whether to use gaia_source_lite
    :return: match table
    """

    if len(_table) >= 2e5:

        # Set up cache directory
        import os, shutil
        cachedir = os.path.join(windows_directories_new.lamodir, "XMATCH_CACHE_GaiaADQL")
        try: # to make the cache directory
            os.mkdir(cachedir)
        except: # shutil remove the directory + then try again. if this fails, exception is thrown and code breaks.
            shutil.rmtree(cachedir, ignore_errors=True)
            os.mkdir(cachedir)

        # Generate subtables of matches (max 2e6 rows for a vizier query.)
        cachefilespaths = []
        list_of_indices_lists = np.array_split(np.arange(0, len(_table), 1),
                                               int(len(_table)/1e6) + 1)
        for num,indices in enumerate(list_of_indices_lists):
            print("Currently matching subtable index ", num)
            filesavepath = os.path.join(cachedir, "XMATCHCACHE_FILE_" + str(num) + ".txt")
            cachefilespaths.append(filesavepath)
            results = _xmatch_gaiadr3_conesearch_gmatch(_table[indices],
                                                        _COLS_INTEREST,
                                                        _MAS_TOLERANCE,
                                                        _G_TOLERANCE,
                                                        _RA_COLNAME,
                                                        _DEC_COLNAME,
                                                        _G_COLNAME,
                                                        _DO_G_TOLERANCE,
                                                        _USE_LITE
                                                        )

            with open(filesavepath, "wb") as f:
                pickle.dump(obj=results, file=f)

        # Load all subfiles/etc
        master_table = []
        for cachefile in cachefilespaths:
            with open(cachefile, "rb") as f:
                master_table.append(pickle.load(f))
        from astropy.table import vstack
        master_table = vstack(master_table)
        shutil.rmtree(cachedir, ignore_errors=True)

    else:

        master_table= _xmatch_gaiadr3_conesearch_gmatch(_table,
                                                        _COLS_INTEREST,
                                                        _MAS_TOLERANCE,
                                                        _G_TOLERANCE,
                                                        _RA_COLNAME,
                                                        _DEC_COLNAME,
                                                        _G_COLNAME,
                                                        _DO_G_TOLERANCE,
                                                        _USE_LITE
                                                        )


    return master_table

def jl_setdiff1dbyindices(
        set1,
        set2
):
    """
    Helper for the Julia setdiff1d.jl functionality. Specifically, return indices of set1 that are not in set2
    as a 1d truefalse. Only tested for |set1|>|set2|. Spawns 8 procs- you should modify the source if you haven't got
    enough free threads.

    Example:
    ------
    >>> a = [1,2,3,5,7,8,9]
    >>> b = [1,2,3,4,5,8,9]
    >>> jl_setdiff1dbyindices(a,b)
    >>> [False, False, False, True, True, False, False]


    :param set1: intarray 1d
    :param set2: intarray 1d
    :return: boolarray 1d
    """

    from juliacall import Main
    import os
    Main.include(os.path.join(windows_directories_new.jl_dir,
                              "setdiff1d.jl"))
    set1_truefalse = Main.setdiff1d(
        set1,
        set2
    )

    return set1_truefalse