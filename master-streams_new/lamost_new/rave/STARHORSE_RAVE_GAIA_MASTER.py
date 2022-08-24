# Use STARHORSE DISTANCES and NOT BDASP
import os
import windows_directories_new

import numpy as np
from astropy.table import Table

ravedir = os.path.join(windows_directories_new.lamodir, "rave")
starhorse_path = os.path.join(ravedir, "rave_starhorse.dat")
starhorse_path_fits = os.path.join(ravedir, "rave_starhorse.fits")
from astropy.io import fits
run_dat_starhorse_fits = False
if run_dat_starhorse_fits == True:

    columns = [
        'rave_obs_id',
        'dr2_source_id',
        'GLON',
        'GLAT',
        'Mass16',
        'Mass50',
        'Mass84',
        'Teff16',
        'Teff50',
        'Teff84',
        'logg16',
        'logg50',
        'logg84',
        'Met16',
        'Met50',
        'Met84',
        'Dist16',
        'Dist50',
        'Dist84',
        'AV16',
        'AV50',
        'AV84'
    ]
    dtypes = [
        str,
        str,
        float,
        float,
        float,
        float,
        float,
        float,
        float,
        float,
        float,
        float,
        float,
        float,
        float,
        float,
        float,
        float,
        float,
        float,
        float,
        float
    ]

    with open(starhorse_path, "r") as f:
        lines = f.readlines()
        lines = [line.split(" ") for line in lines]
        lines = [[d for d in u if d != ''] for u in lines]
        lines = [[d for d in u if d != '\n'] for u in lines]
        lines = [line[0:22] for line in lines]
        table = Table()
        for num,column in enumerate(columns):
            coldata = [line[num] for line in lines]
            table[column] = np.array(coldata, dtypes[num])
        table.rename_columns(['Dist50','Met50','GLON','GLAT'], ['dist','feh','l','b'])
        table['edist'] = 1/2 * np.abs(table['Dist84']-table['Dist16'])
        table['efeh'] = 1/2 * np.abs(table['Met84']-table['Met16'])
        table.write(starhorse_path_fits, overwrite=True)

# Take the main catalogue from rave master and get it into a fits (and rename dirty cols)
# While we're at it, crossmatch to Gaia and get proper motions/etc.
run_master_fits_gaia = False
rave_master_path = os.path.join(ravedir, "III_283_master.dat.gz.fits")
rave_gaia_master_path = os.path.join(ravedir, "rave_gaia.fits")
if run_master_fits_gaia == True:

    import gaia_utils_standalone
    from astroquery.gaia import Gaia
    import numpy as np
    from astropy.table import join

    # Try to load the cache by default! :)
    cache_path = os.path.join(ravedir, "mastergaiacache.fits")
    try:
        with fits.open(cache_path, mode='readonly', memmap=False) as hdul:
            master_data = Table(hdul[1].data)

    except:


        with fits.open(rave_master_path, mode='readonly', memmap=False) as hdul:
            master_data = Table(hdul[1].data)
            master_data['ObsID'] = np.array([d.replace(" ", "") for d in master_data['ObsID']],
                                             str)

        master_data.rename_columns(['HRV','e_HRV','ObsID'], ['vlos', 'evlost', 'rave_obs_id'])

        # Get list of all tables on our Gaia, and prune to useful ones (that match our own upload, specifically) for removal
        tables = Gaia.load_tables(only_names=True, include_shared_tables=False)
        remove_names = []
        for table in tables:
            # split by user (to get rid of common tables, i.e. SDSS or DR3/etc)
            namesplit = table.name.split("user_" + gaia_utils_standalone.username)
            if len(namesplit) >= 2:
                remove_names.append(table.name)
        remove_names = np.array(remove_names, dtype=str)
        gaia_utils_standalone.deltables(remove_names)

        # Get all the jobs and delete them
        remove_jobs = Gaia.list_async_jobs()
        remove_jobs = [d.jobid for d in remove_jobs]
        remove_jobs = np.array(remove_jobs, dtype=str)
        gaia_utils_standalone.deljobs(remove_jobs)

        # Set up rave ID table
        table = Table()
        table['rave_obs_id'] = master_data['rave_obs_id']

        # Upload and match to Gaia Dr3 crossmatch table
        table_id = gaia_utils_standalone.upload_table(table, "ravedr6_master_match")
        query = ("""
            SELECT *
            FROM gaiadr3.ravedr6_best_neighbour
            JOIN {0}
                ON gaiadr3.ravedr6_best_neighbour.original_ext_source_id = {0}.rave_obs_id
        """).format(table_id)
        job = Gaia.launch_job_async(query=query, name="ravedr6_master_match")
        Gaia.upload_table_from_job(job)
        Gaia.remove_jobs([job.jobid])
        gaia_utils_standalone.deltables([table_id])
        job_table_id = "user_" + gaia_utils_standalone.username + ".t" + str(job.jobid)

        # Select from Gaia dr3 (all)
        query = ("""
            SELECT *
            FROM gaiadr3.gaia_source
            JOIN {0}
                ON gaiadr3.gaia_source.source_id = {0}.source_id
        """).format(job_table_id)
        job = Gaia.launch_job_async(query=query, name="ravedr6_gaia_datamatch")
        Gaia.upload_table_from_job(job)
        Gaia.remove_jobs([job.jobid])
        gaia_utils_standalone.deltables([job_table_id])
        job_table_id = "user_" + gaia_utils_standalone.username + ".t" + str(job.jobid)

        # Set columns we want
        default_columns = ["source_id",
                           "original_ext_source_id",
                           "rave_obs_id",
                           "ra", "dec",
                           "parallax", "parallax_error", "pmra",
                           "pmra_error", "pmdec", "pmdec_error",
                           "ruwe"]

        # Set up the string for selection col_1, col_2, col_3... ,col_n
        selection_string = default_columns[0]
        for i in range(1, len(default_columns)):
            selection_string += (",{}").format(default_columns[i])

        # Now select what we want from this
        query = ("""
            SELECT {0}
            FROM {1}
        """).format(selection_string, job_table_id)
        job = Gaia.launch_job_async(query=query, name="download_match_data")
        job_table = job.get_results()
        Gaia.remove_jobs([job.jobid])
        gaia_utils_standalone.deltables([job_table_id])
        job_table_id = "user_" + gaia_utils_standalone.username + ".t" + str(job.jobid)
        gaia_utils_standalone.deltables([job_table_id])

        # Join our ravedr6 data with the job_table data (gaia data)
        master_data = join(master_data, job_table, ['rave_obs_id'])
        try:
            master_data.write(cache_path, overwrite=True)
        except:
            pass

    # Get list of all tables on our Gaia, and prune to useful ones (that match our own upload, specifically) for removal
    tables = Gaia.load_tables(only_names=True, include_shared_tables=False)
    remove_names = []
    for table in tables:
        # split by user (to get rid of common tables, i.e. SDSS or DR3/etc)
        namesplit = table.name.split("user_" + gaia_utils_standalone.username)
        if len(namesplit) >= 2:
            remove_names.append(table.name)
    remove_names = np.array(remove_names, dtype=str)
    gaia_utils_standalone.deltables(remove_names)

    # Get all the jobs and delete them
    remove_jobs = Gaia.list_async_jobs()
    remove_jobs = [d.jobid for d in remove_jobs]
    remove_jobs = np.array(remove_jobs, dtype=str)
    gaia_utils_standalone.deljobs(remove_jobs)

    # Due to matching problems (the stupid ass Starhorse catalogue hasn't got correct rave_source_ids...
    # We need to get source dr2_source_ids! Good fucking luck! :D
    source_ids_table = Table()
    master_data.rename_column('source_id', 'dr3_source_id')
    source_ids_table['dr3_source_id'] = master_data['dr3_source_id']

    # Upload and match to Gaia Dr3 crossmatch table
    table_id = gaia_utils_standalone.upload_table(source_ids_table, "dr3_dr2_ravedr6_match")
    query = ("""
        SELECT *
        FROM {0}
        JOIN gaiadr3.dr2_neighbourhood
            ON gaiadr3.dr2_neighbourhood.dr3_source_id = {0}.dr3_source_id
    """).format(table_id)
    job = Gaia.launch_job_async(query=query, name="dr3_dr2_ravedr6_match")
    job_table = job.get_results()
    Gaia.remove_jobs([job.jobid])
    try:
        gaia_utils_standalone.deltables([job_table_id])
    except:
        pass
    job_table_id = "user_" + gaia_utils_standalone.username + ".t" + str(job.jobid)
    gaia_utils_standalone.deltables([job_table_id])
    master_data['dr3_source_id'] = np.array(master_data['dr3_source_id'], np.int64)
    master_data = join(master_data, job_table, ['dr3_source_id'])

    # Finally save
    master_data.write(rave_gaia_master_path, overwrite=True)

# Match the master fits with Gaia to Starhorse
run_master_starhorse_join = False
master_starhorse_path = os.path.join(ravedir, "master_starhorse.fits")
if run_master_starhorse_join == True:

    with fits.open(rave_gaia_master_path, mode='readonly', memmap=False) as hdul:
        master_data = Table(hdul[1].data)

    with fits.open(starhorse_path_fits, mode='readonly', memmap=False) as hdul:
        starhorse_data = Table(hdul[1].data)

    from astropy.table import join

    import numpy as np
    starhorse_data['dr2_source_id'] = np.array(starhorse_data['dr2_source_id'], np.int64)
    master_data = join(master_data, starhorse_data, ["dr2_source_id"])

    # Specify all the columns of interest (just to keep things crisp.)
    columns_of_interest = [
        'rave_obs_id',
        'vlos',
        'evlos',
        'ra',
        'dec',
        'parallax',
        'parallax_error',
        'pmra',
        'pmra_error',
        'pmdec',
        'pmdec_error',
        'ruwe',
        'dr2_source_id',
        'dr3_source_id',
        'l',
        'b',
        'feh',
        'efeh',
        'dist',
        'edist'
    ]
    colnames = master_data.colnames
    columns_to_remove = []
    for column in colnames:
        if column not in columns_of_interest:
            columns_to_remove.append(column)
    master_data.remove_columns(columns_to_remove)

    master_data.write(master_starhorse_path, overwrite=True)

# Do the standard GCC removal, preliminary plots, etc.
do_gcc_removal = False
gcc_removed_path = os.path.join(ravedir, "gcc_removed.fits")
if do_gcc_removal == True:

    with fits.open(master_starhorse_path, mode='readonly', memmap=True) as hdul:
        # Data (which is an astropy table when loaded.)
        data = Table(hdul[1].data)

    # Get the galactic coordinates from these equatorial ones
    from energistics_new import orbigistics
    from galcentricutils_new import angular

    orbigist = orbigistics()

    # Get GAL from ICRS
    data = orbigist.converter.nowrite_ICRS_to_GAL(data, has_cosfactor=True)

    # Get galactic errors. The cosfactor in the pmra/pmdec still exists- no longer in dmu_l, dmu_b though.
    import lamost_utils

    # Note that in this particular catalogue, l/b is GLON/GLAT
    data = lamost_utils.monte_ICRSGAL_table(data)

    # Remove cosfactor from ICRS to match galconversion convention.
    import numpy as np
    data['pmra'] /= np.cos(np.deg2rad(data['dec']))
    data['pmra_error'] /= np.cos(np.deg2rad(data['dec']))

    # Get Galcent + Momenta
    data = orbigist.converter.nowrite_GAL_to_GALCENT(data)
    data = angular().get_momentum(data)

    # Remove nuisance stars (outside min_radius, but within the max_radius of interest.)
    min_radius, max_radius = 15,150
    data = data[[True if max_radius > r > min_radius else False for r in data['r']]]

    # Grab the table of globular clusters
    with fits.open(windows_directories_new.baumgardt_fits, mode='readonly', memmap=False) as hdul:
        gc_table = Table(hdul[1].data)

    # Add dummy values for proper motions
    gc_table.add_column(0., name='pmra')
    gc_table.add_column(0., name='pmdec')
    gc_table.add_column(0., name='pmra_error')
    gc_table.add_column(0., name='pmdec_error')
    gc_table.add_column(0., name='vlos')
    gc_table.add_column(0., name='evlost')

    # Get Galcent (note that baumgardt catalogue does have a cosfactor attached- cos(dec).)
    gc_table = orbigist.converter.nowrite_ICRS_to_GAL(gc_table, has_cosfactor=True)
    gc_table = orbigist.converter.nowrite_GAL_to_GALCENT(gc_table)

    # Convert tidal radius to kpc
    gc_table['rt'] /= 1000

    # Remove stars within gcs
    import gcc_utils

    gc_to_remove = gcc_utils.remove_gcs(np.array(data['x'], float), np.array(data['y'], float),
                                        np.array(data['z'], float),
                                        np.array(gc_table['x'], float), np.array(gc_table['y'], float),
                                        np.array(gc_table['z'], float),
                                        np.array(gc_table['rt'], float))
    data = data[gc_to_remove]

    # Save fits and ascii.
    data.write(gcc_removed_path, overwrite=True)

    # Do a preliminary plot
    datatab = np.array([data['Lx'], data['Ly'], data['Lz']]).T
    # clustered = cluster3d().listhdbs(data, [1000, 20])
    import graphutils_new
    graphutils_new.threed_graph().kmeans_L_array(datatab, [1 for d in datatab[:, 0]], False, browser=True,
                                                 outliers=True)

with fits.open(gcc_removed_path, mode='readonly', memmap=True) as hdul:
    # Data (which is an astropy table when loaded.)
    data = Table(hdul[1].data)
dists = data['dist']
print(dists)
