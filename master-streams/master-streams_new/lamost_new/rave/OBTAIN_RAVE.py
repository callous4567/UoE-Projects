from pkg_resources import parse_version
import requests
import pyvo as vo
import pandas as pd
import os
import windows_directories_new
import pickle
from astropy.table import hstack

# Run RAVE DR6 queries or not? See the match path for the match. See below for GAIA-related pmra/pmdec match.
run_rave = False
# Run gaia queries or not?
run_gaia = False
# Get galcent/etc and do general quality assurance/etc?
run_preprocessing = False
just_plot_preprocessed = False
# Bdasp distance histogram
run_bdasp_distance_hist = False
ravedir = os.path.join(windows_directories_new.lamodir, "rave")
rave_bdasp_match_path = os.path.join(ravedir, "ravedr6_bdasp.txt")
rave_bdasp_match_gaia_path = os.path.join(ravedir, "ravedr6_bdasp_gaia.fits")
preprocessed_path = os.path.join(ravedir, "preprocessed.fits")
if run_rave == True:

    service_name = 'rave-survey.org' # set up tap
    url = "https://www.rave-survey.org/tap"
    token = "387e0bbae743459bf29744ebf1c436ee6cbc48d9"

    # Setup authorization
    tap_session = requests.Session()
    tap_session.headers['Authorization'] = token
    tap_service = vo.dal.TAPService(url, session=tap_session)


    # Set up query to obtain the BDASP catalogue of interest
    lang = 'PostgreSQL'
    query = ('''
    SELECT rave_obs_id, distance_bdasp, distance_error_bdasp, m_h_bdasp, m_h_error_bdasp
    FROM ravedr6.dr6_bdasp
    WHERE ravedr6.dr6_bdasp.distance_bdasp >= 5
    AND ravedr6.dr6_bdasp.distance_bdasp < 100
    ''')

    job = tap_service.run_sync(query, language=lang)
    bdasp_table = job.to_table()

    query = ('''
    SELECT rave_obs_id, hrv_sparv, hrv_error_sparv
    FROM ravedr6.dr6_sparv
    ''')

    job = tap_service.run_sync(query, language=lang)
    ravedr6_table = job.to_table()

    from APOGEE_standalone import strmatch_general
    matches, found = strmatch_general(bdasp_table['rave_obs_id'], ravedr6_table['rave_obs_id'])
    import numpy as np
    matches = np.array(matches, int)
    found = np.array(found, bool)
    ravedr6_table = ravedr6_table[matches]
    ravedr6_table = ravedr6_table[found]
    bdasp_table = bdasp_table[found]
    rave_bdasp_stack = hstack([ravedr6_table, bdasp_table])
    with open(rave_bdasp_match_path, "wb") as f:
        pickle.dump(obj=rave_bdasp_stack, file=f)
if run_gaia == True:

    import gaia_utils_standalone
    from astroquery.gaia import Gaia
    import numpy as np

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

    with open(rave_bdasp_match_path, "rb") as f:
        rave_bdasp_stack = pickle.load(file=f)

    rave_obs_id = rave_bdasp_stack['rave_obs_id_1']
    from astropy.table import Table

    # Set up rave ID table
    table = Table()
    table['rave_obs_id_1'] = rave_obs_id

    # Upload and match to Gaia Dr3 crossmatch table
    table_id = gaia_utils_standalone.upload_table(table, "ravedr6_bdasp_table")
    query = ("""
        SELECT *
        FROM gaiadr3.ravedr6_best_neighbour
        JOIN {0}
            ON gaiadr3.ravedr6_best_neighbour.original_ext_source_id = {0}.rave_obs_id_1
    """).format(table_id)
    job = Gaia.launch_job_async(query=query, name="ravedr6_dr3_match")
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
                       "rave_obs_id_1",
                       "angular_distance",
                       "number_of_neighbours",
                       "ra", "dec",
                       "parallax", "parallax_error", "pmra",
                       "pmra_error", "pmdec", "pmdec_error", "radial_velocity",
                       "radial_velocity_error",
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
    job_table = hstack([rave_bdasp_stack, job_table])
    job_table.write(rave_bdasp_match_gaia_path, overwrite=True)

if run_preprocessing == True:

    from astropy.io import fits

    from astropy.table import Table

    with fits.open(rave_bdasp_match_gaia_path, mode='readonly') as hdul:
        data = Table(hdul[1].data)

    data.rename_columns(['hrv_sparv', 'hrv_error_sparv',
                         'distance_bdasp', 'distance_error_bdasp'],
                        ['vlos', 'evlost', 'dist', 'edist'])

    # Get the galactic coordinates from these equatorial ones
    from energistics_new import orbigistics
    from galcentricutils_new import angular

    orbigist = orbigistics()

    # Get GAL from ICRS
    data = orbigist.converter.nowrite_ICRS_to_GAL(data, has_cosfactor=True)

    # Get galactic errors. The cosfactor in the pmra/pmdec still exists-
    # no longer in dmu_l, dmu_b though.
    import lamost_utils
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
    data.write(preprocessed_path, overwrite=True)

    # While we're at it, get kmeans_L
    datatab = np.array([data['Lx'], data['Ly'], data['Lz']]).T
    import graphutils_new
    graphutils_new.threed_graph().kmeans_L_array_path(datatab,
                                                      [1 for d in datatab[:, 0]],
                                                      os.path.join(ravedir, "preprocessed_kmeans_test.html"),
                                                      browser=True,
                                                      outliers=True,
                                                      limit_view=True)

if just_plot_preprocessed == True:

    from astropy.io import fits
    import numpy as np
    from astropy.table import Table

    with fits.open(preprocessed_path, mode='readonly') as hdul:
        data = Table(hdul[1].data)

    import matplotlib.pyplot as plt
    dists = data['dist'][[True if d < 100 else False for d in data['dist']]]
    plt.hist(dists, bins=50)
    plt.show()

    """
    data = data[[True if ruwe < 1.4 else False for ruwe in data['ruwe']]]
    data = data[[True if d < 1000 else False for d in np.sqrt(data['vx']**2 +
                                                              data['vy']**2 +
                                                              data['vz']**2)]]

    # While we're at it, get kmeans_L
    datatab = np.array([data['Lx'], data['Ly'], data['Lz']]).T
    import graphutils_new
    graphutils_new.threed_graph().kmeans_L_array_path(datatab,
                                                      [1 for d in datatab[:, 0]],
                                                      os.path.join(ravedir, "preprocessed_kmeans_test.html"),
                                                      browser=True,
                                                      outliers=True,
                                                      limit_view=False)
    """

if run_bdasp_distance_hist == True:

    service_name = 'rave-survey.org' # set up tap
    url = "https://www.rave-survey.org/tap"
    token = "387e0bbae743459bf29744ebf1c436ee6cbc48d9"

    # Setup authorization
    tap_session = requests.Session()
    tap_session.headers['Authorization'] = token
    tap_service = vo.dal.TAPService(url, session=tap_session)


    # Set up query to obtain the BDASP catalogue of interest
    lang = 'PostgreSQL'
    query = ('''
    SELECT distance_bdasp, distance_error_bdasp
    FROM ravedr6.dr6_bdasp
    ''')

    job = tap_service.run_sync(query, language=lang)
    bdasp_table = job.to_table()
    dist = bdasp_table['distance_bdasp']
    import matplotlib.pyplot as plt
    plt.hist(dist, bins=100)
    plt.show()


