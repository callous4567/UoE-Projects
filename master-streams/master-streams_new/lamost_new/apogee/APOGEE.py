"""
Take APOGEE DR17 data and SEGUE data. TODO: SEGUE
Match to GAIA DR3 (using EDR3 source IDs if available, else using standard match ra/dec.)
Note that APOGEE and SEGUE have distant Julian dates- account for this
Once matched, match to our old KGs/BHBs...
- KGs https://arxiv.org/abs/1211.0549 BHBs https://arxiv.org/pdf/0902.1781.pdf
"""

import os
from astropy.io import fits
from astropy.table import Table

import ascii_info
import hdfutils
import windows_directories_new

apogee_dir = os.path.join(windows_directories_new.lamodir, "apogee")
starhorse_path = os.path.join(apogee_dir, "APOGEE_DR17_EDR3_STARHORSE_v2.fits")
starhorse_path_new = os.path.join(apogee_dir, "APOGEE_DR17_EDR3_STARHORSE_v2_processed.fits")
apogee_path = os.path.join(apogee_dir, "allStar-dr17-synspec_rev1.fits")
apogee_starhorse_path = os.path.join(apogee_dir, "allStar_Starhorse.fits")
apogee_starhorse_path_nogcs = os.path.join(apogee_dir, "allStar_Starhorse_nogcs.fits")

do_gaia_starhorse = False
do_apogee_starhorse = False
load = False
do_final_process = False
do_replot = False

if do_gaia_starhorse == True and __name__ == "__main__":

    # Load in starhorse
    import numpy as np
    with fits.open(starhorse_path, mode='readonly', memmap=True) as hdul:
        # Data (which is an astropy table when loaded.)
        data = Table(hdul[1].data)
        data['EDR3_source_id'] = np.array(data['EDR3_source_id'], np.int64)

    # Isolate just the source IDs
    source_id_table = Table()
    source_id_table['EDR3_source_id'] = data['EDR3_source_id']

    # Imports for Gaia
    from astroquery.gaia import Gaia
    import gaia_utils_standalone

    # Get list of all tables on our Gaia, and prune to useful ones (that match our own upload, specifically) for removal
    tables = Gaia.load_tables(only_names=True, include_shared_tables=False)
    remove_names = []
    for table in tables:
        # split by user (to get rid of common tables, i.e. SDSS or DR3/etc)
        namesplit = table.name.split("user_" + gaia_utils_standalone.username)
        if len(namesplit) >= 2:
            remove_names.append(table.name)
    import numpy as np
    remove_names = np.array(remove_names, dtype=str)
    gaia_utils_standalone.deltables(remove_names)

    # Get all the jobs and delete them
    remove_jobs = Gaia.list_async_jobs()
    remove_jobs = [d.jobid for d in remove_jobs]
    remove_jobs = np.array(remove_jobs, dtype=str)
    gaia_utils_standalone.deljobs(remove_jobs)

    # Upload our source ID table
    apogee_tablocal_id = gaia_utils_standalone.upload_table(source_id_table, "apogee_source_ids")

    # Take this table, and retrieve the gaia columns.
    default_columns = ["source_id",
                       "ra", "dec",
                       "parallax", "parallax_error", "pmra",
                       "pmra_error", "pmdec", "pmdec_error", "radial_velocity",
                       "radial_velocity_error",
                       "ruwe"]

    # Set up the string for selection col_1, col_2, col_3... ,col_n
    selection_string = default_columns[0]
    for i in range(1, len(default_columns)):
        selection_string += (",{}").format(default_columns[i])

    # Set up query
    query = ("""
        SELECT *
        FROM {0}
        JOIN gaiadr3.gaia_source
            ON gaiadr3.gaia_source.source_id = {0}.EDR3_source_id 
    """).format(apogee_tablocal_id)

    # Launch query
    job = Gaia.launch_job_async(query=query, name="apogee_DR3_match")
    Gaia.upload_table_from_job(job)
    Gaia.delete_user_table(table_name=apogee_tablocal_id, force_removal=True)
    Gaia.remove_jobs([job.jobid])
    job_table_id = "user_" + gaia_utils_standalone.username + ".t" + str(job.jobid)

    # Take the joined table and select what we want
    query = ("""
        SELECT {0}
        FROM {1}
    """).format(selection_string, job_table_id)

    # Launch query
    job = Gaia.launch_job_async(query=query, name="apogee_DR3_obtain")
    Gaia.upload_table_from_job(job)
    Gaia.delete_user_table(table_name=job_table_id, force_removal=True)
    Gaia.remove_jobs([job.jobid])
    job_table_id = "user_" + gaia_utils_standalone.username + ".t" + str(job.jobid)

    # Grab final results
    job_results = job.get_results()
    Gaia.delete_user_table(table_name=job_table_id)

    # Re-combine with the original table
    from astropy.table import hstack
    data = hstack([job_results, data])

    # Save it
    data.write(starhorse_path_new, overwrite=True)

# Note:
# JOIN ON APOGEE_ID
if do_apogee_starhorse == True and __name__ == "__main__":

    # Load in apogee.
    import numpy as np
    with fits.open(apogee_path, mode='readonly', memmap=True) as hdul:
        # Data (which is an astropy table when loaded.)
        apogee = hdul[1].data
        apogee_ID = np.array(apogee['APOGEE_ID'], str)

    # Load in starhorse (which we will modify.)
    with fits.open(starhorse_path_new, mode='readonly', memmap=True) as hdul:
        # Data (which is an astropy table when loaded.)
        starhorse = Table(hdul[1].data)
        starhorse_ID = np.array(starhorse['APOGEE_ID'], str)

    import pickle

    if load == False:

        # Get indices for splits
        indices = np.arange(0, len(starhorse_ID), 1)
        size = 10000
        split_indices = np.array_split(indices, int(len(indices) / size + 1))

        # Set up for the zips
        zipped = []
        for split_indices_sub in split_indices:
            zipped.append([starhorse_ID[split_indices_sub], apogee_ID])

        # Pool for matches.
        import APOGEE_multi, multiprocessing
        ncores = 10
        pool = multiprocessing.Pool(ncores)
        results = pool.map(APOGEE_multi.do_strmatch, zipped)
        pool.close()
        with open(apogee_dir + "\\results.txt", "wb") as f:
            pickle.dump(obj=results, file=f)

        load = True

    if load == True:

        with open(apogee_dir + "\\results.txt", "rb") as f:
            results = pickle.load(file=f)

        # Trim apogee
        import itertools
        matches = np.array(list(itertools.chain.from_iterable(results)), np.int64)
        apogee_table = apogee[matches]

        # Some bookkeeping on Starhorse... first, get distance/met
        starhorse.rename_columns(['dist50','met50'], ['dist','feh'])

        # Next, get standard deviation (half the percentile range)
        starhorse['edist'] = 1/2 * np.abs((starhorse['dist84'] - starhorse['dist16']))
        starhorse['efeh'] = 1/2 * np.abs((starhorse['met84'] - starhorse['met16']))

        # Verify
        neq = 0
        for i,j in zip(apogee_table['APOGEE_ID'],starhorse['APOGEE_ID']):
            if i != j:
                neq += 1
                print("Not equal ... ", i,j)
        print("Total not equal is ... ", neq)


        # Add columns for vlos/evlost
        starhorse.add_column(9999., name='vlos')
        starhorse.add_column(9999., name='evlost')

        # Extract from apogee
        starhorse['vlos'] = apogee_table['VHELIO_AVG']
        starhorse['evlost'] = apogee_table['VERR']

        # Save the new starhorse
        starhorse.write(apogee_starhorse_path, overwrite=True)

if do_final_process == True and __name__ == "__main__":

    with fits.open(apogee_starhorse_path, mode='readonly', memmap=True) as hdul:
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
    data.write(apogee_starhorse_path_nogcs, overwrite=True)

    # Do a preliminary plot
    datatab = np.array([data['Lx'], data['Ly'], data['Lz']]).T
    # clustered = cluster3d().listhdbs(data, [1000, 20])
    import graphutils_new
    graphutils_new.threed_graph().kmeans_L_array(datatab, [1 for d in datatab[:, 0]], False, browser=True,
                                                 outliers=True)

if do_replot == True and __name__ == "__main__":

    import numpy as np

    with fits.open(apogee_starhorse_path_nogcs, mode='readonly', memmap=True) as hdul:
        # Data (which is an astropy table when loaded.)
        data = Table(hdul[1].data)
        data = data[[True if vv < 1000 else False for vv in np.sqrt(data['vx']**2 + data['vy']**2 + data['vz']**2)]]
        data = data[[True if r < 100 else False for r in data['r']]]
        data['source'] = "APOGEE_sstraszak"


    # QUALITY CUTS ======================================^^^^
    from galcentricutils_quality_cuts import quality_cuts
    qcs = quality_cuts()
    data = data[qcs.zmax_cuts(data, 4,
                              3e9,
                              3000,
                              True,
                              True,
                              "apogee.txt")]
    metmax = -1/2
    data = data[qcs.feh_cut(data, metmax)]
    #data = data[qcs.LSR_cut(data, 220)]
    vphi_plotlims = [-400,400]
    data = data[qcs.bound_cut(data)]
    # QUALITY CUTS ======================================^^^^



    # Grab hold of the other data
    import windows_directories_new
    writer = hdfutils.hdf5_writer(windows_directories_new.datadir, ascii_info.asciiname)
    table = writer.read_table(ascii_info.fullgroup, ascii_info.fullset)
    from astropy.table import vstack
    data = vstack([data, table])



    # QUALITY CUTS ======================================^^^^
    data = data[qcs.fancy_LSR_cut(data, 50)]
    # QUALITY CUTS ======================================^^^^



    # Set up the flat hdbscan run
    from hdbscan import flat
    # Do a preliminary plot
    datatab = np.array([data['Lx'], data['Ly'], data['Lz']]).T
    clusterer = flat.HDBSCAN_flat(X=datatab,
                                  n_clusters=28,
                                  min_cluster_size=12,
                                  min_samples=6,
                                  metric='l2',
                                  algorithm='best',
                                  prediction_data=True)

    # clustered = cluster3d().listhdbs(data, [1000, 20])
    import graphutils_new
    graphutils_new.threed_graph().kmeans_L_array(datatab, clusterer.labels_, False, browser=True,
                                                 outliers=False, limit_view=False)

    # generate for 13
    from graphutils_new import spec_graph, twod_graph
    for clid in range(np.max(clusterer.labels_)):

        spec_graph().clust_radec(data, clusterer.labels_, cluster_id=clid,
                                            savedexdir="lamotest_" + str(clid) + "_clusttest_lb",
                                            lb=True, vasiliev=False,
                                            flatfork=False)
        spec_graph().clust_radec(data, clusterer.labels_, cluster_id=clid,
                                            savedexdir="lamotest_" + str(clid) + "_clusttest_ra",
                                            lb=False, vasiliev=False,
                                            flatfork=False)

        feh = data[[True if d == clid else False for d in clusterer.labels_]]['feh']
        twod_graph().hist_fancy(feh, 20, 'FeH', r'$\rho$', savepath=os.path.join(windows_directories_new.imgdir,
                                                                                 str(clid) + "_apogee_fehhist.png"),
                                lims=[-2,0])