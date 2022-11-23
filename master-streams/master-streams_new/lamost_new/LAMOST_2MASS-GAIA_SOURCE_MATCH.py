import multiprocessing
import os
import pickle
import matplotlib.pyplot as plt
import numba
import numpy as np
from astropy.io import votable, fits
from astropy.table import Table, vstack
from astroquery.gaia import Gaia
from numba import njit
import ascii_info_new
import graphutils_new
import windows_directories_new
from energistics_new import fast_energistics_new

# TODO NOTE: general catalogue FITS broken- use CSV (need to process.)
# TODO NOTE: the value-added-catalogue is a DAT file, too. Double oof.
process_vat = False
raw_dir = os.path.join(windows_directories_new.lamodir, "dr7_lamost_LR_release.dat")
vat_dir = os.path.join(windows_directories_new.lamodir, "dr7_LRS_vat.fits")
dr7_2mass_gaia_dir = os.path.join(windows_directories_new.lamodir, "dr7_2MASS_aidists.fits")
if __name__ == "__main__" and process_vat == True:

    # Path to the Value-Added-Catalogue for DR7 that has the distances + sources/etc.
    # http://dr7.lamost.org/doc/vac
    path = raw_dir
    with open(os.path.join(windows_directories_new.lamodir, path), "r") as f:
        lines = f.readlines()
    lines = [d.replace("\n", "") for d in lines]
    lines = [d.replace("#", "") for d in lines]
    lines = [d.split() for d in lines]
    columns = lines[0]
    columns = [d.replace("(", "_") for d in columns]
    columns = [d.replace(")", "_") for d in columns]
    columns = [d.strip() for d in columns]
    columns = columns[1:]
    name = [line[0] for line in lines[1:]]
    name = np.array(name, str)
    data = [line[1:] for line in lines[1:]]
    data = np.array([np.array(d, float) for d in data]).T
    dr7_LRS_vat = Table()
    dr7_LRS_vat["Name_2MASS"] = name
    for num, column in enumerate(columns):
        dr7_LRS_vat[column] = data[num]

    # Save it as a fits.
    dr7_LRS_vat.write(vat_dir, overwrite=True)

# Should I crossmatch the 2MASS-LAMOST with Gaia dataset, based on sky coordinate?
LAMOST_2MASS = False
if __name__ == "__main__" and LAMOST_2MASS == True:
    import gaia_utils_standalone

    # DO SOME CLEANING UP ON OUR GAIA WORKSPACE (REMOVE ALL TABLES/JOBS THAT WE HAVE DONE.)
    ###################################################################################################################

    ###################################################################################################################
    # Get list of all tables on our Gaia, and prune to useful ones (that match our own upload, specifically) for removal
    tables = Gaia.load_tables(only_names=True, include_shared_tables=False)
    remove_names = []
    for table in tables:
        # split by user (to get rid of common tables, i.e. SDSS or DR3/etc)
        namesplit = table.name.split("user_" + gaia_utils_standalone.username)
        if len(namesplit) >= 2:
            remove_names.append(table.name)
    remove_names = np.array(remove_names, dtype=str)

    # Split the table name list into manageable multicore chunks
    size = 10
    remove_names = np.array_split(remove_names, int(len(remove_names) / size + 1))

    # Pool over each of them and delete all our old tables
    ncores = 1
    pool = multiprocessing.Pool(ncores)
    pool.map(gaia_utils_standalone.deltables, remove_names)
    pool.close()

    # Get all the jobs and delete them
    remove_jobs = Gaia.list_async_jobs()
    remove_jobs = [d.jobid for d in remove_jobs]
    remove_jobs = np.array(remove_jobs, dtype=str)

    # Split the table name list into manageable multicore chunks
    size = 10
    remove_jobs = np.array_split(remove_jobs, int(len(remove_jobs) / size + 1))

    # Pool over each of them
    pool = multiprocessing.Pool(ncores)
    pool.map(gaia_utils_standalone.deljobs, remove_jobs)
    pool.close()
    ###################################################################################################################

    ###################################################################################################################


    with fits.open(vat_dir, mode='readonly', memmap=True) as hdul:

        # Data (which is an astropy table when loaded.)
        # "gaia_source_id" is dr2_source_id inside this.
        dr7_LRS_vat = Table(hdul[1].data)

    # Get stars a minimum of 4 kpc from us
    dr7_vat = dr7_LRS_vat[[True if d > 4 else False for d in dr7_LRS_vat['d_kpc_']]]
    preupload_names = dr7_vat.colnames

    # Upload the dr7_vat and do a match against the catalogue based on 2MASS source id.
    dr7_vat_uploaded_id = gaia_utils_standalone.upload_table(dr7_vat, "test_dr7_vat")
    source_id_col = preupload_names[0]
    query = ("""
        SELECT *
        FROM {0}
        JOIN gaiadr3.tmass_psc_xsc_best_neighbour
            ON gaiadr3.tmass_psc_xsc_best_neighbour.original_ext_source_id = {0}.{1}
    """).format(dr7_vat_uploaded_id, source_id_col)
    job = Gaia.launch_job_async(query=query, name="DR3_2MASS_match")
    Gaia.upload_table_from_job(job)
    Gaia.delete_user_table(table_name=dr7_vat_uploaded_id, force_removal=True)
    Gaia.remove_jobs([job.jobid])
    job_table_id = "user_" + gaia_utils_standalone.username + ".t" + str(job.jobid)
    job_results = job.get_results()
    new_job_name = "DR3_2MASS_match"

    # Take this table, and retrieve the gaia columns.
    default_columns = ["ra", "dec",
                                                  "parallax", "parallax_error", "pmra",
                                                  "pmra_error", "pmdec", "pmdec_error", "radial_velocity",
                                                  "radial_velocity_error",
                       "ruwe"]
    default_columns += list(preupload_names) # append the default names, too- for the join.

    # Set up the string for selection col_1, col_2, col_3... ,col_n
    selection_string = default_columns[0]
    for i in range(1, len(default_columns)):
        selection_string += (",{}").format(default_columns[i])

    # Join our DR3_ID table with the main gaia_source table
    job = Gaia.launch_job_async(query=("""
                                          SELECT * 
                                          FROM {0}
                                          JOIN gaiadr3.gaia_source AS gaiadr3
                                              ON gaiadr3.source_id = {0}.source_id
                                          """).format(job_table_id),
                                name=new_job_name + "_1")
    Gaia.upload_table_from_job(job)
    Gaia.delete_user_table(table_name=job_table_id)
    Gaia.remove_jobs([job.jobid])
    joined_id = "user_" + gaia_utils_standalone.username + ".t" + job.jobid

    # Grab the columns we're interested in
    job = Gaia.launch_job_async(query=("""
                                          SELECT {0}
                                          FROM {1}
                                          """).format(selection_string, joined_id),
                                name=new_job_name + "_2")
    results = job.get_results()
    Gaia.remove_jobs([job.jobid])
    Gaia.delete_user_table(table_name=joined_id)
    results = Table(results)

    # Trim the table down based on RUWE
    ruwe_threshold = 1.4
    results = results[[True if ruwe < 1.4 else False for ruwe in results['ruwe']]]
    results.write(dr7_2mass_gaia_dir, overwrite=True)

# Should I subtract globular clusters & set up galactocentric coordinates?
do_galcent_gcs = False
produce_data = False
plot = False
min_radius = 15
max_radius = 150 # same as orbigistics interpolator for circularity
final_dir = os.path.join(windows_directories_new.lamodir, "LAMOST_master_2MASS.fits")
if __name__ == "__main__" and do_galcent_gcs == True:


    if produce_data:


        with fits.open(dr7_2mass_gaia_dir, mode='readonly', memmap=True) as hdul:
            # Data (which is an astropy table when loaded.)
            data = Table(hdul[1].data)

        data.rename_columns(['rv', 'erv', 'd_kpc_', 'e_d_kpc_'], ['vlos', 'evlost', 'dist', 'edist'])

        # Get the galactic coordinates from these equatorial ones
        from energistics_new import orbigistics
        from galcentricutils_new import angular

        orbigist = orbigistics()

        # Get GAL from ICRS
        data = orbigist.converter.nowrite_ICRS_to_GAL(data, has_cosfactor=True)

        # Get galactic errors. The cosfactor in the pmra/pmdec still exists- no longer in dmu_l, dmu_b though.
        import lamost_utils
        data = lamost_utils.monte_ICRSGAL_table(data)

        # Remove cosfactor from ICRS to match galconversion convention.
        data['pmra'] /= np.cos(np.deg2rad(data['dec']))
        data['pmra_error'] /= np.cos(np.deg2rad(data['dec']))

        # Get Galcent + Momenta
        data = orbigist.converter.nowrite_GAL_to_GALCENT(data)
        data = angular().get_momentum(data)

        # Remove nuisance stars (outside min_radius, but within the max_radius of interest.)
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
        data.write(final_dir, overwrite=True)

    if plot: # do some quality-of-life plots (angular momentum, energy, etc.)

        with fits.open(final_dir, mode='readonly', memmap=True) as hdul:
            # Data (which is an astropy table when loaded.)
            # "gaia_source_id" is dr2_source_id inside this.
            data = Table(hdul[1].data)

        # need to run orbifits.orbiliary(table)/generate the circularity column if you want this.
        # data = data[[True if 1 - abs(circ) > 0.9 else False for circ in data['circ']]]
        #data = data[[True if abs(d) < 300 else False for d in data['vx']]]
        #data = data[[True if abs(d) < 300 else False for d in data['vy']]]
        #data = data[[True if abs(d) < 300 else False for d in data['vz']]]
        #data = fast_energistics_new().default_E_c(data)

        # Quick conversion (not sure why this works- probs forgot process something in CSV/DAT digestion)
        data['x'], data['y'], data['z'], \
        data['vx'], data['vy'], data['vz'] = np.array(data['x'], float), \
                                             np.array(data['y'], float), \
                                             np.array(data['z'], float), \
                                             np.array(data['vx'], float), \
                                             np.array(data['vy'], float), \
                                             np.array(data['vz'], float)

        # Do an L-space preplot
        datatab = np.array([data['Lx'], data['Ly'], data['Lz']]).T
        graphutils_new.threed_graph().kmeans_L_array(datatab, [1 for d in datatab[:,0]],
                                                     False, browser=True, outliers=True)


        #writer = hdfutils.hdf5_writer(windows_directories_new.datadir, "stardata.hdf5")
        #data = writer.read_table("LAMOST_K_FULL_edr3", "astrotable")
        #data = data[[True if r > 15 else False for r in data['r']]]

        # Velocity-space histogram (just another quality check.)
        plt.hist(data['vx'], bins=1000, label='vx')
        plt.hist(data['vy'], bins=1000, label='vy')
        plt.hist(data['vz'], bins=1000, label='vz')
        plt.xlim([-1000,1000])
        plt.legend()
        plt.savefig(os.path.join(windows_directories_new.lamodir, "velspace_histogram_GAIA-GAIA.png"), dpi=300)

        # E-Lz plot, too.
        clustered = [1 for d in data['x']]
        data = fast_energistics_new().default_E_c(data)
        graphutils_new.spec_graph().energy_plot(data,
                                            clustered,
                                            list(set(clustered)),
                                            list(set(clustered)),
                                            [-10000, 10000],
                                            [-200000, 10000],
                                            "lamost_test.png",
                                            False,
                                            False)


# Should I grab metallicities from the main LRS based on object ID for the value-added-catalogue?
do_get_feh = False
new_final_dir = "LAMOST_master_2MASS_feh.fits"
if __name__ == "__main__" and do_get_feh == True:

    # Load in the 2MASS/etc fits
    with fits.open(final_dir, mode='readonly', memmap=True) as hdul:
        data = Table(hdul[1].data)

    # Load in the basic LRS
    with fits.open(windows_directories_new.LRS_path, mode='readonly', memmap=True) as hdul:
        # "gaia_source_id" is dr2_source_id inside this.
        fits_LRS = hdul[1].data

    # Convert the data objid column back to an integer (we turned it into a float when taking it from the dat file.)
    data['objid'] = np.array(data['objid'], int)

    # Match them (inefficient but whatever- expecting a 1-1 match. Some are unmatched- weird.)
    objids_lrs = np.array(fits_LRS['obsid'], int)
    @njit(fastmath=True)
    def match(labs1, labs2):

        """
        Match integers (by index) from labs1 with labs2, returning a list of tuples in which you have

        (index labs1 match, index labs2 match)

        :param labs1: int list array
        :param labs2: int list array
        :return: list tuples matches
        """

        tuples = np.zeros((len(labs1), 2), numba.types.int64)

        for i in range(0, len(labs1)): # if objids match then this gives a valid tuple, else both remain 0/ (0,0).
            for j in range(0, len(labs2)):
                if labs1[i] == labs2[j]:
                    tuples[i][0] = i
                    tuples[i][1] = j

        is_matched = [] # QoL check to ensure that match is had (0-0 match is false, even if true only loss of 1 point)
        for i in range(0, len(labs1)):
            if tuples[i][0] != 0:
                if tuples[i][1] != 0:
                    is_matched.append(True)
            else:
                is_matched.append(False)

        return tuples, is_matched

    tuples, is_matched = match(data['objid'], objids_lrs)
    tuples = tuples[is_matched]

    try: # if column exists this fails- just in case.
        data.add_column(0., name="feh")
        data.add_column(0., name="efeh")
    except:
        pass

    for tuple in tuples: # set radial velocities and fehs. weird column names (same as Jorges ASCII tho.)

        data[tuple[0]]['vlos'], data[tuple[0]]['evlost'], \
        data[tuple[0]]['feh'], data[tuple[0]]['efeh'] \
            = fits_LRS[tuple[1]]['rv'], fits_LRS[tuple[1]]['rv_err'], \
              fits_LRS[tuple[1]]['feh'], fits_LRS[tuple[1]]['feh_err']

    # Update the fits + save
    data.write(new_final_dir, overwrite=True)


# Should I now crossmatch to the old catalogue to remove duplicates, and then stack the tables? (whilst removing slag)
do_finalstack = False
finalstack_load = True # load matches instead, thus avoiding crossmatch comp cost
matchcache = os.path.join(windows_directories_new.lamodir, "finalstack_matchcache.txt")
finalstack_dir = os.path.join(windows_directories_new.lamodir, "LAMOST_final.fits")
finalstack_ascii = os.path.join(windows_directories_new.lamodir, "LAMOST_final.dat")
if __name__ == "__main__" and do_finalstack == True:

    # Load in the 2MASS/etc fits
    with fits.open(new_final_dir, mode='readonly', memmap=True) as hdul:
        # Data (which is an astropy table when loaded.)
        data = Table(hdul[1].data)

    # Add a column to specify the "SOURCE"
    data['source'] = "LAMOST_sstraszak"

    # Open up the HDF that has the LAMOST KGs
    import hdfutils
    from energistics_new import orbigistics
    from galcentricutils_new import angular
    writer = hdfutils.hdf5_writer(windows_directories_new.datadir, ascii_info_new.asciiname)
    lamokgs = writer.read_table(ascii_info_new.lamostk, ascii_info_new.set_raw)

    # convert to ICRS
    orbigist = orbigistics()
    lamokgs = lamokgs[[True if min_radius < r < max_radius else False for r in lamokgs['r']]]


    # Grab ra/decs from both (in radians)!!!!!!!!!
    lamokgs_ra, lamokgs_dec, \
    data_ra, data_dec = lamokgs['ra'], lamokgs['dec'], \
                        data['ra'], data['dec']
    lamokgs_ra, lamokgs_dec, \
    data_ra, data_dec = np.deg2rad(lamokgs_ra), np.deg2rad(lamokgs_dec), \
                        np.deg2rad(data_ra), np.deg2rad(data_dec)

    # Get matches between the older LAMOST KGs and the newer (obtained) LAMOST KGs + remove the rows that match.
    lamokgs_ra = lamokgs_ra
    lamokgs_dec = lamokgs_dec

    try:

        if finalstack_load == True:

            print("Trying to load matchcache...")

            with open(matchcache, "rb") as f:
                matches, truefalse = pickle.load(f)
                print("Successfully loaded the cachematch for finalstack matches.")

        else:

            raise RuntimeError("Will generate match instead...")

    except:

        from APOGEE_standalone import radecmatch_argmin_memlim

        matches, truefalse = radecmatch_argmin_memlim(lamokgs_ra, lamokgs_dec, data_ra, data_dec)
        matches, truefalse = list(matches), list(truefalse)
        matches, truefalse = np.asarray(matches), np.asarray(truefalse)

        with open(matchcache, "wb") as f:
            pickle.dump(obj=[matches, truefalse], file=f)
            print("Successfully generated the cachematch for final stack.")

    matches = matches[truefalse]
    matches = matches.T

    # Generate correlation plots for the matches to check integrity of catalogue source. Note that lamokgs is our
    # old data and "data" is the new data from the value-added-catalogue from LAMOST DR7.
    def do_correlation_plots(old_catalogue, new_catalogue,
                             columns, columns_errs,
                             columns_labels_old, columns_labels_new,
                             savedir):

        """

        Generate correlation plots between an old catalogue and a new catalogue, of same size.

        :param old_catalogue: astropy table or pandas
        :param new_catalogue: astropy table or pandas
        :param columns: list of table labels
        :param columns_errs: list of table err labels
        :param columns_labels_old: labels for matplotlib
        :param columns_labels_new: labels for matplotlib
        :param savedir: directory to save graphs
        :return: -

        """

        # Make savedir
        try:
            os.mkdir(savedir)
        except:
            pass

        # Go over each column and do the graphs/etc
        for column,column_err,\
            column_label_old,column_label_new in zip(columns,columns_errs,
                                                     columns_labels_old,columns_labels_new):

            # Make the save path
            savepath = os.path.join(savedir, column + "_correlation_plot.png")

            # Do the graph
            import graphutils_new
            graphutils_new.twod_graph().corr_plot(old_catalogue[column],
                                                  new_catalogue[column],
                                                  old_catalogue[column_err],
                                                  new_catalogue[column_err],
                                                  column_label_old,
                                                  column_label_new,
                                                  savepath)

    lamokgs_matches, data_matches = matches
    lamokgs_rows, data_rows = lamokgs[lamokgs_matches], data[data_matches]
    razip = np.array([lamokgs_rows['ra'], data_rows['ra']]).T
    columns = ["dmu_l", "dmu_b", "vlos", "dist"]
    columns_errs = ["edmu_l", "edmu_b", "evlost", "edist"]
    labels = [r'$\mu_l$', r'$\mu_b$', r'$\textrm{v}_\textrm{los}$', r'$\textrm{d}$']
    savedir = os.path.join(windows_directories_new.lamodir, "correlation_plots_lamost")
    do_correlation_plots(lamokgs_rows, data_rows,
                         columns, columns_errs,
                         labels, labels,
                         savedir)

    # Remove the rows that matched (from the old lamokgs)
    matches = matches[0]
    lamokgs.remove_rows(matches)

    # Now we can recombine with the rest of our dataset.
    # Load in the rest of the tables (sans LAMOST old)
    tables = [writer.read_table(d, ascii_info_new.set_raw) for d in ascii_info_new.all_groups[0:3]]
    tables += lamokgs # these have had match rows removed
    tables += data # add the new LAMOST data

    # Vstack
    data = vstack(tables)

    # Remove the nuisance columns (including objid/etc.)
    to_remove_columns = ["sgrstrong", "FeH_1", "Sgr Lambda", "Sgr Beta", "Belokurov Flag", "corrcoef",
                         "ra_2", "dec_2", "void", "void_1", "void_2", "name",
                         "plx_mas_", "eplx_mas_", "idup", "objid", "name_2mass", "flag2mass", "snrg",
                         "cvlos", "cdist"]
    data.remove_columns(to_remove_columns)

    # Find out which rows are nuisance rows
    to_remove = [False if np.isnan(d) == True else True for d in data['vx']]
    data = data[to_remove]
    to_remove = [False if np.isnan(d) == True else True for d in data['x']]
    data = data[to_remove]

    print("Saving the finalstack...")

    # Save
    data.write(finalstack_dir, overwrite=True)
    from astropy.io import ascii
    ascii.write(data, finalstack_ascii, overwrite=True)

"""
Load data in and do some cuts on it 
"""
do_cuts = True
cuts_fits_path = os.path.join(windows_directories_new.lamodir, "cut_LAMOST_final.fits")
if __name__ == "__main__" and do_cuts == True:

    # Load in the finalstack
    with fits.open(finalstack_dir, mode='readonly', memmap=True) as hdul:
        # Data (which is an astropy table when loaded.)
        orig_data = Table(hdul[1].data)

        # Specifically retrieve the new LAMOST data subsection
        indices = []
        for i in range(len(orig_data)):
            if orig_data[i]['source'] == "LAMOST_sstraszak":
                if np.abs(orig_data[i]['feh']) != 0.:
                    indices.append(i)

        # Set data
        data = orig_data[indices]

    # QUALITY CUTS ======================================^^^^
    from galcentricutils_quality_cuts import quality_cuts
    qcs = quality_cuts()
    metmax = -1/2
    LSR_VEL = 220
    FANCY_LSR_VEL = 100
    ZMAX = 5/2
    vphi_plotlims = [-550, 550]

    # Get indices for all the cuts necessary (from the original) dataset- for illustration purposes
    FSLR_KEEP, FSLR_REMOVE = qcs.fancy_LSR_cut(data, FANCY_LSR_VEL)
    RLSR_KEEP, RSLR_REMOVE = qcs.LSR_cut(data, LSR_VEL)
    FFEH_KEEP, FFEH_REMOVE = qcs.fancy_feh_cut(data)
    RFEH_KEEP, RFEH_REMOVE = qcs.feh_cut(data, metmax)
    ZMAX_KEEP, ZMAX_REMOVE = qcs.zmax_cuts(data, ZMAX,
                                 3e9,
                                 3000,
                                 True,
                                 True,
                                 "lamost_dr7_cuts_fulldata.txt")
    BOND_KEEP, BOND_REMOVE = qcs.bound_cut(data)
    indices_lists = [FSLR_KEEP, RLSR_KEEP, FFEH_KEEP, RFEH_KEEP, ZMAX_KEEP, BOND_KEEP]
    remove_lists = [FSLR_REMOVE, RSLR_REMOVE, FFEH_REMOVE, RFEH_REMOVE, ZMAX_REMOVE, BOND_REMOVE]
    indices_names = ["fancy_LSR", "regular_LSR", "fancy_feh", "regular_feh", "zmax", "bound"]

    # Get vT/etc anticipating plotting
    from energistics_new import orbigistics
    R, vR, vT, z, vz, phi = orbigistics().get_leftgalpy(data)
    feh = data['feh']

    import matplotlib.pyplot as plt
    for keep, remove, cut_type in zip(indices_lists, remove_lists, indices_names):

        fig, axs = plt.subplots(nrows=1, ncols=1)
        axs.set(xlim=vphi_plotlims,
                ylim=[-5 / 2, metmax])
        axs.grid(True, which='major', color='pink')
        axs.set(xlabel=r'$\textrm{v}_\phi$',
                ylabel=r'$[\textrm{Fe/H}]$')

        # Get the data we aren't cutting out + plot it
        VT_PLOT, FEH_PLOT = vT[keep], feh[keep]
        axs.scatter(VT_PLOT, FEH_PLOT, s=0.1, color='green')

        # Get the data we are cutting + plot it
        VT_PLOT, FEH_PLOT = vT[remove], feh[remove]
        axs.scatter(VT_PLOT, FEH_PLOT, s=0.1, color='red')

        plt.savefig(os.path.join(qcs.savedir, "finalstack_LAMOST_metplot_" + cut_type + ".png"), dpi=300)

    # Get non-lamost indices, recombine to LAMOST_sstraszak indices, to reclaim whole dataset
    no_indices = np.where(orig_data['source'] != "LAMOST_sstraszak")[0]
    orig_data = orig_data[no_indices]
    orig_data = vstack([orig_data, data])

    """
    # While we're at it, also grab the APOGEE data we obtained
    with fits.open(apogee_starhorse_path_nogcs, mode='readonly', memmap=True) as hdul:

        apodata = Table(hdul[1].data)
        apodata = apodata[[True if r < 150 else False for r in apodata['r']]]
        apodata['source'] = "APOGEE_sstraszak"
        orig_data = vstack([orig_data, apodata]) """

    # Generate base plot for debug
    R, vR, vT, z, vz, phi = orbigistics().get_leftgalpy(orig_data)
    orig_data['vT'] = vT
    plt.scatter(vT, orig_data['feh'], s=0.1)
    plt.ylim([-3,1])
    plt.savefig(os.path.join(qcs.savedir, "orig_data_beforecuts.png"), dpi=300)

    # Carry out cuts on orig_data (the full dataset) in anticipation of use
    # noinspection PyRedeclaration
    KEEP, REMOVE = qcs.fancy_LSR_cut(orig_data, FANCY_LSR_VEL)
    orig_data = orig_data[KEEP]
    # noinspection PyRedeclaration
    KEEP, REMOVE = qcs.zmax_cuts(orig_data, ZMAX,
                                 3e9,
                                 3000,
                                 True,
                                 True,
                                 "lamost_dr7_cuts_origdata.txt")
    orig_data = orig_data[KEEP]
    plt.clf()
    plt.close()
    plt.cla()
    plt.scatter(orig_data['vT'], orig_data['feh'], s=0.1)
    plt.xlim([-500,500])
    plt.ylim([-3,1])
    print("Showing")
    plt.show()
    # noinspection PyRedeclaration
    #KEEP, REMOVE = qcs.fancy_feh_cut(orig_data)
    #orig_data = orig_data[KEEP]


    # Save data to fits
    orig_data.write(cuts_fits_path, overwrite=True)


    # Write data to table
    import hdfutils
    writer = hdfutils.hdf5_writer(windows_directories_new.datadir, ascii_info_new.asciiname)
    writer.write_table(ascii_info_new.fullgroup, ascii_info_new.fullset, orig_data)

    # Set up the flat hdbscan run
    data_array = np.array([orig_data['Lx'],orig_data['Ly'],orig_data['Lz']]).T
    from hdbscan import flat
    clusterer = flat.HDBSCAN_flat(X=data_array,
                                  n_clusters=40,
                                  min_cluster_size=10,
                                  max_cluster_size=int(0.4*len(orig_data)),
                                  min_samples=7,
                                  metric='l2',
                                  algorithm='best')#,
                                  #prediction_data=False)
    set_of_labels = list(set(clusterer.labels_))
    sizes = [len(np.where(clusterer.labels_==d)[0]) for d in set_of_labels]
    sorted_sizes = np.sort(sizes, axis=0)
    GSE_SIZE = np.flip(sorted_sizes)[1]

    # Take our cleaned data sample and produce angular plot
    savepath = os.path.join(qcs.savedir, "kmeans_L_LAMOGEEFINAL_test.html")
    graphutils_new.threed_graph().kmeans_L_array_path(data_array,
                                                      clusterer.labels_,
                                                      savepath,
                                                      True,
                                                      True,
                                                      True)

    # Savedir
    savedir = os.path.join(qcs.savedir, "orbifits_tests")
    try:
        os.mkdir(savedir)
    except:
        pass

    # For all our clusters found :3
    run_orbifitplots = False
    if run_orbifitplots == True:

        for clust in list(set(clusterer.labels_)):
            import energistics, graphutils
            orbifit = energistics.orbifitter().finetune_galpy_fitting_nomemb(orig_data,
                                                                             clusterer.labels_,
                                                                             clust,
                                                                             2000,
                                                                             0.6e9,
                                                                             1000,
                                                                             False,
                                                                             False,
                                                                             False,
                                                                             False)

            # Forward Integral
            import copy
            from astropy import units as u

            forward = copy.deepcopy(orbifit)
            backward = copy.deepcopy(orbifit)
            forward.integrate((np.linspace(0, 1e9, 2000) * u.yr), energistics.orbigistics().pot)
            backward.integrate((np.linspace(0, -1e9, 2000) * u.yr), energistics.orbigistics().pot)
            llsf, bbsf, llsb, bbsb = forward.ll((np.linspace(0, 0.7e9, 2000) * u.yr)).value, \
                                     forward.bb((np.linspace(0, 0.7e9, 2000) * u.yr)).value, \
                                     backward.ll((np.linspace(0, -0.15e9, 2000) * u.yr)).value, \
                                     backward.bb((np.linspace(0, -0.15e9, 2000) * u.yr)).value
            lls, bbs = np.concatenate([llsf, llsb]), np.concatenate([bbsf, bbsb])
            lls = [d - 360 if d > 180 else d for d in lls]

            specgrapher = graphutils.spec_graph()
            fig, axs = specgrapher.lb_orbits(orig_data[[True if d == clust
                                                        else False
                                                        for d in clusterer.labels_]],
                                             0.3e9, [-180, 180],
                                             [-90, 90], None,
                                             line=False, points=4000)
            axs.scatter(lls, bbs, color='lime', marker='o', s=3)
            # Misc axis things
            axs.set_facecolor("k")
            axs.grid(True, which='major', alpha=1, linewidth=0.25, color='white')
            axs.grid(True, which='minor', alpha=1, linewidth=0.25, color='white')
            axs.grid(color="white")
            savepath = os.path.join(savedir, str(clust) + ".png")
            import matplotlib.pyplot as plt

            plt.savefig(savepath, dpi=300, transparent=False)
            plt.close()

    # Check mets
    for clust in list(set(clusterer.labels_)):

        # Get metallicity + dispersion
        indices = np.where(clusterer.labels_==clust)[0]
        fehs = orig_data['feh'][indices]
        from astropy.stats import sigma_clipped_stats
        mean, med, std = sigma_clipped_stats(fehs)
        print(clust, mean, med, std)
