import os
from astropy.io import fits
import numpy as np
import windows_directories_new
import pickle

num_match_stars = 10000 # number of stars for dr7-4 comparison
load = True # should I load the match?
dr4dir = os.path.join(windows_directories_new.lamodir, "LAMOST_value_addcat_DR4.fits") # where is VAC for dr4
dr7dir = os.path.join(windows_directories_new.lamodir, "dr7_LRS_vat.fits") # where is VAC for dr7

with fits.open(dr4dir, mode='readonly', memmap=True) as hdul:
    dr4data = hdul[1].data
with fits.open(dr7dir, mode='readonly', memmap=True) as hdul:
    dr7data = hdul[1].data

if load != True: # generate the match

    # Select a random assortment from dr7data (num_match_stars of them.)
    selection = np.random.randint(0, len(dr7data), num_match_stars)
    dr7data_subset = dr7data[selection]

    # Only need the ra/decs anyway.
    ra1, dec1, ra2, dec2 = dr7data_subset['ra'], dr7data_subset['dec'], dr4data['RA'], dr4data['DEC']
    ra1, dec1, ra2, dec2 = np.radians(ra1), np.radians(dec1), np.radians(ra2), np.radians(dec2)
    dr4data, dr7data = None, None

    from APOGEE_standalone import radecmatch_argmin_memlim
    matches, truefalse = radecmatch_argmin_memlim(ra1, dec1, ra2, dec2)
    """
    radecmatch_argmin_memlim uses argmin 1-by-1 to generate matches (slow process but avoids generating an N^2 matrix.)
    if you want a more conclusive match, use radecmatch
    this will use the munkres algorithm at the cost of generating an N^2 matrix (computationally expensive/memory exp.)
    """
    matches, truefalse = list(matches), list(truefalse) # julia types not usable
    matches, truefalse = np.asarray(matches), np.asarray(truefalse)
    found = np.where(truefalse==True)[0]
    with open(os.path.join(windows_directories_new.lamodir, "matchfound.txt"), 'wb') as f:
        pickle.dump(obj=[selection, matches, found], file=f)

    load = True

if load == True: # plot/compare match graphically

    with fits.open(dr4dir, mode='readonly', memmap=True) as hdul:
        dr4data = hdul[1].data
    with fits.open(dr7dir, mode='readonly', memmap=True) as hdul:
        dr7data = hdul[1].data

    with open(os.path.join(windows_directories_new.lamodir, "matchfound.txt"), 'rb') as f:
        selection, matches, found = pickle.load(file=f)

    # Clip the matches/etc/data/etc appropriately
    matches = matches.T[1]
    dr7data_subset = dr7data[selection][found]
    dr4data = dr4data[matches][found]
    truefalse = np.where(dr4data['DIST']>0.1)[0] # remove odd artefact (bad datapoints with DIST ~= 0)
    dr7data_subset=dr7data_subset[truefalse]
    dr4data=dr4data[truefalse]

    dr4data['DIST']/=1000 # convert to kpc
    dr4data['ERR_DIST']/=1000

    # Columns for comparison
    columns = ["d_kpc_", "RV"]
    columns_2 = ["DIST", "VR"]
    columns_errs = ["e_d_kpc_", "eRV"]
    columns_errs_2 = ["ERR_DIST", "ERR_VR"]
    labels = ["d", r'$v_{los}$']
    savedir = os.path.join(windows_directories_new.lamodir, "correlation_plots_lamost_dr4-7")
    def do_correlation_plots(old_catalogue, new_catalogue,
                             columns, columns_errs,
                             columns_2, columns_errs_2,
                             columns_labels_old, columns_labels_new,
                             savedir, lims):

        """

        Generate correlation plots between an old catalogue and a new catalogue, of same size.

        :param old_catalogue: astropy table or pandas
        :param new_catalogue: astropy table or pandas
        :param columns: list of table labels
        :param columns_errs: list of table err labels
        :param columns_2: list of table labels for second table
        :param columns_errs_2: list of table err labels for second table
        :param columns_labels_old: labels for matplotlib
        :param columns_labels_new: labels for matplotlib
        :param savedir: directory to save graphs
        :param lims: list/array of tuples of limits for each (set None to default to sigclip.) See twod_graph.corr_plot
        :return: -

        """

        # Make savedir
        try:
            os.mkdir(savedir)
        except:
            pass

        # Go over each column and do the graphs/etc
        for column,column_err,\
            column_2,column_err_2,\
            column_label_old,column_label_new, lim in zip(columns, columns_errs,
                                                     columns_2, columns_errs_2,
                                                     columns_labels_old, columns_labels_new, lims):

            # Make the save path
            savepath = os.path.join(savedir, column + "_correlation_plot.png")

            # Do the graph
            import graphutils_new
            graphutils_new.twod_graph().corr_plot(old_catalogue[column],
                                                  new_catalogue[column_2],
                                                  old_catalogue[column_err],
                                                  new_catalogue[column_err_2],
                                                  column_label_old,
                                                  column_label_new,
                                                  savepath,
                                                  lim)

    do_correlation_plots(dr7data_subset,dr4data,
                         columns, columns_errs,
                         columns_2, columns_errs_2,
                         labels, labels,
                         savedir,
                         [
                             [[0,5],[0,5]],
                             [[-300,300],[-300,300]]
                         ])