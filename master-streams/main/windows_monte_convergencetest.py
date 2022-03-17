import os
import time
import matplotlib.ticker as pltick
from matplotlib import rc
from matplotlib.ticker import MaxNLocator
from numpy import random
import galcentricutils, hdfutils, windows_directories, ascii_info
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing
plt.rcParams['animation.ffmpeg_path'] = 'C:\\Users\\Callicious\\Documents\\Prog\\pycharm\\venv\\ffmpeg\\bin\\ffmpeg.exe'
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
from matplotlib.patches import Patch
""" 
Monte-carlo Convergence Test
"""
# Get tables (with angular momentum) for our datasets (individually.)
def get_L():
    # Grab Data
    writer = hdfutils.hdf5_writer(windows_directories.datadir, ascii_info.asciiname)
    groups = ascii_info.all_groups
    tables = [writer.read_table(group, ascii_info.set_raw) for group in groups]

    # Grab angular momenta
    tables = [galcentricutils.angular().get_momentum(table) for table in tables]

    return tables

# Select n random rows from each table, return rows (useful for testing Monte-Carlo estimation.
def get_n(tables, n):
    random_integers = [random.default_rng().integers(0, high=len(d), size=n, dtype=np.int64, endpoint=False)
                       for d in tables]
    random_tablerow = [tables[d][random_integers[d]] for d in range(len(tables))]
    return random_tablerow

# For each table, run a convergence test on the Monte-Carlo estimation in errors.
"""
For each row:
- Monte-Carlo for a range of (n) samples from [0,max_samples]
- Ensure same range is used for all tables, with small step in [0->max_samples] via arange
Once mu/var(n) for each row has been obtained, get fractional difference to final average (mu/var(n_max)) 
- Do for each row
- "final average" over some avg_length
Overplot fractional difference in mu/var for the table (i.e. combine row tables into 1 big table of all rows) 
- Repeat for each table. 
See which (n) provides the best fractional convergence/within tolerance. 
"""
# Does the above for one table. Domain of n is split in 2: lower half gets n_step, upper gets 2*n_step.
# avg_length defines the number of final elements to average over to define the "good estimate"
def monte_converge(table, table_name, n_max, n_step, avg_length):
    # Set up the table/environment for plotting, alongside axis limits.
    """
    Left side = Fractional difference to final mean, Right side - stdev
    Top: Lx, Bot: Lz
    """
    # Set up domain of (n) to Monte over and hence x axes
    lower, upper = np.arange(1, int(n_max / 2), int(n_step / 2)), np.arange(int(n_max / 2), n_max, n_step)
    full = np.append(lower, upper)
    ylims = [-1,1]
    numPlotsY = 3
    numPlotsX = 2

    # Set up figure/axes and subplot adjust
    fig, axs = plt.subplots(nrows=numPlotsY, ncols=numPlotsX, sharey='none', squeeze=True, sharex='none')

    # Define the monte object
    monte = galcentricutils.monte_angular()
    strs = 'l','b','dist','dmu_l','dmu_b','vlos','edist','edmu_l','edmu_b','evlost'

    # Monte over the domain for each row and plot data
    for row in table:
        vec = np.array([row[str] for str in strs])
        monte_results = np.array([monte.vec_monte(vec, n) for n in full]) # each monte is [mux,muy,muz,varx,vary,varz]
        avg_montvecs = monte_results[-1-avg_length:-1]
        avg_finalvec = np.average(avg_montvecs, axis=0)
        monte_fracdif = (monte_results - avg_finalvec)/avg_finalvec # fracdif(n)[mux,muy...]
        monte_fracdif = np.split(monte_fracdif.T, [3])
        for rowindex in range(numPlotsY):
            for colindex in range(numPlotsX):
                axs[rowindex,colindex].plot(full, monte_fracdif[colindex][rowindex])

    # Sort out shit ticks: https://stackoverflow.com/questions/37970334/remove-only-overlapping-ticks-in-subplots-grid
    for y_i in range(0, numPlotsY):
        for x_i in range(0, numPlotsX):
            ax = axs[y_i, x_i]

            # If left-most column:  Remove all overlapping y-ticks
            # Else:                 Remove all ticks
            if x_i == 0:
                # If top left subplot:      Remove bottom y-tick
                # If bottom left subplot:   Remove top y-tick
                # Else:                     Remove top and bottom y-ticks
                if y_i == 0:
                    nbins = len(ax.get_yticklabels())
                    ax.yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='lower'))
                elif y_i == numPlotsY - 1:
                    nbins = len(ax.get_yticklabels())
                    ax.yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))
                else:
                    nbins = len(ax.get_yticklabels())
                    ax.yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='both'))
            else:
                ax.yaxis.set_ticks([])

            # If bottom row:    Remove all overlapping x-ticks
            # Else:             Remove all ticks
            if y_i == numPlotsY - 1:
                # If bottom left subplot:   Remove right x-tick
                # If bottom right subplot:  Remove top left x-tick
                # Else:                     Remove left and right x-ticks
                if x_i == 0:
                    nbins = len(ax.get_xticklabels())
                    ax.xaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))
                elif x_i == numPlotsX - 1:
                    nbins = len(ax.get_xticklabels())
                    ax.xaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='lower'))
                else:
                    nbins = len(ax.get_xticklabels())
                    ax.xaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='both'))
            else:
                ax.xaxis.set_ticks([])

    # Set up texts for axis
    left_texts = [r'$\mu_{L_x}$', r'$\mu_{L_y}$', r'$\mu_{L_z}$']
    right_texts = [r'$\sigma_{L_x}$', r'$\sigma_{L_y}$', r'$\sigma_{L_z}$']
    texts = [left_texts, right_texts]


    # Sort out lims and text
    for i in range(3):
        for j in range(2):
            axs[i,j].set(ylim=ylims)
            axs[i,j].text(s=texts[j][i], x=0.9, y=0.15, va="top", ha="left", transform=axs[i,j].transAxes, wrap=True)

    fig.suptitle(("Monte Carlo Fractional Error for " + r'$\mu_{L_i}$' + " \& " + r'$\sigma_{L_i}$' + ("\n{0} data,").format(table_name.replace("_", "\_"))) + (" for {} ").format(len(table)) + "randomly sampled row elements")
    fig.subplots_adjust(top=0.875)
    plt.subplots_adjust(wspace=0, hspace=0)

    # X and Y axis labels
    x,y = ["(n) samples", "frac err"]
    fig.text(0.5, 0.02, x, ha='center')
    fig.text(0.04, 0.5, y, va='center', rotation='vertical')

    # Define filename
    filename = table_name + ".monte_converge.png"
    # Try to save.
    try:
        os.mkdir(windows_directories.imgdir + "\\monte_converge")
        plt.savefig(windows_directories.imgdir + "\\monte_converge\\" + filename, dpi=300)
    except:
        plt.savefig(windows_directories.imgdir + "\\monte_converge\\" + filename, dpi=300)

    plt.show()

# Run a convergence test to see how many steps are roughly needed for the monte-carlo estimate to converge to a final.
def converge_test():
    test_rows = get_n(get_L(),40)
    names = ascii_info.all_groups
    for num, group_rows in enumerate(test_rows):
        test_monte = monte_converge(group_rows, names[num], 1000, 20, 10)

# Run the convergence test
converge_test()