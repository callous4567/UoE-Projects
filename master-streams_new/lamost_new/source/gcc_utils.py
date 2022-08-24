import numba
from astropy.table import Table
from numba import njit

import windows_directories_new
import numpy as np
import os

do_setup = False
def setup_master():

    # Globular Cluster Catalogue (for tidal radius info)
    # https://people.smp.uq.edu.au/HolgerBaumgardt/globular/
    # baumgardt_ascii.txt is name for raw
    filename = "baumgardt_ascii_edited.csv"
    with open(os.path.join(windows_directories_new.lamodir, filename), "r") as f:
        lines = f.readlines()
    lines = [d.replace("\n", "") for d in lines]
    lines = [d.split(";") for d in lines]
    columns = lines[0]
    columns = [d.strip() for d in columns]
    columns = columns[1:]
    name = [line[0] for line in lines[2:]]
    name = np.array(name, str)
    data = [line[1:] for line in lines[2:]]
    data = np.array([np.array(d, float) for d in data]).T
    gc_table = Table()
    gc_table["name"] = name
    for num, column in enumerate(columns):
        gc_table[column] = data[num]

    # Get the orbital parameter table, too.
    filename = "baumgardt_orbits_table_edited.csv"
    with open(os.path.join(windows_directories_new.lamodir, filename), "r") as f:
        lines = f.readlines()
    lines = [d.replace("\n", "") for d in lines]
    lines = [d.split(";") for d in lines]
    columns = lines[0]
    columns = [d.strip() for d in columns]
    columns = columns[1:]
    name = [line[0] for line in lines[2:]]
    name = np.array(name, str)
    data = [line[1:] for line in lines[2:]]
    data = np.array([np.array(d, float) for d in data]).T
    orbi_table = Table()
    orbi_table["name"] = name
    for num, column in enumerate(columns):
        orbi_table[column] = data[num]

    # Assuage which names are in which (and add appropriate columns)
    gc_table_names = gc_table.columns
    orbi_table_nam = orbi_table.columns
    not_in = []
    for colname in orbi_table_nam:

        if colname not in gc_table_names:

            gc_table.add_column(float(0), name=colname)
            not_in.append(colname)

    # Next, for each row in the orbitable, find the name in the gc_table, and add in the values
    for row in orbi_table:

        where_in_gc = np.where(gc_table['name'] == row['name'])[0]
        if len(where_in_gc) == 1:

            for colname in not_in:
                gc_table[where_in_gc[0]][colname] = row[colname]

    # Dirtily set up a flag for if the GC is "usable" and clip to usables
    gc_table['usable_flag'] = [True if d != 0
                               else False for d in gc_table['Rsun']]
    gc_table = gc_table[gc_table['usable_flag']]

    # Rename columns to match our usual format
    ascii_names = ['ra','dec','dist']
    gc_names = ['RA','DEC','R_Sun']
    for asciiname, gcname in zip(ascii_names, gc_names):
        gc_table.rename_column(gcname, asciiname)

    # Save the final table as a FITS (for inspection purposes)
    gc_table.write(windows_directories_new.baumgardt_fits, overwrite=True)

if do_setup == True:
    setup_master()

# Threshold of number of tidal radii within which stars are removed
rt_mult = 3
@njit(fastmath=True)
def remove_gcs(x_s,y_s,z_s,
               x_gc,y_gc,z_gc,
               _rt_gcs):

    _xyz_stars = np.empty((len(x_s), 3), numba.types.float64)
    _xyz_gcs = np.empty((len(x_gc), 3), numba.types.float64)
    _xyz_stars[:,0], _xyz_stars[:,1], _xyz_stars[:,2] = x_s,y_s,z_s
    _xyz_gcs[:,0], _xyz_gcs[:,1], _xyz_gcs[:,2] = x_gc,y_gc,z_gc

    # True/False for whether to retain/keep the star (if it's in a GC, then False. Default True)
    retain = np.ones(len(_xyz_stars), numba.types.boolean)

    # Go through all the globular clusters
    for i in range(1, len(_xyz_gcs)):

        # Go through all the stars
        for j in range(1,len(_xyz_stars)):

            # If within rt_mult tidal radii, then do not retain
            if np.sqrt(np.sum((_xyz_gcs[i] - _xyz_stars[j])**2)) < rt_mult * _rt_gcs[i]:

                retain[j] = False

    # Return
    return retain

