import time
from astropy.table import vstack, Table
import ascii_info
import hdfutils
import numpy as np
import windows_directories

"""
Note that this no longer just has stacking routines but some other combination routines.
Guess that makes the name a bit redundant, hm. 
"""

# Class for stacking up all the individual sets of data according to ascii_info (vertically)
class stacking(object):
    def __init__(self):
        self.null = False

    # Stack all astropy tables
    def tables(self):
        # While we're at it, recreate the full table. This will stack all our astropy tables.
        writer = hdfutils.hdf5_writer(windows_directories.datadir, ascii_info.asciiname)
        tables = [writer.read_table(group, ascii_info.set_raw) for group in ascii_info.all_groups]
        for num,table in enumerate(tables):
            if num != 0:
                tables[0] = vstack([tables[0],table])
        writer.write_table(ascii_info.fullgroup, ascii_info.fullset, tables[0])

    # Stack all dataframes from monte multitables/panda_raw
    def dataframes(self):
        # While we're at it, recreate the full table. This will stack all our astropy tables.
        writer = hdfutils.hdf5_writer(windows_directories.datadir, ascii_info.asciiname)
        dfs = [writer.read_df(group, ascii_info.panda_raw) for group in ascii_info.all_groups]
        for num, df in enumerate(dfs):
            if num != 0:
                dfs[0] = dfs[0].append(df, ignore_index=True, sort=False)
        writer.write_df(ascii_info.fullgroup, ascii_info.fullpanda, dfs[0])


# Class for converting dataframes to astropy tables (quick cleanup after monte multitables)
# Removes covtrix (since we can't nest that into an astropy table.)
class dataframes_to_tables():
    def __init__(self):
        groups = ascii_info.all_groups
        astropysets = [ascii_info.set_raw for group in groups]
        pandasets = [ascii_info.panda_raw for group in groups]
        writer = hdfutils.hdf5_writer(windows_directories.datadir, ascii_info.asciiname)
        # Some housework: convert to Astropy Table, and overwrite the astropy source.
        for group, astropyset, pandaset in zip(groups, astropysets, pandasets):
            df = writer.read_df(group, pandaset)
            # Drop the various arrays we want to keep inside the dataframe but that astropy won't tolerate
            df = df.drop(columns="covtrix")
            df = df.drop(columns="vec_L")
            df = df.drop(columns="covtrix2")
            df = df.drop(columns="vec_4d")
            df = df.drop(columns="covtrix3")
            df = df.drop(columns="vec_4dLE")
            df = df.drop(columns="covtrix6D")
            df = df.drop(columns="vec_6D")
            table = Table.from_pandas(df)
            writer.write_table(group, astropyset, table)

