from astropy.table import vstack
import ascii_info
import hdfutils
import numpy as np
import windows_directories

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

    # Stack all dataframes
    def dataframes(self):
        # While we're at it, recreate the full table. This will stack all our astropy tables.
        writer = hdfutils.hdf5_writer(windows_directories.datadir, ascii_info.asciiname)
        dfs = [writer.read_df(group, ascii_info.panda_raw) for group in ascii_info.all_groups]
        for num, df in enumerate(dfs):
            if num != 0:
                dfs[0] = dfs[0].append(df, ignore_index=True, sort=False)
        writer.write_df(ascii_info.fullgroup, ascii_info.fullpanda, dfs[0])


