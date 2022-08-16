import csv
import pickle
import time

import pandas as pd
from astropy.io import fits
from astroquery.gaia import Gaia

import ascii_info
import hdfutils
import os
import numpy as np
from astropy.table import Table, vstack

# Brief information on format of ASCII files.
import windows_directories

"""
First row is parameter space, semicolon ; delimited, with units in square brackets []. 
- Open ASCII
- Format first row accordingly (string format/title column for astropy table)
- All succeeding rows are rows of table/values. 
asici_port object is the file manipulator object for each individual ascii file, given a directory and filename.
hdf5 directory must bs given to it, too. 
"""

# ASCII .txt file manager: allows import of a single ASCII file to HDF5 format, as an astropy table.
class ascii_port(object):
    def __init__(self, directory, asciiname):
        self.directory = directory
        self.asciiname = asciiname
        self.table = Table()

    # Grab the ASCII file and convert to Astropy Table, in accordance to above. Saves to self.table.
    def astrotabify(self):
        # Various file operations/placeholders
        os.chdir(self.directory)
        file = open(self.asciiname)
        file_read = file.read().splitlines()
        header = False
        header_units = False
        data = []
        # Grab the name of each column, and separately its units.
        for num, line in enumerate(file_read):
            if num == 0:
                line_split = line.split(";")
                line_split_second = []
                for num,item in enumerate(line_split):
                    split = item.split("[")
                    split_second = "[" + split[1].strip()
                    line_split[num] = split[0].strip()
                    # Catch empty column names in ascii: label void.
                    if line_split[num] == "":
                        line_split[num] = "void"
                    line_split_second.append(split_second)
                header = line_split
                header_units = line_split_second
            else:
                line = line.split(" ")
                for num, line_item in enumerate(line):
                    # Catch string error when attempting float conversion
                    try:
                        line[num] = float(line_item)
                    except:
                        pass
                data.append(line)
        # Check for duplicate column names in the header list: duplicates get a _1,2,3, etc etc.
        non_duplicates, indices = [], []
        for num,item in enumerate(header):
            if item in non_duplicates:
                indices[non_duplicates.index(item)] += 1
                header[num] = item + "_" + str(indices[non_duplicates.index(item)])
            else:
                non_duplicates.append(item)
                indices.append(int(0))


        # Write table, one column at a time to ensure we don't hit any strings.
        data_table = Table()
        data_columns = []
        for i in zip(*data):
            data_columns.append(list(i))
        for num, column in enumerate(data_columns):
            data_table[header[num]] = column


        # Add a column to the table for the source/etc (for later vstacking/etc.)
        data_table['source'] = self.asciiname

        # Save table as self.
        self.table = data_table

    # Add an astropy table to this one (vstack). Changes self.table to the new table.
    # If the new table has unknown columns to the old one, "null" is filler.
    def add(self, addition):
        # Grab the table for the addition, and grab the columns that match up with the columns for the current ascii.
        addition_table = addition.table

        # Grab the names in the addition that aren't in the current ascii
        not_in_current = []
        for name in addition_table.colnames:
            if name not in self.table.colnames:
                not_in_current.append(name)

        # Grab the names in the current ascii that aren't in the addition
        not_in_new = []
        for name in self.table.colnames:
            if name not in addition_table.colnames:
                not_in_new.append(name)

        # Add the names in the addition that aren't in the current to the current. Make sure types are compatible.
        for name in not_in_current:
            type_of = type(addition_table[name][0])
            self.table[name] = np.zeros(np.shape(self.table[self.table.colnames[0]]), dtype=type_of)

        # Add the names in the current that aren't in the addition to the addition. Make sure types are compatible.
        for name in not_in_new:
            type_of = type(self.table[name][0])
            addition_table[name] = np.zeros(np.shape(addition_table[addition_table.colnames[0]]), dtype=type_of)

        # Vstack the tables, addition underneath the current ascii.
        self.table = vstack([self.table, addition_table])

    # Save this ascii_port objects table
    def save(self, hdfdir, hdfname, hdfgroup, hdfset):
        # Set up hdf5 file.
        writer = hdfutils.hdf5_writer(hdfdir, hdfname)
        try:
            writer.create()
        except:
            pass
        writer.write_table(hdfgroup, hdfset, self.table)

    # Take an astropy table and write an ascii file (basically, the reverse! :D :D :D)
    def write(self, table, labels):

        # Clip table
        table = table[labels]

        # Write table
        table.write(self.directory + "\\" + self.asciiname, format='ascii', overwrite=True)



"""
# Load Data
table = hdfutils.hdf5_writer(windows_directories.datadir,
                             ascii_info.asciiname).read_table(ascii_info.fullgroup,
                                                              ascii_info.fullset)
asciiwriter = ascii_port(windows_directories.datadir, "full_catalogue.txt")
asciiwriter.write(table, ["x", "y", "z", "vx", "vy", "vz", "l", "b", "dist", "vlos", "dmu_l", "dmu_b", "edist",
                          "evlost", "edmu_l", "edmu_b", "FeH", "Sgr Lambda", "Sgr Beta", "Belokurov Flag", "corrcoef",
                          "source", "Lx", "Ly", "Lz", "L", "r", "theta", "phi", "dLx", "dLy", "dLz", "prelim_clust"]) 
                          """