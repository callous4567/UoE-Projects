import os
import h5py
import nexusformat.nexus as nx
import astropy
import numpy as np
import pandas

# Modified from TGP2020 utils.hdf5_writer. Windows or Unix compatible. Class holds various HDF5 utilities.
class hdf5_writer(object):
    def __init__(self, directory, filename):
        self.directory = directory
        self.filename = filename
        self.owd = os.getcwd()
    # Creates the file if it doesn't exist, does nothing otherwise.
    def create(self):
        os.chdir(self.directory)
        file_created = h5py.File(self.filename, 'a')
        file_created.close()
        os.chdir(self.owd)
    # Returns the ENTIRE FILE TREE of the file
    def info(self):
        os.chdir(self.directory)
        f = nx.nxload(self.filename)
        os.chdir(self.owd)
        return f.tree
    # Returns the groups, only
    def group_info(self):
        os.chdir(self.directory)
        with h5py.File(self.filename, "a") as file:
            os.chdir(self.owd)
            return list(file.keys())
    # Writes ARRAY to the DATASET inside GROUP of FILENAME.
    def write(self, group, dataset, array):
        os.chdir(self.directory)
        with h5py.File(self.filename, 'a') as f:
            # Attempt a replace.
            try:
                dataframe = f[group + "/" + dataset]
                dataframe[...] = array
                os.chdir(self.owd)
            # Replace failed. Try to make group.
            except:
                try:
                    f.create_group(group)
                    os.chdir(self.owd)
                # Group exists (or doesn't... Attempt to delete (if it exists) + (re)make the set
                except:
                    try:
                        del f[group][dataset]
                    except:
                        pass
                    try:
                        f.create_dataset(group + "/" + dataset, shape=np.shape(array), data=array)
                        os.chdir(self.owd)
                    except Exception as e:
                        print("All ARRAY WRITING methods failed... " + str(e) + "..." + group + "_" + dataset)
                        os.chdir(self.owd)
    # Write single strings, duplicate of above. Split to save processing power since usage of strings isn't likely.
    def writestr(self, group, dataset, string):
        os.chdir(self.directory)
        with h5py.File(self.filename, 'a') as f:
            # Attempt a replace.
            stringer = string.encode("ascii", "ignore")
            try:
                dataframe = f[group + "/" + dataset]
                dataframe[...] = stringer
                os.chdir(self.owd)
            # Replace failed. Try to make group.
            except:
                try:
                    f.create_group(group)
                    os.chdir(self.owd)
                # Group exists. Attempt to delete + remake the set
                except:
                    try:
                        del f[group][dataset]
                    except:
                        pass
                    try:
                        f[group].create_dataset(dataset, shape=(1, 1), dtype='S10', data=stringer)
                        os.chdir(self.owd)
                    except Exception as e:
                        print("All STRING writing methods failed... " + str(e) + "..." + group + "," + dataset)
                        os.chdir(self.owd)
    # Delete group, set
    def del_set(self, group, dataset):
        os.chdir(self.directory)
        with h5py.File(self.filename, 'a') as f:
            del f[group, dataset]
        os.chdir(self.owd)
    # Write an astropy table
    def write_table(self, group, dataset, astropytable):
        os.chdir(self.directory)
        astropy.io.misc.hdf5.write_table_hdf5(astropytable, self.filename, path=(group + "/" + dataset), append=True,overwrite=True, serialize_meta=True)
        os.chdir(self.owd)
    # Read an astropy table
    def read_table(self, group, tablename):
        os.chdir(self.directory)
        tab = astropy.io.misc.hdf5.read_table_hdf5(input=self.filename, path=(group + "/" + tablename), character_as_bytes=False)
        os.chdir(self.owd)
        return tab
    # Reads ARRAY from DATASET inside GROUP of FILENAME
    def read(self, group, dataset):
        os.chdir(self.directory)
        data = "null"
        with h5py.File(self.filename, 'r') as f:
            # As a np array.
            data = np.array(f[group][dataset])
        os.chdir(self.owd)
        return data
    # Read a string instead, singular, not tested on arrays.
    def read_string(self, group, dataset):
        os.chdir(self.directory)
        data = "null"
        with h5py.File(self.filename, 'r') as f:
            # As a string (well, bytes.)
            data = str(np.array(f[group][dataset]), "ascii")
        os.chdir(self.owd)
        return data
    # Write a Pandas Dataframe
    def write_df(self, group, dataset, df):
        os.chdir(self.directory)
        df.to_hdf(self.filename, key = group + "/" + dataset, complevel=0, append=False, format="fixed", mode='r+')
        os.chdir(self.owd)
    # Read a Pandas DataFrame
    def read_df(self, group, dataset):
        os.chdir(self.directory)
        df = pandas.read_hdf(self.filename, key=group + "/" + dataset)
        os.chdir(self.owd)
        return df
