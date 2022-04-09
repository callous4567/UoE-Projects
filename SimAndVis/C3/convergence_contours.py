import itertools
import multiprocessing
import os
import warnings
from astropy.table import Table
from matplotlib import pyplot as plt, rc

import hdfutils
from SimAndVis.C3 import windows_multiprocess_functions
from cahn import cahn
import numpy as np
plt.rcParams['animation.ffmpeg_path'] = 'C:\\Users\\Callicious\\Documents\\Prog\\pycharm\\venv\\ffmpeg\\bin\\ffmpeg.exe'
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
run = True
if __name__ == "__main__":

    if run == True:
        # Set up parameter ranges a, k, dx, dt, M, iniphi
        arang, \
        krang, \
        dxrang, \
        dtrang, \
        Mrang, \
        iniphirang = [0.1], [0.1], \
                     np.linspace(-2, 0.5, 40), np.linspace(-2,1, 40), \
                     [0.1], [0.5]
        dxrang = 10**dxrang
        dtrang = 10**dtrang

        # Generate all the unique possible parameters
        zipped = list(itertools.product(arang, krang, dxrang, dtrang, Mrang, iniphirang))
        print("Running for " + str(len(zipped)) + " runs...")

        # Set up the pool for zipped
        pool = multiprocessing.Pool(10)
        results = pool.map(windows_multiprocess_functions.convergence_cahn, zipped)
        pool.close()

        # Generate an array for astropy to work with
        point_array = np.array(zipped)
        sweeps, frees = np.array(results).T
        empty_array = np.empty((len(sweeps), 8))
        empty_array[:,0:6] = point_array
        empty_array[:,6] = sweeps
        empty_array[:,7] = frees

        # Table up
        table = Table(data=empty_array, names=["a","k","dx","dt","M","iniphi","sweep_end","free_end"])

        writer = hdfutils.hdf5_writer(os.getcwd(), "data.hdf5")
        writer.write_table("test_run", "test_table", table)

        run = False

    if run == False:
        writer = hdfutils.hdf5_writer(os.getcwd(), "data.hdf5")
        table = writer.read_table("test_run", "test_table")

        # Grab colour
        X,Y = table['dx'], table['dt']
        Z = table["sweep_end"]

        # Create the plot
        fig, ax1 = plt.subplots(nrows=1, ncols=1)
        ax1.tricontour(X, Y, Z, levels=15, linewidths=0.5, colors='k')
        cntr_fill = ax1.tricontourf(X, Y, Z, levels=15, cmap="RdBu_r")
        fig.colorbar(cntr_fill, ax=ax1, label="convergence sweep")
        ax1.set(xlabel=r'$\delta{x}$', ylabel=r'$\delta{t}$')
        ax1.set(xlim=[np.min(X), np.max(X)], ylim=[np.min(Y), np.max(Y)])
        #ax1.set(xlim=[0.2,0.55], ylim=[0.01, 0.055])
        plt.savefig(os.getcwd() + "dxdtcontourhalf.png", dpi=300)
        plt.show()