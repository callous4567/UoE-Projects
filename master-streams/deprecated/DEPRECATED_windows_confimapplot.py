import math
import os

import numpy as np
import time
import ascii_info
import galcentricutils
import graphutils
import hdfutils
import windows_directories
import matplotlib.pyplot as plt
import matplotlib as mpl
import astropy.io.misc.hdf5
import windows_stack
import windows_confimap

"""
This one actually plots the confidence maps. 
See img for examples.
"""


# For each dtheta/n_phi set up a confidence map and save it
for dtheta in windows_testing.dthetarange:
    for n_phi in windows_testing.nphiarange:
        try:
            print(dtheta, n_phi)
            # For saving within the original group (set within for the saved GCC_TEST_CONFIDENCE)
            string_format = ("{0:.2f}_{1:.0f}").format(dtheta, n_phi)
            writer = hdfutils.hdf5_writer(windows_directories.datadir, ascii_info.asciiname)
            # table_data = writer.read_table("full_raw","astrotable")
            # polar_data = galcentricutils.angular().get_polar(table_data)
            table = writer.read_table(ascii_info.fullgroup, string_format)
            thetas, phis, indices, sigmins = table['theta'], table['phi'], table['index'], table['sigmin']
            #sigmins = sigmins*len(thetas)/2 # bonferroni
            # n_of_phi_region = table['phiregion_n']
            logoversigmins = np.log10(1 / sigmins)
            # logoversigmins = [np.inf if d <= 1 else d for d in logoversigmins]
            # Blank out the infinities with the max of the logoversigmins (without infs)
            # maxi = np.nanmax(logoversigmins[logoversigmins != np.inf])
            # logoversigmins = [maxi if d == np.inf else d for d in logoversigmins]
            # print(logoversigmins, max(logoversigmins))

            # 90 - theta ~ latipolar theta, phi is equivalent for polar/latipolar.
            thetas = 90 - thetas

            fig, axs = plt.subplots(nrows=1, ncols=1, dpi=300, figsize=(8, 4))
            # g = axs.scatter(phis, thetas, marker='x', s=0.1, alpha=0.1,  color="black")
            # cmap = plt.cm.nipy_spectral(np.linspace(0.1,1,1000))
            # cmap = mpl.colors.ListedColormap(cmap[10:,:-1])
            h = axs.scatter(phis, thetas, c=logoversigmins, edgecolors="none", marker='s', s=16, cmap=plt.cm.nipy_spectral)
            fig.colorbar(h, label=r'${\log_{10}}{\left(\frac{1}{p_{gc}}\right)}$')
            axs.grid(False)
            axs.set(xlabel=r'$\phi$',
                    ylabel=r'$\theta$',
                    xlim=[0, 360],
                    ylim=[-60, 60])
            #plt.title("With Bonferroni")
            try:
                os.mkdir(windows_directories.imgdir + "\\" + "confidence_map")
            except:
                pass

            plt.savefig(windows_directories.imgdir + "\\" + "confidence_map" + "\\" + string_format + ".png")
            plt.clf()
            plt.close()

            # Quick number density plot
            # graphutils.twod_graph().thetaphi_twodhist(polar_data, 200)
        except:
            pass

