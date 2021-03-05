import multiprocessing
import os

import tikzplotlib
from matplotlib.patches import Patch
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter

import astrowrapper
import h5py
import numpy as np
import matplotlib.pyplot as plt
utilities = astrowrapper.utils()
#utilities.stat_reduce("NGCU_Uncorrected.fits")
#utilities.stat_reduce("NGCU_Corrected.fits")

#ana = astrowrapper.source_analysis(astrowrapper.rootdir, "data.hdf5")
#ana.ana_fancyplot("M52_PREREDCATALOGUE", "phottable_table_prered", "V_apflux_annuli_manu_atmos_scimag", False, [0.1,200], "M52", 7, 1, [[1,1,0], [0,0,1]])
#ana.ana_fancyplot("NGC7789_PREREDCATALOGUE", "phottable_table_prered", "V_apflux_annuli_manu_atmos_scimag", False, [0.1,200], "NGC7789", 7, 0.7, [[1,1,0], [0,0,1]])

print(astrowrapper.utils().deg_radec([351.46,61.33]))


"""
filer = astrowrapper.hdf5_writer(astrowrapper.rootdir, "data.hdf5")
v_note = "depth_fluxvalues_V.fits"
oneflux, twoflux = filer.read("111773", v_note), filer.read("SA2039", v_note)
measurements = [1,2,3,4]
fig, axs = plt.subplots(1)
fig.set_size_inches(10,5)
axs.scatter(measurements, twoflux, color="blue", marker="x", label="SA 20-39")
axs.scatter(measurements, oneflux, color="red", marker="x", label="111-773")
axs.set(xlim=[0.1, 4.9],
        xlabel="Measurement",
        ylabel=r'$I_\lambda$')

axs.grid(True, which='major', color="pink", alpha=1, linestyle='dotted', lw=0.5)  # Enable grids on subplot

#legend_elements = [Patch(facecolor='red', edgecolor='black', label="SA 20-39"),
#                   Patch(facecolor='blue', edgecolor='black', label="111-773")]
axs.legend(loc='upper right')
plt.savefig(astrowrapper.texdir + "\\" + "magcomparison.png", dpi=300)
tikzplotlib.save(astrowrapper.texdir + "\\" + "magcomparison.tex")
plt.show()"""


"""
def NGC_doer():
        # Load in the information from the filer
        filer = astrowrapper.hdf5_writer(astrowrapper.rootdir, "data.hdf5")
        array = filer.read_table("NGC7789_Cross", "-0.30")
        ages, chisq = array['ages'], array['chisqs']
        chisq = [3*d for d in chisq]
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)

        # Quick plot
        fig, axs = plt.subplots(1)
        argmin = np.argmin(chisq)
        age_min = ages[argmin]
        chisqmin = chisq[argmin]

        # Set up the chisq boundaries
        minchisq = chisqmin
        onesig, twosig, threesig = chisqmin+1, chisqmin+4, chisqmin+9
        onesigcoords, twosigcoords, threesigcoords = [],[],[]
        for num, d in enumerate(chisq):
            if d <= threesig:
                threesigcoords.append(num)
                break
        for num, d in enumerate(chisq):
            if d <= twosig:
                twosigcoords.append(num)
                break
        for num, d in enumerate(chisq):
            if d <= onesig:
                onesigcoords.append(num)
                break
        for num, d in enumerate(chisq):
            one_num = num + onesigcoords[0]
            if chisq[one_num] >= onesig:
                onesigcoords.append(one_num)
                break
        for num, d in enumerate(chisq):
            two_num = num + twosigcoords[0]
            if chisq[two_num] >= twosig:
                twosigcoords.append(two_num)
                break
        for num, d in enumerate(chisq):
            three_num = num + threesigcoords[0]
            if chisq[three_num] >= threesig:
                threesigcoords.append(three_num)
                break

        # Grab the age guess argument for the sigma bounds
        onesigages, twosigages, threesigages = [ages[onesigcoords[0]], ages[onesigcoords[-1]]], [ages[twosigcoords[0]], ages[twosigcoords[-1]]], [ages[threesigcoords[0]], ages[threesigcoords[-1]]]


        axs.axvspan(threesigages[0], threesigages[1], facecolor='red', alpha=1)
        axs.axvspan(twosigages[0], twosigages[1], facecolor='yellow', alpha=1)
        axs.axvspan(onesigages[0], onesigages[1], facecolor='green', alpha=1)
        filer = astrowrapper.hdf5_writer(astrowrapper.rootdir, "data.hdf5")
        sigages = [onesigages, twosigages, threesigages]
        filer.write("NGC7789_Cross", "sigage", sigages)

        axs.plot(ages, chisq, lw=1, color='black')
        axs.set(xlabel="logt / dex",
                ylabel=r'$\chi^2$',
                title=r'$\chi^2$' + " minimum for log(t) = " + ("{0:.3f}").format(age_min))
        #axs.axhspan(chisqmin + 1, chisqmin + 2, facecolor='blue', alpha=0.5)
        #axs.axhspan(chisqmin + 2, chisqmin + 3, facecolor='red', alpha=0.5)
        axs.grid(True, which='major', color="blue", alpha=1, linestyle='dotted', lw=0.5)  # Enable grids on subplot
        axs.grid(True, which='minor', color="pink", alpha=1, linestyle='dotted', lw=0.5)
        legend_elements = [Patch(facecolor='red', edgecolor='black', label=r'$3\sigma$'),
                           Patch(facecolor='yellow', edgecolor='black', label=r'$2\sigma$'),
                           Patch(facecolor='green', edgecolor='black', label=r'$1\sigma$')]
        axs.legend(handles=legend_elements, loc='lower left')
        plt.savefig(astrowrapper.texdir + "\\" + "NGCage.png", dpi=300)
        plt.show()
#NGC_doer()
def M52_doer():
        # Load in the information from the filer
        filer = astrowrapper.hdf5_writer(astrowrapper.rootdir, "data.hdf5")
        array = filer.read_table("M52_Cross", "0.05")
        ages, chisq = array['ages'], array['chisqs']
        ages = ages[50:]
        chisq = chisq[50:]
        chisq = [3*d for d in chisq]
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)
        chisq = savgol_filter(x=chisq, window_length=5, polyorder=0, mode='nearest', deriv=0, delta=0)


        # Quick plot
        fig, axs = plt.subplots(1)
        argmin = np.argmin(chisq)
        age_min = ages[argmin]
        chisqmin = chisq[argmin]

        # Set up the chisq boundaries
        onesig, twosig, threesig = chisqmin+1, chisqmin+4, chisqmin+9
        onesigcoords, twosigcoords, threesigcoords = [],[],[]
        for num, d in enumerate(chisq):
            if d <= threesig:
                threesigcoords.append(num)
                break
        for num, d in enumerate(chisq):
            if d <= twosig:
                twosigcoords.append(num)
                break
        for num, d in enumerate(chisq):
            if d <= onesig:
                onesigcoords.append(num)
                break
        for num, d in enumerate(chisq):
            one_num = num + onesigcoords[0]
            if chisq[one_num] >= onesig:
                onesigcoords.append(one_num)
                break
        for num, d in enumerate(chisq):
            two_num = num + twosigcoords[0]
            if chisq[two_num] >= twosig:
                twosigcoords.append(two_num)
                break
        for num, d in enumerate(chisq):
            three_num = num + threesigcoords[0]
            if chisq[three_num] >= threesig:
                threesigcoords.append(three_num)
                break

        # Grab the age guess argument for the sigma bounds
        onesigages, twosigages, threesigages = [ages[onesigcoords[0]], ages[onesigcoords[-1]]], [ages[twosigcoords[0]], ages[twosigcoords[-1]]], [ages[threesigcoords[0]], ages[threesigcoords[-1]]]
        axs.axvspan(threesigages[0], threesigages[1], facecolor='red', alpha=1)
        axs.axvspan(twosigages[0], twosigages[1], facecolor='yellow', alpha=1)
        axs.axvspan(onesigages[0], onesigages[1], facecolor='green', alpha=1)
        filer = astrowrapper.hdf5_writer(astrowrapper.rootdir, "data.hdf5")
        sigages = [onesigages, twosigages, threesigages]
        filer.write("M52_Cross", "sigage", sigages)
        axs.plot(ages, chisq, lw=1, color='black')
        axs.set(xlabel="logt / dex",
                ylabel=r'$\chi^2$',
                title=r'$\chi^2$' + " minimum for log(t) = " + ("{0:.3f}").format(age_min))

        legend_elements = [Patch(facecolor='red', edgecolor='black', label=r'$3\sigma$'),
                           Patch(facecolor='yellow', edgecolor='black', label=r'$2\sigma$'),
                           Patch(facecolor='green', edgecolor='black', label=r'$1\sigma$')]
        axs.legend(handles=legend_elements, loc='lower left')

        #axs.axhspan(chisqmin + 1, chisqmin + 2, facecolor='blue', alpha=0.5)
        #axs.axhspan(chisqmin + 2, chisqmin + 3, facecolor='red', alpha=0.5)
        axs.grid(True, which='major', color="blue", alpha=1, linestyle='dotted', lw=0.5)  # Enable grids on subplot
        axs.grid(True, which='minor', color="pink", alpha=1, linestyle='dotted', lw=0.5)
        plt.savefig(astrowrapper.texdir + "\\" + "M52.png", dpi=300)
        plt.show()
#M52_doer()

#quickplot2()"""

# M52_PREREDCATALOGUE
# phottable_table_prered
#ana = astrowrapper.source_analysis(astrowrapper.rootdir, "data.hdf5")
#ana.ana_fancyplot("M52_PREREDCATALOGUE", "phottable_table_prered", "V_apflux_annuli_manu_atmos_scimag", 8, [5,200], "M52 Field", 14, 0.5, [[1,1,0], [0,0,1]], "M52_ALL")
#ana.ana_fancyplot("NGC7789_PREREDCATALOGUE", "phottable_table_prered", "V_apflux_annuli_manu_atmos_scimag", 9, [5,200], "NGC 7789 Field", 12, 0.5, [[1,1,0], [0,0,1]], "NGC7789_ALL")
#ana.ana_fancyplot("M52_Cross", "raw_data_second_reduced", "V_apflux_annuli_manu_atmos_scimag", 14, [5,200], "[FSR2007] 0433 Membership", 7, 0.5, [[1,1,0], [0,0,1]], "M52_SECOND")
#ana.ana_fancyplot("M52_Cross", "raw_data_reduced", "V_apflux_annuli_manu_atmos_scimag", 10, [5,200], "M52 Membership", 9, 1, [[1,1,0], [0,0,1]], "M52_MEMBER")
#ana.ana_fancyplot("NGC7789_Cross", "raw_data_reduced", "V_apflux_annuli_manu_atmos_scimag", 10, [5,200], "NGC 7789 Membership", 9, 1, [[1,1,0], [0,0,1]], "NGC_MEMBER")


"""
cmd = astrowrapper.cmd3(astrowrapper.rootdir, "CMD3.hdf5", "CMD3_Test", "isotable")
met = "0.00"
age1, age2, age3 = "7.000", "9.000", "9.500"
iso1, iso2, iso3 = cmd.isochrone(age1, met), cmd.isochrone(age2, met), cmd.isochrone(age3, met)
print(iso1, iso2)
oneb, onev, twob, twov = iso1['Bmag'], iso1['Vmag'], iso2['Bmag'], iso2['Vmag']
threeb, threev = iso3['Bmag'], iso3['Vmag']
onebv, twobv = oneb - onev, twob - twov
threebv = threeb - threev
fig, ax = plt.subplots(1)
ax.scatter(onebv, onev, s=1, label=("log(t = " + str(age1) + ")"))
ax.scatter(twobv, twov, s=1, label=("log(t = " + str(age2) + ")"))
ax.scatter(threebv, threev, s=1, label=("log(t = " + str(age3) + ")"))
ax.set(xlim=[-1, 3],
       ylim=[-10,15],
       xlabel="B - V",
       ylabel=r'$M_v$')
fig.set_size_inches(5,5)
plt.legend(loc='lower right')
plt.gca().invert_yaxis()
plt.savefig(astrowrapper.texdir + "\\isotest.png", dpi=300)
plt.show()"""


"""

#cmd = astrowrapper.cmd3(astrowrapper.rootdir, "CMD3.hdf5", "CMD3_Test", "isotable")
def split(a, n):
    k, m = divmod(len(a), n)
    return (a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n))

cmd.met_age_group("-0.20")
float_ages = cmd.float_ages




def generator():
    mass_range = np.linspace(0.1, 20, 200)
    metallicity = "0.00"
    Ncores = 2
    cmd.met_age_group(metallicity)
    for mass in mass_range:
        cmd.eep_mass_met(mass)

#generator()
print("Done.")

def plotter():
    fig, axs = plt.subplots(1)


    os.chdir(astrowrapper.rootdir + "\\EEEPDump")
    files = os.listdir()
    for file in files:
        try:
            np_array = np.load(file)
            B, V = np_array.T
            BV = B - V
            axs.scatter(BV, V, s=0.5)
        except:
            pass

    # Grab NGC7789 data
    #filer = astrowrapper.hdf5_writer(astrowrapper.rootdir + "\\TEMPFOLDER", "data.hdf5")
    #ngc_table = filer.read_table("NGC7789_Cross", "raw_data_reduced")
    #V_abs, BV_abs = ngc_table['V_abs_dered'], ngc_table['BV_dered']

    #R = 3.1
    #shift = 0
    #V_abs -= shift
    #BV_abs -= shift/R
    #axs.scatter(BV_abs, V_abs, marker="x", color="black", s=2)

    # Also grab an isochrone.
    #isochro = cmd.isochrone("9.230", "-0.20")
    #isoB, isoV = isochro['Bmag'], isochro['Vmag']
    #isoBV = isoB - isoV
    #axs.scatter(isoBV, isoV, s=1, color="red")

    # rootdir + "\\EEEPDump"

    axs.set(xlim=[-0.3, 1.5],
            ylim=[-5, 12.5],
            xlabel="B-V",
            ylabel="V",
            title="EEP Tracks for stars of various masses, [Fe/H] = -0.20")
    axs.grid(True, which='major', color="pink", alpha=1, linestyle='dotted', lw=0.5)  # Enable grids on subplot
    axs.grid(True, which='minor', color="pink", alpha=1, linestyle='dotted', lw=0.5)
    fig.set_size_inches(10,10)
    plt.gca().invert_yaxis()
    plt.savefig("EEPTracks.png", dpi=600)
    plt.show()
#plotter()
#plotter()
print("Plotdone.")

# Planck Plotter
def planck_plotter():
    h, c, kb = 6.626e-34, 3e8, 1.38e-23
    planck_func = lambda lam, T: (2*h*(c**2)/(lam**5))*((np.exp((h*c)/(lam*kb*T)) - 1)**-1)
    lambdas = np.linspace(0, 1200e-9, 1000)
    lambdas_norm = lambdas/(1e-9)
    temp = 12000
    temp = planck_func(lambdas, temp)
    fig, axs = plt.subplots(1)
    axs.set(xlabel=r'$\lambda$' + " / nm",
            ylabel=r'$B / Jm^{-3}s^{-3}$',
            title="Planck Law, example T = 16,000K")
    axs.plot(lambdas_norm, temp)
    axs.grid(True, which='major', color="pink", alpha=1, linestyle='dotted', lw=0.5)  # Enable grids on subplot
    axs.grid(True, which='minor', color="pink", alpha=1, linestyle='dotted', lw=0.5)
    axs.axvspan(530, 570, facecolor='green', alpha=0.5)
    axs.axvspan(425, 465, facecolor='blue', alpha=0.5)
    fig.set_size_inches(10, 5)
    plt.savefig("PlanckPlot.png", dpi=600)
    plt.show()

planck_plotter() """


def M522_doer():
    # Load in the information from the filer
    filer = astrowrapper.hdf5_writer(astrowrapper.rootdir, "data.hdf5")
    array = filer.read_table("M52_Cross", "0.00")
    ages, chisq = array['ages'], array['chisqs']
    chisq = chisq/5


    # Quick plot
    fig, axs = plt.subplots(1)
    argmin = np.argmin(chisq)
    age_min = ages[argmin]
    chisqmin = chisq[argmin]

    # Set up the chisq boundaries
    onesig, twosig, threesig = chisqmin + 1, chisqmin + 4, chisqmin + 9
    onesigcoords, twosigcoords, threesigcoords = [], [], []
    for num, d in enumerate(chisq):
        if d <= threesig:
            threesigcoords.append(num)
            break
    for num, d in enumerate(chisq):
        if d <= twosig:
            twosigcoords.append(num)
            break
    for num, d in enumerate(chisq):
        if d <= onesig:
            onesigcoords.append(num)
            break
    for num, d in enumerate(chisq):
        one_num = num + onesigcoords[0]
        if chisq[one_num] >= onesig:
            onesigcoords.append(one_num)
            break
    for num, d in enumerate(chisq):
        two_num = num + twosigcoords[0]
        if chisq[two_num] >= twosig:
            twosigcoords.append(two_num)
            break
    for num, d in enumerate(chisq):
        three_num = num + threesigcoords[0]
        if chisq[three_num] >= threesig:
            threesigcoords.append(three_num)
            break

    # Grab the age guess argument for the sigma bounds
    onesigages, twosigages, threesigages = [ages[onesigcoords[0]], ages[onesigcoords[-1]]], [ages[twosigcoords[0]],
                                                                                             ages[twosigcoords[-1]]], [
                                               ages[threesigcoords[0]], ages[threesigcoords[-1]]]
    axs.axvspan(threesigages[0], threesigages[1], facecolor='red', alpha=1)
    axs.axvspan(twosigages[0], twosigages[1], facecolor='yellow', alpha=1)
    axs.axvspan(onesigages[0], onesigages[1], facecolor='green', alpha=1)
    filer = astrowrapper.hdf5_writer(astrowrapper.rootdir, "data.hdf5")
    sigages = [onesigages, twosigages, threesigages]
    filer.write("M52_Cross", "sigage", sigages)
    axs.plot(ages, chisq, lw=1, color='black')
    axs.set(xlabel="logt / dex",
            ylabel=r'$\chi^2$',
            title=r'$\chi^2$' + " minimum for log(t) = " + ("{0:.3f}").format(age_min))

    legend_elements = [Patch(facecolor='red', edgecolor='black', label=r'$3\sigma$'),
                       Patch(facecolor='yellow', edgecolor='black', label=r'$2\sigma$'),
                       Patch(facecolor='green', edgecolor='black', label=r'$1\sigma$')]
    axs.legend(handles=legend_elements, loc='lower left')

    # axs.axhspan(chisqmin + 1, chisqmin + 2, facecolor='blue', alpha=0.5)
    # axs.axhspan(chisqmin + 2, chisqmin + 3, facecolor='red', alpha=0.5)
    axs.grid(True, which='major', color="blue", alpha=1, linestyle='dotted', lw=0.5)  # Enable grids on subplot
    axs.grid(True, which='minor', color="pink", alpha=1, linestyle='dotted', lw=0.5)
    plt.savefig(astrowrapper.texdir + "\\" + "M52.png", dpi=300)
    plt.show()
M522_doer()