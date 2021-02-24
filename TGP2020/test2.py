import multiprocessing
import os
import astrowrapper
import h5py
import numpy as np
import matplotlib.pyplot as plt
cmd = astrowrapper.cmd3(astrowrapper.rootdir, "CMD3.hdf5", "CMD3_Test", "isotable")
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


    """
    os.chdir(astrowrapper.rootdir + "\\EEEPDump")
    files = os.listdir()
    for file in files:
        try:
            np_array = np.load(file)
            B, V = np_array.T
            BV = B - V
            axs.scatter(BV, V, s=0.5)
        except:
            pass"""

    # Grab NGC7789 data
    filer = astrowrapper.hdf5_writer(astrowrapper.rootdir + "\\TEMPFOLDER", "data.hdf5")
    ngc_table = filer.read_table("NGC7789_Cross", "raw_data_reduced")
    V_abs, BV_abs = ngc_table['V_abs_dered'], ngc_table['BV_dered']

    R = 3.1
    shift = 0
    V_abs -= shift
    BV_abs -= shift/R
    axs.scatter(BV_abs, V_abs, marker="x", color="black", s=2)

    # Also grab an isochrone.
    isochro = cmd.isochrone("9.230", "-0.20")
    isoB, isoV = isochro['Bmag'], isochro['Vmag']
    isoBV = isoB - isoV
    axs.scatter(isoBV, isoV, s=1, color="red")

    # rootdir + "\\EEEPDump"

    axs.set(xlim=[-0.5, 2.5],
            ylim=[-5, 15],
            xlabel="B-V",
            ylabel="V",
            title="For NGC7789")
    fig.set_size_inches(10,10)
    plt.gca().invert_yaxis()
    plt.savefig("QuickTest.png", dpi=600)
    plt.show()
plotter()
#plotter()
print("Plotdone.")