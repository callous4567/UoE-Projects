import numpy as np
import galcentricutils
import matplotlib.pyplot as plt
import matplotlib.ticker as pltick

# Holds all graphing utilities for project. Notes:
import windows_directories

"""
Units should be as standard.
xyz in kpc
vxvyvz in kms^-1
lxlylz in kpc kms^-1
angles in degrees
"""

# 2D Graphing Tools
class twod_graph(object):
    def __init__(self):
        self.null = "null"

    # Creates LxLyLxLzLyLz Plots. Axis units of 10^3 kpckms^-1
    def tripL(self, table, maxradius, divisor=1e3):
        table = galcentricutils.angular().get_momentum(table)
        fig, axs = plt.subplots(nrows=1, ncols=3, sharey="row", figsize=(15,7), constrained_layout=True)
        plt.subplots_adjust(wspace=0, hspace=0)

        minor_ticks_x = pltick.MultipleLocator(1)
        minor_ticks_y = pltick.MultipleLocator(1)

        texts = ["Lx - Ly", "Lx - Lz", "Ly - Lz"]
        lims = [-maxradius,maxradius]

        for num,ax in enumerate(axs):
            ax.grid(True, which='major', alpha=0, linestyle='dotted')  # Enable grids on subplot
            ax.grid(True, which='minor', alpha=0, linestyle='dotted')
            ax.set(xlim=lims,
                   ylim=lims)
            ax.set(aspect="equal")
            ax.tick_params(axis="x", which="both", direction="in", length=4, bottom=True, left=True, right=True,
                           top=True)
            ax.tick_params(axis="y", which="both", direction="in", length=4, bottom=True, left=True, right=True,
                           top=True)
            ax.xaxis.set_minor_locator(minor_ticks_x)
            ax.yaxis.set_minor_locator(minor_ticks_y)
            ax.text(s=texts[num],x=lims[0],y=lims[1]+0.1)

        axs[0].scatter(table['Lx']/divisor, table['Ly']/divisor, s=1, marker="x",c="black")
        axs[1].scatter(table['Lx']/divisor,table['Lz']/divisor,s=1,marker="x",c="black")
        axs[2].scatter(table['Ly']/divisor, table['Lz']/divisor, s=1, marker="x",c="black")

        axs[1].set(xlabel="(" + r'$10^3$' + " kpc km/s)")
        axs[0].set(ylabel="(" + r'$10^3$' + " kpc km/s)")

        plt.savefig(windows_directories.imgdir + "\\test.png", dpi=300)
        plt.show()

    # Hammer-Aitoff Projection Plots. Need angles in RADIANS for some god damn ridiculous reason.
    # [-180,180],[-90,90] phi,theta
    def haitoff(self, table):
        # First get latitude and azimuth
        table = galcentricutils.angular().get_latipolar(table)
        phis, thetas = np.deg2rad(table['phi']),np.deg2rad(table['theta'])
        for num,phi in enumerate(phis):
            if phi >= np.pi:
                phis[num] = phi - 2*np.pi
        fig, axs = plt.subplots(nrows=1, ncols=1, subplot_kw={'projection': "aitoff"})
        axs.grid(True, which='major', alpha=1, linewidth=0.25, color="black")  # Enable grids on subplot
        axs.grid(True, which='minor', alpha=1, linewidth=0.25, color="black")
        axs.scatter(phis, thetas, s=0.1, alpha=1, color="black")
        plt.show()

    # Create a histogram plot with the variable xname from the table.
    def hist(self, data, nbins, range):
        if range == False:
            plt.hist(data, bins=nbins, density=True)
            plt.show()
        else:
            plt.hist(data, bins=nbins, density=True, range=range)
            plt.show()

