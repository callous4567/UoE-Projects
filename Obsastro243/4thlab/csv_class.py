import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as pltick
from matplotlib.patches import Patch
import matplotlib as mpl
import time
import pandas as pd
import io
import sys
from OOP1 import utilities as utils
# Proudly callicioused by coded.

# Instantiate and use for csv operations, specifically for ThorLabs CSV output files.
class csv_handler(object):
    def __init__(self, dummy): # Dummy object to allow instantation. Idk another way >.<
        self.dummy = dummy

    # Reads in CSV file and returns [x, y] array and the filetension for future operations.
    def csv_reader(self, filetension):
        with open(filetension, 'rb') as f:
            sidebar_delimited = (line.replace(b'[', b'#') for line in f) # b converts to bytes object
            xy = np.genfromtxt(sidebar_delimited, dtype=float, comments='#', delimiter=';') # nm wavelength x vs. intensity y.
            xy = np.array(xy).T # columnstack brings you back to vectors, transposition brings you back to lists.
        return xy, filetension

    # Writes a new CSV file where all elements of "band_filter" are divided by elements of "clear_filter," a transmission curve relative to "clear_filter">
    # Also returns [x,y] array for this curve and the new csv file nametension, for future operations.
    def transmission_curve(self, clear_filter, band_filter):
        clear, transmit = self.csv_reader(clear_filter)[0], self.csv_reader(band_filter)[0]
        transmission_curve = np.array([transmit[1][i]/clear[1][i] for i in range(len(clear[0]))])
        transmission_curve = [float(d) for d in transmission_curve]
        transmission_data = np.asarray([clear[0], transmission_curve]).T
        newcsv = clear_filter + band_filter + ".csv"
        np.savetxt(newcsv, transmission_data, delimiter=";")

        #plt.plot(clear[0], transmission_curve, lw=0.2)
        #plt.ylim([0,1])
        #plt.show()

        return [clear[0], transmission_curve], newcsv

    # Multiplies a transmission curve by another curve. Returns [x,y] (for grapher) and the name of the csv file.
    def multi_curve(self, first_curve, second_curve):
        first, second = self.csv_reader(first_curve)[0], self.csv_reader(second_curve)[0]
        multiplied_curve = np.array([first[1][i]*second[1][i] for i in range(len(first[0]))])
        multiplied_curve = [float(d) for d in multiplied_curve]
        multiplied_data = np.asarray([first[0], multiplied_curve]).T
        newcsv = first_curve + second_curve + ".csv"
        np.savetxt(newcsv, multiplied_data, delimiter=";")

        #plt.plot(first[0], multiplied_curve, lw=0.2)
        #plt.ylim([0,1])
        #plt.show()

        return [first[0], multiplied_curve], newcsv

    # Takes various parameters and makes a graph.
    # is_transmission is "true" or "false" : if true, the y label is set to R(lambda). if false, set to Intensity.
    # xlims and  ylims are limits for wavelength and intensity/transmissivity in array form, i.e. [0,1000] and [0,1]
    # title is the title of the  graph. If you don't want one, just do " "
    # label is the legend label. type "false" if you don't want a legend.
    # colour gives graph colour. 'r' is red, etc
    # linewidth should be 0.2 for pretty graphs.
    # data is the XY data from the CSV file reader/etc.
    # savefigtension is the output and extension. use .png :_:
    def grapher(self, is_transmission, xlims, ylims, title, label, colour, linewidth, data, savefigtension):

        xlabel = '$\lambda$ / nm'
        if is_transmission == "true":   # Customises for transmission curve, I imagine that's desired.
            ylabel = 'R($\lambda$)'
        else:
            ylabel = 'Intensity D, Unknown Units'

        fig = plt.figure()
        ax = fig.add_subplot(111)

        ax.set(xlim=xlims, ylim=ylims, title=title, xlabel=xlabel, ylabel=ylabel)
        ax.plot(data[0], data[1], label=label, lw=linewidth, color=colour)


        if label != "false":
            legend_elements = [Patch(facecolor=colour, edgecolor='black', label=label)]
            ax.legend(handles=legend_elements, loc='upper right', title="Legend")

        mpl.rc('font', size=8)  # Redefine matplotlib rcParams to change font size

        major_ticks_x, major_ticks_y, minor_ticks_x, minor_ticks_y = 100, 0.1, 50, 0.05 # Set your majors and minors. It might be worth modifying this if you intend to clip graphs along the y axis.

        major_ticks_x = pltick.MultipleLocator(major_ticks_x)  # Decide/calculate locators for the major and minor ticks on subplot
        major_ticks_y = pltick.MultipleLocator(major_ticks_y)
        minor_ticks_x = pltick.MultipleLocator(minor_ticks_x)
        minor_ticks_y = pltick.MultipleLocator(minor_ticks_y)

        ax.xaxis.set_major_locator(major_ticks_x)  # Apply locators to ax subplot
        ax.yaxis.set_major_locator(major_ticks_y)
        ax.xaxis.set_minor_locator(minor_ticks_x)
        ax.yaxis.set_minor_locator(minor_ticks_y)

        ax.grid(True, which='major', color="blue", alpha=1, linestyle='dotted', lw=0.5)  # Enable grids on subplot
        ax.grid(True, which='minor', color="pink", alpha=1, linestyle='dotted', lw=0.5)

        plt.plot()
        plt.show()
        fig.savefig(savefigtension, dpi=600)


# GUIDELINES FOR USE

# Instantiate.
csv_handle = csv_handler("dummy")

# Perform operation to get data. EXAMPLE USES OF THE PROGRAM ARE DISPLAYED BELOW! Just copy-paste it and input the CSV files that you want to use.

# REGULAR SINGLE PLOT, NO LABEL
unfiltered = csv_handle.csv_reader("unfiltered bulb.csv")
csv_handle.grapher("false", [0, 1000], [0,1], "Unfiltered Sodium Bulb Spectrum", "false", "r", 0.2, unfiltered[0], "unfiltered.png")

# TRANSMISSION CURVE PLOT WITH A LABEL
tcurve_combo14 = csv_handle.transmission_curve("unfiltered bulb.csv", "combo14.csv")
csv_handle.grapher("true", [0,1000], [0,1], "Transmission Curve of filter wheel combination '1,4'", "FW Combo '1,4'", "b", 0.2, tcurve_combo14[0], "tcurvecombo14.png")

# MULTIPLIED CURVE PLOT WITH A LABEL! In this case, I've first multiplied FW 1 on the default spectrum, then FW4, giving what should be 1,4 in the process! We'll see.
# Anyway it isn't. Whatever. The multiplication method is sound and it works, you get what you pay for.
tcurve_1 = csv_handle.transmission_curve("unfiltered bulb.csv", "filter1.csv")
tcurve_4 = csv_handle.transmission_curve("unfiltered bulb.csv", "filter4fw1.csv")
applied1 = csv_handle.multi_curve("unfiltered bulb.csv", tcurve_1[1])
applied14 = csv_handle.multi_curve(applied1[1], tcurve_4[1])
tcurvepred14 = csv_handle.transmission_curve("unfiltered bulb.csv", applied14[1])
csv_handle.grapher("true", [0,1000], [0,1], "Transmission Curve (predicted) of filter wheel combination '1,4'", "FW Pred Combo '1,4'", "b", 0.2, tcurvepred14[0], "tcurvepred14.png")
