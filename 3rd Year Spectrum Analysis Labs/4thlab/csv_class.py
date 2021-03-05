import numpy as np

import matplotlib.pyplot as plt
import matplotlib.ticker as pltick
from matplotlib.patches import Patch
import matplotlib as mpl
from adjustText import adjust_text

from scipy.optimize import curve_fit
from scipy import constants
from lmfit import minimize, Parameter, Parameters
import lmfit as lmf
from scipy.signal import savgol_filter
from scipy.signal import find_peaks

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
    # locatrix has major_x, major_y, minor_x, minor_y locators.
    def grapher(self, is_transmission, xlims, ylims, title, label, colour, linewidth, data, savefigtension, locatrix):

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

        major_ticks_x, major_ticks_y, minor_ticks_x, minor_ticks_y = locatrix[0], locatrix[1], locatrix[2], locatrix[3] # Set your majors and minors. It might be worth modifying this if you intend to clip graphs along the y axis.

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
        fig.savefig(savefigtension, dpi=1200)
    # Grapher except this one will put peaks on using scipy.peaks
    def peakgrapher(self, is_transmission, xlims, ylims, title, label, colour, linewidth, data, savefigtension, locatrix, scalefactor, savfactor):

        xlabel = '$\lambda$ / nm'
        if is_transmission == "true":   # Customises for transmission curve, I imagine that's desired.
            ylabel = 'R($\lambda$)'
        else:
            ylabel = 'Intensity D, Unknown Units'

        fig = plt.figure()
        ax = fig.add_subplot(111)

        x,y = data[0]*0.1, data[1]
        y_flattened = savgol_filter(y, window_length=savfactor, polyorder=0, deriv=0, delta=0, mode='nearest') # Savgol flatten for noise/etc


        # SCIPYFINDPEAKSPARAMETERS
        peakdata = find_peaks(y_flattened, height=None, distance=20, prominence=max(y_flattened)/scalefactor)
        peakxvalues = [x[d] for d in peakdata[0]]


        ax.set(xlim=xlims, ylim=ylims, title=title, xlabel=xlabel, ylabel=ylabel)

        textdata = [plt.text(x[i], y_flattened[i], x[i], ha='center', va='center') for i in peakdata[0]]
        adjust_text(textdata, arrowprops=dict(arrowstyle='->', color='green'))

        ax.plot(x, y_flattened, label=label, lw=linewidth, color=colour)
        #ax.plot(x, y, label=label, lw=linewidth, color='blue')

        if label != "false":
            legend_elements = [Patch(facecolor=colour, edgecolor='black', label=label)]
            ax.legend(handles=legend_elements, loc='upper right', title="Legend")

        mpl.rc('font', size=8)  # Redefine matplotlib rcParams to change font size

        major_ticks_x, major_ticks_y, minor_ticks_x, minor_ticks_y = locatrix[0], locatrix[1], locatrix[2], locatrix[3] # Set your majors and minors. It might be worth modifying this if you intend to clip graphs along the y axis.

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
        fig.savefig(savefigtension, dpi=1200)

    # Planck Minimise Function. Params is 1x4 array of ckth.
    # Useful if you try scipy optimization, which didn't work well enough with brute force LSQ.
    # Since we don't know the "normalization" of the planck curve here, we'll toss in an extra constant.
    def plancker(self, params, wldata, yvaldata):
        h, c, k, T, norm = constants.h, constants.c, constants.Boltzmann, params['T'].value, params['norm'].value
        resivalues = []
        for d in range(len(wldata)):
            planckercalc = ((2*norm*h*c**2)/((T**4)*wldata[d]**5))/(np.exp((h*c)/(wldata[d]*k*T)) - 1)
            resivalues.append(planckercalc)

        # Various residual fittings
        argsmax = [np.argmax(resivalues), np.argmax(yvaldata)]
        lambdasmax = [wldata[argsmax[0]], wldata[argsmax[1]]]

        return (resivalues - yvaldata)**4 + (1 - max(resivalues)/max(yvaldata))**2 + (lambdasmax[0] - lambdasmax[1])**2

    def plancker_alternate(self, params, wldata, yvaldata):
        h, c, k, T, norm, shift = constants.h, constants.c, constants.Boltzmann, params['T'].value, params['norm'].value, params['shift'].value

        delta = wldata[np.argmax(yvaldata)]

        resivalues = []
        for d in range(len(wldata)):
            f_value = (wldata[d] - shift)
            planckercalc = ((2*norm*h*c**2)/((T**4)*(f_value)**5))/(np.exp((h*c)/((f_value)*k*T)) - 1)
            resivalues.append(planckercalc)

        # Various residual fittings
        argsmax = [np.argmax(resivalues), np.argmax(yvaldata)]
        lambdasmax = [wldata[argsmax[0]], wldata[argsmax[1]]]

        return (resivalues - yvaldata)**2 + (1 - max(resivalues)/max(yvaldata))**2


    # LSQ solve (lmfit)
    def least_solver(self, xydata, label):
        x, y = [float(d) for d in xydata[0]], [float(d) for d in xydata[1]]
        # Correct for nanometers.
        x = [d * 10**-9 for d in x]

        x, y = x[0:3600], y[0:3600] # Clip to first 3600 elements to remove nastyness afterward.
        # Arrayify

        # Smooth and get an estimate for the blackbody temperature using Weins law.
        yn = savgol_filter(y, window_length=101, polyorder=0, deriv=0, delta=0, mode='nearest')
        yn = savgol_filter(yn, window_length=101, polyorder=0, deriv=0, delta=0, mode='nearest')
        yn = savgol_filter(yn, window_length=101, polyorder=0, deriv=0, delta=0, mode='nearest')
        yn = savgol_filter(yn, window_length=101, polyorder=0, deriv=0, delta=0, mode='nearest')
        y = yn

        # Wein temp estimate based on peak wavelength.
        lambda_wein = x[np.argmax(y)]
        bestimate = constants.Wien/lambda_wein
        print(bestimate)

        # Normalisation estimate based on matching the peaks
        peak_flux = max(y)
        peak_calculated_flux = self.planck(x[np.argmax(y)], *[bestimate, 1])
        normstimate = 1/ (peak_calculated_flux/peak_flux)

        # Shift estimate
        shiftstimate = x[np.argmax(y)] - lambda_wein


        x, y = np.array(x), np.array(y)

        plt.plot(x,y, lw=0.1)
        plt.savefig("savefast.png", dpi=1200)

        params = Parameters()
        params.add('T', bestimate, min=10, max=8000)
        params.add('norm', value=normstimate, min=0.0001, max=10000, brute_step=1)
        params.add('shift', value=shiftstimate, min=0, max=100, brute_step=1*10**-9)

        test1 = minimize(self.plancker, params, args=(x, y), method='least_squares')
        params = test1.params
        print(params)

        h, c, k, T, norm, shift = constants.h, constants.c, constants.Boltzmann, params['T'].value, params['norm'].value, params['shift'].value

        # Plancker and Plancker Alternate.
        owo = lambda wl: ((2*norm*h* (c ** 2)) / ((T**4)* ((wl) ** 5))) / (np.exp(h * c / (wl * k * T)) - 1)
        uwu = lambda wl: ((2*norm*h* (c ** 2)) / ((T ** 4) * (((wl - shift)) ** 5))) / (np.exp(h * c / (((wl - shift)) * k * T)) - 1)

        yn = owo(x)

        # Scale to match y
        factor = (max(y)/max(yn))
        yn = [factor*d for d in yn]

        ##############################
        fig = plt.figure()
        ax = fig.add_subplot(111)

        x = x * 1e9
        print(x)

        ax.set(title=(label + " Savgol Filtered vs. Planck Fitted Spectra"), xlabel="$\lambda$ / nm", ylabel="Intensity D, Unknown Units")
        ax.plot(x, y, label="Raw Data", lw=0.5, color='red')
        ax.plot(x, yn, label="Planck Fitted", lw=0.5, color='green')

        if label != "false":
            legend_elements = [Patch(facecolor='red', edgecolor='black', label="Raw Data"),
                               Patch(facecolor='green', edgecolor='black', label="Planck Fitted")]
            ax.legend(handles=legend_elements, loc='upper right', title="Legend")

        mpl.rc('font', size=8)  # Redefine matplotlib rcParams to change font size


        major_ticks_x, major_ticks_y, minor_ticks_x, minor_ticks_y = 100, 0.2, 10, 0.1

        major_ticks_x = pltick.MultipleLocator(
            major_ticks_x)  # Decide/calculate locators for the major and minor ticks on subplot
        major_ticks_y = pltick.MultipleLocator(major_ticks_y)
        minor_ticks_x = pltick.MultipleLocator(minor_ticks_x)
        minor_ticks_y = pltick.MultipleLocator(minor_ticks_y)

        ax.xaxis.set_major_locator(major_ticks_x)  # Apply locators to ax subplot
        ax.yaxis.set_major_locator(major_ticks_y)
        ax.xaxis.set_minor_locator(minor_ticks_x)
        ax.yaxis.set_minor_locator(minor_ticks_y)

        ax.grid(True, which='major', color="blue", alpha=1, linestyle='dotted', lw=0.5)  # Enable grids on subplot
        ax.grid(True, which='minor', color="pink", alpha=1, linestyle='dotted', lw=0.5)

        plt.text(x[np.argmax(yn)], max(yn),"Fitted T of " + str("{0:.3e}").format(params['T'].value) + " K", horizontalalignment='right')



        plt.plot()
        plt.show()
        fig.savefig(label + ".png", dpi=1200)

    # Generic planck function (not minimiser) w/ params having T and Norm
    def planck(self, wavelength, *params):
        h, c, k = constants.h, constants.c, constants.Boltzmann
        return ((2*h*(c**2))/(params[1]*(wavelength**5)))/(np.exp(h*c/(wavelength*k*params[0]))-1)

    # Solves using scipy curve_fit.
    def scipy_solver(self, xydata):
        x, y = [float(d) for d in xydata[0]], [float(d) for d in xydata[1]]
        # Correct for nanometers.
        x = [d * 10 ** -9 for d in x]

        # y = savgol_filter(y, window_length=101, polyorder=8, deriv=0, delta=1, mode='interp')
        x, y = x[0:3600], y[0:3600]  # Clip to first 3600 elements to remove nastyness afterward.
        # Arrayify
        x, y = np.array(x), np.array(y)

        params = [4000, 1] # initial guess for T and norm constant
        popt, cov = curve_fit(self.planck, xdata=x, ydata=y, p0=params)
        print(popt[0])

        yn = self.planck(x, *popt)
        normalisefactor = max(savgol_filter(y, window_length=101, polyorder=0))/max(yn)
        yn = [normalisefactor*d for d in yn]

        plt.plot(x, y)
        plt.plot(x, yn)

        plt.show()


# GUIDELINES FOR USE

# Instantiate.
csv_handle = csv_handler("dummy")

# Perform operation to get data. EXAMPLE USES OF THE PROGRAM ARE DISPLAYED BELOW! Just copy-paste it and input the CSV files that you want to use.

# REGULAR SINGLE PLOT, NO LABEL. FILL PARAMETERS!
FILETOHANDLE1 = "unfiltered bulb.csv"
RAWDATA = csv_handle.csv_reader(FILETOHANDLE1)
RANGE = [279,1000]
TITLE = "Unfiltered Bulb Spectrum"
istransmission = "false"
majorxmajoryminorxminorylocatormatrix = [100, 0.2, 5, 0.1]
label = "false"
peakscale = 32
savgolwindow = 7

#csv_handle.grapher(istransmission, RANGE, [0,1], TITLE, label, "r", 0.5, RAWDATA[0], FILETOHANDLE1 + ".png", majorxmajoryminorxminorylocatormatrix)
csv_handle.peakgrapher(istransmission, RANGE, [0,1], TITLE, label, "r", 0.5, RAWDATA[0], FILETOHANDLE1 + "peak.png", majorxmajoryminorxminorylocatormatrix, peakscale, savgolwindow)

# TRANSMISSION CURVE PLOT WITH A LABEL
#TRANSIRANGE = [200, 1000]
#locatrix = [100, 0.2, 5, 0.1]
FILEUNFILT, FILETRANS = "unfiltered bulb.csv", "filter5fw1.csv"
tcurve_combo14 = csv_handle.transmission_curve(FILEUNFILT, FILETRANS)
csv_handle.grapher("true", TRANSIRANGE, [0,1], "Filter 5 Transmission Curve", "false", "b", 0.2, tcurve_combo14[0], "tcurve5.png", locatrix)

# MULTIPLIED CURVE PLOT WITH A LABEL! In this case, I've first multiplied FW 1 on the default spectrum, then FW4, giving what should be 1,4 in the process! We'll see.
# Anyway it isn't. Whatever. The multiplication method is sound and it works, you get what you pay for.
#tcurve_1 = csv_handle.transmission_curve("bulb.csv", "filter1.csv")
#tcurve_4 = csv_handle.cransmission_curve("unfiltered bulb.csv", "filter4fw1.csv")
#applied1 = csv_handle.multi_curve("unfiltered bulb.csv", tcurve_1[1])
#applied14 = csv_handle.multi_curve(applied1[1], tcurve_4[1])
#tcurvepred14 = csv_handle.transmission_curve("unfiltered bulb.csv", applied14[1])
#csv_handle.grapher("true", [0,1000], [0,1], "Transmission Curve (predicted) of filter wheel combination '1,4'", "FW Pred Combo '1,4'", "b", 0.2, tcurvepred14[0], "tcurvepred14.png", locatrix)

# Planck Fit
#skyspectrum = csv_handle.csv_reader("unfiltered bulb.csv")
#hello = csv_handle.least_solver(skyspectrum[0], "Unfiltered Bulb")
#hello = csv_handle.scipy_solver(skyspectrum[0])