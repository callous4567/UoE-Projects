import matplotlib.pyplot as plt
import matplotlib.ticker as pltick
import numpy as np
import sys
import time
import random

# Raw data
times = np.array([0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120]) # Defines our time period, 0 to 120, step of 10.
temp_list_1 = np.array([11.5, 11.4, 11.4, 11.4, 11.4, 11.4, 11.4, 11.4, 11.5, 12.0, 13.0, 13.4, 14.4])
temp_list_2 = np.array([11.6, 11.6, 11.6, 11.6, 11.6, 11.6, 11.6, 11.6, 11.6, 11.7, 11.7, 12.0, 12.4])
temp_list_3 = np.array([12.0, 12.0, 12.0, 12.0, 12.0, 12.1, 12.1, 12.1, 12.1, 12.1, 12.2, 11.4, 12.7])
temp_list_4 = np.array([15.2, 16.0, 16.0, 16.1, 16.2, 16.3, 16.4, 16.5, 16.6, 16.7, 16.9, 16.9, 17.1])
temp_list_5 = np.array([30.9, 31.1, 31.2, 31.4, 31.4, 31.5, 31.5, 31.6, 31.7, 31.8, 31.7, 31.9, 31.9]) # You could set up temp_list with the import of an XLSX file or a text document, etc. Won't bother.
plume_height = np.array([0, 1.5, 1.5, 1.5, 1.5, 2.0, 3.0, 3.5, 7.0, 9.0, 9.0, 9.0, 9.0]) # All heights in cm
# Pre-measurement temperatures can be found on the sheet, as can equilibrium temperatures.
def list_former(list):
    new_list = [d for d in list]
    return new_list

temp_list_10 = list_former(temp_list_1)
temp_list_20 = list_former(temp_list_2)
temp_list_30 = list_former(temp_list_3)
temp_list_40 = list_former(temp_list_4)
temp_list_50 = list_former(temp_list_5)

print(len(temp_list_1), len(temp_list_2), len(temp_list_3), len(temp_list_4), len(temp_list_5), len(times), len(plume_height))
# Figure
fig1 = plt.figure(111, figsize=(24, 12))
ax1 = plt.subplot(111)
time_deg = "Temperature $^\circle$C"
time = "$\it{t}$"

# Heat flux formula is H_0 = cdeltaz=dt sum of chagne in Ti

constant_c = 4.2 # JC-1cm%-3
d_z = 2 # In centimetres, distance between the thermistors
d_t = 120 # Time between initial and final readings

new_constant = d_z*constant_c*(1/d_t)
difference = lambda d: d[12] - d[0]
summation = difference(temp_list_10) + difference(temp_list_20) + difference(temp_list_30) + difference(temp_list_40) + difference(temp_list_50)
H_value = new_constant*summation
print(H_value)
# About 0.5 Wcm^-2... okidoke. That's about half <3
# Note something that you should remember is that you need to take account for the difference integrated over time... so...
# So yeah. Let's format the graphs a notch.


ax1.set(xlim=[0,120], ylim=[10, 35], title="Graph of Temperature against Time for Lab 5: Convection, with Thermistors in order of 1 -> 5 for ascending height")

ax1.plot(times, temp_list_10, label="Thermistor #1", color="purple")

ax1.plot(times, temp_list_20,label="Thermistor #2", color="blue")

ax1.plot(times,temp_list_30, label="Thermistor #3", color="green")

ax1.plot(times,temp_list_40, label="Thermistor #4", color="orange")

ax1.plot(times,temp_list_50, label="Thermistor #5", color="red")

ax1.text(0, 25, ("$H_0$ value is equal to " + ("{0:.2f} W/cm").format(H_value)))

major_x = pltick.MultipleLocator(10)
major_y = pltick.MultipleLocator(5)
minor_x = pltick.MultipleLocator(10)
minor_y = pltick.MultipleLocator(0.2)
ax1.xaxis.set_major_locator(major_x)
ax1.xaxis.set_minor_locator(minor_x)
ax1.yaxis.set_major_locator(major_y)
ax1.yaxis.set_minor_locator(minor_y)
ax1.grid(True, which="minor")
ax1.legend()
ax1.grid(True, which="major")
ax1.margins(0,0)
fig1.savefig("temperatures.png")
plt.show()



# Okay thats thermometer height and temperature.
# Plume heights!

fig2 = plt.figure(222, figsize=(24, 16))
ax2 = plt.subplot(111)
ax2.plot(times, plume_height, label="Plume height versus time", linestyle="--", color="red")
ax2.set(xlim=[0, 120], ylim=[0, 10], title="Graphing plume height as a function of Time", xlabel="Time (s)", ylabel="Plume Height (cm)")

ax2.grid(True, which='major')
ax2.grid(True, which='minor')
major_x = pltick.MultipleLocator(10)
major_y = pltick.MultipleLocator(1)
minor_x = pltick.MultipleLocator(10)
minor_y = pltick.MultipleLocator(0.2)
ax2.xaxis.set_major_locator(major_x)
ax2.xaxis.set_minor_locator(minor_x)
ax2.yaxis.set_major_locator(major_y)
ax2.yaxis.set_minor_locator(minor_y)


ax2.text(0, 9, "Final resting plume height (equilibrium mixing level) equals 9 (cm)")


plt.savefig("height.png")
plt.show()

