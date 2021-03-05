import matplotlib.pyplot as plt
import matplotlib.ticker as pltick
import numpy as np
# Okidoke let's go.
fig = plt.figure(1) # Defines figure
ax = plt.subplot(111) # Defines the subplot in our figure
ax.set(xlim=[0, 50], ylim=[0, 50], title=("Test Graph {}").format(r'$\alpha$'), xlabel="X", ylabel="Y")
func_e = lambda e: np.cos(2*np.pi*e)*np.exp(np.sin(e))
x = np.arange(0, 50, 0.01)
ax.plot(x, [float(15)*func_e(x) for x in x], label="Cool Graph", color='red', linewidth=0.5, marker='x', markevery=5, linestyle="--")
ax.grid(True, which='major')
ax.grid(True, which='minor')
majorticks = pltick.MultipleLocator(5)
minorticks = pltick.MultipleLocator(1)
ax.minorticks_on()
ax.xaxis.set_major_locator(majorticks)
ax.yaxis.set_major_locator(majorticks)
ax.xaxis.set_minor_locator(minorticks)
ax.yaxis.set_minor_locator(minorticks)
ax.margins(x=0, y=0)
plt.show()
