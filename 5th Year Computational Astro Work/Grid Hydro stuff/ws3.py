import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
import os
from matplotlib.pyplot import *
import matplotlib.colors as colors
from matplotlib import rc, ticker, cm, patches
import pdb
import netCDF4 as netCDF
from pylab import *
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

# First part
"""
Equation to solve is
du/dt = -u
General solution:
u(t) = Ae^(-t)
"""
# Define left-side boundary condition and analytic solution
boundt, boundu = 0,1
A = boundu/np.exp(-1*boundt)
u = lambda t: A*np.exp(-1*t)


# Define range and step and get the analytic solution for this range
dt = 0.045
tlims = [0,12]
n_steps = np.rint((tlims[1]-tlims[0])/dt,out=np.zeros(1,int),casting='unsafe')[0]
tspread = np.linspace(*tlims, n_steps + 1)
analytic_solution = u(tspread)

# Set up various integrators
solutions = [[boundu] for d in range(3)] # forward, backward, centered
implicit_new = [np.exp(-dt),boundu] # before 0 by dt and also the current one at t=0.

# All integration schemes
for step in range(n_steps):
    forward = solutions[0][step]*(1-dt)
    backwrd = solutions[1][step]/(1+dt)
    centerd = solutions[2][step]*(2-dt)/(2+dt)
    solutions[0].append(forward)
    solutions[1].append(backwrd)
    solutions[2].append(centerd)
    implciitnew = -1*(2*dt*implicit_new[step]) + implicit_new[step-1]
    implicit_new.append(implciitnew)

# Get fractional error to true analytic value
solutions = [np.array(d) for d in solutions]
errors = [(solution - analytic_solution)/analytic_solution for solution in solutions]
solerrs = [solutions,errors]

# Get the graphs
fig, axs = plt.subplots(nrows=2,ncols=1, sharex='all',figsize=(7,10))
plt.subplots_adjust(wspace=0, hspace=0.05)
for ax in axs:
    ax.tick_params(axis="x", which="both", direction="in", length=4, bottom=True, left=True, right=True,
                   top=True)
    ax.tick_params(axis="y", which="both", direction="in", length=4, bottom=True, left=True, right=True,
                   top=True)

# Set up the ticker for xaxis
xticks = ticker.MultipleLocator(2)
xticksminor = ticker.MultipleLocator(0.5)
for ax in axs:
    ax.xaxis.set_major_locator(xticks)
    ax.xaxis.set_minor_locator(xticksminor)

# Now do the yaxis
yticks = ticker.MultipleLocator(0.2)
yticksminor = ticker.MultipleLocator(0.1)
axs[0].yaxis.set_major_locator(yticks)
axs[0].yaxis.set_minor_locator(yticksminor)
yticks1 = ticker.MultipleLocator(0.1)
yticksminor1 = ticker.MultipleLocator(0.05)
axs[1].yaxis.set_major_locator(yticks1)
axs[1].yaxis.set_minor_locator(yticksminor1)

# For the solutions
linestyles = ['dashed','dashdot','dotted']

# Plot solutions and set xlims
for axnum,ax in enumerate(axs):
    for num,linestyle in enumerate(linestyles):
        ax.plot(tspread, solerrs[axnum][num], linestyle=linestyle, color='black')
        ax.set(xlim=tlims)

# Set the ylims
axs[0].set(ylim=[-1.5,1])
axs[1].set(ylim=[-1.5,1])

# Do the analytics, too
axs[0].plot(tspread, analytic_solution, color='black')
axs[1].axhline(0, color='black')

# Add some labels
axs[0].set(ylabel=r'$u(t)$')
axs[1].set(ylabel=r'$\epsilon(t)/w(t)$',
           xlabel=r'$t$')

plt.savefig("top_right.png", dpi=300)
plt.show()

plt.clf()
plt.close()
plt.plot(tspread, implicit_new[1:], color='black')
plt.plot(tspread, analytic_solution, color='red')
plt.ylim([-1.5,1])
plt.savefig("implicitnewnew.png", dpi=300)
plt.show()