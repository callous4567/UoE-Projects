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
"""
# Define left-side boundary condition and analytic solution
boundt, boundu = 0,1
A = boundu/np.exp(-1*boundt)
u = lambda t: A*np.exp(-1*t)

# Define range and step and get the analytic solution for this range
dt = 0.2
tlims = [0,2.6]
n_steps = np.rint((tlims[1]-tlims[0])/dt,out=np.zeros(1,int),casting='unsafe')[0]
tspread = np.linspace(*tlims, n_steps + 1)
analytic_solution = u(tspread)

# Set up various integrators
solutions = [[boundu] for d in range(3)] # forward, backward, centered

# All integration schemes
for step in range(n_steps):
    forward = solutions[0][step]*(1-dt)
    backwrd = solutions[1][step]/(1+dt)
    centerd = solutions[2][step]*(2-dt)/(2+dt)
    solutions[0].append(forward)
    solutions[1].append(backwrd)
    solutions[2].append(centerd)

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
xticks = ticker.MultipleLocator(0.4)
xticksminor = ticker.MultipleLocator(0.2)
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
axs[0].set(ylim=[0,1])
axs[1].set(ylim=[-0.2,0.2])

# Do the analytics, too
axs[0].plot(tspread, analytic_solution, color='black')
axs[1].axhline(0, color='black')

# Add some labels
axs[0].set(ylabel=r'$u(t)$')
axs[1].set(ylabel=r'$\epsilon(t)/w(t)$',
           xlabel=r'$t$')

# Create the inset
axins = zoomed_inset_axes(axs[0], 10, loc=1)  # zoom=6
for num, linestyle in enumerate(linestyles):
    axins.plot(tspread, solerrs[0][num], linestyle=linestyle, color='black')
axins.plot(tspread, analytic_solution, color='black')
axins.set(ylabel=r'$u(t)$',
          xlabel=r'$t$')
insxlocator = ticker.MultipleLocator(0.02)
axins.xaxis.set_major_locator(insxlocator)
# Region to zoom
xx,yy = [2, 2.12-1e-6],[0.1, 0.16]
axins.set_xlim(*xx)
axins.set_ylim(*yy)

# Add a rectangular patch to finish it up :)
rectangle = patches.Rectangle((xx[0],yy[0]),0.12, 0.06, edgecolor='black',facecolor='none')
axs[0].add_patch(rectangle)

plt.savefig("figure2.2bodenheimer.png", dpi=300)
plt.show()
"""

# Next part. Using periodic boundary conditions as usual.

# v
v = 100

# Steps and number of 'em
dx = 1.0
nx = 600
tlims=[0,4]
dt = 0.002
stability_dt = dx/v
nt = np.rint(max(tlims)/dt, out=np.zeros(1,int), casting='unsafe')[0]
d = v*dt/dx

# Set up domain
x_domain = np.arange(0, dx*nx, dx) # nx elements
rhox_domain = np.arange(dx/2,(dx*(nx-1))-(dx/2) + 1e-10, dx)
lenrho = nx-1
rho_domain = np.zeros(lenrho) # nx-1 elements (staggered to halves)

# Boundary rho
rho_domain[0] = 0
rho_domain[0:6] = 1
rho_domain[6:51] = 0.5

# Set up list of lists for time series
rho_over_time = [rho_domain]

# Step forward rho_over_time. Constant flow in at x=0
# (for some reason this is the only way I can get this to work.)
# (if that isn't what the guy in the book did, then I have no idea how to make my code work)
# (spent a few hours trying but... yeah... it's 6 AM and I need sleep :_:)
for step in range(1,nt):
    # Set up empty rho_over_time
    rho_step = np.zeros(lenrho)
    # Tile it: same periodic method as previously
    rho_tiled = np.tile(rho_over_time[step-1],3)
    rho_tiled[lenrho-1] = 1
    rho_tiled[2*lenrho-10:3*lenrho] = 0
    # Get all new rho elements
    for j in range(lenrho):
        rho_step[j] = -d*(rho_tiled[j] - rho_tiled[(j-1)+lenrho]) + rho_tiled[j]
    # Append
    rho_over_time.append(rho_step)

# Convert a time to a timestep
timestep = lambda t: np.rint(t/dt, out=np.zeros(1,int), casting='unsafe')[0]

# Set up graphing environent
fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(7,5))
xlocators = [MultipleLocator(100),MultipleLocator(20)]
ylocators = [MultipleLocator(0.2),MultipleLocator(0.05)]
ax.xaxis.set_major_locator(xlocators[0])
ax.xaxis.set_minor_locator(xlocators[1])
ax.yaxis.set_major_locator(ylocators[0])
ax.yaxis.set_minor_locator(ylocators[1])
ax.tick_params(axis="x", which="both", direction="in", length=4, bottom=True, left=True, right=True,
               top=True)
ax.tick_params(axis="y", which="both", direction="in", length=4, bottom=True, left=True, right=True,
               top=True)
ax.set(xlim=[-10,410],
       ylim=[-0.05,1.05])

yticks = ax.yaxis.get_major_ticks()
yticks[-1].set_visible(False)

ax.plot(rhox_domain, rho_over_time[0],linestyle='dashed', color='black')
ax.text(20, 0.7, r'$t=0$')
ax.plot(rhox_domain, rho_over_time[timestep(1)], color='black')
ax.text(140, 0.6, r'$t=1$')
ax.plot(rhox_domain, rho_over_time[timestep(2)], color='black')
ax.text(250, 0.5, r'$t=2$')
ax.plot(rhox_domain, rho_over_time[timestep(3)], color='black')
ax.text(360, 0.4, r'$t=3$')
ax.margins(x=0.05, y=0.1)
ax.set(ylabel=r'$\rho(x)$',
       xlabel=r'$x$')
ax.text(335, 0.9, r'$\Delta{t} = 0.002$')
ax.text(335,0.85,r'$\Delta{x} = 1.0$')
plt.savefig("secondplot.png",dpi=300)
plt.show(dpi=300)

