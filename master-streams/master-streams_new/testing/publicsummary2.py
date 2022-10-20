import pickle
import numpy as np
import plotly.express as px
import plotly.io as plio
from IPython.core.display import display
from galpy import orbit
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation

from energistics_new import orbigistics, orbifitter

plio.renderers.default = 'browser'
import plotly.graph_objects as go
import ascii_info_new
import galcentricutils_new
import graphutils_new
import hdfutils
import windows_directories_new
import chart_studio
import chart_studio.plotly as py
import astropy.units as u

chart_studio.tools.set_credentials_file(username='Callicious', api_key='tm4K20edMOkiH8hfgsO6')
chart_studio.tools.set_config_file(world_readable=True, sharing='public')

table = hdfutils.hdf5_writer(windows_directories_new.datadir,
                             ascii_info_new.asciiname).read_table(ascii_info_new.fullgroup,
                                                              ascii_info_new.fullset)
with open(windows_directories_new.clusterdir + "\\" + "fullgroup.cluster.txt", 'rb') as f:
    clustered = pickle.load(file=f)

# Load in the "mean" percenttable and map
writer = hdfutils.hdf5_writer(windows_directories_new.datadir, ascii_info_new.asciiname)
membership_table = writer.read_table(ascii_info_new.fullgroup, "percent_table_greatfitted")
clustered_final = membership_table['probable_clust']

# Get orbifitter
orbifit = orbifitter()

# Fit orbit (for example case! :D)
fit = orbifit.galpy_final_fitting(table, 1, 2000, 1e9, 1000, True, False, True, True)

# Get just for the orphan
truefalse = [True if d in [1] else False for d in clustered_final]
table = table[truefalse]
orbigist = orbigistics()
orbits = orbigist.orbits(table)
times = np.linspace(0, 4e9, 1000) * u.yr


fit.integrate(t=times, pot=orbigist.pot)
lls, bbs = fit.ll(np.linspace(0, 4e9, 10000)*u.yr), fit.bb(np.linspace(0, 4e9, 10000)*u.yr)

# Plots/etc
orbits.integrate(t=times, pot=orbigist.pot)

# Set axes 
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 8))

# Scatter the "true orbit"
ax.scatter(lls, bbs, color='green', s=10, marker='x', zorder=1)

# Make it pretty
ax.set(xlabel="l / deg",
       ylabel="b / deg",
       xlim=[0,360],
       ylim=[-90,90])
ax.set_facecolor("k")
plt.gca().invert_xaxis()
plt.gcf().set_facecolor('k')

# Set up figure, axes, etc
axis = ax.scatter(orbits.ll(0), orbits.bb(0), s=10, marker='x', color='red', zorder=2)

# Now, we need to generate some false orbits
fake = table[0]
fake['l'] -= 5
fake['b'] -= 5
fake['dist'] *= 1.2
fake['dmu_l'] *= np.cos(np.degrees(fake['b']))
fake['dmu_l'] *= 0.4
fake['dmu_b'] *= 2
vxvv = np.array([fake['l'], fake['b'], fake['dist'], fake['dmu_l'], fake['dmu_b'], fake['vlos']])
newfake = orbit.Orbit(vxvv=vxvv, ro=orbigist.rovo[0],
                   vo=orbigist.rovo[1], zo=orbigist.zo, lb=True)
newfake.integrate(np.linspace(0, 20e9, 10000)*u.yr, orbigist.pot)
lls, bbs = newfake.ll(np.linspace(0, 20e9, 20000)*u.yr), newfake.bb(np.linspace(0, 20e9, 20000)*u.yr)
ax.scatter(lls, bbs, color='yellow', s=3, marker='x', zorder=0)



def animate(i):
    # Run it until we hit an absorbing state.
    print(i)
    ll, bb = orbits.ll(times[i]).value, orbits.bb(times[i]).value
    data = np.array([ll,bb]).T
    axis.set_offsets(data)
    return axis,

anim = FuncAnimation(fig=fig, func=animate, repeat=False, interval=2000 / 30, frames=len(times)-10, blit=True)
plt.show()
anim.save("orphan_test.mp4")