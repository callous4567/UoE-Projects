import ascii_info
import galcentricutils
import hdfutils
import windows_directories
import numpy as np
import pickle
import time
# Plotly
import plotly.express as px
import plotly.io as plio
plio.renderers.default = 'browser'
import plotly.graph_objects as go

# Some useful examples of Plotly 3D plotting- keep this for reference.

"""
data = hdfutils.hdf5_writer(windows_directories.datadir,
                            ascii_info.asciiname).read_df(ascii_info.fullgroup,
                                                          ascii_info.fullpanda)
table = hdfutils.hdf5_writer(windows_directories.datadir,
                            ascii_info.asciiname).read_table(ascii_info.fullgroup,
                                                             ascii_info.fullset)
row = table[10]
l,b,distance,dmul,dmub,vlos,edist,edmul,edmub,evlos = row['l'],row['b'],row['dist'],row['dmu_l'],row['dmu_b'],\
                                                      row['vlos'],row['edist'],row['edmu_l'],row['edmu_b'],\
                                                      row['evlost']
vec = [l,b,distance,dmul,dmub,vlos,edist,edmul,edmub,evlos]
cov, stds, mu = galcentricutils.monte_angular().vec_covmonte(vec, 400) # cov, stdevs, vec_L
rng = np.random.default_rng()
points = rng.multivariate_normal(mean=mu,cov=cov,size=1000) """

# Set up traces for each duplitable to visualize errors
traces = []
for id in ascii_info.duplimonte_saveids:
    with open(windows_directories.duplimontedir + "\\" + ascii_info.bhb + "\\" + id + ".txt", 'rb') as f:
        vec_L = pickle.load(f)
    vec_L = np.array(vec_L)
    vecx, vecy, vecz = vec_L.T
    vectrace = go.Scatter3d(x=vecx,
                            y=vecy,
                            z=vecz,
                            name="Data",
                            marker=dict(color='rgb(34,163,192)'),
                            mode='markers')
    traces.append(vectrace)
fig = go.Figure(data=traces)
fig.show()
