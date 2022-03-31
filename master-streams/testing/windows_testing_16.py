import pickle

import numpy as np

# Plotly
import plotly.express as px
import plotly.io as plio
plio.renderers.default = 'browser'
import plotly.graph_objects as go

import ascii_info
import galcentricutils
import graphutils
import hdfutils
import windows_directories


# Set up data and run the clustering
panda = hdfutils.hdf5_writer(windows_directories.datadir,
                             ascii_info.asciiname).read_df(ascii_info.fullgroup,
                                                           ascii_info.panda_raw)

with open(windows_directories.clusterdir + "\\" + "fullgroup.cluster.txt", 'rb') as f:
    clustered = pickle.load(f)


# Generate a HTML for this clustering under imgdir\\clustered and embed the number of clusters found.
data = list(panda['vec_L'].to_numpy())
graphutils.threed_graph().kmeans_L_array(data, clustered, savedexdir=False, browser=True)

# Take the panda and clustering, sort it, then plot each clustering
panda['cluster'] = clustered
xyz = np.array([panda['x'].to_numpy(),panda['y'].to_numpy(),panda['z'].to_numpy()]).T
groups = panda.groupby(by='cluster')
for group in range(10):
    sub_panda = groups.get_group(group)
    x, y, z = sub_panda['x'].to_numpy(), sub_panda['y'].to_numpy(), sub_panda['z'].to_numpy()
    fig = px.scatter_3d(x=x, y=y, z=z, color=[group for d in x])
    fig.show()