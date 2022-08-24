import pickle
import numpy as np
import plotly.express as px
import plotly.io as plio
plio.renderers.default = 'browser'
import plotly.graph_objects as go
import ascii_info_new
import galcentricutils_new
import graphutils_new
import hdfutils
import windows_directories_new
import chart_studio
import chart_studio.plotly as py

chart_studio.tools.set_credentials_file(username='Callicious', api_key='tm4K20edMOkiH8hfgsO6')
chart_studio.tools.set_config_file(world_readable=True, sharing='public')

table = hdfutils.hdf5_writer(windows_directories_new.datadir,
                             ascii_info_new.asciiname).read_table(ascii_info_new.fullgroup,
                                                              ascii_info_new.fullset)
with open(windows_directories_new.clusterdir + "\\" + "fullgroup.cluster.txt", 'rb') as f:
    clustered = pickle.load(file=f)
    numclust = galcentricutils_new.compclust().nclust_get(clustered)


# Set up data and run the clustering
panda = hdfutils.hdf5_writer(windows_directories_new.datadir,
                             ascii_info_new.asciiname).read_df(ascii_info_new.fullgroup,
                                                           ascii_info_new.panda_raw)
data = np.array(panda['vec_L'])
data = list(data)


# Generate a HTML for this clustering under imgdir\\clustered and embed the number of clusters found.
array = np.array(data)
clusterdata = clustered

# Remove all the noise elements, and also remove the GSE.
remove = [-1, 16]
truefalse = [False if d in remove else True for d in clusterdata]
array, clusterdata = array[truefalse], clusterdata[truefalse]
array = np.array(array)

clusterdata = [-1 for d in clusterdata]
array += np.random.default_rng().normal(0, 1000, np.shape(array))

"""
# Get all clusts
all_clusters = []
for clust in clusterdata:
    if clust not in all_clusters:
        all_clusters.append(clust) 

# Define rng
rng = np.random.default_rng()

# On a per-cluster basis, obfuscate the data a bunch (to prevent people managing to successfully grab it.
for clust in all_clusters:
    truefalse = [False if clust == d else True for d in clusterdata]
    array[truefalse] += rng.normal(loc=0, scale=2000, size=3)
    array[truefalse] += rng.normal(loc=0, scale=50, size=np.shape(array[truefalse])) """

# Need arrays. Make sure.
if type(array) != "numpy.ndarray":
    array = np.array(array)
if type(clusterdata) != "numpy.ndarray":
    clusterdata = np.array(clusterdata)

# Verify if array is only Lx Ly Lz or if it also energy, and if so, then clip to just angular momentum.
if len(array[0]) > 3:
    array = array.T
    array = array[0:3]
    array = array.T
else:
    pass

x, y, z = array.T

fig = px.scatter_3d(x=x, y=y, z=z, color=clusterdata)
py.plot(fig, filename="angular_threed_noise", auto_open=True)