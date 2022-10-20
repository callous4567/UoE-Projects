import galcentricutils_new
import numpy as np
import graphutils_new

# Example of generating random gaussian mixtures using genclust3d in galcentricutils_new.
# Keep this.


# Define cluster properties
k = 8
n_points = [30, 100, 1000, 10000, 300, 500, 200, 60]
covscales = [25, 10, 1200, 4000, 340, 1000, 280, 50]
noisescale = 100
n_mesh = 1000
mubounds = np.array([[-4e3, 1e3],
                     [-5e3, 2.5e3],
                     [-5e3, 1e3]]).T

# Generate mucovs
mucovs = galcentricutils_new.genclust3d().gen_mucov(k, covscales, mubounds)

# Generate data
#data = galcentricutils_new.genclust3d().noisy_data(mucovs, n_points, noisescale)

# Do the plot
graphutils_new.threed_graph().kmeans_L_multinormal_generated(mucovs, n_points, n_mesh)
