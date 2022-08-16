import os
import windows_directories
import numpy as np
clust1 = np.array([1,2,3,4,5,6,6,6,6,7,8,7,7,7])
clust2 = np.array([2,3,4,5,6,7,7,7,7,8,9,8,8,8])

from juliacall import Main, convert
Main.include(os.path.join(windows_directories.jl_dir, "munkres.jl"))
uuu = np.array(Main.compclust_multilabel(clust1, clust2, 9999))
