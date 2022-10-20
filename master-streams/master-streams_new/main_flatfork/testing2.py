
import os
import windows_directories_new
owd = windows_directories_new.sourcedir


import numpy as np


from galcentricutils_new import compclust
cc = compclust()
import pickle
with open(os.path.join(owd, "dump1.txt"),"rb") as f:
    cut_clustered = pickle.load(file=f)
with open(os.path.join(owd, "dump2.txt"),"rb") as f:
    trial_clustered = pickle.load(file=f)

cc.compclust_multilabel_julia(cut_clustered, trial_clustered, 30)
cc.compclust_multilabel(cut_clustered, trial_clustered, 30)
def test_pure():
    cc.compclust_multilabel(cut_clustered, trial_clustered, 30)
def test_julia():
    np.array(cc.compclust_multilabel_julia(cut_clustered, trial_clustered, 30))
import timeit
print(timeit.timeit(test_pure, number=1000))
print(timeit.timeit(test_julia, number=1000))

"""

import os, windows_directories_new
from julia import Main
import numpy as np
Main.include(windows_directories_new.jl_dir + "\\munkres.jl")
vv = Main.solve_hungarian(np.ones((50,25), np.int64))
"""