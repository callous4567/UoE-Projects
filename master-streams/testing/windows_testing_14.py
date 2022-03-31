
# Save the table for the sake of effort
import pickle

from matplotlib import pyplot as plt

import windows_directories
with open(windows_directories.datadir + "\\multistream_dump_onlyhalo.txt", "rb") as f:
    mds, mnfws, cnfws = pickle.load(f)

plt.hist(cnfws)
plt.show()