import numpy as np

import ascii_info
import galcentricutils
import graphutils
import hdfutils
import windows_directories
import matplotlib.pyplot as plt
import astropy.io.misc.hdf5

#n_points = 2.4e4
#coordgrid = galcentricutils.greatcount().fullgridgen(n_points)

writer = hdfutils.hdf5_writer(windows_directories.datadir, ascii_info.asciiname)
table = writer.read_table("full_raw","GCC_TEST_CONFIDENCE")
thetas, phis, indices, sigmins = table['theta'], table['phi'], table['index'], table['sigmin']
fig, axs = plt.subplots(nrows=1, ncols=1, dpi=300)
h = axs.scatter(phis, thetas, c=sigmins, marker='s', s=0.1, cmap='viridis')
fig.colorbar(h, label=r'${\log_{10}}{\frac{1}{p_{gc}}}$')
axs.grid(False)
axs.set(xlabel=r'$\phi$',
        ylabel=r'$\theta$')
plt.show()