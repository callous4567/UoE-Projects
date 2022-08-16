import multiprocessing
import warnings
import windows_multiprocessing
import ascii_info

"""
This file will run clustering for all the individual duplimonte'd angulars we have and produce clustered lists (pickled)
"""

# Suppress Loky Warnings
warnings.filterwarnings("ignore", message="Loky-backed parallel loops cannot be called in a multiprocessing, setting n_jobs=1")

arrayinfominpars = []
group = ascii_info.fullgroup
minpar = ascii_info.fulldata_minpar
for saveid in ascii_info.duplimonte_saveids:
    arrayinfominpars.append([[group,saveid],minpar])

# Regularly map/pool :)
if __name__ == '__main__':
    pool = multiprocessing.Pool(8)
    results = pool.map(windows_multiprocessing.flatfork_do_hdbscan, arrayinfominpars)
    pool.close()



