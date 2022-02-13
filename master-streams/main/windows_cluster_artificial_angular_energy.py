import multiprocessing
import ascii_info
import windows_multiprocessing

# The above, but for the full group/data. More convenient than dealing with stacking/etc.
# TODO: >>>>>>>>>>>>>IMPORTANT NOTE FOR THE FINAL REPORT HERE!!>>>>>>>>>>>>>
"""
We ignore the fact that error correlations/etc may be different for each source of data
Bring this up and ask Jorge specifically about how to report on this. See our notes.
"""

arrayinfominpars = []
group = ascii_info.fullgroup
minpar = ascii_info.fulldata_minpar_LE
for saveid in ascii_info.duplimonte_LE_saveids:
    arrayinfominpars.append([[group,saveid],minpar])

# Regularly map/pool :)
if __name__ == '__main__':
    pool = multiprocessing.Pool(6)
    results = pool.map(windows_multiprocessing.do_hdbscan_LE, arrayinfominpars)
    pool.close()


