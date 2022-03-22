

# The above, but for the full group/data. More convenient than dealing with stacking/etc.
# TODO: >>>>>>>>>>>>>IMPORTANT NOTE FOR THE FINAL REPORT HERE!!>>>>>>>>>>>>>
"""
We ignore the fact that error correlations/etc may be different for each source of data
Bring this up and ask Jorge specifically about how to report on this. See our notes.
"""
import multiprocessing
import ascii_info
import windows_multiprocessing

arrayinfominpars = []
group = ascii_info.fullgroup
minpar = ascii_info.fulldata_minpar_L4D
for saveid in ascii_info.duplimonte_L4D_saveids:
    arrayinfominpars.append([[group,saveid],minpar])

# Regularly map/pool :)
if __name__ == '__main__':
    pool = multiprocessing.Pool(6)
    results = pool.map(windows_multiprocessing.do_hdbscan_L4D, arrayinfominpars)
    pool.close()