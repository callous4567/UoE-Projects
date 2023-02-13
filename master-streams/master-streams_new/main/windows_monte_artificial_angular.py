import multiprocessing
import ascii_info_new
import windows_multiprocessing_new

"""
Just to generate artificial momentum sets for clustering. 
Now also generates the spatial bit, too, and packages LXYZ with XYZ inside a table and saves this, too.
Important note: You need to ensure that panda_raw has covariance matrices and means already in there.
Just run monte_multitable to set this up and you'll be fine (that generates panda_raw anyway.)
"""

"""
m = ascii_info_new.duplimonte_number

# Set up list to iterate over
groupsetm_list = []
for group in ascii_info_new.all_groups:
    groupsetm = [group, ascii_info_new.panda_raw, m]
    groupsetm_list.append(groupsetm)

# Regularly map/pool :)
if __name__ == '__main__':
    pool = multiprocessing.Pool(6)
    pool.map(windows_multiprocessing_new.do_duplimonte_table, groupsetm_list)
    pool.close() """


# The above, but for the fullgroup/set instead.
m = ascii_info_new.duplimonte_number

# Set up list to iterate over
groupsetm_list = []
for group in [ascii_info_new.fullgroup]:
    groupsetm = [group, ascii_info_new.panda_raw, m]
    groupsetm_list.append(groupsetm)

# Regularly map/pool :)
if __name__ == '__main__':
    pool = multiprocessing.Pool(8)
    pool.map(windows_multiprocessing_new.do_duplimonte_table, groupsetm_list)
    pool.close()

