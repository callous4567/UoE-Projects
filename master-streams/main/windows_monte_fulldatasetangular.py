import multiprocessing
import ascii_info
import windows_multiprocessing

# Just to generate artificial momentum sets for clustering.
# Important note: You need to ensure that panda_raw has covariance matrices and means already in there.
# Just run monte_multitable to set this up and you'll be fine (that generates panda_raw anyway.)

"""
m = ascii_info.duplimonte_number

# Set up list to iterate over
groupsetm_list = []
for group in ascii_info.all_groups:
    groupsetm = [group, ascii_info.panda_raw, m]
    groupsetm_list.append(groupsetm)

# Regularly map/pool :)
if __name__ == '__main__':
    pool = multiprocessing.Pool(6)
    pool.map(windows_multiprocessing.do_duplimonte_table, groupsetm_list)
    pool.close() """

# The above, but for the fullgroup/set instead.
m = ascii_info.duplimonte_number

# Set up list to iterate over
groupsetm_list = []
for group in [ascii_info.fullgroup]:
    groupsetm = [group, ascii_info.panda_raw, m]
    groupsetm_list.append(groupsetm)

# Regularly map/pool :)
if __name__ == '__main__':
    pool = multiprocessing.Pool(6)
    pool.map(windows_multiprocessing.do_duplimonte_table, groupsetm_list)
    pool.close()


