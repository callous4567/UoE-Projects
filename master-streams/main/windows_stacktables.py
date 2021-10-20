from astropy.table import vstack
import ascii_info
import hdfutils
import windows_directories

# While we're at it, recreate the full table.
writer = hdfutils.hdf5_writer(windows_directories.datadir, ascii_info.asciiname)
tables = [writer.read_table(group, ascii_info.set_raw) for group in ascii_info.all_groups]
for num,table in enumerate(tables):
    if num != 1:
        tables[0] = vstack([tables[0],table])
writer.write_table(ascii_info.fullgroup, ascii_info.fullset, tables[0])