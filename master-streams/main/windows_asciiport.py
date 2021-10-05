import asciiutils
import os
from windows_directories import datadir, asciidir
from ascii_info import asciiname, asciigroup, asciiset

# Future notes for using columns
"""
To grab the data:
go under group
grab raw via hdf table reader (this is the astropy table)

To grab the unit type for each column:
go under group
use column id as dataset
read string using this dataset
returns original unit string from ascii file
"""

# Navigate to asciidir, grab filenames with/without .txt
os.chdir(asciidir)
ascii_list = os.listdir()

# Set up ascii-port objects for all of them.
asciis = [asciiutils.ascii_port(asciidir, d) for d in ascii_list]

# Grab the name of the ascii without the .txt extension
stringnames = [d.replace(".txt","") for d in ascii_list]

# Get ascii tables and save them (individually for debug)
for n,d in enumerate(asciis):
    d.astrotabify()
    d.save(datadir, asciiname, asciigroup, stringnames[n])

# Add all ascii tables to the first ascii table in the list.
for num, ascii_object in enumerate(asciis):
    if num != 0:
        asciis[0].add(ascii_object)

# Save the first ascii table in full: see ascii_info asciiname
asciis[0].save(datadir, asciiname, asciigroup, asciiset)
