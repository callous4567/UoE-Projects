import graphutils
from galcentricutils import galconversion, angular
from windows_directories import datadir, sourcedir
from astropy.io import ascii
import hdfutils

"""
Class to handle importing and using Vasilievs "Tango for Three" Simulation Snapshot
"""
class vasiliev(object):
    def __init__(self):
        self.simdir = datadir + "\\vasiliev\\Sgr_snapshot"
        self.simfile = "simdata.hdf5"

    # Import and save the regular Sgr_snapshot
    def import_Sgr_snapshot(self):
        sim_table = ascii.read(self.simdir + "\\stars.txt")
        sim_table['dist'] = sim_table['distance']
        writer = hdfutils.hdf5_writer(self.simdir, self.simfile)
        writer.write_table("Sgr_snapshot", "astrotable", sim_table)
        # Instantiate converter for ICRS/Galactic/Galcentric Polar
        galcon = galconversion()
        galcon.solinfo_grab(sourcedir, "solar_info.dat")
        galcon.solgal_set()
        # Convert ICRS to Galactic.
        galcon.ICRS_to_GAL(self.simdir,
                           self.simfile,
                           "Sgr_snapshot",
                           "astrotable")
        # Go ahead and get galactocentric polars while you're at it (NOT LATIPOLAR.)
        sim_table = writer.read_table("Sgr_snapshot", "astrotable")
        sim_table = angular().get_polar(sim_table)
        writer.write_table("Sgr_snapshot", "astrotable", sim_table)

    # Test/display ra-dec plot for the simulation data
    def test_display(self):
        writer = hdfutils.hdf5_writer(self.simdir, self.simfile)
        sim_table = writer.read_table("Sgr_snapshot", "astrotable")
        graphutils.twod_graph().radecplot(sim_table)