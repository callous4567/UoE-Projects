# The index for the entire star catalogue
asciiname = "stardata.hdf5"

# Group/set of the combined catalogue
fullgroup, fullset = "full_raw", "astrotable"

# Set indices (same asciigroup) for the raw data for each individual ascii file we ripped
bhb = "BHB_edr3_metal_superclean"
gcs = "GCs_edr3_name"
kgiants = "KGiant_edr3_metal"
lamostk = "LAMOST_K_FULL_edr3"
set_raw = "astrotable"
all_groups = [bhb, gcs, kgiants, lamostk]

""" IMPORTANT NOTES
JORGES PROPER MOTIONS IN L ALREADY HAVE A COSB FACTOR CHRIST ON A STICK
Code has been adapted to remove cos(b) on import.
Double check, though!!!
"""