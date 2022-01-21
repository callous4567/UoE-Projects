""" IMPORTANT NOTES
JORGES PROPER MOTIONS IN L ALREADY HAVE A COSB FACTOR CHRIST ON A STICK
Code has been adapted to remove cos(b) on import.
Double check, though!!!
"""

# TODO LIST
# TODO See galcentricutils: there's a few todos

# The index for the entire star catalogue
asciiname = "stardata.hdf5"

# Group/set of the combined catalogue
fullgroup, fullset = "full_raw", "astrotable"
fullpanda = "astroframe"

# Set indices (same asciigroup) for the raw data for each individual ascii file we ripped
bhb = "BHB_edr3_metal_superclean"
gcs = "GCs_edr3_name"
kgiants = "KGiant_edr3_metal"
lamostk = "LAMOST_K_FULL_edr3"
set_raw = "astrotable"
all_groups = [bhb, gcs, kgiants, lamostk]
panda_raw = "astroframe"

# Just some information relating to the duplimonte number. Note that duplimontes are stored as .txt pickles.
""" using pickle

with open(windows_directories.duplimontedir + "\\" + ascii_info.bhb + "\\L_0.txt", 'rb') as f:
    owo = pickle.load(f)
    print(owo)
    
add .cluster.txt for cluster list
add .txt for pure [l1,l2,l3] vector list 
"""

duplimonte_number = 200
duplimonte_saveids = [("L_{}").format(d) for d in [str(d) for d in range(duplimonte_number)]]

# Set minpars for each group for cluster_duplimonte. minpar is SIZE,SAMPLES.
bhb_minpar = [9,7]
gcs_minpar = [5,5]
kgiant_minpar = [8,12]
lamostk_minpar = [8,12]
fulldata_minpar = [15,15]
minpars_allgroups = [bhb_minpar, gcs_minpar, kgiant_minpar, lamostk_minpar]


