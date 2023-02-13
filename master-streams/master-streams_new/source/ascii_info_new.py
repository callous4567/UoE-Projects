""" IMPORTANT NOTES
JORGES PROPER MOTIONS IN L ALREADY HAVE A COSB FACTOR CHRIST ON A STICK
Code has been adapted to remove cos(b) on import.
Double check, though!!!
"""

# TODO LIST
# TODO See galcentricutils_new: there's a few todos

# The index for the entire star catalogue
asciiname = "stardata.hdf5"
flatfork_asciiname = "flatfork_stardata.hdf5"

# Group/set of the combined catalogue ***see LAMOST_2MASS_GAIA... this now has the DR7 LAMOST instead.***
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

with open(windows_directories_new.duplimontedir + "\\" + ascii_info_new.bhb + "\\L_0.txt", 'rb') as f:
    owo = pickle.load(f)
    print(owo)
    
add .cluster.txt for cluster list
add .txt for pure [l1,l2,l3] vector list 
"""

duplimonte_number = 1000
duplimonte_saveids = [("L_{}").format(d) for d in [str(d) for d in range(duplimonte_number)]]
duplimonte_L4D_saveids = [("L4D_{}").format(d) for d in [str(d) for d in range(duplimonte_number)]]
duplimonte_LE_saveids = [("LE_{}").format(d) for d in [str(d) for d in range(duplimonte_number)]]
duplimonte_LXYZ_saveids = [("LXYZ_{}").format(d) for d in [str(d) for d in range(duplimonte_number)]]

# Set minpars for each group for cluster_duplimonte. minpar is SIZE,SAMPLES.
bhb_minpar = [8,7]
gcs_minpar = [5,5]
kgiant_minpar = [8,10]
lamostk_minpar = [8,10]
fulldata_oldminpar = [15,15]
fulldata_minpar = [8, 7] # using [8,9] for the initial clustering, then [8,10] for the monte carlo'd ones. # 8,10 is our final borderline acceptance (which works with our data cleaning.). 15,15 is the original. [15,15]
fulldata_minpar_L4D = fulldata_minpar
fulldata_minpar_LE = fulldata_minpar
minpars_allgroups = [bhb_minpar, gcs_minpar, kgiant_minpar, lamostk_minpar]

# ID for the GSE in the flatfork
flatfork_GSE_ID = 23

# Custom-set the range for orbifits/etc
import numpy as np
clusters_to_maindata_orbifit = np.arange(0,24,1)
flatfork_clusters_to_maindata_orbifit = np.arange(0,24,1)

# The save-ids for all the generated orbit fits
n_carlo = 50
orbifit_saveids = [("orbifit_instance_{0:.0f}").format(d) for d in range(n_carlo)]
flatfork_orbifit_saveids = [("flatfork_orbifit_instance_{0:.0f}").format(d) for d in range(n_carlo)]

# For the maindata
n_carlo_maindata = 50
orbifit_maindata_saveids = [("orbifit_instance_maindata_{0:.0f}").format(d) for d in range(n_carlo_maindata)]
flatfork_orbifit_maindata_saveids = [("flatfork_orbifit_instance_maindata_{0:.0f}").format(d) for d in range(n_carlo_maindata)]

# Save the monte results
# monte_fits_address = os.path.join(windows_directories_new.datadir, fullgroup + "_" + fullset +".fits")