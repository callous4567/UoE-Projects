import ascii_info
import os

# TODO: Now updated with os.path.join such that you can use these on linux.
# TODO: DO THE SAME FOR ALL THE INDIVIDUAL FILES THAT HAVE TEXT ADDITION

# Various Directories and placeholders that are useful. Class for building directories.
rootdir = "C:\\Users\\Callicious\\Documents\\Prog\\pycharm\\master-streams"
datadir = os.path.join(rootdir, "data")
sourcedir = os.path.join(rootdir,"source")
jl_dir = os.path.join(sourcedir, "jl_source")
imgdir = os.path.join(rootdir,"img")
asciidir = os.path.join(rootdir,"ascii")
duplimontedir = os.path.join(rootdir,"duplimonte")
clusterdir = os.path.join(rootdir,"clustered")
clusterdir_fine = os.path.join(rootdir,"clustered_finetune")
clusterdir_flat = os.path.join(rootdir,"clustered_flatfork")
orbitsdir = os.path.join(rootdir,"orbits")
orbitsfitdir = os.path.join(rootdir,"orbifitdir")
clusterer_flat_path = os.path.join(clusterdir, "flatfork_clusterer.txt")
clusterer_flat_labels = os.path.join(clusterdir, "flatfork_fullgroup_cluster.txt")
clusterer_fine_path = os.path.join(clusterdir, "finetune_clusterer.txt")
clusterer_fine_labels = os.path.join(clusterdir, "finetune_fullgroup_cluster.txt")

# LAMOST-related directories + paths + gaia stuff
lamodir = os.path.join(rootdir,"lamost_new")

gaia_distantmatch = os.path.join(lamodir, "distant_match.fits")
gaia_distantmatch_galcent = os.path.join(lamodir, "distant_match_galcent.fits")

# LAMOST-related directories + paths (Old Lamodir)
lamodir_old = os.path.join(rootdir,"lamost")

LRS_path = os.path.join(lamodir, "dr7_v2.0_LRS_stellar.fits")
MRS_path = os.path.join(lamodir, "dr7_v2.0_MRS_stellar.fits")

LRS_cleaned_path = os.path.join(lamodir, "dr7_v2.0_LRS_stellar_cleaned.fits")
MRS_cleaned_path = os.path.join(lamodir, "dr7_v2.0_MRS_stellar_cleaned.fits")

LRS_ids = os.path.join(lamodir, "LRS_ids.txt")
MRS_ids = os.path.join(lamodir, "MRS_ids.txt")

LRS_dir = os.path.join(lamodir, "LRS_dir")
MRS_dir = os.path.join(lamodir, "MRS_dir")

master_LRS_fits = os.path.join(lamodir, "master_LRS.fits")
master_MRS_fits = os.path.join(lamodir, "master_MRS.fits")

master_fits = os.path.join(lamodir, "LAMOST_master.fits")
master_fits_galcen = os.path.join(lamodir, "LAMOST_master_galcen.fits")
master_fits_galcen_ascii = os.path.join(lamodir, "LAMOST_master_galcen.txt")

baumgardt_fits = os.path.join(lamodir, "baumgardt_master.fits")