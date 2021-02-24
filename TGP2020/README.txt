############################################################################################################
This is the readme for the data processing pipeline that has been created for TGP2020 project.

If you intend to reuse any significant parts of this code without alteration,
it'd be appreciated if you'd reference the author (Callicious)

Any and all parts of this code (save for a few snips from the Astropy Wiki for Astrometry looping)
were written by Callicious, unless otherwise stated.
############################################################################################################


Generic notes on the code, how to run, and the two main files.
============================================================================================================
astrowrapper.py
------------------------------------------------------------------------------------------------------------
fits_alignment (CLASS)
Holds various utilities for FITS alignment (alignment of FITS files based on sources within a few pixels via
distance minimization, a few hundred pixels based on WCS, rotations, and generic applications of these on
directories/etc in line with the file structures that have been listed below) that may be adapted to work
with a given file  structure/etc.

hdf5_writer (CLASS)
Utility for writing/reading arrays from a predefined HDF5 file within groups + datasets. Overwrites old data
and reads in new data. Both fits_alignment and utils are dependent on this class.

utils (CLASS)
The main data processing module for FITS files that allows all FITS file handling (except for rotation/aligns)
and also all extra analysis of files from photometry to extraction of 2D gaussian intensity plots/etc. fits_
alignment is dependent on this class.

source_analysis (CLASS)
Holds various tools for analysing sources and that kind of thing (photometry) inc. Gliese catalogue text file.

misc_votable_hanlder (CLASS)
Tool for importing GAIA catalogue votables from DR2 and also for import of Hipparcos VIZIER data


############################################################################################################
FILE STRUCTURE
############################################################################################################

Root Directory (Set in astrowrapper for use. main_runtime will rip from astrowrapper.
------------------------------------------------------------------------------------------------------------
rootdir = root directory, i.e. D:\TGP2020


Base Data: Cluster Example, M52 : Offset Example, 0 (INTEGER)
------------------------------------------------------------------------------------------------------------
rootdir -> Sci -> Clusters -> Offsets -> Data
rootdir -> Cali -> U,B,V,BIAS,DARKS -> Data


Calibrat(ed)ion Data
------------------------------------------------------------------------------------------------------------
rootdir -> Sci -> Clusters -> Offsets -> CALIBRATED -> Data
rootdir -> Cali -> U,B,V,BIAS,DARKS -> maserflatU,B,V, masterbias, masterdark


Align(ed) Offset Data
------------------------------------------------------------------------------------------------------------
rootdir -> Sci -> Clusters -> Offsets -> CALIBRATED -> Data


Stack(ed) Offset Data
------------------------------------------------------------------------------------------------------------
rootdir -> Sci -> Clusters -> Offsets -> CALIBRATED -> fully_stacked_CLUSTER_OFFSET_BAND


Stack(ed) Offset Placeholder Folder
------------------------------------------------------------------------------------------------------------
rootdir -> Sci -> Clusterstacks -> Clusters -> U,B,V -> fully_stacked_CLUSTER_OFFSET_BAND


Final Stacked Aligned Folder Structure/etc
------------------------------------------------------------------------------------------------------------
rootdir -> Sci -> ClusterstacksDone -> Clusters -> U,B,V -> WCSDONE (1st), fully_stacked_CLUSTER_OFFSET_BAND
WCSDONE (1st)  -> FinalAligns, WCSDONE (2nd), fully_stacked_CLUSTER_OFFSET_BAND.fitswcs ====================

FinalAligns -> fully_stacked_CLUSTER_OFFSET_BAND.fitswcs.fitsrot_trans_aligned.fits.fitstrimmed.fitsfixed.fi
ts.fitsfitsfitsaligned21 (Final aligned offsets), fully_stacked_images (name unknown)

WCSDONE (2nd) -> Placeholder for WCS of offets used by fits_rot_aligner (holds WCS-fitted offset (1,2) before
os.replace moves it back to WCSDONE (1st) for alignment relative to 0.