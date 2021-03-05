That big mess at the bottom of this document is the old readme. It has some info for how you can get this
code working on your desktop. It's old + hasn't been updated since christmas (Yeah, a long time.) 

astrowrapper.py has various utilities to actually process data. I won't go into the detail, but to sum it up 

class: utils: lots of utilities used throughout, with virtually all atmospheric calibration tools
class: fitser: image alignment routines. See TGP report for information. Will attach later post-grading to avoid turnitin.
class: hdf5_writer: convenience class for handling hdf5 format. 
class: turnoff_analysis: least-squares fitting of isochrones + other such nonsense: HR diagrams + colour-colour + etc
class: misc_file_handler: handles random file imports/etc, i.e. vizier catalogues and other nonsense 
class: cmd3: import PARSEC/COLIBRI isochrone files (given directory or individual) to make one big catalogue used by turnoff_analysis 
class: GAIA: download GAIA DR2 catalogues and a crossmatch query tool (incomplete, need to finish. Now understand queries but have to add.) 
class: source_analysis: field photometry tool tailored to our fields. HR/Colour-Colour/Etc included.
class: member_finder: KDE for estimating cluster centroid + proper motions (error tbd. + also estimate for radius tbd.)

main_runtime.py actually leverages astrowrapper for the data structure we had. 
Aside from manually moving final deep images to a new directory for field photometry + onward (wasn't coded up, lazy.) all file handling is automatic.
As this all works recursively for the given file structure, care is needed in setting it up, if you intend to actually run this code. See below for guidance. 

TO RUN THIS CODE IN ITS COMPLETE FORM YOU WILL NEED A COSMOS LOG-IN TO GAIA DATA 
YOU WILL NEED AN ASTROMETRY.NET API KEY
AND YOU WILL NEED TO MANUALLY ALTER THE NUMBER OF CORES AVAILABLE FOR MULTIPROCESSING (or remove core limiter) 

test2.py and testingfile.py are just a pair of misc files used for testing random stuff/etc. They do not hold any bearing beyond generating some pretty figures.

############################################################################################################
This is the initial readme for the data processing pipeline that has been created for TGP2020 project.

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