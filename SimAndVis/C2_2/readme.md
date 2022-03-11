# 2D SIRS Model, in Python with Numba
Note: Unlike the typical SIRS model, I was instructed to use the Ising near-neighbours (so only 4, not 8.) To my awareness,
this is atypical, but I wouldn't know- I'm no epidemiologist. There is funtionality included to allow for a fraction of souls to remain
permanently immunized, and so forth, including functionality for sequential and parallel update schemes. 

**I have included my checkpoint submission for this particular project. Enjoy.**

## fast_sirs.py
Contains the Numba code. 

## hdfutils.py 
Contains various utilities for writing data and results: see TGP2020, master-streams, and SHons2021, for other usages.
This class uses HDF5 and h5py to save data in a user-interpretable way- nice and easy to use.

## multirun.py 
Example file showing how to leverage multiprocessing to generate many soups simultaneously for "activity histogram"
used to estimate mean equilibration times for soups of grid size N^2. 

## twod_sirs.py 

Main file to use, which contains the bulk of the code. This runs the model, handles the model, graphs the model,
and in general is the model. 

## windows_multiprocess_functions.py 

Has the functions that multirun.py needs to correctly run. See this, if looking at multirun.
