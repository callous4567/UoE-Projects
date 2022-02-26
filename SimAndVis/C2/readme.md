# 2D Game of Life Cellular Automata, in Python with Numba
See my video post [here](https://www.youtube.com/watch?v=nGk46hcFbZc) to see it in action on a 400^2 grid. To run this
code, simply download it and execute "twod_gol." This will run the "Checkpoint" produced for my submission to MVis 
here at UoE. Note that the file structures involved are designed for Windows, and will not function as intended on Unix
or any other OS.

## fast_gol.py
Contains the Numba code. Also includes a method for calculating the centre of mass of Gliders.

## hdfutils.py 
Contains various utilities for writing data and results: see TGP2020, master-streams, and SHons2021, for other usages.
This class uses HDF5 and h5py to save data in a user-interpretable way- nice and easy to use.

## multirun.py 
Example file showing how to leverage multiprocessing to generate many soups simultaneously for "activity histogram"
used to estimate mean equilibration times for soups of grid size N^2. 

## twod_ising.py 

Main file to use, which contains the bulk of the code. This runs the model, handles the model, graphs the model,
and in general is the model. 

## windows_multiprocess_functions.py 

Has the functions that multirun.py needs to correctly run. See this, if looking at multirun.
