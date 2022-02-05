# 2D Ising Model, in Python with Numba
**Which I assume you know about if you're here, so, let's get to it!**

If you want to try this code out, download "Example_Prepared_Code.bz" and try running it via
```
python3 twod_ising.py 
```
after you've extracted it. Keep an eye on the console to answer some questions it'll ask you, defining the simulation run. _This won't produce a save! Only do a nice animation for the run!_
## fast_ising.py
Contains the Numba code. This will allow you to do fast energy calculations
alongside nice and fast Glauber/Kawasaki steps. 

_Note that there there is unexpected behaviour with the method "fast_ising.kawaglauber" which will
not function as intended- the reason is Numba-related and a request will be filed eventually._

## hdfutils.py 
Contains various utilities for writing data and results: see TGP2020, master-streams, and SHons2021, for other usages.
This class uses HDF5 and h5py to save data in a user-interpretable way- nice and easy to use.

## multirun.py 
Example file showing how to leverage multiprocessing to generate results for the Ising model, i.e. to 
produce graphs of averages against temperature, etc.

## twod_ising.py 

Main file to use, which contains the bulk of the code. This runs the model, handles the model, graphs the model,
and in general is the model. 

If you want to do an example run, take a look at the class **checkpoint** within this, specifically the method **run**.


## windows_multiprocess_functions.py 

Has the functions that multirun.py needs to correctly run. See this, if looking at multirun.
