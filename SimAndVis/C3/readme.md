# Cahn-Hilliard Phase Separation, in Python with Numba
The model here follows the notes laid out by Davide in the Computer Modelling and Visualization course
here at the University of Edinburgh. It's not the general form that you find on Wikipedia, for example.

In any case, the phenomenology and such is the same. You can get spinodal decomposition, globule formation, and so forth- i.e.

![alt text](https://github.com/callous4567/UoE-Projects/blob/master/SimAndVis/C3/bridges.PNG)

See [here](https://www.youtube.com/watch?v=ksZ-GeKRyag) for a video showing the code in action. The formation of globules also occurs if you alter the overall value of the composite order field- this is included with a whole sleuth of parameters to modify the model on run. 

## fast_cahn.py
Contains the Numba code. 

## hdfutils.py 
Contains various utilities for writing data and results: see TGP2020, master-streams, and SHons2021, for other usages.
This class uses HDF5 and h5py to save data in a user-interpretable way- nice and easy to use.

## cahn.py 
Contains the main code for running the simulation (plotting, etc) that makes use of fast_cahn.
 
