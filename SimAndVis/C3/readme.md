# Cahn-Hilliard Phase Separation, in Python with Numba
The model here follows the notes laid out by Davide in the Computer Modelling and Visualization course
here at the University of Edinburgh. It's not the general form that you find on Wikipedia, for example.

In any case, the phenomenology and such is the same. You can get spinodal decomposition, globule formation, and so forth- i.e.

![alt text](https://github.com/callous4567/UoE-Projects/blob/master/SimAndVis/C3/bridges.PNG)

See [here](https://www.youtube.com/watch?v=ksZ-GeKRyag) for a video showing the code in action. The formation of globules also occurs if you alter the overall value of the composite order field- this is included with a whole sleuth of parameters to modify the model on run. 

## fast_cahn.py
Contains the Numba code. I've also included some non-numba code for the simulation run, to illustrate how much faster using numba for loops is than using np.roll explicitly with python (I'm positive you could work out np.roll with flattened arrays for this, but numba *does not* support 2D or higher rolling- you're thus limited in that regard, with for loops providing a decent option (faster than default non-compiled python) such as here 

![alt text](https://github.com/callous4567/UoE-Projects/blob/master/SimAndVis/C3/fast_test.png)

You can clearly see that as grid size increases, the ratio of numpy to numba generally levels off, somewhere around a factor of 8-9. By all accounts, in this regime (single-threaded grid sizes upwards of 1000+- I've tested these manually-) for loops with numba are the way to go. 

## hdfutils.py 
Contains various utilities for writing data and results: see TGP2020, master-streams, and SHons2021, for other usages.
This class uses HDF5 and h5py to save data in a user-interpretable way- nice and easy to use.

## cahn.py 
Contains the main code for running the simulation (plotting, etc) that makes use of fast_cahn.
 
