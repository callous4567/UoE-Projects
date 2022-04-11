# Gauss-Seidel with Successive Over Relaxation Poisson Solver, Electric Charge Distribution Edition, in Python with Numba
The model here follows the notes laid out by Davide in the Computer Modelling and Visualization course
here at the University of Edinburgh. It's not the general form that you find on Wikipedia, for example.

You'll find options for the standard Jacobian relaxation, Gauss-Seidel rendition of it, and SOR rendition of that. There's an option
to create random charge distributions (with normals and so forth for each charge) and so on. See the video
[here](https://www.youtube.com/watch?v=ksZ-GeKRyag) for an animation.

You can also create some neat plots to compare the simulation of a point charge to the actual analytical
distribution, i.e.

![alt text](https://github.com/callous4567/UoE-Projects/blob/master/SimAndVis/C3_2/example_fieldplot_better.png)



## fast_relax.py
Contains the Numba code for potential/meshgrids/so forth, alongside Gauss-Seidel/SOR.

## fast_electrorelax.py
Contains the Numba code for the electric field discretisation.

## hdfutils.py 
Contains various utilities for writing data and results: see TGP2020, master-streams, and SHons2021, for other usages.
This class uses HDF5 and h5py to save data in a user-interpretable way- nice and easy to use.

## relax.py 
Contains the main code for running the simulation (plotting, etc) that 
makes use of fast_relax,fast_electrorelax, and multiprocessing_functions.
 
## multiprocessing_functions 
Contains a quick code-snippet for running many SOR runs in parallel to estimate the ideal
value of omega.
