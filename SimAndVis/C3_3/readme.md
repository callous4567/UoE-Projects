# Gauss-Seidel with Successive Over Relaxation Poisson Solver, Magnetic Potential Edition, in Python with Numba
The model here follows the notes laid out by Davide in the Computer Modelling and Visualization course
here at the University of Edinburgh. It's not the general form that you find on Wikipedia, for example.

You'll find options for the Gauss-Seidel rendition of it, and SOR rendition of that. There's an option
to create random z-parallel wire distributions (with normals and so forth for each wire) and so on. See the video
[here](https://www.youtube.com/watch?v=ksZ-GeKRyag) for an animation. *This is for simulating the magnetic field in the
x-y axis for a bunch of wires parallel to z, but you can generalise it for other distributions.*

There's a function to create random charge distributions- pretty fun to play with. You'll always get x-y projections with the plotting I've added for the live-view animation- here's an example snip :) 

![alt text](https://github.com/callous4567/UoE-Projects/blob/master/SimAndVis/C3_3/wireplot.PNG)


You can also create some neat plots to compare the simulation of a single wire to the actual analytical
distribution, i.e.

![alt text](https://github.com/callous4567/UoE-Projects/blob/master/SimAndVis/C3_3/example_wireplot.png)

There's also some curve-fitting tools to evaluate a comparison between the model numerically to the model analytically- this basically rests on the solution
for the A-field you get for a limited-length wire. 


## fast_relax.py
Contains the Numba code for potential/meshgrids/so forth, alongside Gauss-Seidel/SOR.

## fast_magnetorelax.py
Contains the Numba code for the magnetic field discretisation.

## hdfutils.py 
Contains various utilities for writing data and results: see TGP2020, master-streams, and SHons2021, for other usages.
This class uses HDF5 and h5py to save data in a user-interpretable way- nice and easy to use.

## relax.py 
Contains the main code for running the simulation (plotting, etc) that 
makes use of fast_relax,fast_magnetorelax, and multiprocessing_functions.
 
