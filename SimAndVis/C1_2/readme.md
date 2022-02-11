# 2D XY Ising Model, in Python with Numba
See [here](https://www.youtube.com/watch?v=5ZzX12x073k) to get a video of it in action, demonstrating formation of low-temperature vortices (indicative of the KT transition, with a first-order gaussian-smoothed divergence colourmap. Click [here](https://www.youtube.com/watch?v=8lncCFxXaWM) for a VERY long run, showing collapse of vortices. 

Critical slowdown is a problem: I have included an option to use the Wolff Algorithm to get past this. Vortices do not form here, however (unsure if this is correct: my point is, this is ungraded work I did for fun, and thus I do not care if it's correct or not.) 

See [here](https://www.youtube.com/watch?v=r5w7IQqBi7Q) to get a video of it in action with vortices illustrated with a phase colourmap (0,2pi.) 

## fast_xy.py
Contains the Numba code. 

## hdfutils.py 
Contains various utilities for writing data and results: see TGP2020, master-streams, and SHons2021, for other usages.
This class uses HDF5 and h5py to save data in a user-interpretable way- nice and easy to use.

## xy_ising.py 

Contains the code. Minimal working example including fast_xy, used purely to generate the videos above. You can take this, modify it, and make it a full piece like with twod_xy. 
