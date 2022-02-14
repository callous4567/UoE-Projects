# 2D XY Ising Model, in Python with Numba
See [here](https://www.youtube.com/watch?v=5ZzX12x073k) to get a video of it in action, demonstrating formation of low-temperature vortices (indicative of the KT transition, with a first-order gaussian-smoothed divergence colourmap.) Click [here](https://www.youtube.com/watch?v=8lncCFxXaWM) for a VERY long run, showing collapse of vortices. See [here](https://www.youtube.com/watch?v=r5w7IQqBi7Q) to get a video of it in action with vortices illustrated with a phase colourmap (0,2pi.) 

Critical slowdown is a problem: I have included an option to use the Wolff Algorithm to get past this. Vortices do not form here, however (duh- we avoid critical slowdown.) 

**Here's an example plot that you can get using the Wolff implementation I've put in (just demonstrating critical temperature.)** Since this matches the value you find on Wikipedia [here](https://en.wikipedia.org/wiki/Classical_XY_model#/media/File:XY_SpecificHeat.svg) I am guessing that the code I put in works somewhat. I just felt like I needed at least a pseudocode that worked along the same principle before I stopped coding for this little project. Anyhow, with that, this work is complete. [here](https://www.youtube.com/watch?v=jbjJqt6ZBK8) is a video of the Wolff algorithm in action.

![alt text](https://github.com/callous4567/UoE-Projects/blob/master/SimAndVis/C1_2/g_multi.png)



## fast_xy.py
Contains the Numba code. 

## hdfutils.py 
Contains various utilities for writing data and results: see TGP2020, master-streams, and SHons2021, for other usages.
This class uses HDF5 and h5py to save data in a user-interpretable way- nice and easy to use.

## xy_ising.py 

Contains the code. Minimal working example including fast_xy, used purely to generate the videos above. You can take this, modify it, and make it a full piece like with twod_xy. 
