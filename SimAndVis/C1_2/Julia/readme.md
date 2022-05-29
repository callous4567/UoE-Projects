# 2D XY Ising Model, in Julia with Plots 
See [here](https://www.youtube.com/watch?v=34T5GDJ-SoM&t=1s) to get a video of it in action, all the way to closing of vortices.

Note that while one Python step takes 6000 nanoseconds, with Julia this is down to 300 nanoseconds. If you switch to explicitly
2D arrays (instead of 3D) and use angle coordinates, this decreases to 100 nanoseconds (i.e. with contiguous flattened arrays.)

There are many files here: see oned_test for the final working implementation used for the video example. The other
examples included include the standard implementation (3D not 2D) and non-contiguous 2D. 