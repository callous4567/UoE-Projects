# Modelling and Visualization Directory
Poorly labelled (since I mistook it as "Simulation and Visualization! lol. Here you'll find the work I produce(d) for the three checkpoints afforded by this course. 

## C1 
The two-dimensional Ising model, with both Glauber and Kawasaki Dynamics, sped up by Numba, with various graphing (and video-producing) utilities. Video examples are linked. **Not graded or required by course** Wolff Algorithm module for Glauber is included. 

## C1_2
**Not graded or required by course**, to simulate the 2D XY Model (which is a vector-version of the 2D Ising Film.) Video examples linked. Quiver plotting with first-order gaussian-convolved divergence maps demonstrating vorticity. Phase plotting (0,2pi) demonstrating transitions. etc. Wolf Algorithm implementation is present, tested, and produces expected results (seemingly.) 

## C2 
Conway's Game of Life, with a few other miscellaneous metric tools (like calculating the speeds of gliders, or for example, estimating the distribution of lifetimes for soups randomly generated at initialization.

## C2_2
The SIRS Model of Epidemics, albeit not "true SIRS" since it only uses four nearest neighbours instead of eight (though implementing this would only necessitate a few line changes in fast_sirs.) 

## C3
The Cahn-Hilliard phase-separation model in 2D, powered by numba. Spinodal decomposition and the like are seen readily. There's options for free-energy evaluation.

## C3_2
Relaxation solver for the poisson equation- specifically modelled here to allow representation of electric fields and the like for fields of point charges (or arbitrary charge distributions, technically.) Jacobian relaxation, Gauss-Seidel, and Successively-Over-Relaxed implementation. Comparison tools/metrics for the analytical form are also contained. There's also options for staggered grids or distance-relaxed grids when calculating potentials/fields numerically from density distributions.

## C3_3
Same as C3_2, except in this case, we're considering a field of z-orientated wires, and modelling the z-axis magnetic potential (since the x-y components vanish) and subsequently calculating the x-y axis magnetic field. Basically identical except for that change. I've altered the boundary conditions, also, to make convergence a bit faster (given that we expect every x-y plane to be identical this is justified.)  
