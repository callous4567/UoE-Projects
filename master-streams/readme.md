# Frequentist characterization of structure in angular momentum space with hierarchical clustering

This is the GitHub repository for my (poorly named, I know- I don't know what I was thinking!) Masters Thesis.
You'll find all the code I used here. I'll do my best to clean it up a bit *in the future* and perhaps add some documentation,
but I make no promises. I'm a somewhat busy guy, after all! *wink* 

*Note* This repository is still in use- any undocumented material is well... undocumented! :D 

## What the hell does that title even mean?
The title refers to the work taken in the thesis. This broadly fits into the following steps-
- Take kinematic data, process it a bit to be for the outer halo, 
- Cluster data hierarchically in angular momentum space, 
- Characterise clusters found, based mainly on orbit integrations and their appearance on the sky,

Naturally I did more- I always do- but adding "And studying the Minimum Entropy Method in practice" 
to the title would have been a stretch, as would have been making assumptions on anything else. It's a crisp,
clean, neat title, and it says what I did in a way that any decent dynamicist would understand it.

## Neat. Did you find anything out?

56.5% of the dataset were notable substructures. 6 were known. 2 statistically significant ones that we found
were actually brand-spanking new ones. If you don't consider astronomical *cough overestimated* error then
that 2 increases to 6, if you consider their orbits. The Minimum Entropy Method ends up plagued by degeneracy,
and thus only works in practice one-dimensionally. Very good one-dimensionally though!!! 

Overall I'd qualify that as a success. There are a lot of improvements that could be effected- hopefully I'll get
to do these in the future. The thesis covers these- I won't publish that here until the far future, though.

## Summary of the codebase available here 

Here's a brief summary of all the various bits of code. Now, I have to admit, while this was "frequentist"
I did end up making a lot of code in the hope of using Bayesian Gaussian Mixture Models (GMMS.) So there may be a few 
snips around that are useful for this. I think I even made a thing to generate random mixtures to test my GMMs on 

### deprecated 
Has old bits of code not made use of in final work. This includes an (at the time) full implementation of GCCs as made use of by
[Charles King III and pals](https://iopscience.iop.org/article/10.1088/0004-637X/750/1/81) for pole counts, alongside implementations
of Jo Bovy's [XD](https://github.com/jobovy/extreme-deconvolution) package to try and use GMMs to get substructure in L-space. There's generally just lots and lots
of random code snips.

### main 

Makes use of base classes and methods to effect something greater than the sum of its parts. Basically, get the results
using the methods written. Toward the end of the project, I was hasty- consequently, there are some full-blown pieces
of code in here- by and large though it's just multiprocessing implementations of the source material to allow for Monte Carlo
error propagation and that sort of thing.

### source 

Contains all the individual base classes and methods used. There is *a lot* here. This includes *original* methods for...

- Evaluating potentials with Numba as used in the thesis
- Minimum Entropy Method implementation in 3D with quasi-MCMC, Grid Search, and Monte Carlo grid search
- Numerical 3D probability density function generation for arrays via rectangular binning, including Gauss kernel smoothing functionality
- Least-squares-powered orbit fitting, making use of [galpy](https://github.com/jobovy/galpy) for orbit generation
- An overabundance of graphing tools, including for the replication of plots from [Sarah Sofie's thesis](https://fse.studenttheses.ub.rug.nl/24089/) for examining a comparison against the canonical integrals of motion
- The most recent rendition of hdfutils, my class that allows easy use of [h5py](https://www.h5py.org/)
- Various convenience tools for coordinate transformations of [astropy](https://docs.astropy.org/en/stable/index.html) tables astronomically and mathematically
- Implementation of solving the rectangular assignment problem specific to this thesis, and for combination of many frequency tables on a per-label basis
- Tools for euler-angle rotations of tables in the galactocentric frame
- Tools for generating quasi-random Gaussian mixtures with noise 
- Methods for effecting Great Circle Cell counts as [here](https://ui.adsabs.harvard.edu/abs/1996ASPC...92..483J/abstract) with further refinement to allow azimuthal splitting of the GCCs, alongside numerical generation of GCC greatcircles
- Least-squares-powered fitting of GCCs to clusters, with weighting functionality
- Monte Carlo propagation of error from the Galactic to Galactocentric frame in various parameters

and I'm sure there's some more mixed in there. As you can tell, there's lots! :D Hopefully other students who inevitably
find their way here will enjoy it thoroughly.

### testing 

As deadlines closed in, I got lazy and sloppy. I fell from trees, got depressed, and drank- mainly lots of cider. This contains code that was used to dirtily produce plots and last-minute
results. This includes some code for visualizing degeneracy of the Minimum Entropy Method- other than that? Hah.

## Miscellaneous words from myself 

Hopefully you find this repository useful. If you're another student following this kind of field, good luck! 
I hadn't a clue about stellar streams or anything else, and frankly I'm garbage with statistics. Nonetheless, I had fun,
and hopefully my code will help you with whatever you are trying to do. 

# master_streams_new and more... 

Lots of recent updates will result in a new write-up t.b.d. In the mean time, some crossmatch routines, catalogue reduction, etc
and that sort of thing. Will update later. 