Hi! Thanks for opening. 

"Prog_Files" contains the various program files and config files. 

orbit_solver.py is to be used when... y'know, solving orbits. Just run it and it'll solve all the orbital data it has.
Currently optimal for 0.1 day timesteps for 1 million total steps (just above the orbital period of Pluto)
That should take a few minutes to run. Anything less (0.01) and it'll take 4+ hours, no idea about that scaling.

NBody.py is to be ran when you want to generate the datafiles and VMD output files.

parameters.dat holds TIMESTEP, NUMBER OF STEPS, GRAVICONSTANT *I'm unsure about that last one. I didn't write it ;-;
particles.dat holds the particle information/ephemerides/etc.

In terms of output files, consider Jupiter as an example.
"Jupiter.dat" is the raw X Y Z VX VY VZ output.
"Jupiter.dat-Solved-Parametes-.txt" has the solved parameters. The units are (obviously) Astronomical in nature.
(They're astronomical units, for Periapsis, Apoapsis, and Semimajor/Semiminor. Eccentricity is dimensionless)

When solving for the moon, "Luna.dat" and "Luna.dat-Solved-...etc" will have the parameters for the moons orbit.

"Example Data" contains the produced data. Timestep of 0.1 days for 1 million steps. 0.01 days is also doable, though I'll save that for the report since it might be fun to talk about later. Plus it's still generating the text files, currently over a gigabyte in size (go figure.) 

Anyway, I hope this was instructive, and thanks for reading it if you did decide to read it!
If not, then it's a shame since we'll be demerited for not labelling the units in the output. But hey, life is life.
