This is the sample readme for the Stellar Streams MSc Project
Student: Sebastian Straszak
Supervisor: Jorge Peñarrubia
It'll hold general working information for the code involved in this project.
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
¬¬¬¬¬¬
SOURCE
¬¬¬¬¬¬
=ascii_info.py
-various shortkeys for raw ascii file

=asciiutils.py
-various tools to deal with import of raw ascii files

=galacticutils.py
-galactic coordinate to ICRS vis-a-vis

=galcentricutils.py
-galactocentric coordinate to ICRS vis-a-vis
-galactocentric angular momentum, great circle selections, and other fantastical things

=hdfutils.py
-hdf5 file handler utility to make it easier to code with hdf5

=solar_info.dat
-heliocentric frame information relative to galaxy, right-handed galactocentric, Bennett-Body/Gravity/Drimmel-Poggio

=solar_info_jorge.dat
-heliocentric frame information but the one that jorge used in his code: see emails

¬¬¬¬
MAIN
¬¬¬¬
=windows_asciiport.py
-imports raw ascii data into a hdf5 table, concatenating it all while also porting individual lists

=windows_helioport.py
-t.b.d