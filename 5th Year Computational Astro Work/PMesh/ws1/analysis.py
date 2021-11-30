import numpy as np
import yt
from yt.units import *
import astropy.units as u
from dm_overdensity import _overdensity
import matplotlib.pyplot as plt
import yt.units as u

"""
Note that results is quite hefty- don't spam too many results up...
KEEP ALL STEPS IN ORDER and for the love of god almighty, Label. Each. Step. WELL. 
"""

# This has the working directory for all (directly) related Enzo stuff.
resdir = "/home/callicious/Documents/pycharm/CompAstro/pmesh/sup-pm/results"
supdir = "/home/callicious/Documents/pycharm/CompAstro/pmesh/sup-pm"

# WSdir for this workshop (i.e. save images here.)
wsdir = "/home/callicious/Documents/pycharm/CompAstro/pmesh/ws1"

# TODO: FIRST STEPS for Q1.
"""
generate inits
inits.exe -d dmonly.inits

remove
DataDumpName = DD
DataDumpDir = DD
dtDataDump = 2.0
from dmonly_unigrid.enzo

change globaldir to appropriate dump directory

run enzo
enzo.exe -d dmonly_unigrid.enzo > dmonly_unigrid.log &

RDNNNN files have data for certain redshift
amr.out has information on output files, i.e. 
CosmologyOutputRedshift[5] = 1 means RD0005 corresponds to 
z = 1, 
a = 1/ (1 + z) = 0.5
"""

# TODO: Do the spectrum plots for Q2.
# Spectrum names
#spectrum0 = supdir + "/PowerSpectrum_z=0.out"
#spectrum60 = supdir + "/PowerSpectrum_z=60.out"
# Reads in the PowerSpectrum file and generates a basic plot showing the power spectrum
"""
def plotspectrum(powerspectrum0, powerspectrum1):
    with open(powerspectrum0, 'r') as f:
        lines = f.readlines()
    # Clean up lines and split
    for num, line in enumerate(lines):
        line = line.replace('\n', "")
        line = line.split(" ")
        line = [x for x in line if x]
        lines[num] = [float(d) for d in line]
    data0 = np.array(lines).T

    with open(powerspectrum1, 'r') as f:
        lines = f.readlines()
    # Clean up lines and split
    for num, line in enumerate(lines):
        line = line.replace('\n', "")
        line = line.split(" ")
        line = [x for x in line if x]
        lines[num] = [float(d) for d in line]
    data1 = np.array(lines).T

    plt.figure()
    plt.plot(np.log10(data1[0]), np.log10(data1[1]), color="red", label="z=60")
    plt.plot(np.log10(data0[0]), np.log10(data0[1]), color="blue", label="z=0")
    plt.xlabel(r'$\log{k}$')
    plt.ylabel(r'$\log{I(k)}$')
    plt.legend()
    plt.suptitle("Power Spectrum z=0 vs z=60")
    plt.savefig("WS1comparespectrum.png", dpi=300)
    plt.show()

    plt.clf()
    plt.figure()
    ratio = np.log10(data0[1]/data1[1])
    plt.plot(np.log10(data0[0]), ratio, color="red")
    plt.xlabel(r'$\log{k}$')
    plt.ylabel(r'$\log$' +  " " + r'$\frac{I(k)z=0}{I(k)z=60}$')
    plt.suptitle("Power Spectrum z=0 vs z=60, ratio, logarithmic")
    plt.savefig("WS1ratiospectrum.png", dpi=300)
    plt.show()
"""
# Make the spectrum plots we want.
#plotspectrum(spectrum0, spectrum60)

# TODO: Do the steps for Q3.
"""
ds0 = yt.load(resdir + "/RD0000/RD0000")
ds0.add_field(name=("gas",'dm_overdensity'),
              function=_overdensity,
              sampling_type="cell")
print(ds0.time_unit.to(u.yr))
print(ds0.time_unit)
"""

# TODO: Steps for Q4
#ds0 = yt.load(resdir + "/RD0000/RD0000")
#print(ds0.derived_field_list)
#print(ds0.hubble_constant)
#ds0.add_field(name=("gas",'dm_overdensity'),
#              function=_overdensity,
#              sampling_type="cell")

# Slice plot centred on z=0.5, x,y,z = 0.5, with a width of 20 kpc. The units of x,y,z are in box normalized units.
#slc0 = yt.SlicePlot(ds0, 'z', 'all_cic', center=[0.5,0.5,0.5],width= (20,'kpc'))
#slc0.save('densitysliceplot20kpc.png')

"""
# Obtain the grid zone size (Q4)
redshift = 60
boxsize_init = 20 # Mpc/h
h = 0.7 # dimensionless hub
boxsize_init = boxsize_init / h # Mpc
new_a = (1 + redshift)**(-1) # via redshift relation
#h = 0.7
#h = H0/(100 km(sMpc)^-1
#H0 = 70 km(sMpc)^-1
new_boxsize = boxsize_init*new_a
per_cell = new_boxsize/128
print(per_cell)
width_mpc_of_slice = 20*1e-3
cells_for_slice = width_mpc_of_slice/per_cell # divide mpc of slice by mpc per cell for cells for slice
print(cells_for_slice)

# For the entire simulation volume
#slc0 = yt.SlicePlot(ds0, 'z', 'all_cic', center=[0.5,0.5,0.5],width= (1))
#slc0.save('densislicefull.png')

"""

# TODO: Steps for Q5
"""
ds0 = yt.load(resdir + "/RD0000/RD0000")
ds0.add_field(name=("gas",'dm_overdensity'),
              function=_overdensity,
              sampling_type="cell")
# Projection plot instead of a slice, 20 kpc as before.
#prj0 = yt.ProjectionPlot(ds0, 'z', 'all_cic', center=[0.5,0.5,0.5],width= (20,'kpc'))
#prj0.save('densityprojectionplot20kpc.png') """


# TODO: Steps for Q6
"""
ds0 = yt.load(resdir + "/RD0000/RD0000")
ds0.add_field(name=("gas",'dm_overdensity'),
              function=_overdensity,
              sampling_type="cell")
# For the entire simulation volume
slc0 = yt.SlicePlot(ds0, 'z', 'all_cic', center=[0.5,0.5,0.5],width= (1))
slc0.save('densislicefull.png')
# For the entire simulation volume- shifted along x
slc0 = yt.SlicePlot(ds0, 'z', 'all_cic', center=[0.6,0.5,0.5],width= (1))
slc0.save('densislicefull_shifted.png') """


# TODO: Steps for Q7
"""
# Z = 60
ds0 = yt.load(resdir + "/RD0000/RD0000")
ds0.add_field(name=("gas",'dm_overdensity'),
              function=_overdensity,
              sampling_type="cell")
slc0 = yt.SlicePlot(ds0, 'z', 'dm_overdensity', center=[0.5,0.5,0.5],width= (1))
yt.AxisAlignedSlicePlot.set_zlim(slc0,field='dm_overdensity',zmin=0.0001,zmax=200)
slc0.save('z60fullsliceoverdensity.png')


# Z = 4
ds1 = yt.load(resdir + "/RD0002/RD0002")
ds1.add_field(name=("gas",'dm_overdensity'),
              function=_overdensity,
              sampling_type="cell")
slc1 = yt.SlicePlot(ds1, 'z', 'dm_overdensity', center=[0.5,0.5,0.5],width= (1))
yt.AxisAlignedSlicePlot.set_zlim(slc1,field='dm_overdensity',zmin=0.0001,zmax=200)
slc1.save('z4fullsliceoverdensity.png')

# Z = 1
ds2 = yt.load(resdir + "/RD0005/RD0005")
ds2.add_field(name=("gas",'dm_overdensity'),
              function=_overdensity,
              sampling_type="cell")
slc2 = yt.SlicePlot(ds2, 'z', 'dm_overdensity', center=[0.5,0.5,0.5],width= (1))
yt.AxisAlignedSlicePlot.set_zlim(slc2,field='dm_overdensity',zmin=0.0001,zmax=200)
slc2.save('z1fullsliceoverdensity.png')

# Z = 0
ds3 = yt.load(resdir + "/RD0009/RD0009")
ds3.add_field(name=("gas",'dm_overdensity'),
              function=_overdensity,
              sampling_type="cell")
slc3 = yt.SlicePlot(ds3, 'z', 'dm_overdensity', center=[0.5,0.5,0.5],width= (1))
yt.AxisAlignedSlicePlot.set_zlim(slc3,field='dm_overdensity',zmin=0.0001,zmax=200)
slc3.save('z0fullsliceoverdensity.png') """

# TODO: Steps for Q8
"""
# Z = 1
ds2 = yt.load(resdir + "/RD0005/RD0005")
ds2.add_field(name=("gas",'dm_overdensity'),
              function=_overdensity,
              sampling_type="cell")

# Obtain the grid zone size (Q4)
redshift = 1
boxsize_init = 20 # Mpc/h
h = 0.7 # dimensionless hub
boxsize_init = boxsize_init / h # Mpc
new_a = (1 + redshift)**(-1) # via redshift relation
new_boxsize = boxsize_init*new_a
per_cell = new_boxsize/128 # in Mpc
eightcells = 8 * per_cell
print(eightcells)

# z=1 8-cell-slice at centre.
slc2 = yt.SlicePlot(ds2, 'z', 'dm_overdensity', center=[0.5,0.5,0.5],width=(eightcells*u.Mpc))
slc2.save('z1fullsliceoverdensity_8by8.png')

# Set depth = 0.
ppl2=yt.ParticleProjectionPlot(ds2,'z', "particle_mass", center=[0.5,0.5,0.5],width=(8/128),depth=2*per_cell*u.Mpc)
ppl2.save('z1fullsliceoverdensity_8by8particle.png')
"""