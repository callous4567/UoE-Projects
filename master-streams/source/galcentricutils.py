from sympy import *
import math
import hdfutils
import os
import numpy as np
import astropy.coordinates as coord
from astropy.coordinates import galactocentric_frame_defaults
import astropy.units as u

# SOME NOTES
"""
- We're using right-handed galactocentric coordinates
- The spherical representation uses declination-based polar angle: not typical mathematical polar angle. 
"""

# Convert from ICRS to Galactocentric and back with custom solar system parameters (set via .dat file of right format.)
class galconversion(object):
    def __init__(self, sourcedir):
        self.source, self.owd, self.sol_params = sourcedir, os.getcwd(), np.zeros((1,4))

    # Grab hold of the parameters of the solar system from source_solar and SETS THESE FOR ASTROPY.
    # Note that you have to provide source_solar: to deal with conflicting sources. See samples.
    def solinfo_grab(self, source_solar):
        os.chdir(self.source)
        solinfo = open(source_solar).read().splitlines()
        os.chdir(self.owd)
        solinfo_splits = [d.split(";") for d in solinfo]
        sol_params = [False,False,False,False]
        for split in solinfo_splits:
            if split[0] == "sol_pos":
                sol_params[0] = split[1].strip().split(",")
            if split[0] == "sol_pos_err":
                sol_params[1] = split[1].strip().split(",")
            if split[0] == "sol_vel":
                sol_params[2] = split[1].strip().split(",")
            if split[0] == "sol_vel_err":
                sol_params[3] = split[1].strip().split(",")
        sol_params = np.array(sol_params, dtype=float)
        self.sol_params = sol_params

    # Set up/define new galactocentric system with solar parameters: new system called "master-streams"
    def solgal_set(self):
        galstate = galactocentric_frame_defaults.get_from_registry('latest')
        galstate["parameters"]["galcen_distance"] = np.abs(self.sol_params[0][0]) * u.kpc
        galstate["parameters"]["galcen_v_sun"] = self.sol_params[2] * u.km / u.s
        galstate["parameters"]["z_sun"] = self.sol_params[0][2] * u.kpc
        galactocentric_frame_defaults.register(name="master-streams", **galstate)
        galactocentric_frame_defaults.set("master-streams")
        # print(Galactocentric())  (TO SEE FRAME PROPERTIES!!!)

    # Work an astropy table under HDFDIR - HDFFILE - GROUP - SET and add ICRS columns to it, following GAIA formalism
    # Assumes Galactocentric are under x,y,z,vx,vy,vz
    # Does NOT do error conversions- see later.
    def to_ICRS(self, hdfdir, hdfname, group, set):
        # Set up HDF and grab table, and SkyCoord objects for all targets.
        writer = hdfutils.hdf5_writer(hdfdir, hdfname)
        table = writer.read_table(group, set)
        skycoords = coord.SkyCoord(x=table['x'] * u.kpc, y=table['y'] * u.kpc, z=table['z'] * u.kpc,
                                   v_x=table['vx'] * u.km/u.s, v_y=table['vy']*u.km/u.s, v_z=table['vz']*u.km/u.s,
                                   frame="galactocentric")
        # Effect conversion to ICRS, work through objects, collect converted quantities. 
        icrs_skycoords = skycoords.transform_to(coord.ICRS)
        ra_list, dec_list, pmra_list, pmdec_list, distance_list, radial_velocity_list = [],[],[],[],[],[]
        for object in icrs_skycoords:
            ra, dec, pmra_cosdec, pmdec, distance, radial_velocity = object.ra/u.deg, \
                                                                       object.dec/u.deg, \
                                                                       object.pm_ra_cosdec / (u.mas * u.yr), \
                                                                       object.pm_dec / (u.mas * u.yr), \
                                                                       object.distance / u.kpc, \
                                                                       object.radial_velocity / (u.km / u.s)
            # Discard the dimensionless unit.
            ra, dec, pmra_cosdec, pmdec, distance, radial_velocity = ra.value, dec.value, \
                                                                       pmra_cosdec.value, pmdec.value, \
                                                                       distance.value, radial_velocity.value
            # Remove cosdec, append to list
            pmra = pmra_cosdec / np.cos(math.radians(dec))
            ra_list.append(ra), dec_list.append(dec), pmra_list.append(pmra), pmdec_list.append(pmdec), \
            distance_list.append(distance), radial_velocity_list.append(radial_velocity)

        # Modify and save table.
        table['ra'] = ra_list
        table['dec'] = dec_list
        table['pmra'] = pmra_list
        table['pmdec'] = pmdec_list
        table['distance'] = distance_list
        table['vlos'] = radial_velocity_list
        writer.write_table(group, set, table)

    # Converts ICRS to Galactocentric instead.
    def to_GALCENT(self, hdfdir, hdfname, group, set):
        # Set up HDF and grab table, and SkyCoord objects for all targets.
        writer = hdfutils.hdf5_writer(hdfdir, hdfname)
        table = writer.read_table(group, set)
        skycoords = coord.SkyCoord(ra=table['ra']*u.deg,
                                   dec=table['dec']*u.deg,
                                   distance=table['distance']*u.kpc,
                                   pm_ra_cosdec=table['pmra']*np.cos(table['dec'])*u.mas/u.yr,
                                   pm_dec=table['pmdec']*u.mas/u.yr,
                                   radial_velocity=table['vlos']*u.km/u.s,
                                   frame="icrs")
        # Effect conversion to Galactocentric, work through objects, collect converted quantities.
        galcent_skycoords = skycoords.transform_to(coord.Galactocentric)
        x_list, y_list, z_list, vx_list, vy_list, vz_list = [], [], [], [], [], []
        for object in galcent_skycoords:
            x,y,z,vx,vy,vz = object.x/u.kpc, object.y/u.kpc, object.z/u.kpc, \
                             object.v_x/(u.km/u.s), object.v_y/(u.km/u.s), object.v_z/(u.km/u.s)
            # Discard the dimensionless unit.
            x,y,z,vx,vy,vz = x.value,y.value,z.value,vx.value,vy.value,vz.value

            # Remove cosdec, append to list
            x_list.append(x), y_list.append(y), z_list.append(z), \
            vx_list.append(vx), vy_list.append(vy), vz_list.append(vz)

        # Modify and save table.
        table['x'],table['y'],table['z'],table['vx'],table['vy'],table['vz'] = x_list,y_list,z_list,\
                                                                               vx_list,vy_list,vz_list
        writer.write_table(group, set, table)

# Some tools for angular-related things.
class angular(object):
    def __init__(self):
        self.null = "null"

    # Credit to Vincent for fast code.
    # https://stackoverflow.com/questions/4116658/faster-numpy-cartesian-to-spherical-coordinate-conversion
    def asCartesian(self, rthetaphi):
        # takes list rthetaphi (single coord)
        r = rthetaphi[0]
        theta = rthetaphi[1] * pi / 180  # to radian
        phi = rthetaphi[2] * pi / 180
        x = r * sin(theta) * cos(phi)
        y = r * sin(theta) * sin(phi)
        z = r * cos(theta)
        return [x, y, z]
    def asSpherical(self, xyz):
        # takes list xyz (single coord)
        x = xyz[0]
        y = xyz[1]
        z = xyz[2]
        r = sqrt(x * x + y * y + z * z)
        theta = acos(z / r) * 180 / pi  # to degrees
        phi = atan2(y, x) * 180 / pi
        return [r, theta, phi]

    # Grab angular momenta for table and save to/return table. Saves magnitude of angular momentum vectors, too.
    def get_momentum(self, table):
        pos, vel = np.array([table['x'],table['y'],table['z']]).T, np.array([table['vx'],table['vy'],table['vz']]).T
        mom = np.array([np.cross(pos[d], vel[d]) for d in range(len(pos))])
        mom_mag = [np.linalg.norm(z) for z in mom]
        table['Lx'], table['Ly'], table['Lz'] = mom.T
        table['L'] = mom_mag
        return table

    # Generate coordinates for system in galactocentric polar form (note: analogous to galactic coordinates.) RHANDED!
    # In this case the latitude is used in lieu of theta, the mathematical spherical polar form: NGP is 90 deg theta
    # Useful for Hammer-Aitoff Projections https://matplotlib.org/3.1.0/gallery/subplots_axes_and_figures/geo_demo.html
    def get_latipolar(self, table):
        table = self.get_polar(table)
        table['theta'] = 90 - table['theta']
        return table

    # Non-latitude definition (regular mathematical) where NGP is 0 deg theta.
    def get_polar(self, table):
        pos = np.array([table['x'], table['y'], table['z']]).T
        radii = np.array([np.linalg.norm(d) for d in pos])
        table['r'] = radii
        thetas = np.zeros(shape=np.shape(radii), dtype=np.float64)
        phis = np.zeros(shape=np.shape(radii), dtype=np.float64)
        for i in range(len(radii)):
            theta = math.acos(pos[i][2] / radii[i])
            phi = math.atan2(pos[i][1], pos[i][0])
            if phi < 0:
                phi += 2*np.pi
            thetas[i] = math.degrees(theta)
            phis[i] = math.degrees(phi)
        table['theta'] = thetas
        table['phi'] = phis
        return table

    # Just get the radii for the array (if polar angles are unnecessary and you want radial distance to system centre)
    def get_radii(self, table):
        pos = np.array([table['x'], table['y'], table['z']]).T
        radii = np.array([np.linalg.norm(d) for d in pos])
        table['r'] = radii
        return table

# Galactocentric system rotations/etc. x-convention. Recalculates angular momentum, too, in the new system.
# https://mathworld.wolfram.com/EulerAngles.html Euler Convention (Goldstein 1980)
# Give angles in degrees.
class galrotation(object):
    def __init__(self, phi, theta, psi):
        self.phi, self.theta, self.psi = math.radians(phi), \
                                         math.radians(theta), \
                                         math.radians(psi)
        # Set up the rotation matrices
        self.d = np.array([[np.cos(math.radians(phi)), np.sin(math.radians(phi)), 0],
                           [-np.sin(math.radians(phi)), np.cos(math.radians(phi)), 0],
                           [0, 0, 1]])
        self.c = np.array([[1, 0, 0],
                           [0, np.cos(math.radians(theta)), np.sin(math.radians(theta))],
                           [0, -np.sin(math.radians(theta)), np.cos(math.radians(theta))]])
        self.b = np.array([[np.cos(math.radians(psi)), np.sin(math.radians(psi)), 0],
                           [-np.sin(math.radians(psi)), np.cos(math.radians(psi)), 0],
                           [0, 0, 1]])

    # Rotate table to give coordinates inside x-convention coordinate system
    def rotate(self, table):
        pos = np.array([table['x'], table['y'], table['z']]).T
        vel = np.array([table['vx'], table['vy'], table['vz']]).T
        rotrix = np.matmul(np.matmul(self.d,self.c),self.b)
        pos_rot = np.array([np.matmul(rotrix,d) for d in pos]).T
        vel_rot = np.array([np.matmul(rotrix,d) for d in vel]).T
        table['x'],table['y'],table['z'] = pos_rot
        table['vx'], table['vy'], table['vz'] = vel_rot
        table = angular().get_momentum(table)
        return table

    # Derotate the table, undoing the x-convention coordinate transformation
    def derotate(self, table):
        pos = np.array([table['x'], table['y'], table['z']]).T
        vel = np.array([table['vx'], table['vy'], table['vz']]).T
        rotrix = np.linalg.inv(np.matmul(np.matmul(self.d,self.c),self.b))
        pos_rot = np.array([np.matmul(rotrix,d) for d in pos]).T
        vel_rot = np.array([np.matmul(rotrix,d) for d in vel]).T
        table['x'], table['y'], table['z'] = pos_rot
        table['vx'], table['vy'], table['vz'] = vel_rot
        table = angular().get_momentum(table)
        return table

    # Return the unit vectors (in the old system) of the new system.
    def unit(self):
        rotrix = np.matmul(np.matmul(self.d,self.c),self.b)
        unitx = np.matmul(np.array([1,0,0]),rotrix)
        unity = np.matmul(np.array([0,1,0]),rotrix)
        unitz = np.matmul(np.array([0,0,1]),rotrix)
        return np.array([unitx, unity, unitz])

# Great Circle Cell Counts in the Galactocentric System, as defined by 1996 Johnston Paper.
# Note: Designed to work in Standard Polar, not Latipolar. ALL IN DEGREES!
class greatcount(object):
    def __init__(self):
        self.null = "null"

    # Given a table, theta, phi, delta-theta, grab all members within this cell from table, produce new table.
    # NOTE: Poor Performance. If running optimization, switch to arrays and indices. Rebuild table post-calculation.
    def gcc_table(self, table, theta, phi, dtheta, radmin):
        # Set up unit vector for GCC
        theta, phi, dtheta = math.radians(theta), math.radians(phi), math.radians(dtheta)
        polar_unit = np.array([np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)])
        # Get unit vectors for all the targets
        pos = np.array([table['x'], table['y'], table['z']]).T
        radii = np.array([np.linalg.norm(d) for d in pos])
        radcondition = [True if d >= radmin else False for d in radii]
        pos_unit = np.array([pos[d]/radii[d] for d in range(len(pos))])
        # Dot them
        dotted = [np.dot(polar_unit, pos_unit[d]) for d in range(len(pos_unit))]
        condition = np.sin(dtheta)
        indices = []
        # Check for conditions and build new table.
        for num,item in enumerate(dotted):
            if radcondition[num] == True:
                if abs(item) <= condition:
                    indices.append(num)
        # Return the newly-built table
        return table[indices]











