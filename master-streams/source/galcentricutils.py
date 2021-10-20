import time

import astropy
import scipy
from astropy.stats import sigma_clipped_stats
from astropy.table import Table
from numba import jit
from numpy.random import bit_generator
from scipy.stats import norm
from sklearn.mixture import GaussianMixture, BayesianGaussianMixture
from sklearn.preprocessing import normalize
from sympy import *
import math
import graphutils
import hdfutils
import os
import numpy as np
from numpy import random
import astropy.coordinates as coord
from astropy.coordinates import galactocentric_frame_defaults
import astropy.units as u
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans, dbscan, DBSCAN, AgglomerativeClustering, OPTICS, AffinityPropagation
import windows_directories
from astrowrapper import utils
import plotly.express as px
import plotly.io as plio
import hdbscan
plio.renderers.default = "browser"
# SOME NOTES
"""
- We're using right-handed galactocentric coordinates
- The spherical representation uses declination-based polar angle: not typical mathematical polar angle. 
- Galactic Frame is right handed, too. X to the centre, Y to the left, Z upwards to NGP. 
- Work in Galactic or Galactocentric- SCREW ICRS APPARENTLY. 
"""

# Convert from ICRS to Galactocentric and back with custom solar system parameters (set via .dat file of right format.)
class galconversion(object):
    def __init__(self):
        self.owd, self.sol_params = os.getcwd(), np.zeros((1,4))
        self.source = "null"

    # Grab hold of the parameters of the solar system from source_solar and SETS THESE FOR ASTROPY.
    # Note that you have to provide source_solar: to deal with conflicting sources. See samples.
    def solinfo_grab(self, sourcedir, source_solar):
        self.source = sourcedir
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

    # Convert Galactic to ICRS (modified for nowrite)
    def nowrite_GAL_to_ICRS(self, table):
        # Set up HDF and grab table, and SkyCoord objects for all targets.
        #writer = hdfutils.hdf5_writer(hdfdir, hdfname)
        #table = writer.read_table(group, set)
        skycoords = coord.SkyCoord(l=table['l'] * u.deg, b=table['b'] * u.deg, distance=table['dist'] * u.kpc,
                                   pm_l_cosb=table['dmu_l']*np.cos(np.deg2rad(table['b'])) * u.mas/u.yr, pm_b=table['dmu_b']*u.mas/u.yr, radial_velocity=table['vlos']*u.km/u.s,
                                   frame="galactic")
        # Effect conversion to ICRS, work through objects, collect converted quantities.
        icrs_skycoords = skycoords.transform_to(coord.ICRS)
        ra_list, dec_list, pmra_list, pmdec_list, distance_list, radial_velocity_list = [],[],[],[],[],[]
        for object in icrs_skycoords:
            ra, dec, pmra_cosdec, pmdec = object.ra/u.deg, \
                                          object.dec/u.deg, \
                                          object.pm_ra_cosdec / (u.mas * u.yr), \
                                          object.pm_dec / (u.mas * u.yr)

            # Discard the dimensionless unit.
            ra, dec, pmra_cosdec, pmdec = ra.value, dec.value, \
                                          pmra_cosdec.value, pmdec.value

            # Remove cosdec, append to list
            pmra = pmra_cosdec / math.cos(math.radians(dec))
            ra_list.append(ra), dec_list.append(dec), pmra_list.append(pmra), pmdec_list.append(pmdec)

        # Modify and save table.
        table['ra'] = ra_list
        table['dec'] = dec_list
        table['pmra'] = pmra_list
        table['pmdec'] = pmdec_list
        return table
        #writer.write_table(group, set, table)

    # Converts ICRS to Galactic instead.
    def ICRS_to_GAL(self, hdfdir, hdfname, group, set):
        # Set up HDF and grab table, and SkyCoord objects for all targets.
        writer = hdfutils.hdf5_writer(hdfdir, hdfname)
        table = writer.read_table(group, set)
        skycoords = coord.SkyCoord(ra=table['ra']*u.deg,
                                   dec=table['dec']*u.deg,
                                   distance=table['dist']*u.kpc,
                                   pm_ra_cosdec=table['pmra']*np.cos(np.deg2rad(table['dec']))*u.mas/u.yr,
                                   pm_dec=table['pmdec']*u.mas/u.yr,
                                   radial_velocity=table['vlos']*u.km/u.s,
                                   frame="icrs")
        # Effect conversion to Galactocentric, work through objects, collect converted quantities.
        gal_skycoords = skycoords.transform_to(coord.Galactic)
        l_list,b_list,dmu_l_list,dmu_b_list = [],[],[],[]
        for object in gal_skycoords:
            l,b,dmub = object.l/u.deg, object.b/u.deg, object.pm_b/(u.mas/u.yr)
            l,b,dmub = l.value,b.value,dmub.value
            dmul = ((object.pm_l_cosb/(u.mas/u.yr)).value)/math.cos(math.radians(b))

            l_list.append(l),b_list.append(b),dmu_b_list.append(dmub),dmu_l_list.append(dmul)

        # Modify and save table.
        table['l'],table['b'],table['dmu_l'],table['dmu_b'] = l_list,b_list,dmu_l_list,dmu_b_list
        writer.write_table(group, set, table)

    # Galactocentric to Galactic.
    def GALCENT_to_GAL(self, hdfdir, hdfname, group, set):
        # Set up HDF and grab table, and SkyCoord objects for all targets.
        writer = hdfutils.hdf5_writer(hdfdir, hdfname)
        table = writer.read_table(group, set)
        skycoords = coord.SkyCoord(x=table['x'] * u.kpc, y=table['y'] * u.kpc, z=table['z'] * u.kpc,
                                   v_x=table['vx'] * u.km / u.s, v_y=table['vy'] * u.km / u.s,
                                   v_z=table['vz'] * u.km / u.s,
                                   frame="galactocentric")
        # Effect conversion to ICRS, work through objects, collect converted quantities.
        gal_skycoords = skycoords.transform_to(coord.Galactic)
        l_list,b_list,dmu_l_list,dmu_b_list, vlos_list, dist_list = [],[],[],[],[],[]
        for object in gal_skycoords:
            l,b,dmub = object.l/u.deg, object.b/u.deg, object.pm_b/(u.mas/u.yr)
            l,b,dmub = l.value,b.value,dmub.value
            vlos = object.radial_velocity/(u.km/u.s)
            vlos = vlos.value
            dmul = ((object.pm_l_cosb/(u.mas/u.yr)).value)/math.cos(math.radians(b))
            dist_list.append((object.distance/u.kpc).value)
            l_list.append(l),b_list.append(b),dmu_b_list.append(dmub),dmu_l_list.append(dmul)
            vlos_list.append(vlos)

        # Modify and save table.
        table['l'],table['b'],table['dmu_l'],table['dmu_b'] = l_list,b_list,dmu_l_list,dmu_b_list
        table['vlos'] = vlos_list
        table['dist'] = dist_list
        writer.write_table(group, set, table)

    # Galactic to Galactocentric. Second option is table in-table out.
    def GAL_to_GALCENT(self, hdfdir, hdfname, group, set):
        # Set up HDF and grab table, and SkyCoord objects for all targets.
        writer = hdfutils.hdf5_writer(hdfdir, hdfname)
        table = writer.read_table(group, set)
        skycoords = coord.SkyCoord(l=table['l'] * u.deg, b=table['b'] * u.deg, distance=table['dist'] * u.kpc,
                                   pm_l_cosb=table['dmu_l'] * np.cos(np.radians(table['b'])) * u.mas / u.yr,
                                   pm_b=table['dmu_b'] * u.mas / u.yr, radial_velocity=table['vlos'] * u.km / u.s,
                                   frame="galactic")
        # Effect conversion to ICRS, work through objects, collect converted quantities.
        galcent_skycoords = skycoords.transform_to(coord.Galactocentric)
        x_list, y_list, z_list, vx_list, vy_list, vz_list = [], [], [], [], [], []
        for object in galcent_skycoords:
            x,y,z,vx,vy,vz = object.x.to(u.kpc), object.y.to(u.kpc), object.z.to(u.kpc), \
                             object.v_x.to(u.km/u.s), object.v_y.to(u.km/u.s), object.v_z.to(u.km/u.s)
            # Discard the dimensionless unit.
            x,y,z,vx,vy,vz = x.value,y.value,z.value,vx.value,vy.value,vz.value

            # Append to list
            x_list.append(x), y_list.append(y), z_list.append(z), \
            vx_list.append(vx), vy_list.append(vy), vz_list.append(vz)

        # Modify and save table.
        table['x'],table['y'],table['z'],table['vx'],table['vy'],table['vz'] = x_list,y_list,z_list,\
                                                                               vx_list,vy_list,vz_list
        writer.write_table(group, set, table)
    def nowrite_GAL_to_GALCENT(self, table):

        # Set up HDF and grab table, and SkyCoord objects for all targets.
        skycoords = coord.SkyCoord(l=table['l'] * u.deg, b=table['b'] * u.deg, distance=table['dist'] * u.kpc,
                                   pm_l_cosb=table['dmu_l'] * np.cos(np.radians(table['b'])) * u.mas / u.yr,
                                   pm_b=table['dmu_b'] * u.mas / u.yr, radial_velocity=table['vlos'] * u.km / u.s,
                                   frame="galactic")
        # Effect conversion to ICRS, work through objects, collect converted quantities.
        galcent_skycoords = skycoords.transform_to(coord.Galactocentric)
        x_list, y_list, z_list, vx_list, vy_list, vz_list = [], [], [], [], [], []
        for object in galcent_skycoords:
            x,y,z,vx,vy,vz = object.x.to(u.kpc), object.y.to(u.kpc), object.z.to(u.kpc), \
                             object.v_x.to(u.km/u.s), object.v_y.to(u.km/u.s), object.v_z.to(u.km/u.s)
            # Discard the dimensionless unit.
            x,y,z,vx,vy,vz = x.value,y.value,z.value,vx.value,vy.value,vz.value

            # Append to list
            x_list.append(x), y_list.append(y), z_list.append(z), \
            vx_list.append(vx), vy_list.append(vy), vz_list.append(vz)

        # Modify and save table.
        table['x'],table['y'],table['z'],table['vx'],table['vy'],table['vz'] = x_list,y_list,z_list,\
                                                                               vx_list,vy_list,vz_list
        return table

    # Single point GAL to GALCENT. [X,Y,Z,VX,VY,VZ]
    def vec_GAL_to_GALCENT(self, l, b, distance, dmul, dmub, vlos):
        skycoords = coord.SkyCoord(l=l * u.deg, b=b * u.deg, distance=distance * u.kpc,
                                   pm_l_cosb=dmul * np.cos(np.radians(b)) * u.mas / u.yr,
                                   pm_b=dmub * u.mas / u.yr, radial_velocity=vlos * u.km / u.s,
                                   frame="galactic")
        # Effect conversion to ICRS, work through objects, collect converted quantities.
        galcent_skycoords = skycoords.transform_to(coord.Galactocentric)
        x, y, z, vx, vy, vz = galcent_skycoords.x / u.kpc, galcent_skycoords.y / u.kpc, galcent_skycoords.z / u.kpc, \
                              galcent_skycoords.v_x / (u.km / u.s), galcent_skycoords.v_y / (u.km / u.s), \
                              galcent_skycoords.v_z / (u.km / u.s)
        x, y, z, vx, vy, vz = x.value, y.value, z.value, vx.value, vy.value, vz.value
        return [x,y,z,vx,vy,vz]

    # GALCENT to GAL
    def vec_GALCENT_to_GAL(self, x, y, z, vx, vy, vz):
        skycoords = coord.SkyCoord(x=x * u.kpc, y=y * u.kpc, z=z * u.kpc,
                                   v_x=vx * u.km / u.s, v_y=vy * u.km / u.s,
                                   v_z=vz * u.km / u.s,
                                   frame="galactocentric")
        # Effect conversion to ICRS, work through objects, collect converted quantities.
        gal_skycoords = skycoords.transform_to(coord.Galactic)
        l, b, dmub = gal_skycoords.l / u.deg, gal_skycoords.b / u.deg, gal_skycoords.pm_b / (u.mas / u.yr)
        l, b, dmub = l.value, b.value, dmub.value
        distance = (gal_skycoords.distance/u.kpc).value
        vlos = gal_skycoords.radial_velocity / (u.km / u.s)
        vlos = vlos.value
        dmul = ((gal_skycoords.pm_l_cosb / (u.mas / u.yr)).value) / math.cos(math.radians(b))
        return [l, b, distance, dmul, dmub, vlos]

    # Fully-manual GAL to GALCENT
    def manual_vec_GAL_to_GALCENT(self, l, b, distance, dmul, dmub, vlos):
        l,b = math.radians(l),math.radians(b)
        # Sun position: Galactic Cartesian in kpc.
        sun_position = np.array(self.sol_params[0])
        # Proper Motions in mas/year: convert to radians/year. Vlos is fine (kms^-1)
        dmul, dmub = dmul*(1e-3)*(1/3600)*(np.pi/180), dmub*(1e-3)*(1/3600)*(np.pi/180)
        # Now radians/second
        year = ((1*u.yr).to(u.s)).value
        dmul, dmub = dmul/year, dmub/year
        # Distance in kpc: go to metres.
        distance_m = ((distance*u.kpc).to(u.m)).value
        # Calculate speed along l/b in ms^-1
        vl, vb = dmul*distance_m*math.cos(b), dmub*distance_m
        # Make unit vectors for velocity/position
        b_unit = np.array([-math.cos(l)*math.sin(b),
                           -math.sin(l)*math.sin(b),
                           math.cos(b)])
        l_unit = np.array([-math.sin(l),
                           math.cos(l),
                           0])
        r_unit = np.array([math.cos(l)*math.cos(b),
                           math.sin(l)*math.cos(b),
                           math.sin(b)])
        # This is in Galactic Cartesian (from us.) Assume equivalent unit vectors (just a position offset). Vlos/kms^-1
        target_galactic_cartesian_velocity = vl*l_unit + vb*b_unit + r_unit*vlos*1e3
        # Sun velocity [x,y,z] in kms^-1. Go to ms^-1. This is in galactocentric cartesian fine.
        sun_velocity = np.array(self.sol_params[2])*1e3
        # Finally get the target galactocentric velocity
        target_galactocentric_velocity = (target_galactic_cartesian_velocity + sun_velocity)/1e3
        # Get the target position
        target_galactocentric_position = sun_position + (r_unit*distance)
        posvel = np.append(target_galactocentric_position, target_galactocentric_velocity)
        return posvel

    # Fully-manual GALCENT to GAL
    def manual_vec_GALCENT_to_GAL(self, x, y, z, vx, vy, vz):
        # Position transformation, Galactic Cartesian
        solar_position = np.array(self.sol_params[0])
        galactic_target_position = np.array([x,y,z]) - solar_position
        distance = np.linalg.norm(galactic_target_position)
        # That's all in Galactic XYZ. z = distance*sin(b)
        b = math.asin(galactic_target_position[2]/distance)
        l = math.atan2(galactic_target_position[1],galactic_target_position[0])
        if l < 0:
            l += 2*np.pi

        # Make unit vectors for velocity/position
        b_unit = np.array([-math.cos(l)*math.sin(b),
                           -math.sin(l)*math.sin(b),
                           math.cos(b)])
        l_unit = np.array([-math.sin(l),
                           math.cos(l),
                           0])
        r_unit = np.array([math.cos(l)*math.cos(b),
                           math.sin(l)*math.cos(b),
                           math.sin(b)])

        # Velocity transformation
        solar_velocity = np.array(self.sol_params[2]) # in kms^-1
        galactic_target_velocity = np.array([vx,vy,vz]) - solar_velocity
        # Get components
        vlos = np.dot(galactic_target_velocity, r_unit)
        vb,vl = np.dot(galactic_target_velocity, b_unit), np.dot(galactic_target_velocity, l_unit)
        # Get dmub,dmul in radians/second
        dmub, dmul = (((((vb*u.km/u.s).to(u.kpc/u.yr))/(distance*u.kpc))*u.rad).to(u.mas/u.yr)).value,\
                     (((((vl*u.km/u.s).to(u.kpc/u.yr))/(distance*math.cos(b)*u.kpc))*u.rad).to(u.mas/u.yr)).value
        l,b = math.degrees(l),math.degrees(b)

        return l, b, distance, dmul, dmub, vlos

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
        table['r'] = [np.linalg.norm(d) for d in pos]
        return table

    # Just for a single point: magnitude not needed.
    def vec_get_momentum(self,vec):
        pos, vel = vec[0:3],vec[3:6]
        L = np.cross(pos,vel)
        return L

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

    # All the above operations but re-cast to deal with a VEC vector [xyzvxvyvz] instead, vector as a numpy array.
    def vec_momentum(self,vec):
        pos, vel = vec[0:3],vec[3:6]
        mom = np.cross(pos, vel)
        return mom

    def vec_latipolar(self, vec):
        polar_pos = self.vec_polar(vec)
        polar_pos[1] = 90 - polar_pos[1]
        return polar_pos

    def vec_polar(self, vec):
        pos = vec[0:3]
        radius = np.linalg.norm(pos)
        theta = math.acos(pos[2]/radius)
        phi = math.atan2(pos[1],pos[0])
        if phi < 0:
            phi += 2*np.pi
        theta, phi = math.degrees(theta), math.degrees(phi)
        polar_pos = np.array([radius, theta, phi])
        return polar_pos

    def vec_get_radii(self, vec):
        return np.linalg.norm(vec[0:3])



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



# Monte-Carlo Sampling of Angular Momentum Errors. Includes per-value and table-wide sampling. Tidbits...
"""
[l,b,distance,dmul,dmub,vlos,edist,edmul,edmub,evlos]
Galactic Frame to Galactocentric Frame Sampling!!!
Before running anything make sure to use galdefine to make sure your solar coords are correct for astropy 
"""
class monte_angular(object):
    def __init__(self):
        self.converter = galconversion()
    # Define galconversion for self.
    def galdefine(self, source, solinfo):
        self.converter.solinfo_grab(source, solinfo)
        self.converter.solgal_set()
    # Do a monte-carlo for a given [l,b,distance,dmul,dmub,vlos,edist,edmul,edmub,evlos] vector.
    def vec_monte(self, vec, n):
        # Generate spread for vec (errors in parameter space)
        dists = np.random.default_rng().normal(vec[2], vec[6], n)
        dists = np.abs(dists) # negatives happen- bad for calculating.
        dmuls = np.random.default_rng().normal(vec[3], vec[7], n)
        dmubs = np.random.default_rng().normal(vec[4], vec[8], n)
        vloses = np.random.default_rng().normal(vec[5], vec[9], n)
        ls, bs = [vec[0] for d in dists], [vec[1] for d in dists]
        # Create table
        astrotable = Table()
        astrotable['l'], astrotable['b'], astrotable['dist'], astrotable['dmu_l'], astrotable['dmu_b'], \
        astrotable['vlos'] = ls, bs, dists, dmuls, dmubs, vloses
        # Convert to Galactocentric
        astrotable = self.converter.nowrite_GAL_to_GALCENT(astrotable)
        astrotable = angular().get_momentum(astrotable)
        # Get each axis and get the mean/stds.
        Lx,Ly,Lz = astrotable['Lx'],astrotable['Ly'],astrotable['Lz']
        L = [Lx, Ly, Lz]
        means, stds = [],[]
        for num, comp in enumerate(L):
            mean, median, sigma = sigma_clipped_stats(comp)
            means.append(mean), stds.append(sigma)

        return np.append(means, stds)
    # Given an astrotable, monte-carlo the entire darned thing for step n.
    def table_monte(self, table, n):
        table = self.converter.nowrite_GAL_to_GALCENT(table)
        table = angular().get_momentum(table)
        dLx,dLy,dLz = [],[],[]
        for row in table:
            l, b, dist, dmul, dmub, vlos = row['l'], row['b'], row['dist'], row['dmu_l'], row['dmu_b'], row['vlos']
            edist, edmul, edmub, evlos = row['edist'], row['edmu_l'], row['edmu_b'], row['evlost']
            vec = [l,b,dist,dmul,dmub,vlos,edist,edmul,edmub,evlos]
            monte = self.vec_monte(vec, n)
            dLx.append(monte[3]),dLy.append(monte[4]),dLz.append(monte[5])
        table['dLx'],table['dLy'],table['dLz'] = dLx,dLy,dLz
        return table

    # Same as above, but returns Covariance Matrix.
    # [l,b,distance,dmul,dmub,vlos,edist,edmul,edmub,evlos] vector.
    def vec_covmonte(self, vec, n):
        # Get the angular momentum for this vector (judged as the "mean")
        vec_galcent = self.converter.vec_GAL_to_GALCENT(*vec[0:6])
        vec_L = angular().vec_momentum(vec_galcent)

        # Generate spread for vec (errors in parameter space)
        dists = np.random.default_rng().normal(vec[2], vec[6], n)
        dists = np.abs(dists) # negatives happen- bad for calculating.
        dmuls = np.random.default_rng().normal(vec[3], vec[7], n)
        dmubs = np.random.default_rng().normal(vec[4], vec[8], n)
        vloses = np.random.default_rng().normal(vec[5], vec[9], n)
        ls, bs = [vec[0] for d in dists], [vec[1] for d in dists]

        # Create table
        astrotable = Table()
        astrotable['l'], astrotable['b'], astrotable['dist'], astrotable['dmu_l'], astrotable['dmu_b'], \
        astrotable['vlos'] = ls, bs, dists, dmuls, dmubs, vloses

        # Convert to Galactocentric
        astrotable = self.converter.nowrite_GAL_to_GALCENT(astrotable)
        astrotable = angular().get_momentum(astrotable)

        # Generate the Covariance Matrix
        Lx,Ly,Lz = astrotable['Lx'],astrotable['Ly'],astrotable['Lz']
        L = np.array([Lx, Ly, Lz]).T
        L_dev = L - vec_L # difference from mean, vector representation
        L_dev = L_dev.T # list representation

        # Calculate Covariance Matrix
        cov = np.zeros(shape=(3,3))
        for i in range(3):
            for j in range(3):
                cov[i,j] = np.mean(L_dev[i]*L_dev[j], axis=0)

        # Get the standard deviations
        stdevs = np.zeros(shape=(1,3))
        for i in range(3):
            dev_i = math.sqrt(cov[i,i])
            stdevs[0,i] = dev_i

        # Calculate the Pearson Coefficients for this matrix (returning them as an array, too.)
        pears = np.zeros(shape=(3,3))
        for i in range(3):
            for j in range(3):
                if i != j:
                    pears[i,j] = cov[i,j]/(stdevs[0,i]*stdevs[0,j])

        return stdevs, cov, pears

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

    # Generate a hemispheres-worth of equally separated n-points.
    # A simple scheme for generating nearly uniform distribution of antipodally
    # symmetric points on the unit sphere
    # Cheng Guan Koay∗
    def halfgridgen(self, K_top):
        # Initial guess
        n = np.sqrt(K_top*np.pi/8)

        # Basic iterator.
        max_iter = 200
        iter_token = 0
        while True:
            n = (K_top/2) * np.sin(np.pi/(4*n))
            iter_token += 1
            if iter_token >= max_iter:
                break

        # Round
        n = np.rint(n, out=np.zeros(1, int), casting='unsafe')[0]

        # Define i-range
        i_range = np.arange(1, n+1, 1)

        # Generate sinangles
        thetas = np.array([(i - (1/2))*(np.pi/(2*n)) for i in i_range])
        length = 2*np.pi*np.sin(thetas)

        # Generate the K for each latitude element: using the k_i^5 routine.
        #is_odd = lambda num: num & 0x1 # 0 or 1 for not-odd and odd
        k_range = []
        for i in i_range:
            if i < n:
                ki = length[i-1]*K_top/(np.pi/np.sin(np.pi/(4*n)))
                ki = np.rint(ki, out=np.zeros(1, int), casting='unsafe')[0]
                k_range.append(ki)
            if i == n:
                ki = K_top - np.sum(np.array(k_range), axis=0)
                ki = np.rint(ki, out=np.zeros(1, int), casting='unsafe')[0]
                k_range.append(ki)

        # Set up the phis for each latitude element
        phis = []
        for num, theta in enumerate(thetas):
            phi_theta = [(j-(1/2))*(2*np.pi)/k_range[num] for j in np.arange(1, k_range[num]+1,1)]
            phis.append(phi_theta)

        # Unpack the phi lists
        thetaphis = []
        for num, theta in enumerate(thetas):
            thetaphis += list(zip([theta for d in phis[num]], phis[num]))
        thetas_orig, phis_orig = np.array(thetaphis).T

        # Specify an index (to match to antipodal point later) on the list (1->len(thetas_orig): can't have 0.)
        # The reflected index will just be -(this index.)
        indices = np.arange(1, len(thetas_orig) + 1, 1)

        # Return the grid
        return thetas_orig, phis_orig, indices

    # Take the half grid and produce a southern hemisphere too, just a 180 degree reflection of the top.
    # Make sure K_full is EVEN.
    def fullgridgen(self, K_full):
        # Get half-hemispheres worth of points
        thetas_orig, phis_orig, indices = self.halfgridgen(np.rint(K_full/2, out=np.zeros(1, int), casting='unsafe')[0])
        thetaphis = list(zip(thetas_orig,phis_orig))

        # Reflected Phis = Phi + 180 degrees, Reflected Thetas = Pi - Theta
        newthetaphis = []
        new_indices = []
        for index, thetaphi in list(zip(indices,thetaphis)):
            new_indices.append(index*int(-1))
            thetaphi_arr = np.zeros(shape=(2))
            thetaphi_arr[0] = np.pi - thetaphi[0]
            thetaphi_arr[1] = thetaphi[1] + np.pi
            if thetaphi_arr[1] > 2*np.pi:
                thetaphi_arr[1] -= 2*np.pi
            newthetaphis.append(thetaphi_arr)
        # Return thetas, phis, indices
        thetas,phis = np.array(newthetaphis).T
        thetas,phis = np.append(thetas_orig, thetas), np.append(phis_orig, phis)
        indices = np.append(np.array(indices),np.array(new_indices))
        return thetas, phis, indices



# Clustering
class cluster(object):
    def __init__(self):
        self.null = "null"

    # Remove outlying L-values (units of sigma). Assumes r has already been cleaned out (see: gcc_table)
    def clean(self, table, sig_tolerance):

        Ls = [table['Lx'],table['Ly'],table['Lz']]
        mean_L, std_L = np.array([np.mean(d) for d in Ls]),np.array([np.std(d) for d in Ls])
        table_cleaned = table
        for num,row in enumerate(table):
            L = np.array([row['Lx'],row['Ly'],row['Lz']])
            L_dif = L - mean_L
            mag_L_dif = np.array([abs(d) for d in L_dif])
            sig_L_dif = mag_L_dif/std_L
            for i in sig_L_dif:
                if i >= sig_tolerance:
                    table_cleaned.remove_row(num)
                    break
        return table_cleaned


    # kmeans-cluster the given table: returns the inertia of the table, too.
    def kmeans(self, table, k, savedex, browser):
        # Set up vectors/positions/etc
        table = self.clean(table, 5)
        L = np.array([table['Lx'], table['Ly'], table['Lz']]).T
        km = KMeans(n_clusters=k,n_init=10,max_iter=300,algorithm="full")
        kmfit = km.fit_predict(L) # list with indices for cluster
        inertia = km.inertia_
        table['k_index'] = np.array(kmfit)
        graphutils.threed_graph().kmeans_L(table, savedex + ".html", browser)
        return table, inertia

    # DBSCAN. params are "eps" and "min_samples." Browser=True opens in browser (append to others.)
    def dbs(self, table, eps, min_samples, browser):
        table = self.clean(table, 5)
        L = np.array([table['Lx'], table['Ly'], table['Lz']]).T
        dbs = DBSCAN(eps=eps, min_samples=min_samples, metric="l1", leaf_size=5)
        dbsfit = dbs.fit_predict(L)
        table['k_index'] = np.array(dbsfit)
        save_format = ("DBS_TEST_EPS{0}_MINSAMP{1}" + ".html").format(eps, min_samples)
        graphutils.threed_graph().kmeans_L(table, save_format, browser)
        return table, dbs

    # Hierarchical DBS
    def hdbs(self, table, browser):
        table = self.clean(table, 5)
        L = np.array([table['Lx'], table['Ly'], table['Lz']]).T
        hdbs = hdbscan.HDBSCAN(min_cluster_size=25,
                               min_samples=15,
                               metric="l2")
        hdbsfit = hdbs.fit_predict(L)
        table['k_index'] = np.array(hdbsfit)
        #save_format = ("HDBS_TEST_EPS{0}_MINSAMP{1}" + ".html").format(eps, min_samples)
        graphutils.threed_graph().kmeans_L(table, "test_hdbs.html", browser)
        graphutils.threed_graph().xmeans_L(table, "test_hdbs.html", browser)
        return table

    # OPTICS
    def optics(self, table):
        table = self.clean(table, 5)
        L = np.array([table['Lx'], table['Ly'], table['Lz']]).T
        optics = OPTICS(metric="l1", min_cluster_size=20, leaf_size=20, eps=1e4, max_eps=1e6)
        opticsfit = optics.fit_predict(L)
        table['k_index'] = np.array(opticsfit)
        graphutils.threed_graph().kmeans_L(table, "test_plot", False)
        return table

    # Agglomerative
    def aglom(self, table, k):
        table = self.clean(table,5)
        L = np.array([table['Lx'], table['Ly'], table['Lz']]).T
        aglom = AgglomerativeClustering(n_clusters=k, linkage='ward')
        aglomfit = aglom.fit_predict(L)
        table['k_index'] = np.array(aglomfit)
        graphutils.threed_graph().kmeans_L(table, "test_plot", False)
        return table

    # Affinity Propagation
    def afprop(self, table):
        table = self.clean(table, 5)
        L = np.array([table['Lx'], table['Ly'], table['Lz']]).T
        afprop = AffinityPropagation(damping=0.5)
        afpropfit = afprop.fit_predict(L)
        table['k_index'] = np.array(afpropfit)
        graphutils.threed_graph().kmeans_L(table, False, True)
        return table

    # Gaussian Mixture with Variational Bayes
    def bayesgaussmix(self, table, k, savedex, browser):
        # Set up vectors/positions/etc
        table = self.clean(table, 4)
        L = np.array([table['Lx'], table['Ly'], table['Lz']]).T
        gm = BayesianGaussianMixture(n_components=k)
        gmfit = gm.fit_predict(L) # list with indices for cluster
        table['k_index'] = np.array(gmfit)
        graphutils.threed_graph().kmeans_L(table, savedex + ".html", browser)
        return table

    # Gaussian Mixture - Euclidean
    def gaussmix(self, table, k, savedex, browser, graph):
        # Set up vectors/positions/etc
        table = self.clean(table, 5)
        L = np.array([table['Lx'], table['Ly'], table['Lz']]).T
        gm = GaussianMixture(n_components=k)
        gmfit = gm.fit_predict(L) # list with indices for cluster
        table['k_index'] = np.array(gmfit)
        if graph == True:
            graphutils.threed_graph().kmeans_L(table, savedex + ".html", browser)
        bic, aic = gm.bic(L), gm.aic(L)
        return table, bic, aic

    # Grab aics/bics for varioues values of (k) for table.
    def gaussaicsbics(self, table, k_max, savedex):
        y = np.arange(1, k_max, 1)
        bics, aics = [], []
        for i in y:
            gmeansdone = self.gaussmix(table, i, "test", browser=False, graph=False)
            bics.append(gmeansdone[1]), aics.append(gmeansdone[2])
        plt.plot(y, bics, color="red", label="bics")
        plt.plot(y, aics, color="green", label="aics")
        plt.legend()
        try:
            os.mkdir(windows_directories.imgdir + "\\gauss_aics_bics_test")
        except:
            pass

        plt.savefig(windows_directories.imgdir + "\\gauss_aics_bics_test\\" + savedex + ".png")
        plt.clf()
        plt.close()

