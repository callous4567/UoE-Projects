import copy
import pickle
import time
import astropy
import pandas
import scipy
import munkres
from astropy.stats import sigma_clipped_stats
from astropy.table import Table
from matplotlib import pyplot as plt
from numpy.random import bit_generator
from scipy.optimize import linear_sum_assignment
from sklearn.metrics.cluster import contingency_matrix
from sympy import *
import math
import seaborn as sns
import ascii_info
import energistics
import graphutils
import hdfutils
import os
import numpy as np
from numpy import random
import astropy.coordinates as coord
from astropy.coordinates import galactocentric_frame_defaults
import astropy.units as u
import plotly.io as plio
import hdbscan
from scipy.spatial import ConvexHull
import windows_directories
from hdbscan import flat
plio.renderers.default = "browser"
# SOME NOTES
"""
- We're using right-handed galactocentric coordinates
- The spherical representation uses declination-based polar angle: not typical mathematical polar angle. 
- Galactic Frame is right handed, too. X to the centre, Y to the left, Z upwards to NGP. 
- Work in Galactic or Galactocentric- SCREW ICRS APPARENTLY. 
"""

# TODO: 2D clustering like Sofie Groningen Paper (see for final method comparison.)

# TODO: Note that we're trying to make sure code is oriented to work with both pandas and astropy tables.
# TODO: Table Covmonte needs numpy multivariate normal generation instead (faster than generating each set individually)

# Handles Galactocentric Cartesian/Galactic Cartesian transformations (with ICRS optional-see deprected)
class galconversion(object):
    def __init__(self):
        self.owd, self.sol_params = os.getcwd(), np.zeros((1,4))
        self.source = "null"
        # Note about sol_params
        """
        0 is pos, 1 is err, 2 is vel, 3 is err
        x y z respectively.
        """

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
    # cartesian galcen_v_sun and cartesian Z in galcentric
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

# Some tools for spherical-polar angular-related things.
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
        r = math.sqrt(x * x + y * y + z * z)
        theta = math.acos(z / r) * 180 / np.pi  # to degrees
        phi = math.atan2(y, x) * 180 / np.pi
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

    # Convert to cylindrical coordinates: R, phi, vR, vT - necessary for galpy orbits module.
    # vR is radial velocity, vT is transverse velocity: it can be a negative quantity.
    def get_cylindrical(self, table):
        # Evaluate positions (R, phi) and unit vectors
        pos = np.array([table['x'], table['y']]).T
        norms = np.array([np.linalg.norm(d) for d in pos])
        R_unit = np.array([pos[d]/norms[d] for d in range(len(pos))])

        # Get R and Phi
        table['R'] = np.array([np.linalg.norm(d) for d in pos])
        phis = np.zeros(shape=np.shape(table['R']), dtype=np.float64)
        for i in range(len(table['R'])):
            phi = math.atan2(pos[i][1], pos[i][0])
            if phi < 0:
                phi += 2 * np.pi
            phis[i] = math.degrees(phi)
        table['phi'] = phis

        # Evaluate phi tangent unit vectors and get vT/vR
        """
        r = R[cos, sin]
        dr/dphi = R[-sin, cos]
        phi_unit = [-sin, cos]
        """
        phi_unit = np.array([np.array([-np.sin(phi),np.cos(phi)]) for phi in np.deg2rad(table['phi'])])
        v_xy = np.array([table['vx'],table['vy']]).T
        vR = np.array([np.dot(ru, vxy) for ru, vxy in zip(R_unit, v_xy)])
        vT = np.array([np.dot(pu, vxy) for pu, vxy in zip(phi_unit, v_xy)])
        table['vR'], table['vT'] = vR, vT

        return table

# Galactocentric Cartesian system rotations/etc.
# x-convention. Recalculates angular momentum, too, in the new system.
# Passive transformation: gives values in the new system, alongside the unit vectors of the new system (in the old.)
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
        dists = np.abs(dists) # negatives happen- bad for calculating. Introduces bias for low dist (not interested.)
                              # TODO: Even if a bias occurs here, it'll be so few that it shouldn't impact. Check later.
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
    # Given an astrotable, monte-carlo the entire darned thing for step n with the same vector as above.
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
    # TODO: Four other parameters to implicate in covariance estimate: L LZ E CIRCULARITY
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

        # Generate the deviations for each Monte'd point
        Lx,Ly,Lz = astrotable['Lx'],astrotable['Ly'],astrotable['Lz'] # list representation
        L = np.array([Lx, Ly, Lz]).T # get vector representation
        L_dev = L - vec_L # difference from mean, vector representation
        L_dev = L_dev.T # list representation

        # Calculate Covariance Matrix, FOR STANDARD L CLUSTERING
        cov = np.zeros(shape=(3,3))
        for i in range(3):
            for j in range(3):
                cov[i,j] = np.mean(L_dev[i]*L_dev[j], axis=0)

        # Get the standard deviations FOR STANDARD L CLUSTERING
        stdevs = []
        for i in range(3):
            dev_i = math.sqrt(cov[i, i])
            stdevs.append(dev_i)

        # Calculate Energy, L_magnitude, and Circularity (for 4D clustering) and their means
        orbi = energistics.orbigistics()
        table = orbi.orbilarity(astrotable)
        E_mean = np.mean(table['E'])
        circ_mean = np.mean(table['circ'])
        L_mag = np.array([np.linalg.norm(d) for d in L])
        L_mean = np.linalg.norm(vec_L)


        # Generate Deviations
        L_mag_dev = L_mag - L_mean
        E_dev = table['E'] - E_mean
        circ_dev = table['circ'] - circ_mean
        deviation_vector = [L_mag_dev,L_dev[2],E_dev,circ_dev]

        # Generate second Covariance Matrix, FOR 4D ABS(L), LZ, E, CIRC CLUSTERING!
        cov2 = np.zeros(shape=(4, 4))
        for i in range(4):
            for j in range(4):
                cov2[i,j] = np.mean(deviation_vector[i]*deviation_vector[j], axis=0)

        # Since we have vec_L, get a vector with these four variables, too (judged as the mean.)
        vec_4d = np.array([L_mean, vec_L[2], E_mean, circ_mean])

        # Get the standard deviations FOR 4D ABS(L), LZ, E, CIRC CLUSTERING!
        stdevs2 = []
        for i in range(4):
            dev_i = math.sqrt(cov2[i,i])
            stdevs2.append(dev_i)

        # Jorge also wants it for [Lx, Ly, Lz, E]
        cov3 = np.zeros(shape=(4, 4))
        deviation_vector = [L_dev[0], L_dev[1], L_dev[2], E_dev]
        vec_LE = [vec_L[0], vec_L[1], vec_L[2], E_mean]
        for i in range(4):
            for j in range(4):
                cov3[i,j] = np.mean(deviation_vector[i]*deviation_vector[j], axis=0)

        return cov, stdevs, vec_L, cov2, stdevs2, vec_4d, cov3, vec_LE
    # Returns a pandas dataframe instead: useful to avoid write conflicts. Cleanup after.
    # Same monte as above, except with covariance matrices under ['covtrix'] too.
    def table_covmonte(self, table, n):
        table = self.converter.nowrite_GAL_to_GALCENT(table)
        table = angular().get_momentum(table)
        covtrices, stdevslist, mean_angulars = [], [], [] # standard [lx ly lz] covtrices
        covtrices2, stdevslist2, mean_4ds = [], [], [] # sarah sofie lovdals [L Lz E circ]
        covtrices3, mean_4dLEs = [],[] # jorges request [lx ly lz E]
        for row in table:
            l, b, dist, dmul, dmub, vlos = row['l'], row['b'], row['dist'], row['dmu_l'], row['dmu_b'], row['vlos']
            edist, edmul, edmub, evlos = row['edist'], row['edmu_l'], row['edmu_b'], row['evlost']
            vec = [l, b, dist, dmul, dmub, vlos, edist, edmul, edmub, evlos]
            covtrix, stdevs, vec_L, covtrix2, stdevs2, vec_4d, covtrix3, vec_LE = self.vec_covmonte(vec, n)
            covtrices.append(covtrix), stdevslist.append(stdevs), mean_angulars.append(vec_L)
            covtrices2.append(covtrix2), stdevslist2.append(stdevs2), mean_4ds.append(vec_4d)
            covtrices3.append(covtrix3), mean_4dLEs.append(vec_LE)

        df = table.to_pandas()

        df['covtrix'] = covtrices # lx ly lz covariance matrices
        df['vec_L'] = mean_angulars # array of angular momenta [l1 l2 l3 l4...]
        df['dLx'], df['dLy'], df['dLz'] = np.array(stdevslist).T

        df['vec_4d'] = mean_4ds # array of 4d SSLovdal vectors [L Lz E circ]
        df['covtrix2'] = covtrices2 # covariance matrices for sarah sofie lovdal
        df['L'], df['Lz'], df['E'], df['circ'] = np.array(mean_4ds).T # stored instead as individual values

        df['dL'], df['dLz'], df['dE'], df['dcirc'] = np.array(stdevslist2).T # in order
        df['covtrix3'] = covtrices3 # covariance matrices for [lx ly lz E] format jorge wanted
        df['vec_4dLE'] = mean_4dLEs # list of vectors for the above covariance matrices

        return df
    # Given table with n rows, will generate m multivariate-normal distributed momentum vectors for each row in table.
    # Will return a pandas, with n rows, and m columns: each column is one unique set of momentum vectors.
    # TODO: If you go down and do 2D analysis, you will need to scatter in energy space, too. Hence position.
    # TODO: Update 19/1/22 - 4D Analysis. Need to also monte-up L, Lz, E, Circularity.
    # TODO: Update 27/1/22 - Done. Not necessary, apparently- Jorge doesn't want us to do so.
    def panda_duplimonte(self, panda, m):
        # Create example datasets with LxLyLz
        covtrices = panda['covtrix']
        L_vectors = panda['vec_L']
        list_of_rows = []
        # Generate L components for each row.
        for num,L_vector,covtrix in zip(range(len(covtrices)),L_vectors,covtrices):
            L_generated = random.default_rng().multivariate_normal(mean=L_vector, cov=covtrix, size=m)
            list_of_rows.append(list(L_generated))
        list_of_columns = list(map(list, zip(*list_of_rows)))

        # Create example datasets with L LZ E CIRC
        covtrices2 = panda['covtrix2']
        fourD_vectors = panda['vec_4d']
        list_of_rows_2 = []
        # Generate L components for each row.
        for num,L_vector,covtrix in zip(range(len(covtrices2)),fourD_vectors,covtrices2):
            L_generated = random.default_rng().multivariate_normal(mean=L_vector, cov=covtrix, size=m)
            list_of_rows_2.append(list(L_generated))
        list_of_columns_2 = list(map(list, zip(*list_of_rows)))

        # Create example datasets with LxLyLzE
        covtrices3 = panda['covtrix3']
        LE_vectors = panda['vec_4dLE']
        list_of_rows_3 = []
        # Generate L components for each row.
        for num,L_vector,covtrix in zip(range(len(covtrices3)),LE_vectors,covtrices3):
            L_generated = random.default_rng().multivariate_normal(mean=L_vector, cov=covtrix, size=m)
            list_of_rows_3.append(list(L_generated))
        list_of_columns_3 = list(map(list, zip(*list_of_rows_3)))



        return list_of_columns, list_of_columns_2, list_of_columns_3

# Great Circle Cell Counts in the Galactocentric System, as defined by 1996 Johnston Paper.
# Note: Designed to work in Standard Polar, not Latipolar. ALL IN DEGREES!
# Astropy Tables
class greatcount(object):
    def __init__(self):
        self.null = "null"

    """
    Sagittarius within 273,-13 degrees, maybe 1 degree off (thus 103,273 in theta/phi) 
    dtheta is the HALF WIDTH of the GCC count. HALF WIDTH MATE!!! 
    """
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
        # Also grab the number of stars that satisfy the radcondition while we're at it, not just the ones in the GCC.
        radnumber = 0
        for d in radcondition:
            if d == True:
                radnumber += int(1)
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
        #print(theta, phi, len(indices))
        # Return the newly-built table, alongside the radius conditions number
        return table[indices], radnumber

    # ^ but returns only the indices, no fancier computation.
    def gcc_table_indices(self, table, theta, phi, dtheta, radmin):
        # Set up unit vector for GCC
        theta, phi, dtheta = math.radians(theta), math.radians(phi), math.radians(dtheta)
        polar_unit = np.array([np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)])
        # Get unit vectors for all the targets
        pos = np.array([table['x'], table['y'], table['z']]).T
        radii = np.array([np.linalg.norm(d) for d in pos])
        radcondition = [True if d >= radmin else False for d in radii]
        # Also grab the number of stars that satisfy the radcondition while we're at it, not just the ones in the GCC.
        radnumber = 0
        for d in radcondition:
            if d == True:
                radnumber += int(1)
        pos_unit = np.array([pos[d] / radii[d] for d in range(len(pos))])
        # Dot them
        dotted = [np.dot(polar_unit, pos_unit[d]) for d in range(len(pos_unit))]
        condition = np.sin(dtheta)
        indices = []
        # Check for conditions and build new table.
        for num, item in enumerate(dotted):
            if radcondition[num] == True:
                if abs(item) <= condition:
                    indices.append(num)
        # print(theta, phi, len(indices))
        return indices

    # Get and split a gcc_table into n_phi zones, i.e. [0, 360/n, 2*360/n,... 360] areas.
    def gcc_splittable(self, table, theta, phi, dtheta, radmin, n_phi):
        # Get table (no split)
        table, radnumber = self.gcc_table(table,theta,phi,dtheta,radmin)

        # Generate coordinate system for the great circle for coordinate decomposition
        theta, phi = math.radians(theta), math.radians(phi)
        polar_z = np.array([np.sin(theta) * np.cos(phi),
                            np.sin(theta) * np.sin(phi),
                            np.cos(theta)])
        theta_x, phi_x = theta + (np.pi/2), phi
        if theta_x > np.pi:
            theta_x = 2*np.pi - theta_x
            phi_x += np.pi
        polar_x = np.array([np.sin(theta_x) * np.cos(phi_x),
                            np.sin(theta_x) * np.sin(phi_x),
                            np.cos(theta_x)])
        polar_y = np.cross(polar_z, polar_x)

        # Get xx yy zz coordinates IN THE GREAT CIRCLE SYSTEM for the table, XX,YY,ZZ!
        v = np.array([table['x'],
                      table['y'],
                      table['z']])
        polx, poly = np.matmul(polar_x, v), np.matmul(polar_y, v)

        # Get the value of phi (counter-clockwise) to the new x axis.
        phis = np.arctan2(poly, polx)
        phis = np.array([d + (2*np.pi) if d < 0 else d for d in phis])

        # Get full/half-width of each phi zone
        fw = 2*np.pi/(n_phi)

        # Set up the regions for splitting the phi's up
        phi_regions = fw*np.arange(0, n_phi+1, 1) # ex: nphi = 3, [0, 1, 2, 3], 3 areas: 0->1, 1->2, 2->3.
        phi_regions_indices = [[] for d in range(n_phi)]

        # Get the region that each phi belongs to and save it to an index of indices!
        for num, phi in enumerate(phis):
            for region_index in range(n_phi):
                if phi_regions[region_index] <= phi < phi_regions[region_index + 1]:
                    phi_regions_indices[region_index].append(num)

        # Clip the gcc_table into all the individual subsections
        subtables = [table[indices] for indices in phi_regions_indices]

        # Return the subtables
        return subtables, radnumber

    # Generate a hemispheres-worth of equally separated n-points. RETURNS IN RADIANS.
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
    # Make sure K_full is EVEN. Returns are in RADIANS.
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

    # Within an angular distance (degrees) select K_partial points from theta,phi
    """
    Generate a fullgrid with the correct number of points: 
    - first get point density for the partial grid by getting its K_partial/solid-angle
    - get total number for the entire 4pi-steradian-sphere and use to generate full-grid
    Dot each point in the fullgrid with the vector for this "partial pole" and expect abs(costheta) <= cosdtheta
    """
    def partialgen(self, K_partial, dtheta, theta, phi):
        # Get in radians.
        dtheta, theta, phi = math.radians(dtheta), math.radians(theta), math.radians(phi)

        # Get the fullgrid with the correct number required
        solidangle = 2*np.pi*(1 - math.cos(dtheta))
        sadensity = K_partial/solidangle
        K_full = 4*np.pi*sadensity
        K_full = np.rint(K_full, out=np.zeros(1, int), casting='unsafe')[0]
        if (K_full % 2) != 0:
            K_full += 1
        fullgrid = self.fullgridgen(K_full)

        # Dot and select the coordinates satisfying lying within dtheta of our point.
        unit = np.array([math.sin(theta) * math.cos(phi),
                         math.sin(theta) * math.sin(phi),
                         math.cos(theta)])
        thetaphiindices = list(zip(fullgrid[0], fullgrid[1], fullgrid[2]))
        acceptphetaphiindices = []
        # The dot-product lim
        lim = math.cos(dtheta)
        for dex in thetaphiindices:
            thetaa, phii, indexx = dex
            index_unit = np.array([math.sin(thetaa) * math.cos(phii),
                                   math.sin(thetaa) * math.sin(phii),
                                   math.cos(thetaa)])
            dotted = abs(np.dot(unit, index_unit))
            if dotted >= lim:
                acceptphetaphiindices.append(np.array(dex))

        # Get back the list format and return the list
        thetas, phis, indices = np.array(acceptphetaphiindices).T
        return thetas, phis, indices

    # Generate a set of theta, phi (in galcentricpolar) for a given GCC for scatterplot
    def gcc_gen(self, n_points, theta, phi):
        # Generate coordinate system for the great circle for coordinate decomposition
        theta, phi = math.radians(theta), math.radians(phi)
        polar_z = np.array([np.sin(theta) * np.cos(phi),
                            np.sin(theta) * np.sin(phi),
                            np.cos(theta)])
        theta_x, phi_x = theta + (np.pi/2), phi
        if theta_x > np.pi:
            theta_x = 2*np.pi - theta_x
            phi_x += np.pi
        polar_x = np.array([np.sin(theta_x) * np.cos(phi_x),
                            np.sin(theta_x) * np.sin(phi_x),
                            np.cos(theta_x)])
        polar_y = np.cross(polar_z, polar_x)

        # Set up range of phi (within the GCC coordinate frame)
        phirange = np.linspace(0, 2*np.pi, n_points)

        # Set up all the vectors for the periphery of the GCC (unit vectors.)
        gcc_univects = np.array([polar_x*math.cos(d) + polar_y*math.sin(d) for d in phirange])

        # Get theta and phi for gcc_univects in the non-GCC frame
        # polar_pos = np.array([radius, theta, phi])
        polars = np.array([angular().vec_polar(d) for d in gcc_univects]).T
        thetas, phis = polars[1],polars[2]

        return thetas, phis

# Carry out a clustering (HDBSCAN-oriented).
# Needs array [L1,L2,L3...]
class cluster3d(object):
    def __init__(self):
        self.null = "null"
        self.bounds = np.array([[-8e3, 6e3],
                               [-10e3, 10e3],
                               [-7e3, 7e3]])

    # Remove outlying L-values (units of sigma). Takes either Panda or an Astropy Table via our format.
    """
    This assumes r has already been cleaned: see r_clean
    Note that this will check each axis (x,y,z) for the sig_tolerance requirement
    If even one of them do not satisfy it, it removes that value from table ("strong" requirement.) 
    """
    def L_clean(self, table, sig_tolerance):
        # Grab momenta magnitudes (i.e. absolute values.)
        Ls = np.abs(np.array([(table['Lx']),(table['Ly']),(table['Lz'])]).T) # as an array [[lx ly lz], [l2]...]
        # Grab errors
        dLs = np.array([table['dLx'],table['dLy'],table['dLz']]).T # same format as above.
        # Grab differences (absolute)
        difs = np.abs(Ls - dLs)
        # Ratio them
        difratios = difs / dLs
        # If a dataframe, use table.drop (pandas)
        if type(table) == pandas.DataFrame:
            for num,row in enumerate(table):
                for i in difratios[num]:
                    if i >= sig_tolerance:
                        table.drop(num)
                        break
        # Else use table.remove_row (astropy)
        else:
            for num,row in enumerate(table):
                for i in difratios[num]:
                    if i >= sig_tolerance:
                        table.remove_row(num)
                        break
        return table

    # Remove stars within radius of r
    def r_clean(self, table, minimum_radius):
        # Clean Radius
        for num, row in enumerate(table):
            if row['r'] <= minimum_radius:
                table.remove_row(num)

        return table

    # Return probabilistic selection for bounding box as True or False (a generator object.) Replaces clean.
    # Returns [True, False, True...]: do data[getCuboid(data)] for the cube clip you need. Useful for preprocess.
    def getCuboid(self,L_data):
        coords = L_data
        box_limits = self.bounds
        # Get probs: 1 or 0.
        return (coords[:, 0] > box_limits[0, 0]) & (coords[:, 0] < box_limits[0, 1]) \
               & (coords[:, 1] > box_limits[1, 0]) & (coords[:, 1] < box_limits[1, 1]) \
               & (coords[:, 2] > box_limits[2, 0]) & (coords[:, 2] < box_limits[2, 1])

    # Hierarchical DBS on an array of L-Vectors [1,2,...] provided [minclustsize,minsamples]
    # Assumes data already preprocessed (see: getCuboid.)
    def listhdbs(self, array, minpar):
        # Flat-clustering example (select n_clusters instead.)
        #clusterer = flat.HDBSCAN_flat(min_cluster_size=minpar[0],
        #                              min_samples=minpar[1],
        #                              metric='l2',
        #                              algorithm='best',prediction_data=True,n_clusters=10,X=array)
        #clustered = flat.approximate_predict_flat(clusterer, array, 10)
        # return np.array(clustered[0])
        hdbs = hdbscan.HDBSCAN(min_cluster_size=minpar[0],
                               min_samples=minpar[1],
                               metric='l2',
                               algorithm='best')
        hdbsfit = hdbs.fit_predict(array)
        return np.array(hdbsfit)


    # Hierarchical DBS: returns memberships. Astropy Tables.
    # Deprecated basically (designed for kmeans_L/xmeans_L graph)
    # Keeping it around in case it needs to be revived for use. 
    def hdbs(self, table, browser, minpar):
        table = self.clean(table, 5)
        L = np.array([table['Lx'], table['Ly'], table['Lz']]).T
        hdbs = hdbscan.HDBSCAN(min_cluster_size=minpar[0],
                               min_samples=minpar[1],
                               metric='l2',
                               algorithm='best')
        hdbsfit = hdbs.fit_predict(L)
        table['k_index'] = np.array(hdbsfit)
        #save_format = ("HDBS_TEST_EPS{0}_MINSAMP{1}" + ".html").format(eps, min_samples)
        graphutils.threed_graph().kmeans_L(table, False, browser)
        #graphutils.threed_graph().xmeans_L(table, "test_hdbs.html", browser)
        return table


# Create pseudo-random 3D multivariates for cluster testing,
# bounded by cluster3d() bounding cuboid,
# reliant on normal dists.
class genclust3d(object):
    def __init__(self):
        self.null = "null"
        self.rng = np.random.default_rng() #[[minx,miny,minz],[max...]]

    # Randomly generate (n) mu vectors/symmetric covariance matrices in 3 dimensions.
    """
    All elements generated using uniform distributions with maxmin of covbounds=sigma for standard deviation/covtrix.
    Optionally specify covscale on a per-cluster basis as a list [cs1, cs2, cs3...]
    mubounds describes points spatially of mu: just like with cluster3d() bounding box: centroids in here
    [[x0,x1],
     [y0,y1],
     [z0,z1]]
    """
    def gen_mucov(self, n, covscale, mubounds):
        # Get mu's
        mus = self.rng.integers(low=mubounds[0],high=mubounds[1],size=(n,1,3))
        mus = [d[0] for d in mus]
        mus = np.array(mus)

        # Get cov's
        covtrices = []
        for m in range(n):
            scale = covscale
            if type(covscale) == list:
                scale = covscale[m]
            # Create and populate diagonals randomly
            cov_diag = self.rng.uniform(low=-1*scale,high=scale, size=3)
            diagtrix = np.zeros(shape=(3,3))
            for i in range(3):
                diagtrix[i,i] = cov_diag[i]*cov_diag[i]
            # Now create/populate off-diagonals.
            offdiagtrix = np.zeros(shape=(3,3))
            for i in range(1,3):
                for j in range(0,2):
                    if i != j:
                        randi, randj = self.rng.uniform(low=-1*scale,high=scale, size=2)
                        offdiagtrix[i,j] = randi*randj
                        offdiagtrix[j,i] = randi*randj
            # Full covariance matrix
            covtrix = diagtrix + offdiagtrix
            covtrices.append(covtrix)
        covtrices=np.array(covtrices)
        # Return
        return mus, covtrices

    # Visualize the clusters we have created, with N (static)-points per cluster for a convex hull about data.
    def mucov_visualize(self, mucovs,N):
        graphutils.threed_graph().kmeans_L_multinormal_generated(mucovs,N)

    # Generate data using gen_mucov mu/cov
    """
    with n_points a list/tuple for specificity to each cluster.
    covscale is a static value for normal noise to all datapoints/creating per-point covariance matrices
    (points are more likely to have smaller than larger noise, though heavily noised points exist.)
    Set "individual=True" to return list individual data arrays. 
    """
    def noisy_data(self, mucov, n_points, noisescale, individual=False):
        # Define list placeholders, etc
        mus, covtrices = mucov
        pointslist = []
        covslist = []

        # Set up data for each mu/covtrix (we're not trying to be fancy here- just somewhat dirty.
        for num, mu, cov, npoint in zip(range(len(mus)),mus,covtrices,n_points):
            points = self.rng.multivariate_normal(mean=mu,cov=cov,size=npoint)
            covs = []
            for m in range(npoint):
                # Create and populate diagonals randomly
                cov_diag = self.rng.normal(loc=0, scale=noisescale, size=3)
                diagtrix = np.zeros(shape=(3, 3))
                for i in range(3):
                    diagtrix[i, i] = cov_diag[i] * cov_diag[i]
                # Now create/populate off-diagonals.
                offdiagtrix = np.zeros(shape=(3, 3))
                for i in range(1, 3):
                    for j in range(0, 2):
                        if i != j:
                            randi, randj = self.rng.normal(loc=0, scale=noisescale, size=2)
                            offdiagtrix[i, j] = randi * randj
                            offdiagtrix[j, i] = randi * randj
                # Full covariance matrix
                covtrix = diagtrix + offdiagtrix
                covs.append(covtrix)
            pointslist.append(points)
            covslist.append(covs)

        # Decide returns
        if individual == False:
            return np.concatenate(pointslist), np.concatenate(covslist)
        else:
            return pointslist, covslist

# Comparison Metrics for Clusters, specific to our data structure.
# Uses lists of cluster indices.
class compclust(object):
    def __init__(self):
        null = "null"

    # Get number of clusts (including outlier clustering, -1)
    # This assumes that clusters are linearly ordered [-1,0,1,2,3,4...]
    def nclust_get(self, clust):
        return np.max(clust) + 2

    # More complete nclust_get, albeit slower.
    def nclust_get_complete(self, clust):
        clust_list = []
        nclust = 0
        for clust_id in clust:
            if clust_id in clust_list:
                pass
            else:
                clust_list.append(clust_id)
                nclust += 1
        return nclust, clust_list

    # For a group, get the entirety of the set and get the number of clusters per clustering.
    # Also generates a PDF for the set of a given clustering having n_clust.
    def graph_nclust(self, group, energy=False, sofie=False):
        # Decide if using LxLyLzE or just LxLyLz
        if energy == True:
            # Load in the lists as a list of lists!
            arrays = [windows_directories.duplimontedir + "\\" + group + "\\" + d + ".cluster.txt" \
                      for d in ascii_info.duplimonte_LE_saveids]
        elif sofie == True:
            # Load in the lists as a list of lists!
            arrays = [windows_directories.duplimontedir + "\\" + group + "\\" + d + ".cluster.txt" \
                      for d in ascii_info.duplimonte_L4D_saveids]
        else:
            # Load in the lists as a list of lists!
            arrays = [windows_directories.duplimontedir + "\\" + group + "\\" + d + ".cluster.txt" \
                      for d in ascii_info.duplimonte_saveids]

        # Load in the arrays
        for num, array in enumerate(arrays):
            with open(array, 'rb') as f:
                arrays[num] = pickle.load(f)
        # Get the num
        arrays = [self.nclust_get(d) for d in arrays]
        # Pass this to graphutils plotting to generate a plot
        graphutils.twod_graph().nclust_n(arrays, group)

    # Given two clusterings (clust1,clust2), map clust2 onto clust1 and return remapped clust 2 (same original format.)
    # This only works for NxN cardinality (equal set size) and leverages the hungarian algorithm.
    """
    Thank Pallie/Aaron!
    https://stackoverflow.com/questions/55258457/
    find-mapping-that-translates-one-list-of-clusters-
    to-another-in-python
    """
    def compclust(self, clust1, clust2):
        contMatrix = contingency_matrix(clust1, clust2)
        labelMatcher = munkres.Munkres()
        labelTranlater = labelMatcher.compute(contMatrix.max() - contMatrix)
        uniqueLabels1 = list(set(clust1))
        uniqueLabels2 = list(set(clust2))
        tranlatorDict = {}
        for thisPair in labelTranlater:
            print(thisPair)
            tranlatorDict[uniqueLabels2[thisPair[1]]] = uniqueLabels1[thisPair[0]]

        return np.array([tranlatorDict[label] for label in clust2])

    # Given two clusterings, do the above, but in rectangular fashion (i.e. multiple matches to each.)
    # Thanks Sascha for verbatim solution (we're modifying it, though.)
    # https://stackoverflow.com/questions/50305614/finding-an-optimal-selection-in-a-2d-matrix-with-given-constrains
    def compclust_multi(self, clust1, clust2):# Generate cost/etc matrices and do the map for lowest set cardinality
        contmat = contingency_matrix(clust1, clust2)
        max_cost = np.amax(contmat)
        harvest_profit = max_cost - contmat
        row_ind, col_ind = linear_sum_assignment(harvest_profit)
        sol_map = np.zeros(contmat.shape, dtype=bool)
        sol_map[row_ind, col_ind] = True

        # Visualize plots (thanks Sascha!)
        f, ax = plt.subplots(2, figsize=(9, 6))
        sns.heatmap(contmat, annot=True, linewidths=.5, ax=ax[0], cbar=False,
                    linecolor='black', cmap="YlGnBu")
        sns.heatmap(contmat, annot=True, mask=~sol_map, linewidths=.5, ax=ax[1],
                    linecolor='black', cbar=False, cmap="YlGnBu")
        plt.tight_layout()
        plt.show()

        # Remap clust 2
        row_ind, col_ind = row_ind - 1, col_ind - 1
        newclust = np.zeros(shape=np.shape(clust2), dtype=int)
        for thisPair in zip(row_ind, col_ind):
            for num,clust_id in enumerate(clust2):
                if clust_id == thisPair[1]:
                    newclust[num] = thisPair[0]
        return newclust


    # The pre-final rendition:
    """
    Compute cluster match via hungarian method
    If the 2nd clustering has say 9 clusterings while the first has 8, take the "least quality" match and relabel
    The new label will be subject to an index limit (i.e. a "minimum value" to avoid mixing labels from multiple.) 
    """
    def compclust_multilabel(self, clust1, clust2, minimum_excess_index):
        # Hungarian Algorithm Linear Sum Assignment (thanks Stack Exchange)
        contmat = contingency_matrix(clust1, clust2)
        max_cost = np.amax(contmat)
        harvest_profit = max_cost - contmat
        row_ind, col_ind = linear_sum_assignment(harvest_profit)
        sol_map = np.zeros(contmat.shape, dtype=bool)
        sol_map[row_ind, col_ind] = True

        # Visualize plots (thanks Sascha!)
        f, ax = plt.subplots(2, figsize=(9, 6))
        sns.heatmap(contmat, annot=True, linewidths=.5, ax=ax[0], cbar=False,
                    linecolor='black', cmap="YlGnBu")
        sns.heatmap(contmat, annot=True, mask=~sol_map, linewidths=.5, ax=ax[1],
                    linecolor='black', cbar=False, cmap="YlGnBu")
        plt.tight_layout()
        plt.show()

        # Row is "clust1" and col is "clust2" and adjust for the fact that -1 is now 0
        # This is our mapping for the columns (clust2) to the rows (clust1)
        row_ind, col_ind = row_ind - 1, col_ind - 1

        # Generate the new relabelled clustering
        newclust = np.zeros(shape=np.shape(clust2), dtype=int)
        for thisPair in zip(row_ind, col_ind):
            for num,clust_id in enumerate(clust2):
                if clust_id == thisPair[1]:
                    newclust[num] = thisPair[0]

        # Identify which clustering is the one with a higher number of clusters
        """
        If the 2nd cluster is np.max higher than the first, then we have a problem
        This is "col_ind"
        Find out if the max of this is larger for 2 than 1
        If it is, then find the elements in this list not mapped over to the first
        """
        if np.max(clust2) > np.max(clust1):
            # Generate a list of all unique indices for the original clust2
            uniques = []
            for cluster in clust2:
                if cluster not in uniques:
                    uniques.append(cluster)

            # Find the indices that weren't paired by the Hungarian Algorithm (i.e. not remapped)
            not_remapped = []
            for unique in uniques:
                if unique not in col_ind:
                    not_remapped.append(unique)

            # Generate replacement indices (subject to a "minimum_excess_index")
            """
            this index exists to let all "excess clusters" remain unique
            minimum_excess_index is the minimum replacement!!!
            for example, if clust2 is 8 and clust1 is 7, then a wise minimum would be 8:
            the 8th cluster is labelled as an "excess cluster"
            If you have two clusterings to compare to the base clustering, and both have 8
            then comparing 1 to 0 should have a minimum of 8
            comparing 2 to 0 should have a minimum of 9
            this keeps the excess of either of them independent (i.e. the excess in 1-0 is 8, 2-0 is 9)
            then you can individually compare cluster 8 and 9 to see if they are the same or different
            anyway, just note that you should carefully pick the minimum_excess
            and make sure that no excess clusterings between different Monte-carlo's are labelled the same
            """
            not_in_indices = []
            for num,not_index in enumerate(not_remapped):
                not_in_indices.append(minimum_excess_index + num)

            # Go through this relabelled clustering and set all "not_in" clusters to min_new_label or higher
            # This means that if clust1 is np.max of 8, and clust2 is 10, you want min_new_label set to be set to 9.
            for mun,not_in_index in enumerate(not_remapped):
                for num,clust_id in enumerate(newclust):
                    if clust_id == not_in_index:
                        newclust[num] = not_in_indices[mun]

        return newclust