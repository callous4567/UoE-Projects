import math
import hdfutils
import os
import numpy as np
import astropy.coordinates as coord
import astropy.units as u
# See galcentricutils.py for analogies. This is just for Galactic Coordinates, instead.

# To Galactic and Back (via ICRS exclusively.)
class galconversion(object):
    def __init__(self, sourcedir):
        self.source, self.owd, self.sol_params = sourcedir, os.getcwd(), np.zeros((1,4))

    # Note: requires l,b,distance,dmu_l,dmu_b,vlos. Assumes conventions of original ASCII data colnames.
    def to_ICRS(self, hdfdir, hdfname, group, set):
        # Set up HDF and grab table, and SkyCoord objects for all targets.
        writer = hdfutils.hdf5_writer(hdfdir, hdfname)
        table = writer.read_table(group, set)
        skycoords = coord.SkyCoord(l=table['l'] * u.deg, b=table['b'] * u.deg, distance=table['distance'] * u.kpc,
                                   pm_l_cosb=table['dmu_l']*math.cos(math.radians(table['b'])) * u.mas/u.yr, pm_b=table['dmu_b']*u.mas/u.yr, radial_velocity=table['vlos']*u.km/u.s,
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
        writer.write_table(group, set, table)

    # Converts ICRS to Galactic instead.
    def to_GALACTIC(self, hdfdir, hdfname, group, set):
        # Set up HDF and grab table, and SkyCoord objects for all targets.
        writer = hdfutils.hdf5_writer(hdfdir, hdfname)
        table = writer.read_table(group, set)
        skycoords = coord.SkyCoord(ra=table['ra']*u.deg,
                                   dec=table['dec']*u.deg,
                                   distance=table['distance']*u.kpc,
                                   pm_ra_cosdec=table['pmra']*math.cos(table['dec'])*u.mas/u.yr,
                                   pm_dec=table['pmdec']*u.mas/u.yr,
                                   radial_velocity=table['vlos']*u.km/u.s,
                                   frame="icrs")
        # Effect conversion to Galactocentric, work through objects, collect converted quantities.
        galcent_skycoords = skycoords.transform_to(coord.Galactic)
        l_list,b_list,dmu_l_list,dmu_b_list = [],[],[],[]
        for object in galcent_skycoords:
            l,b,dmub = object.l/u.deg, object.b/u.deg, object.pm_b/(u.mas/u.yr)
            l,b,dmub = l.value,b.value,dmub.value
            dmul = ((object.pm_l_cosb/(u.mas/u.yr)).value)/math.cos(math.radians(b))

            l_list.append(l),b_list.append(b),dmu_b_list.append(dmub),dmu_l_list.append(dmul)

        # Modify and save table.
        table['l'],table['b'],table['dmu_l'],table['dmu_b'] = l_list,b_list,dmu_l_list,dmu_b_list
        writer.write_table(group, set, table)
