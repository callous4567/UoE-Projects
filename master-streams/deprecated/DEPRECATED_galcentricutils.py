# Various old disused codes from galcentricutils.
def GALCENT_to_ICRS(self, hdfdir, hdfname, group, set):
    # Set up utilities environment for radec conversion
    utilities = utils()
    # Set up HDF and grab table, and SkyCoord objects for all targets.
    writer = hdfutils.hdf5_writer(hdfdir, hdfname)
    table = writer.read_table(group, set)
    skycoords = coord.SkyCoord(x=table['x'] * u.kpc, y=table['y'] * u.kpc, z=table['z'] * u.kpc,
                               v_x=table['vx'] * u.km/u.s, v_y=table['vy']*u.km/u.s, v_z=table['vz']*u.km/u.s,
                               frame="galactocentric")
    # Effect conversion to ICRS, work through objects, collect converted quantities.
    icrs_skycoords = skycoords.transform_to(coord.ICRS)
    ra_list, dec_list, pmra_list, pmdec_list, distance_list, radial_velocity_list = [],[],[],[],[],[]
    textradeclist = []
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
        pmra = pmra_cosdec / math.cos(math.radians(dec))
        # Grab text values
        radec = [ra,dec]
        textradec = utilities.deg_radec(radec)
        textradeclist.append(textradec)
        ra_list.append(ra), dec_list.append(dec), pmra_list.append(pmra), pmdec_list.append(pmdec), \
        distance_list.append(distance), radial_velocity_list.append(radial_velocity)

    # Modify and save table.
    table['textradec'] = np.array(textradeclist)
    table['ra'] = ra_list
    table['dec'] = dec_list
    table['pmra'] = pmra_list
    table['pmdec'] = pmdec_list
    table['dist'] = distance_list
    table['vlos'] = radial_velocity_list
    writer.write_table(group, set, table)
def ICRS_to_GALCENT(self, hdfdir, hdfname, group, set):
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
    galcent_skycoords = skycoords.transform_to(coord.Galactocentric)
    x_list, y_list, z_list, vx_list, vy_list, vz_list = [], [], [], [], [], []
    for object in galcent_skycoords:
        x,y,z,vx,vy,vz = object.x/u.kpc, object.y/u.kpc, object.z/u.kpc, \
                         object.v_x/(u.km/u.s), object.v_y/(u.km/u.s), object.v_z/(u.km/u.s)
        # Discard the dimensionless unit.
        x,y,z,vx,vy,vz = x.value,y.value,z.value,vx.value,vy.value,vz.value

        # Append to list
        x_list.append(x), y_list.append(y), z_list.append(z), \
        vx_list.append(vx), vy_list.append(vy), vz_list.append(vz)

    # Modify and save table.
    table['x'],table['y'],table['z'],table['vx'],table['vy'],table['vz'] = x_list,y_list,z_list,\
                                                                           vx_list,vy_list,vz_list
    writer.write_table(group, set, table)
def vec_ICRS_to_GALCENT(self, ra,dec,distance,dmura,dmudec,vlos):
    vecskycoord = coord.SkyCoord(ra=ra*u.deg,
                                 dec=dec*u.deg,
                                 distance=distance*u.kpc,
                                 pm_ra_cosdec=dmura*np.cos(np.radians(dec))*u.mas/u.yr,
                                 pm_dec=dmudec*u.mas/u.yr,
                                 radial_velocity=vlos*u.km/u.s,
                                 frame="icrs")
    vecskycoord = vecskycoord.transform_to(coord.Galactocentric)
    x,y,z,vx,vy,vz = vecskycoord.x/u.kpc, vecskycoord.y/u.kpc, vecskycoord.z/u.kpc, \
                     vecskycoord.v_x/(u.km/u.s), vecskycoord.v_y/(u.km/u.s), vecskycoord.v_z/(u.km/u.s)
    x,y,z,vx,vy,vz = x.value,y.value,z.value,vx.value,vy.value,vz.value
    return [x,y,z,vx,vy,vz]


