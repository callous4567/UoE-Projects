# Various old disused codes from galcentricutils_new.
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

# Clustering
class cluster3d(object):
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
        graphutils_new.threed_graph().kmeans_L(table, savedex + ".html", browser)
        return table, inertia

    # DBSCAN. params are "eps" and "min_samples." Browser=True opens in browser (append to others.)
    def dbs(self, table, eps, min_samples, browser):
        table = self.clean(table, 5)
        L = np.array([table['Lx'], table['Ly'], table['Lz']]).T
        dbs = DBSCAN(eps=eps, min_samples=min_samples, metric="l1", leaf_size=5)
        dbsfit = dbs.fit_predict(L)
        table['k_index'] = np.array(dbsfit)
        save_format = ("DBS_TEST_EPS{0}_MINSAMP{1}" + ".html").format(eps, min_samples)
        graphutils_new.threed_graph().kmeans_L(table, save_format, browser)
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
        graphutils_new.threed_graph().kmeans_L(table, "test_hdbs.html", browser)
        graphutils_new.threed_graph().xmeans_L(table, "test_hdbs.html", browser)
        return table

    # OPTICS
    def optics(self, table):
        table = self.clean(table, 5)
        L = np.array([table['Lx'], table['Ly'], table['Lz']]).T
        optics = OPTICS(metric="l1", min_cluster_size=20, leaf_size=20, eps=1e4, max_eps=1e6)
        opticsfit = optics.fit_predict(L)
        table['k_index'] = np.array(opticsfit)
        graphutils_new.threed_graph().kmeans_L(table, "test_plot", False)
        return table

    # Agglomerative
    def aglom(self, table, k):
        table = self.clean(table,5)
        L = np.array([table['Lx'], table['Ly'], table['Lz']]).T
        aglom = AgglomerativeClustering(n_clusters=k, linkage='ward')
        aglomfit = aglom.fit_predict(L)
        table['k_index'] = np.array(aglomfit)
        graphutils_new.threed_graph().kmeans_L(table, "test_plot", False)
        return table

    # Affinity Propagation
    def afprop(self, table):
        table = self.clean(table, 5)
        L = np.array([table['Lx'], table['Ly'], table['Lz']]).T
        afprop = AffinityPropagation(damping=0.5)
        afpropfit = afprop.fit_predict(L)
        table['k_index'] = np.array(afpropfit)
        graphutils_new.threed_graph().kmeans_L(table, False, True)
        return table

    # Gaussian Mixture with Variational Bayes
    def bayesgaussmix(self, table, k, savedex, browser):
        # Set up vectors/positions/etc
        table = self.clean(table, 4)
        L = np.array([table['Lx'], table['Ly'], table['Lz']]).T
        gm = BayesianGaussianMixture(n_components=k)
        gmfit = gm.fit_predict(L) # list with indices for cluster
        table['k_index'] = np.array(gmfit)
        graphutils_new.threed_graph().kmeans_L(table, savedex + ".html", browser)
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
            graphutils_new.threed_graph().kmeans_L(table, savedex + ".html", browser)
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
            os.mkdir(windows_directories_new.imgdir + "\\gauss_aics_bics_test")
        except:
            pass

        plt.savefig(windows_directories_new.imgdir + "\\gauss_aics_bics_test\\" + savedex + ".png")
        plt.clf()
        plt.close()
