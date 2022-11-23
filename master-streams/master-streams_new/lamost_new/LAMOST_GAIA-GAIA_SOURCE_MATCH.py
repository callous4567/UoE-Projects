import os
import pickle
import numba
from numba import njit
import numpy as np
from astropy.io import fits
from astropy.table import Table, vstack
import ascii_info_new
import windows_directories_new

# TODO NOTE: general catalogue FITS broken- use CSV (need to process.)
# TODO NOTE: dr7_v2.0_LRS_stellar.fits works fine as a FITS straight out of their website!
# TODO NOTE: the value-added-catalogue is a DAT file, too. Double oof.
process_vat = False
raw_dir = os.path.join(windows_directories_new.lamodir, "dr7_lamost_LR_release.dat")
stellar_dir = os.path.join(windows_directories_new.lamodir, "dr7_v2.0_LRS_stellar.fits")
vat_dir = os.path.join(windows_directories_new.lamodir, "dr7_LRS_vat.fits")
vatstellarjoin_dir = os.path.join(windows_directories_new.lamodir, "dr7_LRS_vat-stellar-join.fits")
dr7_gaia_gaia_dir = os.path.join(windows_directories_new.lamodir, "mag0.3_mas3000_dr7_GAIA_aidists.fits")
dr7_2mass_gaia_dir = dr7_gaia_gaia_dir  # for now compatibility
if __name__ == "__main__" and process_vat == True:

    # Path to the Value-Added-Catalogue for DR7 that has the distances + sources/etc.
    # http://dr7.lamost.org/doc/vac
    path = raw_dir
    with open(os.path.join(windows_directories_new.lamodir, path), "r") as f:
        lines = f.readlines()
    lines = [d.replace("\n", "") for d in lines]
    lines = [d.replace("#", "") for d in lines]
    lines = [d.split() for d in lines]
    columns = lines[0]
    columns = [d.replace("(", "_") for d in columns]
    columns = [d.replace(")", "_") for d in columns]
    columns = [d.strip() for d in columns]
    columns = columns[1:]
    name = [line[0] for line in lines[1:]]
    name = np.array(name, str)
    data = [line[1:] for line in lines[1:]]
    data = np.array([np.array(d, float) for d in data]).T
    dr7_LRS_vat = Table()
    dr7_LRS_vat["Name_2MASS"] = name
    for num, column in enumerate(columns):
        dr7_LRS_vat[column] = data[num]

    # Save it as a fits.
    dr7_LRS_vat.write(vat_dir, overwrite=True)

    # Generate the join to stellar too.
    with fits.open(stellar_dir, mode='readonly', memmap=True) as hdul:
        # "gaia_source_id" is dr2_source_id inside this.
        dr7_LRS_stellar = Table(hdul[1].data)
        dr7_LRS_stellar.rename_column("obsid", "objID")

    # Create the join
    from astropy.table import join
    dr7_LRS_stellar = join(dr7_LRS_vat, dr7_LRS_stellar, keys=["objID"])
    dr7_LRS_stellar.write(vatstellarjoin_dir, overwrite=True)

# Should I crossmatch the LAMOST dataset with Gaia DR3 using Vizier, obtaining all the columns of interest?
LAMOST_GAIA_GAIA = True
vizier_comparison = False
debug = True
if __name__ == "__main__" and LAMOST_GAIA_GAIA == True:

    # Gaia columns that we are fishing for. Note that objID is turned to objid thanks to ADQL HATING CAPITAL LETTERS!!!
    gaia_COLS_OF_INTEREST = [
        'objid',
        'source_id',
        'ra',
        'ra_obs',
        'dec',
        'dec_obs',
        'phot_g_mean_mag',
        'g_obs',
        'ra_error',
        'dec_error',
        'parallax',
        'parallax_error',
        'pmra',
        'pmra_error',
        'pmdec',
        'pmdec_error',
        'ra_dec_corr',
        'ra_pmra_corr',
        'ra_pmdec_corr',
        'dec_pmra_corr',
        'dec_pmdec_corr',
        'pmra_pmdec_corr',
        'ruwe',
    ]

    if vizier_comparison == True: # get a vizier match for comparison- note that it may fail if Vizier is busy a lot.

        # On top of all that, do a flat match using Vizier for comparison just using ra/dec to compare
        print("Now running an xMatch to Vizier for comparison...")

        with fits.open(vatstellarjoin_dir, mode='readonly', memmap=True) as hdul:
            # "gaia_source_id" is dr2_source_id inside this.
            dr7_LRS_vatstellar = Table(hdul[1].data)[['Ra', 'Dec', 'objID', 'mag3']]

        import gaia_utils_standalone
        vizier_xmatch = gaia_utils_standalone.xmatch_gaiadr3(dr7_LRS_vatstellar,
                                                             'Ra',
                                                             'Dec',
                                                             'objID',
                                                             3000,
                                                             ['angDist', 'RAdeg', 'DEdeg', 'Ra', 'Dec', 'objID',
                                                              'Source', 'Gmag', 'mag3'])

        del dr7_LRS_vatstellar

        vizier_xmatch.rename_columns(['angDist', 'RAdeg', 'DEdeg', 'Ra', 'Dec', 'Source', 'Gmag'],
                                     ['xmatch_angdist', 'ra', 'dec', 'ra_obs', 'dec_obs', 'gaiadr3_source_id',
                                      'phot_g_mean_mag'])
        vizier_xmatch['xmatch_magdist'] = np.abs(vizier_xmatch['mag3']-vizier_xmatch['phot_g_mean_mag'])
        vizier_xmatch.write(os.path.join(windows_directories_new.lamodir, "vizier_test_xmatch.fits"), overwrite=True)
        print("The length of the Vizier equivalent is... ", len(vizier_xmatch))
        del vizier_xmatch

    if debug == False:
        # Now do the regular match with all our extra stuff
        with fits.open(vatstellarjoin_dir, mode='readonly', memmap=True) as hdul:
            # "gaia_source_id" is dr2_source_id inside this.
            dr7_LRS_vatstellar = Table(hdul[1].data)

        # Condense catalogue for the sake of matching
        dr7_crosscat = dr7_LRS_vatstellar[["Ra", "Dec", "objID", "mag3"]]

        del dr7_LRS_vatstellar

        import gaia_utils_standalone

        valid_g_indices = np.where(np.abs(dr7_crosscat['mag3']) < 25)[0] # indices w/ a mag3/lamost-r/~Gaia G magnitude
        nog_indices = np.where(np.abs(dr7_crosscat['mag3']) >= 25)[0] # indices w/o

        # Get the best possible quality matches (to within 3 arcseconds, 0.3 magnitudes)
        matched_by_g = gaia_utils_standalone.xmatch_gaiadr3_conesearch_gmatch(dr7_crosscat[valid_g_indices],
                                                                              gaia_COLS_OF_INTEREST,
                                                                              3000,
                                                                              0.3,
                                                                              'Ra', 'Dec', 'mag3',
                                                                              True,
                                                                              _USE_LITE=False)
        matched_by_g['gmatch_flag'] = 1
        matched_by_g.write("QUICKCACHE_MATCHBYG.fits", overwrite=True)

    else: # debug is true

        # Now do the regular match with all our extra stuff
        with fits.open(vatstellarjoin_dir, mode='readonly', memmap=True) as hdul:
            # "gaia_source_id" is dr2_source_id inside this.
            dr7_LRS_vatstellar = Table(hdul[1].data)

        # Condense catalogue for the sake of matching
        dr7_crosscat = dr7_LRS_vatstellar[["Ra", "Dec", "objID", "mag3"]]

        valid_g_indices = np.where(np.abs(dr7_crosscat['mag3']) < 25)[0] # indices w/ a mag3/lamost-r/~Gaia G magnitude
        nog_indices = np.where(np.abs(dr7_crosscat['mag3']) >= 25)[0] # indices w/o

        del dr7_LRS_vatstellar

        with fits.open("QUICKCACHE_MATCHBYG.fits", mode='readonly', memmap=True) as hdul:
            # "gaia_source_id" is dr2_source_id inside this.
            matched_by_g = Table(hdul[1].data)

    print("Gmatches done #1")

    # Deprecated/too slow.
    #@njit(fastmath=True)
    #def find_unmatched_by_g(base_ids, matched_ids):
    #    """
    #    Return INDICES of elements from base_ids not in matched_ids: comparison by integer value not by index.
    #    :param base_ids: array 1d
    #    :param matched_ids: array 1d
    #    :return: array 1d integers
    #    """
    #    where_in_base_unmatched_bool = np.zeros(np.shape(base_ids), dtype=numba.types.bool_)
    #    for num,i in enumerate(base_ids):
    #        if i not in matched_ids:
    #            where_in_base_unmatched_bool[num] = True
    #    return np.where(where_in_base_unmatched_bool==True)[0]

    print("Now finding unmatched_by_g...")

    import gaia_utils_standalone
    unmatched_by_g = gaia_utils_standalone.jl_setdiff1dbyindices(
        np.array(dr7_crosscat['objID'][valid_g_indices], np.int64),
        np.array(matched_by_g['objid'], np.int64)
    )

    unmatched_by_g = [np.where(d == dr7_crosscat['objID'][valid_g_indices])[0][0] for d in unmatched_by_g]
    print(unmatched_by_g)
    unmatched_by_g = dr7_crosscat[valid_g_indices[unmatched_by_g]]

    print("Found unmatched by g!")

    unmatched_by_g = gaia_utils_standalone.xmatch_gaiadr3_conesearch_gmatch(unmatched_by_g,
                                                                            gaia_COLS_OF_INTEREST,
                                                                            3000,
                                                                            _RA_COLNAME='Ra',
                                                                            _DEC_COLNAME='Dec',
                                                                            _G_COLNAME='mag3',
                                                                            _USE_LITE=False)
    unmatched_by_g['gmatch_flag'] = 0

    # match not considering G-band magnitude similarity
    print("Now matching without g...")
    match_with_nog = gaia_utils_standalone.xmatch_gaiadr3_conesearch_gmatch(dr7_crosscat[nog_indices],
                                                                            gaia_COLS_OF_INTEREST,
                                                                            3000,
                                                                            _RA_COLNAME='Ra',
                                                                            _DEC_COLNAME='Dec',
                                                                            _G_COLNAME='mag3',
                                                                            _USE_LITE=False)
    match_with_nog['gmatch_flag'] = np.NaN

    # stack them up

    print("Vstacking!")
    dr7_crosscat = vstack([matched_by_g, unmatched_by_g, match_with_nog])

    dr7_crosscat.write("QUICKSAVE_TEST.fits", overwrite=True)
    dr7_crosscat.rename_column('source_id', 'gaiadr3_source_id')

    # evaluate distance for the match using small-separation formulae https://en.wikipedia.org/wiki/Haversine_formula
    @njit(fastmath=True)
    def haversine(lon1, lon2, lat1, lat2):

        """
        Compute the haversine angular distance in latitude-longitude coordinates. Input must be in units of radians.
        Note: the haversine distance formula breaks down for points that are antipodal to each-other... be warned!
        :param lon1:
        :param lon2:
        :param lat1:
        :param lat2:
        :return:
        """

        return 2*np.arcsin(
            np.sqrt(
                (np.sin((lat2 - lat1)/2))**2 + (np.cos(lat2))*(np.cos(lat1))*((np.sin((lon2-lon1)/2))**2)
            )
        )
    print("Haversining!")
    dr7_crosscat['xmatch_angdist'] = haversine(
        np.deg2rad(dr7_crosscat['ra']),
        np.deg2rad(dr7_crosscat['ra_obs']),
        np.deg2rad(dr7_crosscat['dec']),
        np.deg2rad(dr7_crosscat['dec_obs']),
    )
    mas = (2*np.pi / 360)*(1/3600)*(1/1000)
    dr7_crosscat['xmatch_angdist'] /= mas
    dr7_crosscat['xmatch_magdist'] = np.abs(dr7_crosscat['phot_g_mean_mag'] - dr7_crosscat['g_obs'])

    # Rename back objID
    dr7_crosscat.rename_column('objid', 'objID')

    del matched_by_g
    del unmatched_by_g
    del match_with_nog

    # All the queries are done and we have xmatch angdists. Now doing all the joins.
    print("All the queries are done and we have xmatch angdists. Now doing all the joins.")

    # Select columns of interest from vatstellar that we can make use of (including for the join.)
    vatstellar_COLS_OF_INTEREST = [
        "d_kpc_",
        "e_d_kpc_",
        "rv",
        "rv_err",
        "feh",
        "feh_err",
        "objID"
    ]
    vatstellar_RENAMED_C_O_I = [
        'dist',
        'edist',
        'vlos',
        'evlost',
        'feh',
        'efeh',
        'objID'
    ]

    print("Now running the join to vatstellar...")
    # get back vatstellar
    with fits.open(vatstellarjoin_dir, mode='readonly', memmap=True) as hdul:
        # "gaia_source_id" is dr2_source_id inside this.
        dr7_LRS_vatstellar = Table(hdul[1].data)[vatstellar_COLS_OF_INTEREST]
        dr7_LRS_vatstellar.rename_columns(vatstellar_COLS_OF_INTEREST, vatstellar_RENAMED_C_O_I)

    # Join it to the master on "objID"
    from astropy.table import join
    dr7_crosscat = join(dr7_crosscat, dr7_LRS_vatstellar, keys="objID")
    del dr7_LRS_vatstellar

    # Write it to disk
    total_rows_found = len(dr7_crosscat)
    print("The total number of rows found in the crosscat from Gaia is... ", total_rows_found)
    dr7_crosscat.write(dr7_gaia_gaia_dir, overwrite=True)

# Section I of the TODO_LIST.txt
SECTION_I = False
if __name__ == "__main__" and SECTION_I == True:


    with fits.open(dr7_gaia_gaia_dir, mode='readonly', memmap=True) as hdul:
        data = Table(hdul[1].data)  # LAMOST-GAIA catalogue to 2 arcsecond crossmatch
        data['source'] = "LAMOKGS_straszaks" # just for the sake of keeping track

    # Get only the K-type stars as a subset... remember each entry is "TYPE" + "0->9" heat/temperature
    # TODO: No matches within 3 arcseconds. Expand to all? Doesn't help.
    #LAMOKGs_truefalse = [True if d == "K" else False for d in [d[0] for d in data['subclass']]]
    #LAMOKGs = data[LAMOKGs_truefalse]


    # Get the stars from the original Jorge dataset (easy function provided)
    jorgeta_cachepath = os.path.join(windows_directories_new.lamodir, "jorgetacache.txt")
    try:
        with open(jorgeta_cachepath, "rb") as f:
            jorgeta = pickle.load(f)
    except:
        import windows_asciihelioport_standalone

        jorgeta = windows_asciihelioport_standalone.GAIA_GAIA_stack()
        with open(jorgeta_cachepath, "wb") as f:
            pickle.dump(obj=jorgeta, file=f)

    # Alternatively use the Astropy match_to_catalogue_sky functionality (uses a KDTree, much faster.)
    import astropy.units as u

    """
    jkgss = jorgeta[np.where(jorgeta['source'] == 'KGiant_edr3_metal')[0]]
    from astropy.coordinates import SkyCoord
    # We get no results for matching to 3 arcseconds. Match quality = zip after that anyway. :/ :/ :/
    JKGSS = SkyCoord(ra=jkgss['ra'] * u.deg, dec=jkgss['dec'] * u.deg, frame='icrs')
    LKGSS = SkyCoord(ra=LAMOKGs['ra'] * u.deg, dec=LAMOKGs['dec'] * u.deg, frame='icrs')"""
    """
    idx, d2d, d3d = JKGSS.match_to_catalog_sky(LKGSS)
    d2d = d2d.to(u.mas).value / 1000
    where_matched = np.where(d2d <= 3)[0]  # to 3 arcseconds
    idx, d2d, jkgss = idx[where_matched], d2d[where_matched], jkgss[where_matched]
    LAMOKGs = LAMOKGs[idx]

    from graphutils_new import spec_graph
    spec_graph().four_correlation_catalogue(jkgss, LAMOKGs, "SKGs", "LKGs")
    """
    # Alright. Bad. We need to try propagating the epoch?
    # TODO: Result: Didn't help at all.
    def propagate_astrometry(phi, theta, parallax, muphistar, mutheta, vrad, t0, t1):
        """
        Same as the PyGaia documentation. Propagate the astrometric parameters of a source from the reference epoch t0 to
        the new epoch t1. Note coordinates are ra-dec/lat-long style, not theta-phi style.

        Parameters
        ----------
        phi : float
            Longitude at reference epoch (radians).
        theta : float
            Latitude at reference epoch (radians).
        parallax : float
            Parallax at the reference epoch (mas).
        muphistar : float
            Proper motion in longitude (including np.cos(latitude) term) at reference
            epoch (mas/yr).
        mutheta : float
            Proper motion in latitude at reference epoch (mas/yr).
        vrad : float
            Radial velocity at reference epoch (km/s).
        t0 : float
            Reference epoch (Julian years).
        t1 : float
            New epoch (Julian years).

        Returns
        -------
        phi1, theta1, parallax1, muphistar1, mutheta1, murad1 : float or array
            Astrometric parameters, including the "radial proper motion" (NOT the radial
            velocity), at the new epoch.
        """
        from juliacall import Main
        Main.include(os.path.join(windows_directories_new.jl_dir, "propagate_epoch.jl"))
        phi1, theta1, parallax1, \
        muphistar1, mutheta1, murad1 = Main.propagate_astrometry(phi, theta, parallax,
                                                                 muphistar, mutheta, vrad,
                                                                 t0, t1)
        return phi1, theta1, parallax1, muphistar1, mutheta1, murad1
    """
    jkgss['parallax'] = 1/jkgss['dist']
    ra1, dec1, parallax1, pmra_cosdec1, pmdec1, pmrad1 = propagate_astrometry(np.deg2rad(jkgss['ra']),
                                                                              np.deg2rad(jkgss['dec']),
                                                                              jkgss['parallax'],
                                                                              jkgss['pmra']*np.cos(np.radians(jkgss['dec'])),
                                                                              jkgss['pmdec'],
                                                                              jkgss['vlos'],
                                                                              2008,2016)

    jkgss['ra'],jkgss['dec'] = ra1, dec1
    JKGSS = SkyCoord(ra=jkgss['ra'] * u.deg, dec=jkgss['dec'] * u.deg, frame='icrs')
    idx, d2d, d3d = JKGSS.match_to_catalog_sky(LKGSS)
    d2d = d2d.to(u.mas).value / 1000
    where_matched = np.where(d2d <= 60)[0]  # to 3 arcseconds
    print(len(where_matched))
    """
    """
    # Another test... take entire LRS catalogue + match to SEGUE KGs... to 3 arcseconds.
    # Get stellar subclass dataset
    with fits.open(stellar_dir, mode='readonly') as hdul:
        steldata = Table(hdul[1].data)
    stelcoords = SkyCoord(ra=steldata['ra_obs']*u.deg,
                          dec=steldata['dec_obs']*u.deg,
                          frame='icrs')
    idx, d2d, d3d = JKGSS.match_to_catalog_sky(stelcoords)
    d2d = d2d.to(u.mas).value / 1000
    where_matched = np.where(d2d <= 3)[0]  # to 3 arcseconds
    idx, d2d, jkgss = idx[where_matched], d2d[where_matched], jkgss[where_matched]
    print("Matching to the entire LRS catalogue... ", len(where_matched)) # only 4...? Come on!
    """

    """
    # Alright, Load the old DR4 catalogue and give it a shot using the observed ra-dec.
    dr4dir = os.path.join(windows_directories_new.lamodir, "dr4_v2_stellar.fits")  # where is VAC for dr4
    with fits.open(dr4dir, mode='readonly', memmap=True) as hdul:
        dr4data = Table(hdul[1].data)
    OBAFGKM = np.array([d[0] for d in dr4data['subclass']], dtype=str)
    is_ktype = np.where(OBAFGKM == "K")
    dr4data = dr4data[is_ktype]
    ra4,dec4 = dr4data['ra_obs'],dr4data['dec_obs']
    dr4_skycoord = SkyCoord(ra=ra4*u.deg, dec=dec4*u.deg, frame='icrs')
    idx, d2d, d3d = JKGSS.match_to_catalog_sky(dr4_skycoord)
    d2d = d2d.to(u.mas).value / 1000
    print(d2d)
    where_matched = np.where(d2d <= 3)[0]  # to 3 arcseconds
    idx, d2d, jkgss = idx[where_matched], d2d[where_matched], jkgss[where_matched]
    dr4data = dr4data[idx]
    print(dr4data) # 0 matches! Mother..."""

    """
    # Alright. DR4 doesn't match at all to 3 arseconds. Neither does DR7. :sob: Let's try an XMatch of our KGs to Gaia
    from gaia_utils_standalone import xmatch_gaiadr3

    print(len(jkgss), " is the length of the original jkgss.")
    jkgss['jkgss_id'] = np.arange(0, len(jkgss), 1)
    jkgss.write("SEGUE_KGS.fits", overwrite=True)
    _cols_of_interest = ['angDist', 'Source', 'jkgss_id',
                         'ra', 'dec', 'RAdeg', 'DEdeg',
                         'Plx', 'e_Plx',
                         'pmRA', 'e_pmRA',
                         'pmDE', 'e_pmDE',
                         'RADEcor', 'RAPlxcor', 'RApmRAcor', 'RApmDEcor',
                         'DEPlxcor', 'DEpmRAcor', 'DEpmDEcor', 'PlxpmRAcor',
                         'PlxpmDEcor', 'pmRApmDEcor',
                         'RUWE',
                         'RV', 'e_RV']
    jkgss_matchedtogaia = xmatch_gaiadr3(jkgss, "ra", "dec", "lkgss_id", 3000, _cols_of_interest)
    print(jkgss_matchedtogaia) # angDist is in arcseconds. Only 103 matches... To GAIA DR3??? ... Damn. Alright.

    # Propagate epoch...?
    ra1, dec1, parallax1, pmra_cosdec1, pmdec1, pmrad1 = propagate_astrometry(np.deg2rad(jkgss['ra']),
                                                                              np.deg2rad(jkgss['dec']),
                                                                              1/jkgss['dist'],
                                                                              jkgss['pmra']*np.cos(np.radians(jkgss['dec'])),
                                                                              jkgss['pmdec'],
                                                                              jkgss['vlos'],
                                                                              2008,2016)
    jkgss['ra'],jkgss['dec'],jkgss['parallax'],jkgss['pmra'],jkgss['pmdec'] = ra1, dec1, parallax1, pmra_cosdec1, pmdec1
    jkgss.write("SEGUE_KGS.fits", overwrite=True)
    jkgss_matchedtogaia = xmatch_gaiadr3(jkgss, "ra", "dec", "lkgss_id", 3000, _cols_of_interest)
    print(jkgss_matchedtogaia) # epoch propagation again didn't help at all. What is going wrong here...? 
    # """

    # Alright. Take the new LAMOST DR7 catalogue, match to Jorges LAMOST stars, and do correlations.
    jLkgss = jorgeta[np.where(jorgeta['source'] == 'LAMOST_K_FULL_edr3')[0]]
    from astropy.coordinates import SkyCoord
    # We get no results for matching to 3 arcseconds. Match quality = zip after that anyway. :/ :/ :/
    JLKGSS = SkyCoord(ra=jLkgss['ra'] * u.deg, dec=jLkgss['dec'] * u.deg, frame='icrs')
    LKGSS = SkyCoord(ra=data['ra'] * u.deg, dec=data['dec'] * u.deg, frame='icrs')
    idx, d2d, d3d = JLKGSS.match_to_catalog_sky(LKGSS)
    d2d = d2d.to(u.mas).value / 1000
    where_matched = np.where(d2d <= 3)[0]  # to 3 arcseconds
    jLkgss = jLkgss[where_matched]
    idx, d2d = idx[where_matched], d2d[where_matched]
    data = data[idx]
    import graphutils_new
    graphutils_new.spec_graph().four_correlation_catalogue(jLkgss,
                                                           data,
                                                           "Jorge LAMOST",
                                                           "DR7 LAMOST",
                                                           os.path.join(windows_directories_new.lamodir,
                                                                        "OLD_LAMOST-NEW_LAMOST-4corr.png"))

# Section II of the TODO_LIST.txt
SECTION_II = False
if __name__ == "__main__" and SECTION_II == True:

    with fits.open(dr7_gaia_gaia_dir, mode='readonly', memmap=True) as hdul:
        data = Table(hdul[1].data)  # LAMOST-GAIA catalogue to 2 arcsecond crossmatch
        data['source'] = "LAMOKGS_straszaks" # just for the sake of keeping track

    # Require parallax/error > 50
    data['parallax_over_error'] = data['parallax']/data['parallax_error']
    data = data[np.where(data['parallax_over_error'] >= 9)[0]]
    data = data[np.where(data['ruwe'] < 1.4)[0]]

    # Require distance/error > 5
    data['dist_over_edist'] = data['dist']/data['edist']
    data = data[np.where(data['dist_over_edist'] >= 3)[0]]

    # Select the most distant 100 by parallax
    pardist = 1 / data['parallax']
    #indices = np.argsort(pardist) # ascending order
    #data = data[indices[0:300000]]
    indices = np.where(pardist > 0.5)[0]
    data = data[indices]
    pardist = pardist[indices]
    indices = np.where(pardist < 3)[0]
    data = data[indices]
    import random
    indices = random.sample(range(len(data)), 300000)
    data = data[indices]

    # Compare the two on a plot (parallax distance vs. LAMOST distance.)
    from graphutils_new import spec_graph
    spec_graph().parallax_versus_distance(data,
                                          os.path.join(windows_directories_new.lamodir,
                                                       "LAMOST_GAIA_DISTPARISON.png"))#_firstfew.png"))

# Pre-process data + generate covariance matrices in galactocentric space whilst removing GCs.
# Hereon, we enact a 4 kpc minimum distance requirement.
do_gcs = False
produce_data = False
plot = False
min_radius = 15
max_radius = 150  # same as orbigistics interpolator for circularity calculation
final_dir = os.path.join(windows_directories_new.lamodir, "LAMOST_master_2MASS.fits")
if __name__ == "__main__" and do_gcs == True:

    if produce_data:
        # with fits.open(dr7_gaia_gaia_dir, mode='readonly', memmap=True) as hdul:
        #    data = Table(hdul[1].data)

        # Convert data to pandas dataframe anticipating covtrices
        # data = data.to_pandas()

        # Set up covariance matrices for the ICRS frame of reference.
        """
        When we Monte, we are doing it in 6D via spherical coordinates using...
        RA, DEC, PARALLAX, PMRA, PMDEC, VLOS 
        We assume no error in RA or DEC (we have none through LAMOST anyway) so reduce covtrix to 4^2 not 6^2
        The remaining covtrices thus have PARALLAX, PMRA, PMDEC, VLOS
        Which is built up via...
        [[eplx, plx-pmra-corr, plx-pmdec-corr, plx-vlos-corr,],
         [pmra-plx-corr, epmra, pmra-pmdec-corr, pmra-vlos-corr],
         [pmdec-plx-corr, pmdec-pmra-corr, epmdec, pmdec-vlos-corr],
         [vlos-plx-corr, vlos-pmra-corr, vlos-pmdec-corr, evlos]]
        Generate this as a column in the pandas dataframe and then we propagate this to galactic coordinates using 
        Monte, and then to galactocentric using Monte- that will be done later only for distant stars, though. 
        """

        with fits.open(dr7_gaia_gaia_dir, mode='readonly', memmap=True) as hdul:
            data = Table(hdul[1].data)

        # Get the galactic coordinates from these equatorial ones
        from energistics_new import orbigistics
        from galcentricutils_new import angular

        orbigist = orbigistics()

        # Get GAL from ICRS
        data = orbigist.converter.nowrite_ICRS_to_GAL(data, has_cosfactor=True)

        # Get galactic errors. The cosfactor in the pmra/pmdec still exists- no longer in dmu_l, dmu_b though.
        import lamost_utils

        data = lamost_utils.monte_ICRSGAL_table(data)  # TODO: Take advantage of correlations from Gaia

        # Remove cosfactor from ICRS to match galconversion convention.
        data['pmra'] /= np.cos(np.deg2rad(data['dec']))
        data['pmra_error'] /= np.cos(np.deg2rad(data['dec']))

        # Get Galcent + Momenta
        data = orbigist.converter.nowrite_GAL_to_GALCENT(data)
        data = angular().get_momentum(data)

        # Remove nuisance stars (outside min_radius, but within the max_radius of interest.)
        data = data[[True if max_radius > r > min_radius else False for r in data['r']]]

        # Grab the table of globular clusters
        with fits.open(windows_directories_new.baumgardt_fits, mode='readonly', memmap=False) as hdul:
            gc_table = Table(hdul[1].data)

        # Add dummy values for proper motions
        gc_table.add_column(0., name='pmra')
        gc_table.add_column(0., name='pmdec')
        gc_table.add_column(0., name='pmra_error')
        gc_table.add_column(0., name='pmdec_error')
        gc_table.add_column(0., name='vlos')
        gc_table.add_column(0., name='evlost')

        # Get Galcent (note that baumgardt catalogue does have a cosfactor attached- cos(dec).)
        gc_table = orbigist.converter.nowrite_ICRS_to_GAL(gc_table, has_cosfactor=True)
        gc_table = orbigist.converter.nowrite_GAL_to_GALCENT(gc_table)

        # Convert tidal radius to kpc
        gc_table['rt'] /= 1000

        # Remove stars within gcs
        import gcc_utils

        gc_to_remove, which_gc = gcc_utils.remove_gcs(np.array(data['x'], float), np.array(data['y'], float),
                                                      np.array(data['z'], float),
                                                      np.array(gc_table['x'], float), np.array(gc_table['y'], float),
                                                      np.array(gc_table['z'], float),
                                                      np.array(gc_table['rt'], float))
        data["in_GGC"] = gc_to_remove
        data["which_GGC"] = which_gc

        # Save fits
        data.write(final_dir, overwrite=True)

    if plot:
        with fits.open(final_dir, mode='readonly', memmap=True) as hdul:
            # "gaia_source_id" is dr2_source_id inside this.
            data = Table(hdul[1].data)

        # data = data[[True if 1 - abs(circ) > 0.9 else False for circ in data['circ']]]
        # data = data[[True if abs(d) < 300 else False for d in data['vx']]]
        # data = data[[True if abs(d) < 300 else False for d in data['vy']]]
        # data = data[[True if abs(d) < 300 else False for d in data['vz']]]
        # data = fast_energistics_new().default_E_c(data)

        # Won't work without this step for some damn reason- type(data['x'][0]) = float anyway :/
        data['x'], data['y'], data['z'], \
        data['vx'], data['vy'], data['vz'] = np.array(data['x'], float), \
                                             np.array(data['y'], float), \
                                             np.array(data['z'], float), \
                                             np.array(data['vx'], float), \
                                             np.array(data['vy'], float), \
                                             np.array(data['vz'], float)

        # plt.scatter(master_fits['Lz'], master_fits['E'], color='blue', s=1)
        # plt.scatter(table['Lz'], table['E'], color='red', s=0.5)
        # plt.hist(table['vlos'], bins=100)
        # plt.hist(data['vlos'], bins=100)
        # plt.show()

        # L-space plotly plot
        datatab = np.array([data['Lx'], data['Ly'], data['Lz']]).T
        import graphutils_new
        import matplotlib.pyplot as plt
        graphutils_new.threed_graph().kmeans_L_array(datatab, [1 for d in datatab[:, 0]],
                                                     False, browser=True, outliers=True)

        # Velocity-space histogram (just another quality check.)
        plt.hist(data['vx'], bins=1000, label='vx')
        plt.hist(data['vy'], bins=1000, label='vy')
        plt.hist(data['vz'], bins=1000, label='vz')
        plt.xlim([-1000, 1000])
        plt.legend()
        plt.savefig(os.path.join(windows_directories_new.lamodir, "velspace_histogram_GAIA-GAIA.png"), dpi=300)

        # E-Lz plot.
        clustered = [1 for d in data['x']]
        import fast_energistics_new
        data = fast_energistics_new().default_E_c(data)
        graphutils_new.spec_graph().energy_plot(data,
                                                clustered,
                                                list(set(clustered)),
                                                list(set(clustered)),
                                                [-10000, 10000],
                                                [-200000, 10000],
                                                "E-Lz_test_GAIA-GAIA.png",
                                                False,
                                                False)

# Should I grab metallicities from the main LRS based on object ID for the value-added-catalogue?
do_get_feh = False
new_final_dir = "LAMOST_master_2MASS_feh.fits"
if __name__ == "__main__" and do_get_feh == True:

    # Load in the 2MASS/etc fits
    with fits.open(final_dir, mode='readonly', memmap=True) as hdul:
        data = Table(hdul[1].data)

    # Load in the basic LRS
    with fits.open(windows_directories_new.LRS_path, mode='readonly', memmap=True) as hdul:

        # "gaia_source_id" is dr2_source_id inside this.
        fits_LRS = hdul[1].data

    # Convert the data objid column back to an integer (we turned it into a float when taking it from the dat file.)
    data['objid'] = np.array(data['objid'], int)

    objids_lrs = np.array(fits_LRS['obsid'], int)


    @njit(fastmath=True)
    def match(labs1, labs2):

        """
        Match integers (by index) from labs1 with labs2, returning a list of tuples in which you have. FINDS MATCHES
        FOR LABS1! The return is an array len(labs1) of form (index labs1 match, index labs2 match)

        :param labs1: int list array
        :param labs2: int list array
        :return: matches, truefalse
        """

        tuples = np.zeros((len(labs1), 2), numba.types.int64)

        for i in range(0, len(labs1)):  # do the match (return (0,0) if no match for a given objid in labs1)
            for j in range(0, len(labs2)):
                if labs1[i] == labs2[j]:
                    tuples[i][0] = i
                    tuples[i][1] = j

        is_matched = []
        for i in range(0, len(labs1)):  # check to ensure that matches exist
            if tuples[i][0] != 0:
                if tuples[i][1] != 0:
                    is_matched.append(True)
            else:
                is_matched.append(False)

        return tuples, is_matched


    tuples, is_matched = match(data['objid'], objids_lrs)
    tuples = tuples[is_matched]

    try:
        data.add_column(0., name="feh")
        data.add_column(0., name="efeh")
    except:
        pass

    # Set the radial velocities/metallicities
    for tuple in tuples:
        data[tuple[0]]['vlos'], data[tuple[0]]['evlost'], \
        data[tuple[0]]['feh'], data[tuple[0]]['efeh'] \
            = fits_LRS[tuple[1]]['rv'], fits_LRS[tuple[1]]['rv_err'], \
              fits_LRS[tuple[1]]['feh'], fits_LRS[tuple[1]]['feh_err']

    # Update the fits + save
    data.write(new_final_dir, overwrite=True)

# Should I now crossmatch to the old catalogue to remove duplicates, and then stack the tables? (whilst removing slag)
do_finalstack = False
finalstack_load = False  # just load the crossmatch and continue (match is expensive process.)
matchcache = os.path.join(windows_directories_new.lamodir, "finalstack_matchcache.txt")
finalstack_dir = os.path.join(windows_directories_new.lamodir, "LAMOST_final.fits")
finalstack_ascii = os.path.join(windows_directories_new.lamodir, "LAMOST_final.dat")
if __name__ == "__main__" and do_finalstack == True:

    # Load in the 2MASS/etc fits
    with fits.open(new_final_dir, mode='readonly', memmap=True) as hdul:

        data = Table(hdul[1].data)

    # Add a column to specify the "SOURCE"
    data['source'] = "LAMOST_sstraszak"

    # Open up the HDF that has the LAMOST KGs
    import hdfutils
    from energistics_new import orbigistics
    from galcentricutils_new import angular

    writer = hdfutils.hdf5_writer(windows_directories_new.datadir, ascii_info_new.asciiname)
    lamokgs = writer.read_table(ascii_info_new.lamostk, ascii_info_new.set_raw)

    # convert to ICRS
    orbigist = orbigistics()
    lamokgs = lamokgs[[True if min_radius < r < max_radius else False for r in lamokgs['r']]]

    # Grab ra/decs from both (in radians)!!!!!!!!!
    lamokgs_ra, lamokgs_dec, \
    data_ra, data_dec = lamokgs['ra'], lamokgs['dec'], \
                        data['ra'], data['dec']
    lamokgs_ra, lamokgs_dec, \
    data_ra, data_dec = np.deg2rad(lamokgs_ra), np.deg2rad(lamokgs_dec), \
                        np.deg2rad(data_ra), np.deg2rad(data_dec)

    # Get matches between the older LAMOST KGs and the newer (obtained) LAMOST KGs + remove the rows that match.
    lamokgs_ra = lamokgs_ra
    lamokgs_dec = lamokgs_dec

    try:

        if finalstack_load == True:

            print("Trying to load matchcache...")

            with open(matchcache, "rb") as f:
                matches, truefalse = pickle.load(f)
                print("Successfully loaded the cachematch for finalstack matches.")

        else:

            raise RuntimeError("Will generate match instead...")

    except:

        from APOGEE_standalone import radecmatch_argmin_memlim

        matches, truefalse = radecmatch_argmin_memlim(lamokgs_ra, lamokgs_dec, data_ra, data_dec)
        matches, truefalse = list(matches), list(truefalse)
        matches, truefalse = np.asarray(matches), np.asarray(truefalse)

        with open(matchcache, "wb") as f:
            pickle.dump(obj=[matches, truefalse], file=f)
            print("Successfully generated the cachematch for final stack.")

    matches = matches[truefalse]
    matches = matches.T


    # Generate correlation plots for the matches to check integrity of catalogue source. Note that lamokgs is our
    # old data and "data" is the new data from the value-added-catalogue from LAMOST DR7.
    def do_correlation_plots(old_catalogue, new_catalogue,
                             columns, columns_errs,
                             columns_labels_old, columns_labels_new,
                             savedir):

        """

        Generate correlation plots between an old catalogue and a new catalogue, of same size.

        :param old_catalogue: astropy table or pandas
        :param new_catalogue: astropy table or pandas
        :param columns: list of table labels
        :param columns_errs: list of table err labels
        :param columns_labels_old: labels for matplotlib
        :param columns_labels_new: labels for matplotlib
        :param savedir: directory to save graphs
        :return: -

        """

        # Make savedir
        try:
            os.mkdir(savedir)
        except:
            pass

        # Go over each column and do the graphs/etc
        for column, column_err, \
            column_label_old, column_label_new in zip(columns, columns_errs,
                                                      columns_labels_old, columns_labels_new):
            # Make the save path
            savepath = os.path.join(savedir, column + "_correlation_plot.png")

            # Do the graph
            import graphutils_new
            graphutils_new.twod_graph().corr_plot(old_catalogue[column],
                                                  new_catalogue[column],
                                                  old_catalogue[column_err],
                                                  new_catalogue[column_err],
                                                  column_label_old,
                                                  column_label_new,
                                                  savepath)


    lamokgs_matches, data_matches = matches
    lamokgs_rows, data_rows = lamokgs[lamokgs_matches], data[data_matches]
    razip = np.array([lamokgs_rows['ra'], data_rows['ra']]).T
    columns = ["dmu_l", "dmu_b", "vlos", "dist"]
    columns_errs = ["edmu_l", "edmu_b", "evlost", "edist"]
    labels = [r'$\mu_l$', r'$\mu_b$', r'$\textrm{v}_\textrm{los}$', r'$\textrm{d}$']
    savedir = os.path.join(windows_directories_new.lamodir, "correlation_plots_lamost")
    do_correlation_plots(lamokgs_rows, data_rows,
                         columns, columns_errs,
                         labels, labels,
                         savedir)

    # Remove the rows that matched (from the old lamokgs)
    matches = matches[0]
    lamokgs.remove_rows(matches)

    # Now we can recombine with the rest of our dataset.
    # Load in the rest of the tables (sans LAMOST old)
    tables = [writer.read_table(d, ascii_info_new.set_raw) for d in ascii_info_new.all_groups[0:3]]
    tables += lamokgs  # these have had match rows removed
    tables += data  # add the new LAMOST data

    # Vstack
    data = vstack(tables)

    # Remove the nuisance columns (including objid/etc.)
    to_remove_columns = ["sgrstrong", "FeH_1", "Sgr Lambda", "Sgr Beta", "Belokurov Flag", "corrcoef",
                         "ra_2", "dec_2", "void", "void_1", "void_2", "name",
                         "plx_mas_", "eplx_mas_", "idup", "objid", "name_2mass", "flag2mass", "snrg",
                         "cvlos", "cdist"]
    data.remove_columns(to_remove_columns)

    # Find out which rows are nuisance rows
    to_remove = [False if np.isnan(d) == True else True for d in data['vx']]
    data = data[to_remove]
    to_remove = [False if np.isnan(d) == True else True for d in data['x']]
    data = data[to_remove]

    print("Saving the finalstack...")

    # Save
    data.write(finalstack_dir, overwrite=True)
    from astropy.io import ascii

    ascii.write(data, finalstack_ascii, overwrite=True)

# Should I do the stuff that Mike wanted (i.e. remove LAMOKGs/SEGUE KGs duplicates + get biases between them?)


"""
Load data in and do some cuts on it 
"""
do_cuts = False
cuts_fits_path = os.path.join(windows_directories_new.lamodir, "cut_LAMOST_final.fits")
if __name__ == "__main__" and do_cuts == True:

    # Load in the finalstack
    with fits.open(finalstack_dir, mode='readonly', memmap=True) as hdul:

        orig_data = Table(hdul[1].data)

        # Specifically retrieve the new LAMOST data subsection
        indices = []
        for i in range(len(orig_data)):
            if orig_data[i]['source'] == "LAMOST_sstraszak":
                if np.abs(orig_data[i]['feh']) != 0.:
                    indices.append(i)

        # Set data
        data = orig_data[indices]

    # QUALITY CUTS ======================================^^^^
    from galcentricutils_quality_cuts import quality_cuts

    qcs = quality_cuts()
    metmax = -1 / 2
    LSR_VEL = 220
    FANCY_LSR_VEL = 100
    ZMAX = 5 / 2
    vphi_plotlims = [-550, 550]

    # Get indices for all the cuts necessary (from the original) dataset- for illustration purposes
    FSLR_KEEP, FSLR_REMOVE = qcs.fancy_LSR_cut(data, FANCY_LSR_VEL)
    RLSR_KEEP, RSLR_REMOVE = qcs.LSR_cut(data, LSR_VEL)
    FFEH_KEEP, FFEH_REMOVE = qcs.fancy_feh_cut(data)
    RFEH_KEEP, RFEH_REMOVE = qcs.feh_cut(data, metmax)
    ZMAX_KEEP, ZMAX_REMOVE = qcs.zmax_cuts(data, ZMAX,
                                           3e9,
                                           3000,
                                           True,
                                           True,
                                           "lamost_dr7_cuts_fulldata.txt")
    BOND_KEEP, BOND_REMOVE = qcs.bound_cut(data)
    indices_lists = [FSLR_KEEP, RLSR_KEEP, FFEH_KEEP, RFEH_KEEP, ZMAX_KEEP, BOND_KEEP]
    remove_lists = [FSLR_REMOVE, RSLR_REMOVE, FFEH_REMOVE, RFEH_REMOVE, ZMAX_REMOVE, BOND_REMOVE]
    indices_names = ["fancy_LSR", "regular_LSR", "fancy_feh", "regular_feh", "zmax", "bound"]

    # Get vT/etc anticipating plotting
    from energistics_new import orbigistics

    R, vR, vT, z, vz, phi = orbigistics().get_leftgalpy(data)
    feh = data['feh']

    import matplotlib.pyplot as plt

    for keep, remove, cut_type in zip(indices_lists, remove_lists, indices_names):
        fig, axs = plt.subplots(nrows=1, ncols=1)
        axs.set(xlim=vphi_plotlims,
                ylim=[-5 / 2, metmax])
        axs.grid(True, which='major', color='pink')
        axs.set(xlabel=r'$\textrm{v}_\phi$',
                ylabel=r'$[\textrm{Fe/H}]$')

        # Get the data we aren't cutting out + plot it
        VT_PLOT, FEH_PLOT = vT[keep], feh[keep]
        axs.scatter(VT_PLOT, FEH_PLOT, s=0.1, color='green')

        # Get the data we are cutting + plot it
        VT_PLOT, FEH_PLOT = vT[remove], feh[remove]
        axs.scatter(VT_PLOT, FEH_PLOT, s=0.1, color='red')

        plt.savefig(os.path.join(qcs.savedir, "finalstack_LAMOST_metplot_" + cut_type + ".png"), dpi=300)

    # Get non-lamost indices, recombine to LAMOST_sstraszak indices, to reclaim whole dataset
    no_indices = np.where(orig_data['source'] != "LAMOST_sstraszak")[0]
    orig_data = orig_data[no_indices]
    orig_data = vstack([orig_data, data])

    """
    # While we're at it, also grab the APOGEE data we obtained
    with fits.open(apogee_starhorse_path_nogcs, mode='readonly', memmap=True) as hdul:

        apodata = Table(hdul[1].data)
        apodata = apodata[[True if r < 150 else False for r in apodata['r']]]
        apodata['source'] = "APOGEE_sstraszak"
        orig_data = vstack([orig_data, apodata]) """

    # Generate base plot for debug
    R, vR, vT, z, vz, phi = orbigistics().get_leftgalpy(orig_data)
    orig_data['vT'] = vT
    plt.scatter(vT, orig_data['feh'], s=0.1)
    plt.ylim([-3, 1])
    plt.savefig(os.path.join(qcs.savedir, "orig_data_beforecuts.png"), dpi=300)

    # Carry out cuts on orig_data (the full dataset) in anticipation of use
    # noinspection PyRedeclaration
    KEEP, REMOVE = qcs.fancy_LSR_cut(orig_data, FANCY_LSR_VEL)
    orig_data = orig_data[KEEP]
    # noinspection PyRedeclaration
    KEEP, REMOVE = qcs.zmax_cuts(orig_data, ZMAX,
                                 3e9,
                                 3000,
                                 True,
                                 True,
                                 "lamost_dr7_cuts_origdata.txt")
    orig_data = orig_data[KEEP]
    plt.clf()
    plt.close()
    plt.cla()
    plt.scatter(orig_data['vT'], orig_data['feh'], s=0.1)
    plt.xlim([-500, 500])
    plt.ylim([-3, 1])
    print("Showing")
    plt.show()
    # noinspection PyRedeclaration
    # KEEP, REMOVE = qcs.fancy_feh_cut(orig_data)
    # orig_data = orig_data[KEEP]

    # Save data to fits
    orig_data.write(cuts_fits_path, overwrite=True)

    # Write data to table
    import hdfutils

    writer = hdfutils.hdf5_writer(windows_directories_new.datadir, ascii_info_new.asciiname)
    writer.write_table(ascii_info_new.fullgroup, ascii_info_new.fullset, orig_data)

    # Set up the flat hdbscan run
    data_array = np.array([orig_data['Lx'], orig_data['Ly'], orig_data['Lz']]).T
    from hdbscan import flat

    clusterer = flat.HDBSCAN_flat(X=data_array,
                                  n_clusters=40,
                                  min_cluster_size=10,
                                  max_cluster_size=int(0.4 * len(orig_data)),
                                  min_samples=7,
                                  metric='l2',
                                  algorithm='best')  # ,
    # prediction_data=False)
    set_of_labels = list(set(clusterer.labels_))
    sizes = [len(np.where(clusterer.labels_ == d)[0]) for d in set_of_labels]
    sorted_sizes = np.sort(sizes, axis=0)
    GSE_SIZE = np.flip(sorted_sizes)[1]

    # Take our cleaned data sample and produce angular plot
    savepath = os.path.join(qcs.savedir, "kmeans_L_LAMOGEEFINAL_test.html")
    graphutils_new.threed_graph().kmeans_L_array_path(data_array,
                                                      clusterer.labels_,
                                                      savepath,
                                                      True,
                                                      True,
                                                      True)

    # Savedir
    savedir = os.path.join(qcs.savedir, "orbifits_tests")
    try:
        os.mkdir(savedir)
    except:
        pass

    # For all our clusters found :3
    run_orbifitplots = False
    if run_orbifitplots == True:

        for clust in list(set(clusterer.labels_)):
            import energistics, graphutils

            orbifit = energistics.orbifitter().finetune_galpy_fitting_nomemb(orig_data,
                                                                             clusterer.labels_,
                                                                             clust,
                                                                             2000,
                                                                             0.6e9,
                                                                             1000,
                                                                             False,
                                                                             False,
                                                                             False,
                                                                             False)

            # Forward Integral
            import copy
            from astropy import units as u

            forward = copy.deepcopy(orbifit)
            backward = copy.deepcopy(orbifit)
            forward.integrate((np.linspace(0, 1e9, 2000) * u.yr), energistics.orbigistics().pot)
            backward.integrate((np.linspace(0, -1e9, 2000) * u.yr), energistics.orbigistics().pot)
            llsf, bbsf, llsb, bbsb = forward.ll((np.linspace(0, 0.7e9, 2000) * u.yr)).value, \
                                     forward.bb((np.linspace(0, 0.7e9, 2000) * u.yr)).value, \
                                     backward.ll((np.linspace(0, -0.15e9, 2000) * u.yr)).value, \
                                     backward.bb((np.linspace(0, -0.15e9, 2000) * u.yr)).value
            lls, bbs = np.concatenate([llsf, llsb]), np.concatenate([bbsf, bbsb])
            lls = [d - 360 if d > 180 else d for d in lls]

            specgrapher = graphutils.spec_graph()
            fig, axs = specgrapher.lb_orbits(orig_data[[True if d == clust
                                                        else False
                                                        for d in clusterer.labels_]],
                                             0.3e9, [-180, 180],
                                             [-90, 90], None,
                                             line=False, points=4000)
            axs.scatter(lls, bbs, color='lime', marker='o', s=3)
            # Misc axis things
            axs.set_facecolor("k")
            axs.grid(True, which='major', alpha=1, linewidth=0.25, color='white')
            axs.grid(True, which='minor', alpha=1, linewidth=0.25, color='white')
            axs.grid(color="white")
            savepath = os.path.join(savedir, str(clust) + ".png")
            import matplotlib.pyplot as plt

            plt.savefig(savepath, dpi=300, transparent=False)
            plt.close()

    # Check mets
    for clust in list(set(clusterer.labels_)):
        # Get metallicity + dispersion
        indices = np.where(clusterer.labels_ == clust)[0]
        fehs = orig_data['feh'][indices]
        from astropy.stats import sigma_clipped_stats

        mean, med, std = sigma_clipped_stats(fehs)
        print(clust, mean, med, std)
