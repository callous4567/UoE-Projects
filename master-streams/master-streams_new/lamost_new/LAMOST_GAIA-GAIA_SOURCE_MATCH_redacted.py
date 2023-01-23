import os
import pickle
import matplotlib.pyplot as plt
import numba
import numpy as np
from astropy.io import fits
from astropy.table import Table, vstack
from numba import njit
import ascii_info_new
import graphutils_new
import windows_directories_new
from energistics_new import fast_energistics_new

# TODO NOTE: general catalogue FITS broken- use CSV (need to process.)
# TODO NOTE: dr7_v2.0_LRS_stellar.fits works fine as a FITS straight out of their website!
# TODO NOTE: the value-added-catalogue is a DAT file, too. Double oof.
process_vat = False
raw_dir = os.path.join(windows_directories_new.lamodir, "dr7_lamost_LR_release.dat")
stellar_dir = os.path.join(windows_directories_new.lamodir, "dr7_v2.0_LRS_stellar.fits")
vat_dir = os.path.join(windows_directories_new.lamodir, "dr7_LRS_vat.fits")
dr7_gaia_gaia_dir = os.path.join(windows_directories_new.lamodir, "dr7_GAIA_aidists.fits")
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

# Should I crossmatch the LAMOST dataset with Gaia DR3 using Vizier, obtaining all the columns of interest?
LAMOST_GAIA_GAIA = False
LAMOST_GAIA_GAIA_STELLARJOIN = True  # join to the LRS stellar catalogue on objID to get spectral class/etc
if __name__ == "__main__" and LAMOST_GAIA_GAIA == True:
    import gaia_utils_standalone

    with fits.open(vat_dir, mode='readonly', memmap=True) as hdul:
        # "gaia_source_id" is dr2_source_id inside this.
        dr7_LRS_vat = Table(hdul[1].data)

    # Get stars a minimum of 4 kpc from us (hashed out for now- let's avoid clipping any for now- need KGs correlation.)
    # dr7_vat = dr7_LRS_vat[[True if d > 4 else False for d in dr7_LRS_vat['d_kpc_']]]
    # dr7_crossmatch_vat = dr7_LRS_vat
    # preupload_names = dr7_vat.colnames

    # Condense dr7_LRS_vat catalogue to include only the RA, DEC, and objID
    dr7_crossmatch_vat = dr7_LRS_vat[["Ra", "Dec", "objID", 'RV', 'eRV', 'd_kpc_', 'e_d_kpc_']]

    # Rename for uniqueness the ra/dec (gaia/vizier/CDS doesn't like capitals, apparently.)
    dr7_crossmatch_vat.rename_columns(["Ra", "Dec", "objID",
                                       'RV', 'eRV',
                                       'd_kpc_', 'e_d_kpc_'],
                                      ["ra_dr7", "dec_dr7", "objid",
                                       'LAMOST_VLOS', 'LAMOST_EVLOST',
                                       'LAMOST_DIST', 'LAMOST_EDIST'])

    # Specify columns to take from the return catalogue (including the distance/radial velocity we put in from LAMOST.)
    """
    NOTE THE USUAL APPLIES! pmra is MULTIPLIED BY COS(DEC) and parallax IS IN MAS/YR and RUWE > 1.4 IS NEEDED LATER.
    """
    cols_of_interest_ = ['angDist', 'objid', 'Source',
                         'ra_dr7', 'dec_dr7', 'RAdeg', 'DEdeg',
                         'Plx', 'e_Plx',
                         'pmRA', 'e_pmRA',
                         'pmDE', 'e_pmDE',
                         'RADEcor', 'RAPlxcor', 'RApmRAcor', 'RApmDEcor',
                         'DEPlxcor', 'DEpmRAcor', 'DEpmDEcor', 'PlxpmRAcor',
                         'PlxpmDEcor', 'pmRApmDEcor',
                         'RUWE',
                         'RV', 'e_RV',
                         'LAMOST_VLOS', 'LAMOST_EVLOST',
                         'LAMOST_DIST', 'LAMOST_EDIST']

    # Run the xmatch. Remember that tolerance is in mas. We use 1 arcsecond tolerance.
    master_table = gaia_utils_standalone.xmatch_gaiadr3(dr7_crossmatch_vat,
                                                        "ra_dr7", "dec_dr7", "objid",
                                                        2000,
                                                        cols_of_interest_)
    cols_of_renamed_ = ['GDR3_ANG_DIST', 'LAMOST_objID', 'dr3_source_id',
                        'ra', 'dec', 'ra_dr7', 'dec_dr7',
                        'parallax', 'parallax_error',
                        'pmra', 'pmra_error',
                        'pmdec', 'pmdec_error',
                        'radec_corr', 'raparallax_corr', 'rapmra_corr', 'rapmdec_corr',
                        'decparallax_corr', 'decpmra_corr', 'decpmdec_corr', 'parallaxpmra_corr',
                        'parallaxpmdec_corr', 'pmrapmdec_corr',
                        'ruwe',
                        'vlos_gaia', 'evlost_gaia',
                        'vlos', 'evlost',
                        'dist', 'edist']
    master_table.rename_columns(cols_of_interest_, cols_of_renamed_)

    # Get inverted parallax/etc
    master_table['dist_gaia'], \
    master_table['edist_gaia'] = 1 / master_table['parallax'], \
                                 master_table['parallax_error'] / (master_table['parallax'] ** 2)

    # Get stellar subclass dataset
    with fits.open(stellar_dir, mode='readonly') as hdul:
        steldata = Table(hdul[1].data)['obsid', 'objtype', 'subclass']
        steldata.rename_column('obsid', 'LAMOST_objID')

    from astropy.table import join

    # Join master_table w the LRS stellar catalogue
    master_table = join(master_table, steldata, keys='LAMOST_objID')

    # Write it to disk
    master_table.write(dr7_gaia_gaia_dir, overwrite=True)

# Section I of the TODO_LIST.txt
SECTION_I = True
if __name__ == "__main__" and SECTION_I == True:

    """
    SECTION I
    ========================================================================================================================
    Take Yanny KGs
    Take LAMOST KGs
    Crossmatch them + extract matches
    Generate correlation plots in vlos-space and examine distance relationships too.
    
    Result:
    - can't get any matches
    - tried matching to LAMOST DR4 and DR7 to no avail- little matches to within 3 arseconds
    - tried epoch propagation to no avail
    - tried matching Yanny KGs to Gaia to 3 arcseconds to no avail- only 104 matches with 3 arseconds.
    - tried epoch propagation to no avail
    - tried matching using the CDS XMatch service and also via Topcat, alongside using Astropy
    - ???
    
    =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    """

    with fits.open(dr7_gaia_gaia_dir, mode='readonly', memmap=True) as hdul:
        data = Table(hdul[1].data)  # LAMOST-GAIA catalogue to 2 arcsecond crossmatch
        data['source'] = "LAMOKGS_straszaks" # just for the sake of keeping track

    # Get only the K-type stars as a subset... remember each entry is "TYPE" + "0->9" heat/temperature
    # TODO: No matches within 3 arcseconds. Expand to all? Doesn't help.
    LAMOKGs_truefalse = [True if d == "K" else False for d in [d[0] for d in data['subclass']]]
    LAMOKGs = data[LAMOKGs_truefalse]


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
    jkgss = jorgeta[np.where(jorgeta['source'] == 'KGiant_edr3_metal')[0]]

    # Alternatively use the Astropy match_to_catalogue_sky functionality (uses a KDTree, much faster.)
    import astropy.units as u

    """
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
    jkgss.write("SEGUE_KGS.fits", overwrite=True) # manually matching with Topcat only gives the same results... wth?
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