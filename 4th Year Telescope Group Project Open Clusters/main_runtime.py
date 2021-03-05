import os, time, multiprocessing, numpy as np
from TGP2020 import astrowrapper # Change 'TGP2020' to the WINDOWS FOLDERNAME of ROOTDIR.
utils = astrowrapper.utils
hdf5_writer = astrowrapper.hdf5_writer
fits_alignment = astrowrapper.fits_alignment
rootdir = astrowrapper.rootdir
source_analysis = astrowrapper.source_analysis


# "main" class. Holds "bulk operations" that make use of utilities.
class main_runtime(object):
    # Init Null
    def __init__(self):
        null = "null"
    # Misc test Funct
    def test_function(self):
        print("SEIKO DANA! Main runtime is fine, babe.")
    # Final Function
    def complete(self):
            print("All operations complete, onee-chan.")
    # Generic function to extract background histogram from a FITS file in directory DIRECTORY (relative to rootdir) (straight outta utils)
    def stat_reducer(self, directory, name):
        # Enable utils and navigate
        utility = utils()
        os.chdir(rootdir + directory)
        utility.stat_reduce(name)
    # Trims, bias/darks/flats, calibration, and astrometry.
    def data_calibrator(self, bands):
        # Navigate
        os.chdir(rootdir)
        # Load utilities and hdf5 writer.
        utility = utils()
        # Execute all scripts in order.
        utility.fittrim_recursive('Sci')
        utility.fittrim_recursive('Cali')
        utility.calib_gen(bands)
        utility.recur_cali('Sci')
    # Standalone astrometry recursifier in case of errors above (do not re-run above in the case of astrometry crash.)
    # see astro_metrifier to alter star detection threshold.
    def astro_recursifier(self, directory):
        # Navigate
        os.chdir(rootdir)
        # Load utilities
        utility = utils()
        # Astrometrify
        utility.astro_recursifier(directory)
    # Generate optical depth data for the standard directories provided within 'Sci\STD'.
    # Standards should be the directories, coordinates the "RA,DEC" in normal literature form, aprads/anrads as a list.
    # ClipPER and frac as with utility.airmass_stats
    # Will automatically try to determine the aperture radii and annular radii for each star.
    def airmass_generator(self, standards, coordinates, clip, threshold, clipper, frac, maxdev, aprad, anradin, anradout):
        # Navigate
        os.chdir(rootdir)
        # Load utilities alongside hdf5.
        utility = utils()
        filer = hdf5_writer(rootdir, 'data.hdf5')
        # Make sure the hdf5 is written.
        filer.create()
        # Generate necessary data and write for each entry in standards.
        for num,u in enumerate(standards):
            # Note down the coordinate for this standard star. COORDINATES SHOULD BE IN DEGREES
            coordsave = utility.radec_deg(coordinates[num])
            filer.write(u, 'radec', coordsave)
            # Note down various variables (clip, threshold)
            filer.write(u, 'clip', clip)
            filer.write(u, 'threshold', threshold)
            # Generate endpoints.
            endpoints = utility.dir_sub_lister('Sci\\STD\\' + u)
            # FWHM holder, since we also want to generate FWHM's for statistical analyses.
            # Generate airmass data for each endpoint.
            bands, fwhms, zeniths = [],[],[]
            for v in endpoints:
                # Initial Airmass Solver Data
                indi_band_data = utility.airmass_solver(v, [coordinates[num]], clip, threshold, aprad, anradin, anradout) # depth, depth_err, alt_dif, band, FWHMS, ZENITHS (for this band)
                bands.append(indi_band_data[3]), fwhms.append(indi_band_data[5]), zeniths.append(indi_band_data[6])
                # Record the (unstatistified) depths and their errors.
                filer.write(u, ("depth_values_{}").format(indi_band_data[3]), indi_band_data[0])
                filer.write(u, ("depth_errors{}").format(indi_band_data[3]), indi_band_data[1])
                # Run statistical reduction on the depth data to remove the non-viables
                depth_statistified = utility.airmass_stats(indi_band_data[0], indi_band_data[1], clipper, frac, maxdev)
                # Write the final value of optical depth alongside its error as array [VALUE, ERROR].
                filer.write(u, ("optical_depth_{}").format(indi_band_data[3]), [depth_statistified[3], depth_statistified[5]])
                # Also write down the new "corrected" array values.
                filer.write(u, ("depth_stat_values_{}").format(indi_band_data[3]), depth_statistified[0])
                filer.write(u, ("depth_stat_errors{}").format(indi_band_data[3]), depth_statistified[1])
                # Also, we need to write down the fluxes used for diagnostic purposes (for this endpoint directory)
                filer.write(u, ("depth_fluxvalues_{}.fits").format(indi_band_data[3]), indi_band_data[4])
            # Run fwhmstats to collect statistics on fwhm/etc
            utility.standard_fwhmstats(bands, fwhms, zeniths, standards[num])
    # Generates (and writes to file) all the corrected instrumental magnitudes, errors, zeropoints, errors, for the standards provided (and their true_mags)
    # The mag arrays should be of form [U,B,V]. Each true_mag is stored as 'true_mag_u,b,v' for each STD star.
    def mag_generator(self, standards, true_mag_arrays, aprad, anradin, anradout):
        # Instantiate Utility and file suite
        utility = utils()
        filer = hdf5_writer(rootdir, 'data.hdf5')
        # Note mag types
        mag_types = ['U','B','V']
        # Extract the "true depth" from the chosen standard. This is just hard-coded at the moment for all-intensive-purposes, but could be altered later. Just set "true_depth" values to the values you wish to use in calculations to alter it.
        utility.true_depth("111773")
        # Write the magnitudes to file.
        for num,u in enumerate(standards):
            for bum,v in enumerate(mag_types):
                filer.write(u, 'true_mag_' + v, true_mag_arrays[num][bum])
            # Generate endpoints.
            endpoints = utility.dir_sub_lister('Sci\\STD\\' + u)
            # Run for each endpoint directory the production of various data.
            for gamma in endpoints:
                os.chdir(gamma)
                utility.fits_zp_solver_recursive(u, gamma, aprad, anradin, anradout)
                os.chdir(rootdir)

        # Sort out ZPs
        utility.true_zpers("111773")



    # Various stacking arguments that are currently used.
    def final_stacker(self, clusters):
        # Instantiate utilities and fits handler.
        utilities = utils()
        fitser = fits_alignment(3)

        # This will run for all cluster subdirectories. It will use the "zeroth" image as a reference. All subsequent images will be aligned to it. Aligns images in anticipation of final averaging.
        for cluster in clusters:
            endpoints = utilities.dir_sub_lister(rootdir + '\\Sci\\' + cluster)
            if __name__ == '__main__':
                pool = multiprocessing.Pool(12)
                pool.map(fitser.fits_dir_aligner, endpoints)

        # Iterate for all subdirectories under the clusters (all endpoints except U,B,V_FINALS). Pure average stack for all cluster subdirs.
        """
        for cluster in clusters:
            endpoints = utilities.dir_sub_lister(rootdir + '\\Sci\\' + cluster)
            for point in endpoints:
                os.chdir(point)
                utilities.average_stacker("average_stacked")
        """
        # Iterate for clusters and bands in "cluster\U,B,V_FINALS" directory (assumes all are registered/stacked/etc correctly.
        """
        for cluster in clusters:
            for band in ['U','B','V']:
                os.chdir(rootdir + '\\Sci\\' + cluster + '\\' + band + '_FINALS')
                utilities.average_stacker(cluster + '_' + band)
        """
    # Aligns aligned images.
    def cluster_stacker(self, clusters):
        # Instantiate utils.
        utilities = utils()
        fitser = fits_alignment(3)

        # This will average only the "aligned" files in endpoint directories.
        for cluster in clusters:
            endpoints = utilities.dir_sub_lister(rootdir + '\\Sci\\' + cluster)
            for point in endpoints:
                os.chdir(point)
                fitser.identifier_stacker()
    # Gathers up offsets FROM ORIGINAL CALIBRATED DIRECTORS.(placeholder, add to the offset_stacked in eventuality when confident enough...)
    def offset_gatherer(self, clusters):
        fitser = fits_alignment(3)
        fitser.offset_stackwcs(clusters)
    # Aligns up and stacks up all the offset files.
    def offset_final_aligner(self, clusters):
        fitser = fits_alignment(3)

        # Generates the files
        bands = ['U', 'B', 'V']
        # Offsets
        offsets = np.arange(0,3,1)
        offsets = [str(d) for d in offsets]
        print(offsets)

        end_dirfile_concats = []
        # Generate all endpoint directories to use.
        for cluster in clusters:
            for band in bands:
                directory = rootdir + "\\Sci\\ClusterstacksDone\\" + cluster + "\\" + band + "\\WCSDONE"
                file_list = []
                dir_files = [directory, file_list]
                end_dirfile_concats.append(dir_files)
                for offset in offsets:
                    file = "fully_stacked_" + cluster + "_" + offset + "_" + band + ".fitswcs.fits"
                    file_list.append(file)
        if __name__ == '__main__':
            print("Running.")
            pool = multiprocessing.Pool(6)
            pool.map(fitser.fit_rot_aligner_all, end_dirfile_concats)
    # Navigates to all the final directories for cluster stacks, WCS fits (upload-entire-image-method) and renames file to appropriate name.
    def offset_final_WCS(self, clusters):
        # Utility gen
        utilities = utils()
        # Generate final file structure.
        bands = ["V","B","U"]
        file_identifier = "average_stack_fully_stacked" # The generic identifier for the stack.
        for cluster in clusters:
            for band in bands:
                band_cluster_directory = rootdir + "\\Sci\\ClusterstacksDone\\" + cluster + "\\" + band + "\\WCSDONE\\FinalAligns"

                # We need to also add on the "idi1,idi2" identifiers, see astrowrapper. In my case, it's the cluster/band respectively.
                file_identifier_current = file_identifier + "_" + cluster + "_" + band + ".fits"

                os.chdir(band_cluster_directory)
                utilities.astro_metrifier_entire_file(file_identifier_current)



    # Generate file structure prior to analysis
    def source_tree(self, clusters, bands):
        # Set up class
        anal = source_analysis(rootdir, "data.hdf5")
        # Sort out file tree
        anal.anatree_maker(clusters, bands)
    # Carry out source analysis, including:
    # - identifying sources from the last band provided, i.e. you should provide bands in order of bright to faint
    # - using sources found for photometry in all three bands, same source RA/DEC as first image
    # - correcting for atmospheric absorption effects in all three bands, using the true_depth argument in the data file
    def source_photometry(self, clusters, bands, threshold, array_of_squarecoords, fwhm_guess, minsep_fwhm, sigma_radius, aprad, anradin, anradout, region_clip, savefigs, sharplo, roundlo, sharphi, roundhi, minfluxes):
        # Set up class
        anal = source_analysis(rootdir, "data.hdf5")
        utilities = utils()

        # Run for each cluster.
        for num, cluster in enumerate(clusters):

            # Navigate and find file for initial source extraction (brightest file)
            os.chdir(rootdir + "\\Sci\\SEX\\" + cluster + "\\" + bands[0])
            file_bright = utilities.fits_identifier()[0] # it's alone, or should be.

            # Get sources for this file
            data_group_targets = cluster + "_targets"
            print("Sources extracted.")
            anal.anasource_extraction(file_bright, data_group_targets, array_of_squarecoords[num], threshold, fwhm_guess, minsep_fwhm, sigma_radius, sharplo, roundlo, sharphi, roundhi, minfluxes[num]) # sharplo, roundlo, sharphi, roundhi


            #(for linear/non parallel)
            # Use obtained sources for photometry in all other three bands
            for band in bands:
                # Navigate and find file for photometry
                os.chdir(rootdir + "\\Sci\\SEX\\" + cluster + "\\" + band)
                file_fits = utilities.fits_identifier()[0]  # it's alone, or should be.

                # Also specify the new identifier (for this band) for final photometry data
                data_group_band = cluster + "_" + band

                # Run photometry
                cwd = os.getcwd()

                # For simple apertures (no fancy sigclip bullshit)
                anal.anasource_simpaertures(file_fits, data_group_band, data_group_targets, aprad, anradin, anradout)

                # For manu apertures (fancy sigclip bullshit)... alters the phottable from above.
                os.chdir(cwd)
                anal.anasource_manuapertures_redo(file_fits, data_group_band, data_group_targets, aprad, anradin, anradout, region_clip, savefigs)

                # Atmos Corr

                # For simple apertures (no fancy sigclip bullshit)
                anal.ana_atmos_correction(data_group_band)
    # Parallel atmos correction (test)
    def source_atmoscorr(self, clusters, bands):
        for cluster in clusters:
            for band in bands:
                datagroup_band = cluster + "_" + band 
                anal = source_analysis(rootdir, "data.hdf5")
                anal.ana_atmos_correction(datagroup_band)
    # Quick cleanup
    def source_cleanup(self, clusters, bands):
        anal = source_analysis(rootdir, "data.hdf5")
        for cluster in clusters:
            anal.ana_mag_cleanup_1(cluster, bands)

# FROM HEREON ALL THE CODE WILL BE DIRTY AS FUCK AND WON'T BE INTENDED TO BE REPRODUCEABLE. SORRY NOT SORRY!

    # Generate a menagerie of "colourful" diagrams (This isn't exactly "automated" and is instead fine-tuned. You'll need to edit values as the crow flies unless this becomes pertinent to the final report, in which case you'll want to go ahead and configurability.
    # Some parts may be dependent on previous imports, i.e. GAIA/Hipparcos/Pleiades/etc catalogues.
    # Use misc_file_handler/GAIA to obtain these before attempting to run them.
    def source_diagrams(self):
        sourcer = source_analysis(rootdir, "data.hdf5")
        # sourcer.ana_tricolour("PLEIADES", "phottable", "", ["V", "B", "U"], False, False, False, [-0.25,2.1], [2.7,-1], False, "Sample Data for Pleiades, UBV Colour-Colour Plot",4)

        sourcer.ana_bicolour("GAIA", "test_query", "", ["g","bp"], False, [-0.5, 2], [10, 20], "", False, 1)
        #sourcer.ana_bicolour("M52_PREREDCATALOGUE", "phottable_table_prered", "_apflux_annuli_manu_atmos_scimag", ["V","B"], False, [0.2, 2.2], [10, 18], "", False, 1)
        #sourcer.ana_tricolour("M52_PREREDCATALOGUE", "phottable_table_prered", "_apflux_annuli_manu_atmos_scimag", ["V", "B", "U"], False, False, [0,2], [-6,3], False, "",1)
        # sourcer.ana_tricolour("M52_PREREDCATALOGUE", "phottable_table_prered", "_apflux_annuli_atmos_scimag", ["V", "B", "U"], False, False, False, [0,2], [3,-6], False, "M52 UBV CC Plot Sample",4)
        # sourcer.ana_tricolour("M52_PREREDCATALOGUE", "phottable_table_prered", "_apflux_sig_atmos_scimag", ["V", "B", "U"], False, False, False, [0,2], [3,-6], False, "M52 UBV CC Plot Sample",4)
        # sourcer.ana_tricolour("GAIA", "test_query", "", ["rp", "g", "bp"], False, False, False, [0,2], [3,-6], False, "M52 GAIA 20' test",1)

        sourcer.ana_bicolour("GAIA", "test_query_2", "", ["g","bp"], False, [-0.5, 2], [10, 20], "", False, 1)
        #sourcer.ana_bicolour("NGC7789_PREREDCATALOGUE", "phottable_table_prered", "_apflux_annuli_manu_atmos_scimag", ["V","B"], False, [0, 2], [10, 20], "", False, 1)
        #sourcer.ana_tricolour("NGC7789_PREREDCATALOGUE", "phottable_table_prered", "_apflux_annuli_manu_atmos_scimag", ["V", "B", "U"], False, False, [0,2], [-6,3], False, "",1)
        # sourcer.ana_tricolour("NGC7789_PREREDCATALOGUE", "phottable_table_prered", "_apflux_annuli_atmos_scimag", ["V", "B", "U"], False, False, False, [0,2], [3,-6], False, "NGC7789 UBV CC Plot Sample",1)
        # sourcer.ana_tricolour("NGC7789_PREREDCATALOGUE", "phottable_table_prered", "_apflux_sig_atmos_scimag", ["V", "B", "U"], False, False, False, [0,2], [3,-6], False, "NGC7789 UBV CC Plot Sample",1)
        # sourcer.ana_tricolour("GAIA", "test_query_2", "", ["rp", "g", "bp"], False, False, False, [0,2], [3,-6], False, "NGC7789 GAIA 20' test",1)

        #sourcer.ana_bicolour("HIPPARCOS", "mag_table", "", ["vt", "bt"], False, [-0.5, 2.5], [-5, 15], "", "vtabs", 1)
        #sourcer.ana_bicolour("HIPPARCOS", "mag_table_errored", "", ["vt", "bt"], False, [-0.5, 2.5], [-5, 15], "", "vtabs", 1)
    # For turnoff point reduction/etc.
    # I'm a bit impatient so from hereon it won't be customizable, and you won't be able to alter the names/etc. It's all hard coded.
    # Everything in astrowrapper is generalized but, this won't be.
    # Import all of Syazza's cluster files (contains colour + etc) and all of Meghan/Olivias membership files (has parallax/etc)
    def turn_import(self):
        # Set up all utilities
        handler = astrowrapper.misc_file_handler(rootdir, "data.hdf5")
        filer = hdf5_writer(rootdir, "data.hdf5")
        utilities = utils()
        turnoff_handler = astrowrapper.turnoff_analysis(rootdir, "data.hdf5")

        # Sort out tree
        turnoff_handler.turn_treemaker()

        # Get all the files we need.
        syaz_dir = rootdir + "\\Sci\\TURNOFF\\VOT"

        # Columns to keep. Differs dependent on whose table you're using. Sebi uses ra,dec/ Vinny/Meg use ra_x, dec_x
        cols = ["parallax", "parallax_error", "V_apflux_annuli_manu_atmos_scimag", "B_apflux_annuli_manu_atmos_scimag", "U_apflux_annuli_manu_atmos_scimag"]
        cols_own = ["ra", "dec", "pmra", "pmdec", "parallax", "parallax_error", "V_apflux_annuli_manu_atmos_scimag", "B_apflux_annuli_manu_atmos_scimag", "U_apflux_annuli_manu_atmos_scimag", "phot_g_mean_mag", "phot_bp_mean_mag", "V_apflux_annuli_manu_atmos_scimag_err", "B_apflux_annuli_manu_atmos_scimag_err"]

        # AV + Clusts from Syazzmatazz
        av_value = [1.8955172413793102, 1.1758620689655173]
        clust = ["M52", "NGC7789"]

        # Save the AV values to file.
        filer.write("M52", "AV", av_value[0])
        filer.write("NGC7789", "AV", av_value[1])
        set = "raw_data" # Dataset name for imported tables

        # The entire match.
        files = ["M52_Cross.vot", "NGC7789_Cross.vot"]
        groups_own_more = [d.split(".vot")[0] for d in files]
        for i in files:
            handler.votable_import(syaz_dir, i, i.split(".vot")[0], set, cols_own)

        """

        # The entire match.
        files = ["M52_Cross.vot", "NGC7789_Cross.vot"]
        groups_own_more = [d.split(".vot")[0] for d in files]
        for i in files:
            handler.votable_import(syaz_dir, i, i.split(".vot")[0], set, cols_own)

        # Our stuff
        files = ["M52FullMatch.vot", "NGC7789FullMatch.vot"]
        groups = [d.split(".vot")[0] for d in files]
        for i in files:
            handler.votable_import(syaz_dir, i, i.split(".vot")[0], set, cols)
            
            
        # Classmate stuff
        files_mates = ["M52ClusterMembersMatched13arcminrad.vot", "NGC7789EverythingMatched16arcminrad.vot"]
        groups_mates = [d.split(".vot")[0] for d in files_mates]
        for i in files_mates:
            handler.votable_import(syaz_dir, i, i.split(".vot")[0], set, cols)

        # Larger-radius classmate stuff
        files_mates_wideclust = ["M52ClusterMembersMatched13arcminrad.vot", "NGC7789EverythingMatched16arcminrad.vot"]
        groups_mates_wideclust = [d.split(".vot")[0] for d in files_mates_wideclust]
        for i in files_mates_wideclust:
            handler.votable_import(syaz_dir, i, i.split(".vot")[0], set, cols) """

        # Deredden sets using the given AV values we have.
        av_value = np.array(av_value)
        av_error = np.array([0.04,0.02]) # Error provided by shazmatazz
        EBminusV = av_value/3.1
        EBminusV_error = av_error/3.1

        # Set up groups
        #M52_groups = [groups[0], groups_mates[0], groups_mates_wideclust[0], groups_own_more[0]]
        #NGC7789_groups = [groups[1], groups_mates[1], groups_mates_wideclust[1], groups_own_more[1]]

        M52_groups = [groups_own_more[0]]
        NGC7789_groups = [groups_own_more[1]]

        # For M52
        for group in M52_groups:
            table = filer.read_table(group, set)
            bv = table['B_apflux_annuli_manu_atmos_scimag'] - table['V_apflux_annuli_manu_atmos_scimag']
            table['bv'] = bv
            bv_dereddened = bv - EBminusV[0]
            V_dereddenned = table['V_apflux_annuli_manu_atmos_scimag'] - av_value[0]
            table['V_dered'], table['BV_dered'] = V_dereddenned, bv_dereddened
            V_abs, V_abs_dered = [], []
            BVDEREDERRORS, VABSDEREDERRORS = [],[]
            VABSERRORS = []
            CMDERRS = []
            for i in table:

                # Grab BV/V errors
                bv_error = np.sqrt((i['V_apflux_annuli_manu_atmos_scimag_err'])**2 + (i['B_apflux_annuli_manu_atmos_scimag_err'])**2)
                v_error = i['V_apflux_annuli_manu_atmos_scimag_err']

                # Set up the error in the dereddened apparent V
                V_dered_error = np.sqrt(v_error**2 + (av_error[0])**2)
                BV_dered_error = np.sqrt(bv_error**2 + (EBminusV_error[0])**2)
                BVDEREDERRORS.append(BV_dered_error)

                # Get V_Abs_Dered & V_Abs
                dist, disterr = utilities.par_to_dist(i['parallax'], i['parallax_error'])
                vabs, vabserr = utilities.abs_to_app(i['V_apflux_annuli_manu_atmos_scimag'], v_error, dist, disterr)
                vabsdered, vabsderederr = utilities.abs_to_app(i['V_dered'], V_dered_error, dist, disterr)
                VABSDEREDERRORS.append(vabsderederr)
                V_abs.append(vabs)
                V_abs_dered.append(vabsdered)
                VABSERRORS.append(vabserr)

                # We'll also estimate a distance error for the CMD plot.
                cmderr = np.sqrt(vabsderederr**2 + BV_dered_error**2)
                CMDERRS.append(cmderr)

            table['V_abs'] = V_abs
            table['V_abs_dered'] = V_abs_dered
            table['BV_dered_err'] = BVDEREDERRORS
            table['V_abs_dered_err'] = VABSDEREDERRORS
            table['V_abs_err'] = VABSERRORS
            table['CMDERR'] = CMDERRS
            filer.write_table(group, set, table)



        # For NGC7789
        for group in NGC7789_groups:
            table = filer.read_table(group, set)
            bv = table['B_apflux_annuli_manu_atmos_scimag'] - table['V_apflux_annuli_manu_atmos_scimag']
            table['bv'] = bv
            bv_dereddened = bv - EBminusV[1]
            V_dereddenned = table['V_apflux_annuli_manu_atmos_scimag'] - av_value[1]
            table['V_dered'], table['BV_dered'] = V_dereddenned, bv_dereddened
            V_abs, V_abs_dered = [], []
            BVDEREDERRORS, VABSDEREDERRORS = [], []
            VABSERRORS = []
            CMDERRS = []
            for i in table:
                # Grab BV/V errors
                bv_error = np.sqrt((i['V_apflux_annuli_manu_atmos_scimag_err']) ** 2 + (i['B_apflux_annuli_manu_atmos_scimag_err']) ** 2)
                v_error = i['V_apflux_annuli_manu_atmos_scimag_err']

                # Set up the error in the dereddened apparent V
                V_dered_error = np.sqrt(v_error ** 2 + (av_error[1]) ** 2)
                BV_dered_error = np.sqrt(bv_error ** 2 + (EBminusV_error[1]) ** 2)
                BVDEREDERRORS.append(BV_dered_error)

                # Get V_Abs_Dered & V_Abs
                dist, disterr = utilities.par_to_dist(i['parallax'], i['parallax_error'])
                vabs, vabserr = utilities.abs_to_app(i['V_apflux_annuli_manu_atmos_scimag'], v_error, dist, disterr)
                vabsdered, vabsderederr = utilities.abs_to_app(i['V_dered'], V_dered_error, dist, disterr)
                VABSDEREDERRORS.append(vabsderederr)
                V_abs.append(vabs)
                V_abs_dered.append(vabsdered)
                VABSERRORS.append(vabserr)

                # We'll also estimate a distance error for the CMD plot.
                cmderr = np.sqrt(vabsderederr ** 2 + BV_dered_error ** 2)
                CMDERRS.append(cmderr)

            table['V_abs'] = V_abs
            table['V_abs_dered'] = V_abs_dered
            table['BV_dered_err'] = BVDEREDERRORS
            table['V_abs_dered_err'] = VABSDEREDERRORS
            table['V_abs_err'] = VABSERRORS
            table['CMDERR'] = CMDERRS
            filer.write_table(group, set, table)

    # Generate bicolour images for V-BV.
    # Will also overplot against Hipparcos.
    def turn_bicolours(self):
        turnanalysis = astrowrapper.turnoff_analysis(rootdir, "data.hdf5")
        anal = astrowrapper.source_analysis(rootdir, "data.hdf5")
        clusts = ["M52", "NGC7789"]

        # This is for our custom done field stuff
        groups = ["M52_Cross", "NGC7789_Cross"]  # groups = ["M52_MEMBERS", "NGC7789_MEMBERS"] #
        sets = "raw_data_reduced"  # "corrected_ubv" #
        for num, group in enumerate(groups):
            turnanalysis.turn_bicolour_hipparcos(group, sets, "V_abs_dered", "BV_dered", True, False, False, "", 4, clusts[num])

            """
                        turnanalysis.turn_bicolour(group, sets, "V_abs_dered", "BV_dered", True, "null", "null",
                                       "Absolute Bicolour for matched data, " + clusts[num], 2)
            turnanalysis.turn_bicolour_hipparcos(group, sets, "V_abs_dered", "BV_dered", True, "null", "null",
                                                 "Absolute Bicolour for matched data " + clusts[num], 2, clusts[num])
            """
        """
        # This is for Syaza's stuff
        groups = ["M52_MEMBERS", "NGC7789_MEMBERS"]
        sets = "corrected_ubv"
        for num, group in enumerate(groups):
            turnanalysis.turn_bicolour(group, sets, "V_abs", "BV", True, "null", "null", "Absolute Bicolour for AV deduction data, " + clusts[num], 2)
            turnanalysis.turn_bicolour_hipparcos(group, sets, "V_abs", "BV", True, "null", "null", "Absolute Bicolour for AV deduction of " + clusts[num], 2, clusts[num])




        # This is for the entire field (our own crossmatch)
        groups = ["M52FullMatch", "NGC7789FullMatch"]
        sets = "raw_data"
        for num,group in enumerate(groups):
            turnanalysis.turn_bicolour(group, sets, "V_abs", "bv", True, "null", "null", "Absolute Bicolour for field image of " + clusts[num], 2)
            turnanalysis.turn_bicolour_hipparcos(group, sets, "V_abs", "bv", True, "null", "null", "Absolute Bicolour for field image of " + clusts[num], 2, clusts[num])



        # This is for the wider field data (larger radius, dereddened via our method)
        groups = ["M52ClusterMembersMatched13arcminrad", "NGC7789EverythingMatched16arcminrad"]
        sets = "raw_data"
        for num,group in enumerate(groups):
            turnanalysis.turn_bicolour(group, sets, "V_abs_dered", "BV_dered", True, "null", "null", "Absolute Bicolour for mag-non-limited " + clusts[num], 2)
            turnanalysis.turn_bicolour_hipparcos(group, sets, "V_abs_dered", "BV_dered", True, "null", "null", "Absolute Bicolour for mag-non-limited " + clusts[num], 2, clusts[num])
        """
    # Generate a (modelled) logL/logT graph
    def turn_loglog(self):
        # Grab turnanalysis
        turnanal = astrowrapper.turnoff_analysis(rootdir, "data.hdf5")

        # Groups and ID's for OUR CROSSMATCH + REDUCED DATA
        groups = ["M52_Cross", "NGC7789_Cross"]  # groups = ["M52_MEMBERS", "NGC7789_MEMBERS"] #
        set = "raw_data_reduced"  # "corrected_ubv" #
        absV, BV = "V_abs_dered", "BV_dered"  # BV
        lims = [[[3, 4.5], [-2, 6]],[[3, 4.5], [-2, 6]]]
        lims_ubv = [[[-1, 3], [-4, 8]],[[-1, 3], [-4, 8]]]
        age_guess = "7.80", "9.14"
        #metallicity_guess = "0.00", "-0.30"
        metallicities = ["0.05", "-0.30"] # ["0.05", "0.00", "-0.05", "-0.20", "-0.25",
        maxshift = [0, 0]
        shift_increment = 0.02
        age_ranges = [[7.6, 8.20], [9.0, 9.30]]
        age_increment = [0.001, 0.001]
        for num, group in enumerate(groups):
            #turnanal.turn_mspectype_estimator(group, set, absV, BV, lims[num], age_ranges[num], age_increment[num], maxshift[num], shift_increment)
            #turnanal.turn_mspectype_estimator_UBV(group, set, absV, BV, lims_ubv[num], age_ranges[num], age_increment[num], maxshift[num], shift_increment)
            if num != 1:
                turnanal.turn_mspectype_estimator_UBV_CMD(group, set, absV, BV, lims_ubv[num], age_ranges[num], age_increment[num], metallicities[num])
        """
        # Groups and ID's for FULLMATCH
        groups = ["M52FullMatch", "NGC7789FullMatch"] #groups = ["M52_MEMBERS", "NGC7789_MEMBERS"] #
        set = "raw_data" # "corrected_ubv" #
        absV, BV = "V_abs_dered", "BV_dered" # BV
        for group in groups:
            turnanal.turn_mspectype_estimator(group, set, absV, BV)

        # Groups and ID's for Syazamatch
        groups = ["M52_MEMBERS", "NGC7789_MEMBERS"] #groups = ["M52_MEMBERS", "NGC7789_MEMBERS"] #
        set = "corrected_ubv"
        absV, BV = "V_abs", "BV" # BV
        for group in groups:
            turnanal.turn_mspectype_estimator(group, set, absV, BV)"""


    # Membership stuff
    def member(self):
        finder = astrowrapper.member_finder(rootdir, "data.hdf5")
        filer = astrowrapper.hdf5_writer(rootdir, "data.hdf5")
        # Plot PM plots


        groups = ["M52_Cross", "NGC7789_Cross"]
        sets = "raw_data"
        """
        titles = ["M52", "NGC7789"]
        pmlims = [[[-2.7,-1.25],[-0.6,-1.6]],[[-1.5,-0.4],[-1.5,-2.35]]]
        for num, group in enumerate(groups):
            finder.pmscat(group, sets, 'pmra', 'pmdec', pmlims[num], titles[num])
            #finder.radecscat(group, sets, 'ra', 'dec', False,titles[num])
            #finder.radcalc(group, sets, 12.5, 0.25, 0.4)
        """
        # This is just for some fun. Trialling the second cluster in the image for M52. Misc stuff.
        # Obtain the GAIA data for the region
        #gaia = astrowrapper.GAIA(rootdir, "data.hdf5")
        #filer = hdf5_writer(rootdir, "data.hdf5")
        #racent, deccent = filer.read(groups[0], "centroid_second")
        #gaia.obtain_data(60, racent, deccent, "SECOND_TEST")

        # Generate proper motion/ra/dec data using the GAIA set instead of ours.
        #finder.pmscat3("GAIA", "SECOND_TEST", "pmra", "pmdec")
        #finder.radecscat3("GAIA", "SECOND_TEST", "ra", "dec")
        #finder.m52third("GAIA", "SECOND_TEST", 7, 0.3, 0.4)


        finder.m52second(groups[0], sets, 5, 0.2, 1)
        finder.pmscat2(groups[0], sets, "pmra", "pmdec")
        finder.radecscat2(groups[0], sets, "ra", "dec")
        """
        newsets = "raw_data" + "_second_reduced"
        anal = astrowrapper.source_analysis(rootdir, "data.hdf5")
        #anal.ana_bicolour(groups[0], newsets, "_apflux_annuli_manu_atmos_scimag", ["V", "B"], True, False, False,"Test for " + groups[0] + "_second", False, 8)
        anal.ana_tricolour(groups[0], newsets, "_apflux_annuli_manu_atmos_scimag", ["V", "B", "U"], True, False, False, False, False, "Test for " + groups[0] + "_second", 8)


        # Av is 2.232 (ish) so we do a quick deredden
        utilities = utils()
        m52secondtable = filer.read_table("M52_Cross", "raw_data_second_reduced")
        V, B = m52secondtable['V_apflux_annuli_manu_atmos_scimag'], m52secondtable['V_apflux_annuli_manu_atmos_scimag']
        BV = B - V
        Av = 1.9817241379310346
        V_dered = V - Av
        BV_dered = BV - (Av/3.1)
        V_abs_dered = []
        for num,item in enumerate(V_dered):
            vabs = utilities.abs_to_app(item, 0, m52secondtable[num]['dist'], 0)
            V_abs_dered.append(vabs[0])
        m52secondtable['BV_dered'], m52secondtable['V_abs_dered'] = BV_dered, V_abs_dered


        """
        #turner = astrowrapper.turnoff_analysis(rootdir, "data.hdf5")
        #turner.turn_bicolour_hipparcos("M52_Cross", "raw_data_second_reduced", "V_abs_dered", "BV_dered", True, True, True, "Test for 2nd clust in M52", 16, "M52_Second")
        #metallicity_tests = ["0.00"]
        #for met in metallicity_tests:
        #    turner.turn_mspectype_estimator_UBV_CMD("M52_Cross", "raw_data_second_reduced", "V_abs_dered", "BV_dered", [[-0.5,0.5],[-2,3]], [6.6,9.3], 0.01,  met)


    def isofitter(self):
        groups = ["M52_Cross", "NGC7789_Cross"]
        sets = "raw_data_reduced"
        turnanal = astrowrapper.turnoff_analysis(rootdir, "data.hdf5")

        lims = [[0, 1.5], [-2, 8]]
        turnanal.isofit(groups[1], sets, "V_abs_dered", "BV_dered", 9.17,lims)

        lims = [[-1, 2], [-5, 5]]
        turnanal.isofit(groups[0], sets, "V_abs_dered", "BV_dered", 7.7,lims)
        #turnanal.isofit(groups[0], sets, "V_abs_dered", "BV_dered", 8.3, lims)
        #turnanal.isofit(groups[0], sets, "V_abs_dered", "BV_dered", 8.25, lims)


    # Misc testing
    def misc_tester(self):
        groups = ["M52_Cross", "NGC7789_Cross"]
        sets = "raw_data_reduced"
        ana = astrowrapper.source_analysis(rootdir, "data.hdf5")
        for group in groups:
            ana.ana_tricolour(group, sets, "_apflux_annuli_manu_atmos_scimag", ["V", "B", "U"], False, True, False, False, False, False, "Test  for " + group, 1)
    # Various deprecated tools used for stacking and the like.
    """
        # Runs wcs_finish.
        def offset_wcsfinish(self, clusters):
            fitser = fits_alignment(3)
            fitser.offset_wcsfinish(clusters)
        """
    """
        # Stacks all final offset stacks. Don't recommend running this.
        def offset_stacker(self, clusters):
            # Utils, writer
            utilities = utils()
            fitser = fits_alignment(3)

            # Run recursive
            fitser.offset_recursifier(clusters)
        """

    # Generate IMF/etc for the stuff
    def imf(self):
        turnoff = astrowrapper.turnoff_analysis(rootdir, "data.hdf5")
        groups = ["M52_Cross", "NGC7789_Cross"]
        nbins = [200, 200]
        set = "raw_data_reduced"
        alphas = [2.65, 2.35]
        mass_vals = [0.5,0.5]
        for num, group in enumerate(groups):
            turnoff.turn_imf(group, set, "V_abs_dered", [0,20],nbins[num], alphas[num], mass_vals[num])

    # Generate IMF/etc for the stuff
    def ilf(self):
        turnoff = astrowrapper.turnoff_analysis(rootdir, "data.hdf5")
        groups = ["M52_Cross", "NGC7789_Cross"]
        set = "raw_data_reduced"
        for group in groups:
            turnoff.turn_ilf(group, set, "V_abs_dered", [0,4],100)

# Set up main runtime environment
mainer = main_runtime()


# DATA CALIBRATION AND AIRMASS CALIBRATION GENERATOR/ETC.
# Note: APRAD ANRADIN ANRADOUT HAVE BEEN DEFINED AS THE RADIUS OF THE GAUSSIAN KERNEL, NOT THE DIAMETER!
#mainer.test_function()
#mainer.data_calibrator(['U', 'B', 'V'])
#mainer.astro_recursifier(rootdir + "\\Sci\\STD")
#mainer.airmass_generator(['111773','SA2039'], ["19:37:15.83, +00:10:59.5", "00:45:34.14, +45:36:48.0"], 40, 80, [0.1,0.7], 0.1, 0.1, 3.5/2, 3.5/2, 4.5/2)
#mainer.mag_generator(['111773','SA2039'],[[8.963 + 0.206 + -0.21, 8.963 + 0.206, 8.963],[9.353 + 0.472 + 0.092, 9.353 + 0.472, 9.353]], 3.5/2, 3.5/2, 4.5/2)

# GENERATES ALL THE OFFSET STACKED IMAGES (RAWS)
#mainer.final_stacker(['NGC7789','M52'])
#mainer.cluster_stacker(['NGC7789','M52'])
#mainer.offset_gatherer(['NGC7789','M52'])
#mainer.offset_final_aligner(['NGC7789','M52'])
#mainer.offset_final_WCS(['NGC7789','M52'])

# PHOTOMETRY AND SOURCE HANDLING
"""
The roundness statistic is defined as follows roundness = (hx - hy) / (hx + hy) where and hx and hy are the amplitudes of the best fitting 1-D gaussian
of fwhmpsf equal to the current value of the fwhmpsf parameter in pixels
at the position of a detected object. For completely round stars hx ~= hy
and roundness is ~= 0.0. Objects elongated in x have roundness < 0.0 and
objects elongated in y have roundness > 0.0. This statistic is good for
eliminating objects which are elongated along the image axes, e.g. bad
rows and columns. The sharpness parameter is defined as sharpness = (data(max) - <data>nei) / h where h is the amplitude of the best fitting 2-D gaussian of fwhmpsf
equal to the current value of the parameter fwhmpsf, data(max) is the data
maximum of the detected objects and <data>nei is the mean data value
of the surrounding neighbours. For completely round gaussian stars
and nsigma (which along with fwhmpsf determines the amount of
data used to compute all these things) equal to its default value of 1.5,
sharpness is ~=.66. "Hot" pixels have sharpness values > 1, whereas
cold pixels have low sharpness values ~=0.0. Hope this answers your questions. I have improved the manual
pages in the new version of daophot, and am in the process of
completing a reference manual for the package which discusses
some of these things in more detail. You should also look at the original
daophot paper (Stetson, 1987, PASP).
Lindsey Davis"""
#mainer.source_tree(["M52","NGC7789"], ["V","B","U"]) # minsep_fwhm, sigma_radius, aprad, anradin, anradout):
#mainer.source_photometry(["M52","NGC7789"], ["V","B","U"], 20, [["1031, 80", "3955, 3826"], ["120, 113", "3825, 3899"]], 3, 3, 3, 2.5, 3, 4, 20, True, 0, 0, 9999, 0.3, [200,200,200]) # sharplo, roundlo, sharphi, roundhi, minfluxes
#mainer.source_cleanup(["M52","NGC7789"], ["V","B","U"])
#mainer.source_diagrams()


#mainer.turn_import()
#mainer.turn_bicolours()
mainer.member()
#mainer.turn_loglog()
#mainer.isofitter()
#mainer.misc_tester()
#mainer.imf()
#mainer.ilf()

mainer.complete()

