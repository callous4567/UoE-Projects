import ascii_info_new, windows_directories_new, galcentricutils_new, hdfutils



# This is deprecated: we are dealing with Galactic-Galactocentric conversions, not ICRS-Galactocentric conversions :/

# Set up handler for conversion
galcent = galcentricutils_new.galconversion(windows_directories_new.sourcedir)

# First remove Jorge's conversion from the data, moving back to ICRS
galcent.solinfo_grab("solar_info_jorge.dat")
galcent.solgal_set()
# For each original ascii
groups = ascii_info_new.all_groups
for group in groups:
    galcent.GALCENT_to_ICRS(windows_directories_new.datadir,
                            ascii_info_new.asciiname,
                            group,
                            ascii_info_new.set_raw)
# For the "full" ascii
galcent.GALCENT_to_ICRS(windows_directories_new.datadir,
                        ascii_info_new.asciiname,
                        ascii_info_new.fullgroup,
                        ascii_info_new.fullset)

# This gives us back the original ICRS data. Next: update Galactocentric using our own conversion (updated values.)
galcent.solinfo_grab("solar_info.dat")
galcent.solgal_set()
# For each original ascii
groups = ascii_info_new.all_groups
for group in groups:
    galcent.ICRS_to_GALCENT(windows_directories_new.datadir,
                            ascii_info_new.asciiname,
                            group,
                            ascii_info_new.set_raw)
# For the "full" ascii
galcent.ICRS_to_GALCENT(windows_directories_new.datadir,
                        ascii_info_new.asciiname,
                        ascii_info_new.fullgroup,
                        ascii_info_new.fullset)


