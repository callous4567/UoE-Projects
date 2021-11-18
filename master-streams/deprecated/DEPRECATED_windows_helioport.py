import ascii_info, windows_directories, galcentricutils, hdfutils



# This is deprecated: we are dealing with Galactic-Galactocentric conversions, not ICRS-Galactocentric conversions :/

# Set up handler for conversion
galcent = galcentricutils.galconversion(windows_directories.sourcedir)

# First remove Jorge's conversion from the data, moving back to ICRS
galcent.solinfo_grab("solar_info_jorge.dat")
galcent.solgal_set()
# For each original ascii
groups = ascii_info.all_groups
for group in groups:
    galcent.GALCENT_to_ICRS(windows_directories.datadir,
                            ascii_info.asciiname,
                            group,
                            ascii_info.set_raw)
# For the "full" ascii
galcent.GALCENT_to_ICRS(windows_directories.datadir,
                        ascii_info.asciiname,
                        ascii_info.fullgroup,
                        ascii_info.fullset)

# This gives us back the original ICRS data. Next: update Galactocentric using our own conversion (updated values.)
galcent.solinfo_grab("solar_info.dat")
galcent.solgal_set()
# For each original ascii
groups = ascii_info.all_groups
for group in groups:
    galcent.ICRS_to_GALCENT(windows_directories.datadir,
                            ascii_info.asciiname,
                            group,
                            ascii_info.set_raw)
# For the "full" ascii
galcent.ICRS_to_GALCENT(windows_directories.datadir,
                        ascii_info.asciiname,
                        ascii_info.fullgroup,
                        ascii_info.fullset)


