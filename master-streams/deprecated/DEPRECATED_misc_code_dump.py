"""

writer = hdfutils.hdf5_writer(windows_directories.datadir, ascii_info.asciiname)
bhb_data = writer.read_table(ascii_info.bhb, ascii_info.set_raw)
# Grab just the ra/dec/etc for one row
row = bhb_data[1000]
vlos, dist = row['vlos', 'dist']
l,b,dmul,dmub = row['l'],row['b'],row['dmu_l'],row['dmu_b']
x,y,z,vx,vy,vz = row['x','y','z','vx','vy','vz']

uwu = galconversion(sourcedir)
uwu.solinfo_grab("solar_info.dat")
uwu.solgal_set()

original_galactic_coords = [l,b,dist,dmul,dmub,vlos]
original_galactic_coords[3] = original_galactic_coords[3]/math.cos(math.radians(original_galactic_coords[1]))
original_cartesian_coords = [x,y,z,vx,vy,vz]
manu_galactic_coords = uwu.manual_vec_GALCENT_to_GAL(*original_cartesian_coords)
auto_galactic_coords = uwu.vec_GALCENT_to_GAL(*original_cartesian_coords)
print(original_galactic_coords)
print(manu_galactic_coords)
print(auto_galactic_coords)


"""

"""
# NGC 288 from Vizier. ICRS.
#radec = [13.1885, -26.5826111]
#pmradec = [4.22, -5.65]
vlos = -44.45
distance = 8.988
lb = [151.285,-89.3804]
dmulb = [6.35972,-3.06293]
uwa = galconversion(sourcedir)
uwa.solinfo_grab("solar_info_jorge.dat")
uwa.solgal_set()
print(*lb, distance, *dmulb, vlos)
owo = uwa.vec_GAL_to_GALCENT(*lb, distance, *dmulb, vlos)
print(owo)
owo = uwa.vec_GALCENT_to_GAL(*owo)
print(owo)"""
