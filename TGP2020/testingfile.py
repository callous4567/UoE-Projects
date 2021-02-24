from TGP2020 import astrowrapper
import numpy as np
import matplotlib.pyplot as plt
import shapely.geometry as geom
from scipy.interpolate import interp1d
from astropy.table import Table
import copy
filer = astrowrapper.hdf5_writer(astrowrapper.rootdir, "data.hdf5")

# Grab Pleiades + Meanstars
pleiades = filer.read_table("PLEIADES", "phottable")
meanstars = filer.read_table("MEANSTARS", "spec_table")

# Grab our Av Data (just for M52, as example)
groups = ["M52_Cross", "NGC7789_Cross"]
#magclip = 100
# Bands
bands = ["U", "B", "V"]
band_delimiter = "_apflux_annuli_manu_atmos_scimag"
band_in_table = [d + band_delimiter for d in bands]
tab1, tab2 = [filer.read_table(d, "raw_data_reduced") for d in groups]
tab1.sort("U_apflux_annuli_manu_atmos_scimag", reverse=True), tab2.sort("U_apflux_annuli_manu_atmos_scimag", reverse=True)
#tab1, tab2 = tab1[np.arange(0, magclip, 1)], tab2[np.arange(0, magclip, 1)]
print(tab1['U_apflux_annuli_manu_atmos_scimag'])
tables = [tab1, tab2]
UBVs = [[np.array(table[d]) for d in band_in_table] for table in tables]
UBBVs = [[d[0] - d[1], d[1] - d[2]] for d in UBVs]

# U-B, B-V
# plei_UB, plei_BV = pleiades['U'] - pleiades['B'], pleiades['B'] - pleiades['V']
plei_UB, plei_BV = UBBVs[0]


# FOR SHAZZLES
from MeanStars import MeanStars
MSdata = MeanStars().data
meanUB, meanBV = [],[]
for i in MSdata:
    if i['U-B'] != np.nan:
        if i['B-V'] != np.nan:
            meanUB.append(i['U-B']), meanBV.append(i['B-V'])

"""

plt.scatter(plei_BV, plei_UB, color="black", label="Pleiades", s=0.5)
plt.plot(meanBV, meanUB, color="red", label="Meanstars")
plt.xlabel("B - V")
plt.ylabel("U - B")
plt.gca().invert_yaxis()
plt.legend()
plt.title("Pleiades vs. Rochester Colour-Colour")
plt.show() """

# Chi-Sq fit to get E(B-V) and Av
R = 3.1 # (Cardelli et al.)
meanbvub = np.array([meanBV, meanUB]).T # ([x,y], [x,y]....etc)
string = geom.LineString(meanbvub) # The linestring for the meanstars catalogue.

# U B V in nanometres.
U, B, V = 370, 430, 550 # Notes
UBdivBV = (V*((B/U) - 1)) / (B*((V/B) - 1)) # E(U-B)/E(B-V)

shift_range = [0, 1]
shift_nums = 1000
shift_vals = np.linspace(shift_range[0], shift_range[1], shift_nums)

chisq_list = []
for shift in shift_vals:
    plei_BV_new = copy.deepcopy(plei_BV) - shift
    plei_UB_new = copy.deepcopy(plei_UB) - UBdivBV*shift
    coords = np.array([plei_BV_new, plei_UB_new]).T
    summation = 0
    for coord in coords:
        point = geom.Point(coord)
        dist = point.distance(string)

        # DISTANCE WEIGHTINGS! Dist <= 0.1, 0.5, 2, 4
        dist_weights = 0, 0.1, 12000, 0
        dist_lims = 0, 0.2, 0.3, 10, 100
        lens = len(dist_lims) - 1
        for uwa in range(lens):
            if dist >= dist_lims[uwa]:
                pass
            elif dist <= dist_lims[uwa + 1]:
                if dist <= np.inf:
                    summation += (dist * dist_weights[uwa])
                else:
                    summation += 32

        #if dist <= 10:
        #    summation += dist**2
    chisq_list.append(summation)

argmin = np.argmin(chisq_list)
shift_min = shift_vals[argmin]


plt.scatter(plei_BV, plei_UB, color="black", label="Pleiades", s=0.25)
plt.scatter(plei_BV - shift_min, plei_UB - UBdivBV*shift_min, color="blue", marker="x", label="Pleiades shifted", s=2)
plt.plot(meanBV, meanUB, color="red", label="Meanstars")
plt.xlabel("B - V")
plt.ylabel("U - B")
plt.gca().invert_yaxis()
plt.legend()
plt.title("Pleiades vs. Rochester Colour-Colour")
plt.show()

av = shift_min*R

plt.plot(shift_vals, chisq_list)
plt.xlabel("E(B-V)")
plt.ylabel("Chisq")
plt.title("For Pleiades, Av = " + str(av))
plt.show()