from matplotlib import pyplot as plt, rc

import ascii_info
import hdfutils
import windows_directories
# Enable TeX
plt.rcParams['animation.ffmpeg_path'] = 'C:\\Users\\Callicious\\Documents\\Prog\\pycharm\\venv\\ffmpeg\\bin\\ffmpeg.exe'
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

table = hdfutils.hdf5_writer(windows_directories.datadir,
                             ascii_info.asciiname).read_table(ascii_info.fullgroup,
                                                              ascii_info.set_raw)

# Create axis
fig, axs = plt.subplots(nrows=1,ncols=1,figsize=(7,4))
axs.set(xlim=[-180,180])
axs.set(ylim=[-90,90])


table['l'] = [d - 360 if d > 180 else d for d in table['l']]
axs.scatter(table['l'], table['b'], color='red', marker='s', s=5)

# Grid-up
# plt.legend(handles=legend_elements, loc='upper right')
axs.grid(which='major', color='pink')
axs.set_facecolor("k")
axs.set(xlabel="l / deg",
        ylabel="b / deg")
plt.gca().invert_xaxis()
savepath = windows_directories.imgdir + "\\all_sky_coverage.png"
if savepath != None:
    try:
        plt.savefig(savepath, dpi=300)
    except:
        plt.savefig(savepath)
plt.show(dpi=200)