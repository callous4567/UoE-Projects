import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as const
from matplotlib.patches import Patch

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 12}
import matplotlib as mpl
mpl.rc('font', **font)
# Set up constants
G = const.G
M = 5.972e24
R = 6371e3

# Set up plot to work off
minmax = [R, 10*R]
x = np.linspace(*minmax,1000)
y = -G*M/x

# Pick example points you're going to be working with
random_r = np.sort(np.random.default_rng().integers(*minmax, size=5))
random_V = -G*M/random_r

"""
V = GM/r: ratio V_1/V_2 = r_2/r_1
We should be able to "predict" r_2 based on V_1/V_2 via
r_2 = (V_1/V_2)*r_1
Try this: if it works, it most certainly is a 1/r plot.
"""
# Predict r_2 based on V(r), assuming V(r) proportionate to 1/r. Redo V for this, too.
random_r2 = np.array([random_V[0]*random_r[0]/random_V[d] for d in range(len(random_V))])
random_V2 = -G*M/random_r2

# Set up graphing environment/etc
fig, ax = plt.subplots(nrows=1,ncols=1,squeeze=True)
ax.set(xlabel=r'$r$',
       ylabel=r'$V(r)$')
ax.plot(x,y,color="black",lw=1)
ax.scatter(random_r,random_V,marker="o",s=50,color="blue")
ax.scatter(random_r2,random_V2, marker="x",s=25, color="yellow")
legend_elements = [Patch(edgecolor="black",facecolor="blue",label="True"),
                   Patch(edgecolor="black",facecolor="yellow", label="Guess")]
plt.legend(loc="lower right",handles=legend_elements)
plt.savefig("TSRtest.png", dpi=300)
plt.show()