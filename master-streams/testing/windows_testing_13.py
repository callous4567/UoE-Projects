from galpy.potential import plotRotcurve, MWPotential2014
from galpy.potential.Irrgang13 import Irrgang13I
from galpy.potential.McMillan17 import McMillan17
from matplotlib import pyplot as plt
from matplotlib.pyplot import legend

plotRotcurve(MWPotential2014,label=r'$\mathrm{MWPotential2014}$',ro=8.,vo=220.) # need to set ro and vo explicitly, because MWPotential2014 has units turned off
plotRotcurve(McMillan17,overplot=True,label=r'$\mathrm{McMillan\, (2017)}$')
plotRotcurve(Irrgang13I,overplot=True,label=r'$\mathrm{Irrgang\ et\ al.\, (2017), model\ I}$')
plt.show()