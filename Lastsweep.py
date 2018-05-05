#%%
import numpy as np
import scipy.constants as const
import matplotlib as matlib
import matplotlib.pyplot as mplot
from pylab import *
import Laser
import SpaceSolver
import Plasmaunit as punit
import Rampfunctions
import Ionization as Ion
import matplotlib.ticker as ticker
mplot.rcParams['mathtext.fontset'] = 'cm'

#%%

c = const.speed_of_light
epsilon = const.epsilon_0
q_e = const.elementary_charge
m_e = const.electron_mass
bohr = const.value("Bohr radius")
eV_to_joule = const.value("electron volt-joule relationship")
hbar = const.hbar

#%%

fig, ax=mplot.subplots()
r = np.arange(-6,6,1e-1)
Coulomb = -(1/(r)**2)+1
mplot.axis([r[0],r[-1],-0.5,2])
mplot.plot(r,Coulomb,'-.',linewidth=2,label=r'$V_C$')

x = np.array([r[0],r[-1]])
y = np.array([0,0])
ax.plot(x,y,':',linewidth=2)

linje = 0.9-0.2*r
ax.plot(r,linje+0.1,'--',linewidth=2,label=r'$V_L$')

linje2 = linje*Coulomb
ax.plot(r,linje2,'black',linewidth=2)

mplot.xlabel(r'$r$', fontsize=20)
#mplot.ylabel(r'$E$', fontsize=20)
mplot.title(r'$\mathrm{Potential\ barrier}$',fontsize=20)
legend = ax.legend(loc='upper right', shadow=True, fontsize=20)

mplot.text(-7,0,r'$E_A$',fontsize=20,rotation=90)
mplot.text(0.05,0.8,r'$(\mathrm{a})$',fontsize=20,transform=ax.transAxes)

mplot.xticks([])
mplot.yticks([])

mplot.savefig("Tunneling_ionization.pdf",format="pdf",bbox_inches="tight")

