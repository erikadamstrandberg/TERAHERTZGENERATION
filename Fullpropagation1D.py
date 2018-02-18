import numpy as np
import scipy.constants as const
import matplotlib as matlib
import matplotlib.pyplot as mplot
from pylab import *
import Laser
import Ionization as Ion
import SpaceSolver
import Plasmaunit

#%% Full propagation! With W1 and W2

TIME = 3000
SIZE = 1500
dim = [TIME,SIZE]
E = np.zeros(SIZE)
Etemp = np.zeros(SIZE)
B = np.zeros(SIZE)
J = np.zeros(SIZE)
n = np.zeros(SIZE)
ntemp = np.zeros(SIZE)
W1 = np.zeros(SIZE)
W2 = np.zeros(SIZE)
Nat = np.ones(SIZE)
Ni1 = np.zeros(SIZE)
Ni1tot = np.zeros(dim)
Ni1temp = np.zeros(SIZE)
Ni0 = np.ones(SIZE)
Ni0tot = np.zeros(dim)
Ni0temp = np.zeros(SIZE)
ne = np.zeros(SIZE)
netemp = np.zeros(SIZE)
netot = np.zeros(dim)
Etot = np.zeros(dim)
Btot = np.zeros(dim)
Jtot = np.zeros(dim)
ntot = np.zeros(dim)
W1tot = np.zeros(dim)
k = np.arange(SIZE)

n[150:300] = 100
dt = 0.01
dz = 0.01
nu = 0

PULSELENGTH = 500
PULSESTART = 0
f = 5e15
OMEGAREAL = 2*np.pi*f
OMEGAPRIM = 25                  # this is the plasma omega, use this everywhere in the code
OMEGA_0 = OMEGAREAL/OMEGAPRIM   # this is the arbitrary omega, use this as argument in punits
VAR = 0.1
E0 = 0.063
Laser.Gauss_forward(E,B,E0,PULSELENGTH,PULSESTART,OMEGAPRIM,VAR,dt)

PLASMASTART = 500
PLASMASTOPP = 1000

for i in range(1,TIME):
    
    E = SpaceSolver.E(E,B,J,dt,dz)
    B = SpaceSolver.B(E,B,dt,dz)
    
    for z in range(len(E)):
        if abs(E[z]) == 0:
            W1[z] = 0
            W2[z] = 0
        else:
            if z > PLASMASTART and z < PLASMASTOPP:
                W1[z] = Ion.Landau(E[z],OMEGA_0,1,dt)
                W2[z] = Ion.Landau(E[z],OMEGA_0,2,dt)
            else:
                W1[z] = 0
                W2[z] = 0

    
    for z in range(len(E)):
        Ni0[z] = Nat[z]-Ni1[z]
        Ni1[z] = (Ni1[z]*(1-(dt/2)*W2[z])+(dt/2)*W1[z]*(Ni0[z]+Ni0temp[z]))/(1+(dt/2)*W2[z])
        ne[z] = 1*Ni1[z]
        Ni0temp[z] = Ni0[z]

    J = SpaceSolver.J(E,J,ne,ntemp,nu,dt,dz)
    netemp = ne
    
    
    Etot[i] = E
    Jtot[i] = J
    W1tot[i] = W1
    netot[i] = ne
    Ni0tot[i] = Ni0
    Ni1tot[i] = Ni1
    


#%%

plotz = np.arange(SIZE)
plott = np.arange(TIME)

t = 500
z = 100
mplot.plot(plotz,Etot[t]*1,'b')
#mplot.plot(k,W1tot[t]*1e-2)
mplot.plot(plotz,netot[t]*1,'r')
mplot.plot(plotz,Jtot[t]*300,'y')
#mplot.plot(plott,Etot[:,300])
#mplot.plot(plott,W1tot[:,20])
#mplot.plot(plott,Ni1tot[:,z])

#%%

Eefter = Etot[0:1500,2*PULSELENGTH+10]
#plotz = np.arange(len(Eefter))
#mplot.plot(plotz,Eefter)
Efft = np.fft.fft(Eefter)
mplot.plot(plotz,np.abs(Efft))
Efore = Etot[0:1500,PULSELENGTH-1]
plotz = np.arange(len(Efore))
Efft = np.fft.fft(Efore)
mplot.plot(plotz,np.abs(Efft))

#%% Real units

nreal = Plasmaunit.nreal(n[200],OMEGAUNIT)
print(nreal)
omegareal = Plasmaunit.omegareal(OMEGA,OMEGAUNIT)
print(omegareal)