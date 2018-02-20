# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 09:15:23 2018

@author: Hanna
"""
import numpy as np
import matplotlib.pyplot as mplot
import Ionization as ion
import Pulse 
import Plasmaunit as punit
import scipy.constants as const
#%% REAL PARAMETERS
c = const.speed_of_light                # [m/s]
epsilon = const.epsilon_0               # [F/m]
LAMBDA = 800e-9                         # [m] (800 nm)
f = c/LAMBDA
OMEGAREAL = 2*np.pi*f       
t0REAL = 50e-15                         # [s] (50 fs)
I0 = 4e18                               # [W/m^2] (4e14 W/cm^2)
E0REAL = np.sqrt(2*I0/(epsilon*c))      # [V/m]
NatREAL = 3e25                          # [1/m^3] (3e19 [1/cm^3])

# PARAMERERS IN PLASMA UNITS
tmin = 0
tmax = 1000
dt = 0.1

T = np.arange(tmin,tmax,dt)

SIZE = len(T)
#error = np.zeros(SIZE)

OMEGAPRIM = 1                        # this is the plasma omega, use this everywhere in the code
OMEGA = OMEGAREAL/OMEGAPRIM          # this is the omvandlings omega, use this as argument in punits
E0 = punit.Eplasma(E0REAL,OMEGA)
t0 = punit.tplasma(t0REAL,OMEGA)
    
E = E0*Pulse.Gauss(SIZE,OMEGAPRIM,t0,dt)
Eenv = E0*Pulse.Envelope(SIZE,t0,dt)
    
mplot.plot(T,E)

ne = np.zeros(SIZE)   # elektroner
Ni0 = np.zeros(SIZE)  # joner +0
Ni1 = np.zeros(SIZE)  # joner +1
Ni2 = np.zeros(SIZE)  # joner +2
Ni3 = np.zeros(SIZE)  # joner +3
Nat = punit.nplasma(NatREAL,OMEGA)
Ni0[0] = Nat
Ni0[1] = Nat

W1 = np.zeros(SIZE)
W2 = np.zeros(SIZE)
W3 = np.zeros(SIZE)
W4 = np.zeros(SIZE)

for t in range(len(T)):
    W1[t] = ion.Landau(E[t],OMEGA,1)
    W2[t] = ion.Landau(E[t],OMEGA,2)
    W3[t] = ion.Landau(E[t],OMEGA,3)
    
for t in range(1,len(T)):
    Ni1[t] = (Ni1[t-1]*(1-(dt/2)*W2[t-1])+(dt/2)*W1[t-1]*(Ni0[t]+Ni0[t-1]))/(1+(dt/2)*W2[t-1])
    Ni2[t] = (Ni2[t-1]*(1-(dt/2)*W3[t-1])+(dt/2)*W2[t-1]*(Ni1[t]+Ni1[t-1]))/(1+(dt/2)*W3[t-1])
    Ni3[t] = (Ni3[t-1]*(1-(dt/2)*W4[t-1])+(dt/2)*W3[t-1]*(Ni2[t]+Ni2[t-1]))/(1+(dt/2)*W4[t-1])
    Ni0[t] = Nat-Ni1[t]-Ni2[t]-Ni3[t]
    ne[t] = 1*Ni1[t]+2*Ni2[t]+3*Ni3[t]
#%%
start = 0
stop = SIZE
mplot.plot(T[start:stop],W1[start:stop])
#%%
start = 2500
stop = 7500
mplot.plot((T[start:stop]/OMEGA)*1e15-100,ne[start:stop]/Nat,'r')
#mplot.plot(T[start:stop],E[start:stop])
mplot.plot((T[start:stop]/OMEGA)*1e15-100,(Eenv[start:stop]/E0)**2,'g')

mplot.ylabel('E^2,ne (units of E0, Nat)')
mplot.xlabel('t (fs)')

tspanplasma = T[stop]-T[start]
tspan = punit.treal(tspanplasma,OMEGA)
print(tspan)

#%% Ne in real units

neREAL = punit.nreal(ne,OMEGA)
TREAL = punit.treal(T,OMEGA)
EREALenv = punit.Ereal(Eenv,OMEGA)
IREAL = epsilon*c*EREALenv**2/2

mplot.plot(TREAL,neREAL/NatREAL)
mplot.plot(TREAL,IREAL/I0)








