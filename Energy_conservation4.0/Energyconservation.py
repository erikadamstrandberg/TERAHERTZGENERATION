# -*- coding: utf-8 -*-
"""
Created on Mon Mar  5 18:42:14 2018

NOTE: program for checking stuff for the run on Illias

@author: Hanna
"""
#%%
import numpy as np
import scipy.constants as const
import Plasmaunit as punit
import Laser
import matplotlib.pyplot as mplot
import Rampfunctions
import SpaceSolver
from datetime import datetime
#%%
TIMEINTERVAL = 400  # timeinterval in plasmaseconds
SPACEINTERVAL = 400 # spaceinterval in plasmameter

dt = 0.5
dz = 0.5

TIME = int(TIMEINTERVAL/dt)
SIZE = int(SPACEINTERVAL/dz)
dim = [TIME,SIZE]

T = np.arange(TIME)*dt
Z = np.arange(SIZE)*dz

# Create matrices
E = np.zeros(SIZE)
Etemp = np.zeros(SIZE)
B = np.zeros(SIZE)
J = np.zeros(SIZE)
W1 = np.zeros(SIZE)
W2 = np.zeros(SIZE)
W3 = np.zeros(SIZE)
Ni2 = np.zeros(SIZE)
Ni2temp = np.zeros(SIZE)
Ni1 = np.zeros(SIZE)
Ni1temp = np.zeros(SIZE)
Ni0 = np.ones(SIZE)
Ni0temp = np.zeros(SIZE)
ne = np.zeros(SIZE)
netemp = np.zeros(SIZE)
#Ni0tot = np.zeros(dim)
#Ni1tot = np.zeros(dim)
netot = np.zeros(dim)
#Etot = np.zeros(dim)
#Btot = np.zeros(dim)
#Jtot = np.zeros(dim)
#ntot = np.zeros(dim)
#W1tot = np.zeros(dim)

# Constants 
nu = 0   
c = const.speed_of_light
epsilon = const.epsilon_0 
LAMBDA = 800e-9 
f = c/LAMBDA
OMEGAREAL = 2*np.pi*f
OMEGAPRIM = 1                       # this is the plasma omega, use this everywhere in the code
OMEGA_0 = OMEGAREAL/OMEGAPRIM       # this is the arbitrary omega, use this as argument in punits
t0REAL = 5e-15 
I0 = 4e18                           # real peak intensity of pulse
NatREAL = 3e25                      # real atomdensity of the argon gas
E0REAL = np.sqrt(2*I0/(epsilon*c))
E0 = punit.Eplasma(E0REAL,OMEGA_0)
t0 = punit.tplasma(t0REAL,OMEGA_0)
Natpunit = punit.nplasma(NatREAL,OMEGA_0)
Nat = np.zeros(SIZE)*Natpunit        # OBS should be zeros!?

PULSESTART = int(0/dz)              # at what index the pulse starts
PULSELENGTH = int(100/dz)           # pulselength in 'index units'
PLASMASTART = PULSELENGTH+PULSESTART
PLASMASTOP = PLASMASTART+4*PULSELENGTH
RAMP_DAMP = 0.1

Laser.Gauss_forward(E,B,E0,PULSELENGTH,PULSESTART,OMEGAPRIM,t0,dz)
Rampfunctions.Ramp_exp(PLASMASTART,PLASMASTOP,RAMP_DAMP,Natpunit,Nat,SIZE,dt)

start = int(0/dz)
stop = int(SIZE/dz)
mplot.plot(Z[start:stop],E[start:stop])
mplot.plot(Z[start:stop],Nat[start:stop])

mplot.xlabel('Z (plasmaunits)')
mplot.ylabel('E (plasmau.), Nat (plasmau.)')
#mplot.savefig('start_setting_run2.pdf')

E_Efield = np.zeros(TIME)
E_Bfield = np.zeros(TIME)
E_kin = np.zeros(TIME)
U_th = np.zeros(TIME)
Total_energy = np.zeros(TIME)
kvot = np.zeros(SIZE)
kvot2 = np.zeros(SIZE)

Energy_start = (1/2)*np.sum(E**2)+(1/2)*np.sum(B**2)
#%%    
print(str(datetime.now())+': Beginning simulation.')
for i in range(1,TIME):
    E = SpaceSolver.E(E,B,J,dt,dz)
    B = SpaceSolver.B(E,B,dt,dz)   
    ne = SpaceSolver.N(E,Nat,Ni0,Ni1,Ni2,Ni0temp,Ni1temp,ne,W1,W2,W3,OMEGA_0,dt)
    Ni0temp = Ni0
    Ni1temp = Ni1
    J = SpaceSolver.J(E,J,ne,netemp,nu,dt,dz)

    #Etot[i] = E
    #Btot[i] = B
    #Jtot[i] = J
    netot[i] = ne
    for z in range(SIZE):
        if ne[z]**2 == 0:
            kvot[z] = 0
        else:
            kvot[z] = J[z]**2/ne[z]
            kvot2[z] = (1/dt)*(ne[z]-netemp[z]) #*(J[z]**2)/(ne[z]**2) OBS att (netot[i,z]-netot[i-1,z]) funkar
    netemp = ne
    U_th[i] = U_th[i-1] + np.sum(kvot2)         
    E_kin[i] = (1/2)*np.sum(kvot)
    E_Efield[i] = (1/2)*np.sum(E**2)
    E_Bfield[i] = (1/2)*np.sum(B**2)
    Total_energy[i] = (E_Efield[i]+E_Bfield[i]+E_kin[i]+U_th[i])/Energy_start

print(str(datetime.now())+': Simulation complete.')


#%%
start = int(0/dt)
stop = int(TIME-1)
mplot.plot(T[start:stop],Total_energy[start:stop])

print('Total energy at stop is:')
print(Total_energy[TIME-1])

#np.savetxt("Total_energy",Total_energy, delimiter='.')
#mplot.savefig('Total_energy_04.pdf')