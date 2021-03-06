# -*- coding: utf-8 -*-
"""
Created on Tue Mar 27 12:40:30 2018

@author: Hanna
"""
#%%
import numpy as np
import scipy.constants as const
import Plasmaunit as punit
import Laser
import matplotlib.pyplot as mplot
import Rampfunctions 
#%%
TIMEINTERVAL = 600  # timeinterval in plasmaseconds
SPACEINTERVAL = 600 # spaceinterval in plasmameter

dt = 0.1
dz = 0.1

TIME = int(TIMEINTERVAL/dt)
SIZE = int(SPACEINTERVAL/dz)
dim = [TIME,SIZE]

T = np.arange(TIME)*dt
Z = np.arange(SIZE)*dz

E = np.zeros(SIZE)
B = np.zeros(SIZE)
J = np.zeros(SIZE)

E_tot = np.zeros(dim)
B_tot = np.zeros(dim)
J_tot = np.zeros(dim)

#E_timemean = np.zeros(SIZE)
#E_temp = np.zeros(SIZE)
ne  = np.zeros(SIZE)
Nat = np.ones(SIZE)

c = const.speed_of_light
epsilon = const.epsilon_0 
LAMBDA = 800e-9 
f = c/LAMBDA
OMEGAREAL = 2*np.pi*f
OMEGAPRIM = 1                       # this is the plasma omega, use this everywhere in the code
OMEGA_0 = OMEGAREAL/OMEGAPRIM       # this is the arbitrary omega, use this as argument in punits
t0REAL = 5e-15 
I0 = 4e18    
E0REAL = np.sqrt(2*I0/(epsilon*c))
E0 = punit.Eplasma(E0REAL,OMEGA_0)
t0 = punit.tplasma(t0REAL,OMEGA_0)
NatREAL = 3e25                      # real atomdensity of the argon gas
Natpunit = punit.nplasma(NatREAL,OMEGA_0)
nu = 0

PULSESTART = int(100/dz)              # at what index the pulse starts
PULSELENGTH = int(400/dz)           # pulselength in 'index units'
PLASMASTART = int(200/dz)#PULSELENGTH+PULSESTART
PLASMASTOP = PLASMASTART+6*PULSELENGTH
RAMP_DAMP = 0.001
RAMPLENGTH = 500

Laser.Gauss_forward(E,B,E0,PULSELENGTH,PULSESTART,OMEGAPRIM,t0,dz)
#Rampfunctions.Ramp_exp(PLASMASTART,PLASMASTOP,RAMP_DAMP,Natpunit,Nat,SIZE,dt)
#Rampfunctions.Ramp_linear_Erik(PLASMASTART,PLASMASTOP,RAMPLENGTH,Natpunit,Nat,SIZE,dz)
Nat = np.ones(SIZE)*Natpunit
ne = Nat

start = int(0/dz)
stop = int(600/dz)
mplot.plot(Z[start:stop],ne[start:stop])
mplot.plot(Z[start:stop],E[start:stop])

#U_E = np.zeros(TIME)
#U_B = np.zeros(TIME)
#U_kin = np.zeros(TIME)
#Total_energy = np.zeros(TIME)
#Energy_start = (1/2)*np.sum(E**2)+(1/2)*np.sum(B**2)

#kvot_kin = np.zeros(SIZE)

#%%

for t in range(1,TIME-1):
    E = E-(dt/dz)*(B-np.roll(B,1))-dt*J     # OBS!! this is wrong should be J[z]   
    B = B-(dt/dz)*(np.roll(E,-1)-E)
    J = J + dt*ne*E
    
    #kvot_kin[z] = J[z]**2/(2*(ne[z]+1e-100))
  
    #E_timemean = (1/2)*(E+E_temp)
    #E_temp = E
    #U_kin[t] = np.sum(kvot_kin)
    #Total_energy[t] = (U_E[t]+U_B[t]+U_kin[t])/Energy_star
    
    #U_E[t] = (1/2)*np.sum(E_timemean**2)
    #U_B[t] = (1/2)*np.sum(B**2)t
        
    E_tot[t] = E
    B_tot[t] = B
    #J_tot[t] = J
        
#%%
J_E = J_tot*E_tot

t = 250
time = int(t/dt)
start = int(0/dz)
stop = int(600/dz)

#mplot.plot(Z[start:stop],Nat[start:stop])
mplot.plot(Z[start:stop],E_tot[time,start:stop],'b')
#mplot.plot(Z[start:stop],B_tot[time,start:stop],'r')
#mplot.plot(Z[start:stop],J_tot[time,start:stop],'g')
mplot.plot(Z[start:stop],J_E[time,start:stop]*1e2)

mplot.grid()
mplot.xlabel('z (plasmaunits)')
mplot.ylabel('E, Ne (plasmaunits)')
#%%
start = int(1/dt)
stop = int(160/dt)
#mplot.plot(T[start:stop],Total_energy[start:stop])

mplot.xlabel('t (plasmaunits)')
mplot.ylabel('E/E_0')



