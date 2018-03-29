# -*- coding: utf-8 -*-
"""
Created on Wed Mar 28 16:42:08 2018

@author: Hanna
"""
#%%
import numpy as np
import scipy.constants as const
import Plasmaunit as punit
import Laser
import matplotlib.pyplot as mplot
import Rampfunctions 
import Ionization as Ion

TIMEINTERVAL = 600  # timeinterval in plasmaseconds
SPACEINTERVAL = 800 # spaceinterval in plasmameter

dz = 0.20
dt = 0.9*dz

TIME = int(TIMEINTERVAL/dt)
SIZE = int(SPACEINTERVAL/dz)
dim = [TIME,SIZE]

T = np.arange(TIME)*dt
Z = np.arange(SIZE)*dz

E = np.zeros(SIZE)
B = np.zeros(SIZE)
J = np.zeros(SIZE)

#E_tot = np.zeros(dim)
#B_tot = np.zeros(dim)
#J_tot = np.zeros(dim)

ne  = np.zeros(SIZE)
ne_temp  = np.zeros(SIZE)
Nat = np.ones(SIZE)

Ni0 = np.zeros(SIZE)
Ni1 = np.zeros(SIZE)
Ni2 = np.zeros(SIZE)
Ni0_temp = np.zeros(SIZE)
Ni1_temp = np.zeros(SIZE)
W1 = np.ones(SIZE)
W2 = np.ones(SIZE)
W3 = np.ones(SIZE)

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

PULSESTART = int(50/dz)              # at what index the pulse starts
PULSELENGTH = int(200/dz)            # pulselength in 'index units'
PLASMASTART = PULSELENGTH+PULSESTART
PLASMASTOP = SIZE-10
RAMP_DAMP = 0.001
RAMPLENGTH = int(50/dz)

Laser.Gauss_forward_Hanna(E,B,E0,PULSELENGTH,PULSESTART,OMEGAPRIM,t0,dt,dz)
Rampfunctions.Ramp_linear_Erik(PLASMASTART,PLASMASTOP,RAMPLENGTH,Natpunit,Nat,SIZE,dz)
Ni0 = Nat


start = int(0/dz)
stop = int(1000/dz)
mplot.plot(Z[start:stop],Nat[start:stop])
mplot.plot(Z[start:stop],E[start:stop])

########### Energy stuff ################################
U_E = np.zeros(TIME)
U_B = np.zeros(TIME)
U_kin = np.zeros(TIME)
U_th = np.zeros(TIME)
Energy_start = (1/2)*np.sum(E**2)+(1/2)*np.sum(B**2)

E_timemean = np.zeros(SIZE)
E_temp = np.zeros(SIZE)
B_spacemean = np.zeros(SIZE)
kvot_kin = np.zeros(SIZE)
kvot_th = np.zeros(SIZE)
#########################################################

c1 = dt/dz
c2 = dt/2
#%%

for t in range(1,TIME-1):
    E = E-c1*(B-np.roll(B,1))-dt*J
    
    ########## Energy stuff ##########################################
    kvot_kin = J**2/(2*(ne+1e-100))
    kvot_th = (W1*Ni0_temp+W2*Ni1_temp)*(J**2)/(2*(ne+1e-100)**2)
    E_timemean = (1/2)*(E+E_temp)
    E_temp = E
    B_spacemean = (1/2)*(B+np.roll(B,1))
    
    U_kin[t] = np.sum(kvot_kin)
    U_th[t] = dt*(U_th[t-1]+np.sum(kvot_th))
    U_E[t] = (1/2)*np.sum(E_timemean**2)
    U_B[t] = (1/2)*np.sum(B_spacemean**2)
    ##################################################################
    
    B = B-c1*(np.roll(E,-1)-E)
    J = ((1-nu*c2)*J+c2*(ne+ne_temp)*E)/(1+nu*c2)
    ne_temp = ne
    
    W1 = Ion.Landau_array_roll(E,OMEGA_0,1)
    W2 = Ion.Landau_array_roll(E,OMEGA_0,2)
    #W3 = Ion.Landau_array_roll(E,OMEGA_0,3)
    
    Ni0 = Nat-Ni1-Ni2
    Ni1 = (Ni1*(1-(dt/2)*W2)+(dt/2)*W1*(Ni0+Ni0_temp))/(1+(dt/2)*W2)
    Ni2 = (Ni2*(1-(dt/2)*W3)+(dt/2)*W2*(Ni1+Ni1_temp))/(1+(dt/2)*W3)
    ne = 1*Ni1 + 2*Ni2
    Ni0_temp = Ni0
    Ni1_temp = Ni1
        
    #E_tot[t] = E
    #B_tot[t] = B
    #J_tot[t] = J
        

#%%
t = 550
time = int(t/dt)
start = int(0/dz)
stop = int(800/dz)

#mplot.plot(Z[start:stop],Nat[start:stop])
#mplot.plot(Z[start:stop],E_tot[time,start:stop],'b')
#mplot.plot(Z[start:stop],B_tot[time,start:stop],'r')
#mplot.plot(Z[start:stop],J_tot[time,start:stop]*1e1,'g')
#mplot.plot(Z[start:stop],J_E[time,start:stop]*1e2)

mplot.xlabel('z (plasmaunits)')
mplot.ylabel('E, Ne (plasmaunits)')

#mplot.savefig('pulse_t400.pdf')
#%%
Total_energy = U_E + U_B + U_kin + U_th

start = int(1/dt)
stop = int(t/dt)
mplot.plot(T[start:stop],Total_energy[start:stop]/Energy_start)
#mplot.plot(T[start:stop],U_th[start:stop])

mplot.xlabel(r'$t$ (plasmaunits)',fontsize=20)
mplot.ylabel(r'$E/E_0$',fontsize=20)

#mplot.savefig('Total_energy_t400.pdf')
#np.savetxt("Total_energy_dz20",Total_energy, delimiter='.')


