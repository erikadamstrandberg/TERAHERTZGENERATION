# -*- coding: utf-8 -*-
"""
Created on Thu Feb  8 13:42:57 2018

@author: Hanna
"""

import numpy as np
import matplotlib.pyplot as mplot
import Ionization as ion
#%%
zmin = 0
zmax = 600
dz = 1
tmin = 0
tmax = 1000
dt = 1

Z = np.arange(zmin,zmax,dz)
T = np.arange(tmin,tmax,dt)

ne = np.zeros((len(T),len(Z)))      # matris med elektrontäthet vid varje tid
E = np.zeros((len(T),len(Z)))       # E-fält vid varje tid

Q = 18                              # max antal joniserade elektroner från en atom
nat = 1                             # homogen atomtäthet
N = np.zeros((Q,len(Z)))            # matris med jontätheterna

# sätt startjontäthet och E-fält
for z in range(len(z)):
    N(0,z) = nat
    E(0,z) = 1
    
for t in range(1,len(T)):
    
    for q in range(1,Q):  
        for z in range(len(Z)):
            N(0,z) = N(0,z)-N(q,z)
       
    for q in range(1,Q):
        for z in range(len(Z)):
            N(q,z) = # f(N(q,z), E(t,z), N(q-1,z), N(q-1,z,t-1))...??
            
#%% testar Q=1
zmin = 0
zmax = 600
dz = 1
tmin = 0
tmax = 700
dt = 1

Z = np.arange(zmin,zmax,dz)
T = np.arange(tmin,tmax,dt)

ne = np.zeros((len(T),len(Z)))      # matris med elektrontäthet vid varje tid
W1 = np.zeros((len(T),len(Z)))      # E-fält vid varje tid
W2 = np.zeros((len(T),len(Z))) 
Q = 1
N = np.zeros((Q*2+2,len(Z)))
ne = np.zeros((len(T),len(Z)))
nat = 1

# sätter startjontäthet och E-fält
for z in range(len(Z)):
    N[0,z] = nat
    N[1,z] = nat
    W1[0,z] = 1
    W2[0,z] = 0.5
    
for t in range(len(T)):
    for z in range(len(Z)):
        N[3,z] = (N[2,z]*(1-(dt/dz)*W2[t,z])+(dt/2)*W1[t,z]*(N[1,z]+N[0,z]))/(1+(dt/2)*W2[t,z])
        
        
#%%
zmin = 0
zmax = 1
dz = 1
tmin = 0
tmax = 100
dt = 0.01
FREK = 1e14

Z = np.arange(zmin,zmax,dz)
T = np.arange(tmin,tmax,dt)

ne = np.zeros((len(T),len(Z)))  #elektroner
Ni0 = np.zeros((len(T),len(Z))) #joner laddning +0
Ni1 = np.zeros((len(T),len(Z))) #joner laddning +1
Ni2 = np.zeros((len(T),len(Z))) #joner laddning +2
E = np.zeros((len(T),len(Z)))

W = np.zeros((len(T),len(Z)))
Nat = 1
K = 0.01

for z in range(len(Z)):
    Ni0[0,z] = Nat
    
for t in range(len(T)):
    E[t,0] = 0.1*np.cos(0.3*t*dt)
for t in range(len(T)):
    W[t,0] = ion.Landau(E[t,z],FREK)

for t in range(1,len(T)):
    for z in range(len(Z)):
        Ni0[t,z] = Nat-Ni1[t-1,z]-Ni2[t-1,z]
        
        Ni1[t,z] = (Ni1[t-1,z]*(1-(dt/dz)*W[t-1,z])+(dt/2)*W[t-1,z]*(Ni0[t,z]+Ni0[t-1,z]))/(1+(dt/2)*W[t-1,z])
        Ni2[t,z] = (Ni2[t-1,z]*(1-(dt/dz)*W[t-1,z])+(dt/2)*W[t-1,z]*(Ni1[t,z]+Ni1[t-1,z]))/(1+(dt/2)*W[t-1,z])
        
        ne[t,z] = 1*Ni1[t,z]+2*Ni2[t,z]

#%%
mplot.plot(T,ne[:,0])
mplot.plot(T,1e-83*E[:,0])
#%%
T = np.arange(100)
N = (Nat/2)*(1-np.exp(-2*W*T))
mplot.plot(T,N)