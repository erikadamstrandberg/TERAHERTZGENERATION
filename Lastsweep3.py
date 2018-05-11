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
import Micro_i as Mic

dt = 0.0495#495 
dz = 0.05#5
nu = 0

TIME = int(4000/dt)
SIZE = int(5000/dz)
PULSELENGTH = int(1250/dz)
PULSESTART = int(1000/dz)

diff = 1e18
start = 1e18
stopp = 1e20+diff
I0 = np.arange(start,stopp,diff)

Laser_input = np.zeros(len(I0))
ne_nat = np.zeros(len(I0))
THz_power = np.zeros(len(I0))

for P in range(len(I0)):
    dim = [TIME,SIZE]
    plott = np.arange(TIME)
    plotz = np.arange(SIZE)
    E = np.zeros(SIZE)
    B = np.zeros(SIZE)
    J = np.zeros(SIZE)
    W1 = np.zeros(SIZE)
    W2 = np.zeros(SIZE)
    W3 = np.zeros(SIZE)
    W4 = np.zeros(SIZE)
    W5 = np.zeros(SIZE)
    W6 = np.zeros(SIZE)
    W7 = np.zeros(SIZE)
    W8 = np.zeros(SIZE)
    W9 = np.zeros(SIZE)
    W10 = np.zeros(SIZE)
    W11 = np.zeros(SIZE)
    W12 = np.zeros(SIZE)
    W13 = np.zeros(SIZE)
    W14 = np.zeros(SIZE)
    W15 = np.zeros(SIZE)
    W16 = np.zeros(SIZE)
    W17 = np.zeros(SIZE)
    W18 = np.zeros(SIZE)
    Ni18 = (np.zeros(SIZE),np.zeros(SIZE))
    Ni17 = (np.zeros(SIZE),np.zeros(SIZE))
    Ni16 = (np.zeros(SIZE),np.zeros(SIZE))
    Ni15 = (np.zeros(SIZE),np.zeros(SIZE))
    Ni14 = (np.zeros(SIZE),np.zeros(SIZE))
    Ni13 = (np.zeros(SIZE),np.zeros(SIZE))
    Ni12 = (np.zeros(SIZE),np.zeros(SIZE))
    Ni11 = (np.zeros(SIZE),np.zeros(SIZE))
    Ni10 = (np.zeros(SIZE),np.zeros(SIZE))
    Ni9 = (np.zeros(SIZE),np.zeros(SIZE))
    Ni8 = (np.zeros(SIZE),np.zeros(SIZE))
    Ni7 = (np.zeros(SIZE),np.zeros(SIZE))
    Ni6 = (np.zeros(SIZE),np.zeros(SIZE))
    Ni5 = (np.zeros(SIZE),np.zeros(SIZE))
    Ni4 = (np.zeros(SIZE),np.zeros(SIZE))
    Ni3 = (np.zeros(SIZE),np.zeros(SIZE))
    Ni2 = (np.zeros(SIZE),np.zeros(SIZE))
    Ni1 = (np.zeros(SIZE),np.zeros(SIZE))
    Ni0 = (np.ones(SIZE),np.zeros(SIZE))
    ne = (np.zeros(SIZE),np.zeros(SIZE))
    Z = np.arange(SIZE)*dz
    T = np.arange(TIME)*dt
    Etera = np.zeros(TIME)
    
    c = const.speed_of_light
    epsilon = const.epsilon_0 
    LAMBDA = 800e-9 
    f = c/LAMBDA
    OMEGAREAL = 2*np.pi*f
    OMEGAPRIM = 1                       # this is the plasma omega, use this everywhere in the code
    OMEGA_0 = OMEGAREAL/OMEGAPRIM       # this is the arbitrary omega, use this as argument in punits
    t0REAL = 15e-15
    NatREAL = 2.7e25
    E0REAL = np.sqrt(2*I0[P]/(epsilon*c))
    
    E0 = punit.Eplasma(E0REAL,OMEGA_0)
    t0 = punit.tplasma(t0REAL,OMEGA_0)
    
    xi = 0.3
    phi = 0
    
    El = np.zeros(PULSELENGTH)
    Bl = np.zeros(PULSELENGTH)
    START = -PULSELENGTH/2  
    STOPP = PULSELENGTH/2
    t = dz*np.arange(START,STOPP,1)
    for i in range(len(El)):
        El[i] = E0*(np.sqrt(1-xi)*np.sin(OMEGAPRIM*t[i])*np.exp(-(t[i]**2/(2*t0**2)))+np.sqrt(xi)*np.sin(2*OMEGAPRIM*t[i]+phi)*np.exp(-(t[i]**2/(t0**2))))
    t = dz*np.arange(START,STOPP,1)+(dz-dt)/2
    for i in range(len(Bl)):
        Bl[i] = E0*(np.sqrt(1-xi)*np.sin(OMEGAPRIM*t[i])*np.exp(-(t[i]**2/(2*t0**2)))+np.sqrt(xi)*np.sin(2*OMEGAPRIM*t[i]+phi)*np.exp(-(t[i]**2/(t0**2))))
    E[PULSESTART:PULSESTART+PULSELENGTH] = El
    B[PULSESTART:PULSESTART+PULSELENGTH] = Bl
    
    
    Natpunit = punit.nplasma(NatREAL,OMEGA_0)
    Nat = np.ones(SIZE)*Natpunit
    
    PLASMASTART = PULSELENGTH+PULSESTART
    PLASMASTOPP = PLASMASTART+Mic.Micro_i(7e-6,Z,OMEGA_0)
    RAMP_DAMP = 2
    
    Rampfunctions.Ramp_exp(PLASMASTART,PLASMASTOPP,RAMP_DAMP,Natpunit,Nat,SIZE,dz)
    Nkritisk = punit.nreal(1,OMEGA_0)
    
    collect_i = PLASMASTOPP+Mic.Micro_i(3e-6,Z,OMEGA_0)#PLASMASTOPP+Mic.Micro_i(1e-6,Z,OMEGA_0)
    
    nuREAL = 100e-15
    nu = punit.omegaplasma(nuREAL,OMEGA_0)
    
    Eft_L = np.fft.fft(E)*2*dt
    Laser_input[P] = sum(np.abs(Eft_L)**2)
    
    Total_Nat = sum(Nat)
    
    SAVE_EVERY_T = 10000
    SAVE_dim = int(TIME/SAVE_EVERY_T)
    dim = [SAVE_dim,SIZE]
    Etot = np.zeros(dim)
    netot = np.zeros(dim)
    
    print(TIME)
    SAVE_EVERY_count = 0
        
    for i in range(1,TIME):
        E = SpaceSolver.E(E,B,J,dt,dz)
        B = SpaceSolver.B(E,B,dt,dz)
        ne = SpaceSolver.N_1e20(E,Nat,Ni0,Ni1,Ni2,Ni3,Ni4,Ni5,Ni6,Ni7,Ni8,ne,W1,W2,W3,W4,W5,W6,W7,W8,OMEGA_0,dt)
        J = SpaceSolver.J(E,J,ne,nu,dt,dz)
        E[0] = 0
        B[SIZE-1] = 0
        
        Etera[i-1] = E[collect_i]
        
        if (i%SAVE_EVERY_T) == 0: 
            Etot[SAVE_EVERY_count] = E
            netot[SAVE_EVERY_count] = ne[0] 
            SAVE_EVERY_count = SAVE_EVERY_count + 1
            
    N = len(Etera)
    Efrekreal = np.fft.fftfreq(N)*OMEGA_0/dt
    THzmax = 30e12
    THz_cut = np.nonzero(Efrekreal>THzmax)
    THz_cut_i = int(THz_cut[0][0])
    Eft = np.fft.fft(Etera)*dt*2
    THz_power[P] = sum(np.abs(Eft[0:THz_cut_i])**2)
    
    Total_ne = sum(ne[0])
    ne_nat[P] = Total_ne/Total_Nat
    ne_max = np.amax(ne[0])#ne[0][collect_i]
    
np.savetxt("THz_power_1e18_1e20_7",THz_power,delimiter=',')
np.savetxt("I0_1e18_1e20_7",I0,delimiter=',')
np.savetxt("Laser_input_1e18_1e20_7",Laser_input,delimiter=',')
np.savetxt("ne_nat__1e18_1e20_7",ne_nat,delimiter=',')
