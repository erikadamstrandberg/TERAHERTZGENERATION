#%%
import numpy as np
from pylab import *
import scipy.constants as const
import matplotlib as matlib
matlib.use('Agg')                    # This allows the standalone application to plot and save figs.
import matplotlib.pyplot as mplot
import SpaceSolver
import Plasmaunit as punit
import Rampfunctions

dt = 0.0495#099
dz = 0.05#1
nu = 0

TIME = int(25000/dt)
SIZE = int(50000/dz)
PULSELENGTH = int(1250/dz)
PULSESTART = int(25000/dz)

E = np.zeros(SIZE)
B = np.zeros(SIZE)
J = np.zeros(SIZE)
W1 = np.zeros(SIZE)
W2 = np.zeros(SIZE)
W3 = np.zeros(SIZE)
Ni3 = (np.zeros(SIZE),np.zeros(SIZE))
Ni2 = (np.zeros(SIZE),np.zeros(SIZE))
Ni1 = (np.zeros(SIZE),np.zeros(SIZE))  # The second array is used for saving values one time loop
Ni0 = (np.ones(SIZE),np.zeros(SIZE))
ne = (np.zeros(SIZE),np.zeros(SIZE))
T = np.arange(TIME)*dt
Z = np.arange(SIZE)*dz

c = const.speed_of_light
epsilon = const.epsilon_0 
LAMBDA = 1e-6
f = c/LAMBDA
OMEGAREAL = 2*np.pi*f
OMEGAPRIM = 1                     # this is the plasma omega, use this everywhere in the code
OMEGA_0 = OMEGAREAL/OMEGAPRIM       # this is the arbitrary omega, use this as argument in punits
t0REAL = 50e-15
I0 = 150e16
NatREAL = 2.7e25
E0REAL = np.sqrt(2*I0/(epsilon*c))
E0 = punit.Eplasma(E0REAL,OMEGA_0)
t0 = punit.tplasma(t0REAL,OMEGA_0)

xi = 0.1
phi = np.pi/2

El = np.zeros(PULSELENGTH)
Bl = np.zeros(PULSELENGTH)
START = -PULSELENGTH/2  
STOPP = PULSELENGTH/2
t = dz*np.arange(START,STOPP,1)
for i in range(len(El)):
    El[i] = E0*(np.sqrt(1-xi)*np.cos(OMEGAPRIM*t[i])*np.exp(-2*np.log(2)*(t[i]/t0)**2)+np.sqrt(xi)*np.cos(2*OMEGAPRIM*t[i]+phi)*np.exp(-8*np.log(2)*(t[i]/t0)**2))
t = dz*np.arange(START,STOPP,1)+(dz-dt)/2
for i in range(len(Bl)):
    Bl[i] = E0*(np.sqrt(1-xi)*np.cos(OMEGAPRIM*t[i])*np.exp(-2*np.log(2)*(t[i]/t0)**2)+np.sqrt(xi)*np.cos(2*OMEGAPRIM*t[i]+phi)*np.exp(-8*np.log(2)*(t[i]/t0)**2))
E[PULSESTART:PULSESTART+PULSELENGTH] = El
B[PULSESTART:PULSESTART+PULSELENGTH] = Bl

print(E[PULSESTART])
print(E[PULSESTART+PULSELENGTH-1])
    
Natpunit = punit.nplasma(NatREAL,OMEGA_0)
Nat = np.ones(SIZE)*Natpunit

PLASMASTART = PULSELENGTH+PULSESTART
PLASMASTOPP = SIZE
RAMP_DAMP = 2

threemicroREAL = 3e-6
threemicro = punit.splasma(threemicroREAL,OMEGA_0)
threemicro_cut = np.nonzero(Z > threemicro)
threemicro_i = int(threemicro_cut[0][0])

#tenmicroREAL = 10e-6
#tenmicro = punit.splasma(tenmicroREAL,OMEGA_0)
#tenmicro_cut = np.nonzero(Z > tenmicro)
#tenmicro_i = int(tenmicro_cut[0][0])
#
#hunderedmicroREAL = 100e-6
#hunderedmicro = punit.splasma(hunderedmicroREAL,OMEGA_0)
#hunderedmicro_cut = np.nonzero(Z > hunderedmicro)
#hunderedmicro_i = int(hunderedmicro_cut[0][0])
#
thousandmicroREAL = 1000e-6
thousandmicro = punit.splasma(thousandmicroREAL,OMEGA_0)
thousandmicro_cut = np.nonzero(Z > thousandmicro)
thousandmicro_i = int(thousandmicro_cut[0][0])

PLASMASTOPP = SIZE
collect_i = int(PLASMASTART+thousandmicro_i)

Rampfunctions.Ramp_exp(PLASMASTART,PLASMASTOPP,RAMP_DAMP,Natpunit,Nat,SIZE,dz)

#Etera_three = np.zeros(TIME)
#Etera_ten = np.zeros(TIME)
#Etera_hundred = np.zeros(TIME)
Etera_thousand = np.zeros(TIME)

#mplot.plot(np.arange(len(E)),E)
#mplot.plot(np.arange(len(Nat)),Nat)
#
#x = np.array([collect_i,collect_i])
#y = np.array([0,Natpunit])
#mplot.plot(x,y)

print(TIME)
for i in range(1,TIME):
    if (i % 10000) == 0:
        print(i)
    # Calculate all fields for current time
    E = SpaceSolver.E(E,B,J,dt,dz)
    E[0] = 0
    B = SpaceSolver.B(E,B,dt,dz) 
    B[SIZE-1] = 0
    ne = SpaceSolver.N(E,Nat,Ni0,Ni1,Ni2,Ni3,ne,W1,W2,W3,OMEGA_0,dt)
    J = SpaceSolver.J(E,J,ne,nu,dt,dz)
    
    Etera_thousand[i-1] = E[collect_i]

    # Save current time

ne_e = ne[0][collect_i]

#Start = 0
#Stopp = 3000/dt
#mplot.plot(T[Start:Stopp],Etera_three[Start:Stopp])
N = len(Etera_thousand)

vikt = N*1e-4
Eft = np.fft.ifft(Etera_thousand)*vikt*2
Efrek = np.fft.fftfreq(N)*OMEGA_0/dt
#mplot.plot(Efrek,np.abs(Eft))
#
#mplot.axis([0,90e12,5e-7,1e-4])
#mplot.yscale('log')

THzmax = 30e12
THz_cut = np.nonzero(Efrek > THzmax)

THz_cut_i = int(THz_cut[0][0])
#x = np.array([Efrek[THz_cut_i],Efrek[THz_cut_i]])
#y = np.array([0,4])
#mplot.plot(x,y,'black',label=r'$30$ THz')
#
fp = punit.omegareal(np.sqrt(ne_e)/(2*np.pi),OMEGA_0)
#x = np.array([fp,fp])
#y = np.array([0,2])
#mplot.plot(x,y,'green',label=r'$\nu_p$ THz')
#
#x = np.array([Efrek[THz_cut_i],Efrek[THz_cut_i]])
#y = np.array([0,4])
#mplot.plot(x,y,'black',label=r'$30$ THz')
#
#mplot.axis([0,10e13,5e-7,1e-4])
#mplot.yscale('log')

np.savetxt("Per100micro_spectrum",Eft,delimiter=',')
np.savetxt("Per100micro_frek",Efrek,delimiter=',')
np.savetxt("Per100micro_ne_thz30,",np.array([fp,Efrek[THz_cut_i]]),delimiter=",")
