import numpy as np
import scipy.constants as const
import matplotlib as matlib
matlib.use('Agg')                    # This allows the standalone application to plot and save figs.
import matplotlib.pyplot as mplot
from pylab import *
import Laser
import Ionization as Ion
import SpaceSolver
import Plasmaunit as punit
import Rampfunctions
from plotnsave import plotnsave
from progress.bar import ChargingBar
from datetime import datetime
from copy import deepcopy

def main():
    dt = 0.1
    dz = 0.1
    nu = 0
    
    TIME = 14000
    SIZE = 22000
    PULSELENGTH = 10000
    
    dim = [TIME,SIZE]
    E = np.zeros(SIZE)
    B = np.zeros(SIZE)
    J = np.zeros(SIZE)
    W1 = np.zeros(SIZE)
    W2 = np.zeros(SIZE)
    W3 = np.zeros(SIZE)
    Ni2 = (np.zeros(SIZE),np.zeros(SIZE))
    Ni1 = (np.zeros(SIZE),np.zeros(SIZE))  # The second array is used for saving values one time loop
    Ni0 = (np.ones(SIZE),np.zeros(SIZE))
    ne = (np.zeros(SIZE),np.zeros(SIZE))
    Ni0tot = np.zeros(dim)
    Ni1tot = np.zeros(dim)
    netot = np.zeros(dim)
    Etot = np.zeros(dim)
    Btot = np.zeros(dim)
    Jtot = np.zeros(dim)
    ntot = np.zeros(dim)
    
    c = const.speed_of_light
    epsilon = const.epsilon_0 
    LAMBDA = 800e-9 
    f = c/LAMBDA
    PULSESTART = 4000
    OMEGAREAL = 2*np.pi*f
    OMEGAPRIM = 1                       # this is the plasma omega, use this everywhere in the code
    OMEGA_0 = OMEGAREAL/OMEGAPRIM       # this is the arbitrary omega, use this as argument in punits
    t0REAL = 50e-15
    I0 = 150e16 
    NatREAL = 3e25
    E0REAL = np.sqrt(2*I0/(epsilon*c))
    
    E0 = punit.Eplasma(E0REAL,OMEGA_0)
    t0 = punit.tplasma(t0REAL,OMEGA_0)
    
    xi = 0.1
    phi = np.pi/2
    Laser.TwoC_forward(E,B,E0,xi,phi,PULSELENGTH,PULSESTART,OMEGAPRIM,t0,dt)
    
    Natpunit = punit.nplasma(NatREAL,OMEGA_0)
    Nat = np.ones(SIZE)*Natpunit
    
    PLASMASTART = PULSELENGTH+PULSESTART
    PLASMASTOPP = SIZE
    RAMP_DAMP = 0.1
    
    Rampfunctions.Ramp_exp(PLASMASTART,PLASMASTOPP,RAMP_DAMP,Natpunit,Nat,SIZE,dz)
    
    Sample1real = 3e-6
    Sample1 = int(punit.splasma(Sample1real,OMEGA_0))
    
    Sample2real = 10e-6
    Sample2 = int(punit.splasma(Sample2real,OMEGA_0))
    
    Etera1 = np.zeros(TIME)
    Etera2 = np.zeros(TIME)

    bar = ChargingBar('Simulation running', max = TIME)
    print(str(datetime.now())+': Beginning simulation.')

    for i in range(1,TIME):
        # Calculate all fields for current time
        E = SpaceSolver.E(E,B,J,dt,dz)
        B = SpaceSolver.B(E,B,dt,dz)   
        ne = SpaceSolver.N(E,Nat,Ni0,Ni1,Ni2,Ni0tot[i-1],Ni1tot[i-1],ne,W1,W2,W3,OMEGA_0,dt)
        J = SpaceSolver.J(E,J,ne,netot[i-1],nu,dt,dz)

        # Save current time
        Etera1[i-1] = E[PLASMASTART+Sample1/dt]
        Etera2[i-1] = E[PLASMASTART+Sample2/dt]
        Ni0tot[i] = Ni0
        Ni1tot[i] = Ni1
        netot[i] = ne
        bar.next()
        
    bar.next()
    bar.finish()
    
    print(str(datetime.now())+': Simulation complete.')
    
    z = np.arange(len(Etera1))
    plotnsave.plotnsave(z, Etera1, filename = 'Sample1')
    plotnsave.plotnsave(z, Etera2, filename = 'Sample1')

def plog(msg):
    print(str(datetime.now()) + ': ' + msg)
