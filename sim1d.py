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

def runsim(
        dt = 1,               # Time step
        dz = 1,               # Space step
        time = 200,           # Temporal length of simulation window in units of dt
        size = 400,           # Spatial length of simulation window in units of dz
        pulsestart = 100,     # Start position of pulse
        pulselength = 100,    # Spatial length of pulse
        nu = 0,               # Collision rate
        wavelength = 800e-9,  # Laser wavelength
        ramplength = 1/10,    # Electron density initial ramp length as a fraction of size
        rampdamp = 1/10,      # Electron density ramp parameter, (1-exp(-rampdamp*z))
        plasmastopp = deepcopy(size),
        plottime = int(deepcopy(time)/2)): 
    
    dim = [time,size]
    
    # Initialise the various fields to be calculated
    E = np.zeros(size)
    B = np.zeros(size)
    J = np.zeros(size)
    W1 = np.zeros(size)
    W2 = np.zeros(size)
    W3 = np.zeros(size)
    Ni2 = np.zeros(size)
    Ni2temp = np.zeros(size)
    Ni1 = np.zeros(size)
    Ni1temp = np.zeros(size)
    Ni0 = np.ones(size)
    Ni0temp = np.zeros(size)
    ne = np.zeros(size)
    netemp = np.zeros(size)
    
    Ni0tot = np.zeros(dim)
    Ni1tot = np.zeros(dim)
    netot = np.zeros(dim)
    Etot = np.zeros(dim)
    Btot = np.zeros(dim)
    Jtot = np.zeros(dim)
    ntot = np.zeros(dim)
    W1tot = np.zeros(dim)

    LIGHTSPEED = const.speed_of_light
    EPSILON = const.epsilon_0 

    f = LIGHTSPEED/wavelength
    omegareal = 2*np.pi*f
    OMEGAPRIM = 1                       # this is the plasma omega, use this everywhere in the code
    omega_0 = omegareal/OMEGAPRIM       # this is the arbitrary omega, use this as argument in punits
    T0REAL = 50e-16 
    I0 = 4e18 
    NatREAL = 7e26
    E0REAL = np.sqrt(2*I0/(EPSILON*LIGHTSPEED))
    plasmastart = pulselength+pulsestart    # where the electron density starts
    plasmastopp = size                      # Where it stops. Note: for PLASMASTOPP = SIZE the plasma extents all the way to the end of the simlation window 
    E0 = punit.Eplasma(E0REAL,omega_0)
    t0 = punit.tplasma(T0REAL,omega_0)
    Laser.Gauss_forward(E,B,E0,pulselength,pulsestart,OMEGAPRIM,t0,dt) # Sets up the laser pulse in the window
    
    Natpunit = punit.nplasma(NatREAL,omega_0)
    Nat = np.ones(size)*Natpunit
    Rampfunctions.Ramp_exp(plasmastart,plasmastopp,rampdamp,Natpunit,Nat,size,dt) # Creates a ramp for the electron density

    bar = ChargingBar('Simulation running', max = time)
    print(str(datetime.now())+': Beginning simulation.')

    for i in range(1,time):    
        E = SpaceSolver.E(E,B,J,dt,dz)
        B = SpaceSolver.B(E,B,dt,dz)   
        ne = SpaceSolver.N(E,Nat,Ni0,Ni1,Ni2,Ni0tot[i-1],Ni1tot[i-1],ne,W1,W2,W3,omega_0,dt)
        J = SpaceSolver.J(E,J,ne,netot[i-1],nu,dt,dz)

        Etot[i] = E
        Btot[i] = B
        Jtot[i] = J
        Ni0tot[i] = Ni0
        Ni1tot[i] = Ni1
        netot[i] = ne
        bar.next()
        
    bar.next()
    bar.finish()
    
    print(str(datetime.now())+': Simulation complete.')

    z = np.arange(len(Etot[0]))
    plotnsave(z, Etot[plottime],)
    mplot.clf()
